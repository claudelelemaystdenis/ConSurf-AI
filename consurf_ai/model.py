"""ConSurf-AI network (inference-only copy of the training code's AttentionNetwork) and weight loading."""
import hashlib
import os
import shutil
import urllib.request
from pathlib import Path

import torch
import torch.nn as nn
import torch.nn.functional as F

MAX_LEN = 2000  # the rotary-embedding table covers 2000 positions

WEIGHTS_FILE = "consurf_ai_3B_weights.pth"
WEIGHTS_SHA256 = "99c06aa4494a361ecfc1ec306d32109bef1d9b80d1b78be70c64d046ff15205f"
WEIGHTS_URL = f"https://github.com/claudelelemaystdenis/ConSurf-AI/raw/main/weights/{WEIGHTS_FILE}"


# ── Network ──────────────────────────────────────────────────────────────────
class RotaryPositionalEmbedding(nn.Module):
    def __init__(self, head_dim, max_len=MAX_LEN):
        super().__init__()
        inv_freq = 1.0 / (10000 ** (torch.arange(0, head_dim, 2).float() / head_dim))
        freqs = torch.einsum("i,j->ij", torch.arange(max_len).float(), inv_freq)
        emb = torch.cat((freqs, freqs), dim=-1)
        self.register_buffer("cos", emb.cos()[None, None, :, :])
        self.register_buffer("sin", emb.sin()[None, None, :, :])

    def forward(self, q, k):
        L = q.size(-2)
        cos, sin = self.cos[:, :, :L, :], self.sin[:, :, :L, :]

        def rotate_half(x):
            x1, x2 = x[..., : x.shape[-1] // 2], x[..., x.shape[-1] // 2:]
            return torch.cat((-x2, x1), dim=-1)

        return (q * cos) + (rotate_half(q) * sin), (k * cos) + (rotate_half(k) * sin)


class MLP(nn.Module):
    def __init__(self, cfg):
        super().__init__()
        self.c_fc = nn.Linear(cfg["reduced_dim"], cfg["hidden_dimension"], bias=cfg["bias"])
        self.gelu = nn.GELU(approximate="none")
        self.c_proj = nn.Linear(cfg["hidden_dimension"], cfg["reduced_dim"], bias=cfg["bias"])

    def forward(self, x):
        return self.c_proj(self.gelu(self.c_fc(x)))


class AttentionLayer(nn.Module):
    def __init__(self, cfg):
        super().__init__()
        D = cfg["reduced_dim"]
        self.num_heads = cfg["num_heads"]
        self.mlp = MLP(cfg)
        self.ln_1 = nn.LayerNorm(D, bias=cfg["bias"])
        self.ln_2 = nn.LayerNorm(D, bias=cfg["bias"])
        self.rope = RotaryPositionalEmbedding(D // self.num_heads, max_len=MAX_LEN)
        self.q_proj = nn.Linear(D, D, bias=cfg["bias"])
        self.k_proj = nn.Linear(D, D, bias=cfg["bias"])
        self.v_proj = nn.Linear(D, D, bias=cfg["bias"])
        self.out_proj = nn.Linear(D, D, bias=cfg["bias"])

    def forward(self, x, mask):
        h = self.ln_1(x)
        B, L, D = h.shape
        H = self.num_heads
        q = self.q_proj(h).view(B, L, H, D // H).transpose(1, 2)
        k = self.k_proj(h).view(B, L, H, D // H).transpose(1, 2)
        v = self.v_proj(h).view(B, L, H, D // H).transpose(1, 2)
        q, k = self.rope(q, k)
        attn_mask = mask.unsqueeze(1).unsqueeze(2) if mask is not None else None
        out = F.scaled_dot_product_attention(q, k, v, attn_mask=attn_mask, dropout_p=0.0, is_causal=False)
        x = x + self.out_proj(out.transpose(1, 2).contiguous().view(B, L, D))
        return x + self.mlp(self.ln_2(x))


class AttentionNetwork(nn.Module):
    def __init__(self, cfg):
        super().__init__()
        self.input_reduction = nn.Linear(cfg["inferred_input_dim"], cfg["reduced_dim"], bias=cfg["bias"])
        self.layers = nn.ModuleList([AttentionLayer(cfg) for _ in range(cfg["num_layers"])])
        self.ln_f = nn.LayerNorm(cfg["reduced_dim"], bias=cfg["bias"])
        self.lm_head = nn.Linear(cfg["reduced_dim"], 9, bias=cfg["bias"])
        self.conf_head = nn.Linear(cfg["reduced_dim"], 1, bias=cfg["bias"])

    def forward(self, x, mask):
        """x: (B, L, ESM2 dim) embeddings. Returns grade logits (B, L, 9) and uncertainty (B, L) in [0, 8]."""
        x = self.input_reduction(x)
        for layer in self.layers:
            x = layer(x, mask)
        x = self.ln_f(x)
        return self.lm_head(x), 8.0 * torch.sigmoid(self.conf_head(x).squeeze(-1))


# ── Weights ──────────────────────────────────────────────────────────────────
def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def find_weights(path=None):
    """Locate the ConSurf-AI weights, in this order:
    1. `path` (or the CONSURF_AI_WEIGHTS environment variable),
    2. weights/ in a cloned copy of the repository,
    3. the user cache (~/.cache/consurf_ai), downloading the file from GitHub the first time."""
    path = path or os.environ.get("CONSURF_AI_WEIGHTS")
    if path:
        if not os.path.exists(path):
            raise FileNotFoundError(f"Weights file not found: {path}")
        return str(path)
    repo_copy = Path(__file__).resolve().parents[1] / "weights" / WEIGHTS_FILE
    if repo_copy.exists():
        return str(repo_copy)
    cache = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "consurf_ai" / WEIGHTS_FILE
    if not cache.exists():
        if "<" in WEIGHTS_URL:
            raise FileNotFoundError("ConSurf-AI weights not found and no download URL is configured. "
                                    "Pass --weights /path/to/" + WEIGHTS_FILE)
        cache.parent.mkdir(parents=True, exist_ok=True)
        print(f"Downloading ConSurf-AI weights to {cache} ...")
        tmp = cache.with_suffix(".part")
        req = urllib.request.Request(WEIGHTS_URL, headers={"User-Agent": "ConSurf-AI"})
        with urllib.request.urlopen(req, timeout=300) as r, open(tmp, "wb") as f:
            shutil.copyfileobj(r, f)
        if _sha256(tmp) != WEIGHTS_SHA256:
            tmp.unlink()
            raise RuntimeError(f"Checksum mismatch for the file downloaded from {WEIGHTS_URL}.")
        tmp.replace(cache)
    return str(cache)


def load_network(path=None, device="cpu"):
    ckpt = torch.load(find_weights(path), map_location="cpu", weights_only=True)
    net = AttentionNetwork(ckpt["config"])
    state = {k.removeprefix("_orig_mod.").removeprefix("module."): v for k, v in ckpt["model"].items()}
    net.load_state_dict(state, strict=True)
    return net.to(device).eval(), ckpt.get("esm2_name") or "facebook/esm2_t36_3B_UR50D"
