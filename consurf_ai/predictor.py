"""ESM2-3B embeddings + ConSurf-AI network -> per-residue conservation grades."""
import contextlib
import io
import re

import torch
import torch.nn.functional as F

from .model import MAX_LEN, load_network


def default_device():
    if torch.cuda.is_available():
        return torch.device("cuda")
    if getattr(torch.backends, "mps", None) is not None and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


class ConSurfAI:
    """Loads ESM2-3B and the ConSurf-AI network once; `predict(sequence)` then takes ~1 s per chain on a GPU.

    >>> model = ConSurfAI()
    >>> out = model.predict("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMF")
    >>> out["grades"]          # 1 (variable) ... 9 (conserved)
    """

    def __init__(self, weights=None, device=None, esm_dtype=None, verbose=True):
        import transformers
        from transformers import AutoTokenizer, EsmModel
        from transformers.utils import logging as hf_logging

        self.device = torch.device(device) if device else default_device()
        self.net, self.esm_name = load_network(weights, self.device)

        if esm_dtype is None:   # fp32 ESM2-3B needs ~12 GB; use fp16 on smaller GPUs
            esm_dtype = torch.float32
            if self.device.type == "cuda" and torch.cuda.get_device_properties(self.device).total_memory / 2**30 < 14:
                esm_dtype = torch.float16
        self.esm_dtype = esm_dtype
        if verbose:
            print(f"Loading {self.esm_name} on {self.device} ({str(esm_dtype).replace('torch.', '')}); "
                  "the first run downloads ~11 GB from Hugging Face...")
        tv = tuple(int(x) for x in re.findall(r"\d+", transformers.__version__)[:2])
        kw = {"add_pooling_layer": False, ("dtype" if tv >= (4, 56) else "torch_dtype"): esm_dtype}
        if tv < (5, 0):
            try:
                import accelerate  # noqa: F401  (needed for low_cpu_mem_usage: avoids holding two copies in RAM)
                kw["low_cpu_mem_usage"] = True
            except ImportError:
                pass
        # Only the encoder is used, so the masked-LM head stored in the checkpoint is skipped on purpose;
        # silence transformers' report about those unused weights.
        prev = hf_logging.get_verbosity()
        hf_logging.set_verbosity_error()
        try:
            with contextlib.redirect_stderr(io.StringIO()):
                self.tokenizer = AutoTokenizer.from_pretrained(self.esm_name)
                self.esm = EsmModel.from_pretrained(self.esm_name, **kw).to(self.device).eval()
        finally:
            hf_logging.set_verbosity(prev)

    @torch.no_grad()
    def predict(self, sequence):
        """Returns a dict of numpy arrays (one value per residue):
        grades (1-9), probs (L x 9), uncertainty (raw, 0-8, lower = better),
        uncertainty_rounded, confident (rounded uncertainty < 2)."""
        sequence = sequence.strip().upper()
        if not sequence:
            raise ValueError("Empty sequence.")
        if len(sequence) > MAX_LEN:
            raise ValueError(f"Sequence has {len(sequence)} residues; ConSurf-AI supports up to {MAX_LEN}.")
        tokens = self.tokenizer(sequence, return_tensors="pt").to(self.device)
        hidden = self.esm(**tokens).last_hidden_state[0, 1:1 + len(sequence)].float()   # drop BOS/EOS
        if not torch.isfinite(hidden).all():
            raise RuntimeError("ESM2 produced non-finite embeddings (fp16 overflow?). Try esm_dtype=torch.float32.")
        mask = torch.ones(1, len(sequence), dtype=torch.bool, device=self.device)
        logits, unc = self.net(hidden.unsqueeze(0), mask)
        probs = F.softmax(logits[0].float(), dim=-1)
        unc = unc[0].float()
        unc_rounded = torch.floor(unc + 0.5).clamp(0, 8).int()
        return {"grades": (probs.argmax(dim=-1) + 1).cpu().numpy(),
                "probs": probs.cpu().numpy(),
                "uncertainty": unc.cpu().numpy(),
                "uncertainty_rounded": unc_rounded.cpu().numpy(),
                "confident": (unc_rounded < 2).cpu().numpy()}
