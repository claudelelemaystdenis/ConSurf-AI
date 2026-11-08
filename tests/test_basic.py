"""Fast CPU tests (no ESM2 download, no network): weights load, network runs, B-factor files keep everything else."""
import os

import torch

from consurf_ai.model import MAX_LEN, WEIGHTS_SHA256, _sha256, find_weights, load_network
from consurf_ai.structure import chain_residues, load_structure, pdb_with_bfactors, protein_chains

PDB = """\
HEADER    TEST
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 10.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 10.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 10.00           C
ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00 10.00           O
ATOM      5  N   GLY A   2       3.332   1.536   0.000  1.00 10.00           N
ATOM      6  CA  GLY A   2       3.970   2.846   0.000  1.00 10.00           C
ATOM      7  C   GLY A   2       5.486   2.705   0.000  1.00 10.00           C
ATOM      8  O   GLY A   2       6.009   1.593   0.000  1.00 10.00           O
HETATM    9  O   HOH A 101      10.000  10.000  10.000  1.00 20.00           O
CONECT    1    2
END
"""


def test_weights_checksum():
    assert _sha256(find_weights()) == WEIGHTS_SHA256


def test_network_forward():
    net, esm_name = load_network(device="cpu")
    assert esm_name == "facebook/esm2_t36_3B_UR50D"
    x = torch.randn(1, 30, net.input_reduction.in_features)
    logits, unc = net(x, torch.ones(1, 30, dtype=torch.bool))
    assert logits.shape == (1, 30, 9) and unc.shape == (1, 30)
    assert float(unc.min()) >= 0 and float(unc.max()) <= 8
    assert net.layers[0].rope.cos.shape[2] == MAX_LEN


def test_pdb_bfactors_only_change_bfactor(tmp_path):
    src, out = tmp_path / "in.pdb", tmp_path / "out.pdb"
    src.write_text(PDB)
    st = load_structure(str(src))
    assert protein_chains(st) == ["A"]
    assert [r["key"] for r in chain_residues(st, "A")] == ["1", "2"]
    pdb_with_bfactors(str(src), str(out), "A", {"1": 9, "2": 3})
    a, b = src.read_text().splitlines(), out.read_text().splitlines()
    assert len(a) == len(b)
    for la, lb in zip(a, b):
        assert la[:60] == lb[:60] and la[66:] == lb[66:]      # everything but columns 61-66 unchanged
    bf = [float(l[60:66]) for l in b if l.startswith(("ATOM", "HETATM"))]
    assert bf == [9, 9, 9, 9, 3, 3, 3, 3, 0]                    # water: 0
