"""Check a ConSurf-AI installation against reference predictions.

Runs ConSurf-AI on the two test inputs and compares the results with the expected outputs in tests/expected/:
  1. structure path: tests/input/4HHB.cif, chain A        -> tests/expected/4HHB_A/
  2. sequence path:  tests/input/dihydrofolate_reductases.fasta -> tests/expected/dihydrofolate_reductases/

Usage (from the repository folder, in the consurf-ai environment):
    python tests/run_tests.py                 # best available device (GPU if there is one)
    python tests/run_tests.py --device cpu

The first run downloads ESM2-3B (~11 GB). Grades must be identical; uncertainties and probabilities may differ
slightly between machines (GPU vs CPU arithmetic), so differences up to --tolerance (default 0.01) are accepted.
"""
import argparse
import filecmp
import os
import shutil
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))       # use this repository's consurf_ai package

INPUT = os.path.join(HERE, "input")
EXPECTED = os.path.join(HERE, "expected")
STRUCTURE_FILE = os.path.join(INPUT, "4HHB.cif")
STRUCTURE_CHAIN = "A"
FASTA_FILE = os.path.join(INPUT, "dihydrofolate_reductases.fasta")
FASTA_JOB = "dihydrofolate_reductases"
NUMERIC = ["UNCERTAINTY_RAW"] + [f"P_GRADE_{g}" for g in range(1, 10)]


def compare_tables(expected_path, produced_path, keys, tol, label):
    """Returns (ok, list of report lines)."""
    exp = pd.read_csv(expected_path, sep="\t", dtype={k: str for k in keys})
    got = pd.read_csv(produced_path, sep="\t", dtype={k: str for k in keys})
    lines = []
    if list(exp[keys].itertuples(index=False)) != list(got[keys].itertuples(index=False)):
        return False, [f"  ✗ {label}: the residues differ (expected {len(exp)}, got {len(got)}) – "
                       "is the input file the one in tests/input/?"]
    diff = exp.GRADE.values != got.GRADE.values
    numeric = (exp[NUMERIC] - got[NUMERIC]).abs()
    max_num = float(numeric.values.max())
    ok = True
    if diff.any():
        ok = False
        lines.append(f"  ✗ {label}: {diff.sum()} of {len(exp)} residues have a different GRADE:")
        for i in list(diff.nonzero()[0])[:15]:
            e, g = exp.iloc[i], got.iloc[i]
            where = " ".join(f"{k}={e[k]}" for k in keys)
            lines.append(f"      {where} {e.AA}: expected grade {e.GRADE} (p={e[f'P_GRADE_{e.GRADE}']:.3f}), "
                         f"got {g.GRADE} (p={g[f'P_GRADE_{g.GRADE}']:.3f})")
        if diff.sum() > 15:
            lines.append(f"      ... and {diff.sum() - 15} more")
    else:
        lines.append(f"  ✓ {label}: all {len(exp)} grades identical")
    if exp.CONFIDENT.ne(got.CONFIDENT).any():
        ok = False
        lines.append(f"  ✗ {label}: the confidence flag differs for {exp.CONFIDENT.ne(got.CONFIDENT).sum()} residue(s)")
    if max_num > tol:
        ok = False
        worst = numeric.max().sort_values(ascending=False)
        lines.append(f"  ✗ {label}: uncertainties/probabilities differ by up to {max_num:.4f} (tolerance {tol}); "
                     f"largest in {', '.join(worst.index[:3])}")
    else:
        lines.append(f"  ✓ {label}: uncertainties and probabilities within {tol} (max difference {max_num:.4f})")
    return ok, lines


def compare_bfactors(expected_cif, produced_cif, label):
    import gemmi
    e, g = gemmi.read_structure(expected_cif), gemmi.read_structure(produced_cif)
    eb = [(ch.name, str(r.seqid), a.name, a.altloc, round(a.b_iso, 2)) for ch in e[0] for r in ch for a in r]
    gb = [(ch.name, str(r.seqid), a.name, a.altloc, round(a.b_iso, 2)) for ch in g[0] for r in ch for a in r]
    if eb == gb:
        return True, [f"  ✓ {label}: structure file identical (all {len(eb)} atoms, grades in the B-factor column)"]
    n = sum(x != y for x, y in zip(eb, gb)) + abs(len(eb) - len(gb))
    return False, [f"  ✗ {label}: B-factor (grade) differs for {n} atom(s)"]


def compare_texts(expected_dir, produced_dir, names, label):
    bad = [n for n in names if not filecmp.cmp(os.path.join(expected_dir, n), os.path.join(produced_dir, n), shallow=False)]
    if bad:
        return False, [f"  ✗ {label}: differs from the expected file(s): {', '.join(bad)}"]
    return True, [f"  ✓ {label}: identical ({', '.join(names)})"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--device", help="cuda, mps or cpu (default: best available)")
    ap.add_argument("--tolerance", type=float, default=0.01,
                    help="accepted difference in uncertainties and probabilities (default 0.01)")
    ap.add_argument("--out", default=os.path.join(HERE, "output"), help="where to write the produced results")
    ap.add_argument("--update-expected", action="store_true",
                    help="(maintainers) overwrite tests/expected/ with the results of this run")
    args = ap.parse_args()

    from consurf_ai.cli import run_sequences, run_structure
    from consurf_ai.predictor import ConSurfAI
    from consurf_ai.sequences import read_fasta

    shutil.rmtree(args.out, ignore_errors=True)
    model = ConSurfAI(device=args.device)
    import torch
    print(f"Device: {model.device} (ESM2 in {str(model.esm_dtype).replace('torch.', '')}), torch {torch.__version__}\n")

    print("Running the structure path (4HHB, chain A)...")
    s_dir = run_structure(model, STRUCTURE_FILE, [STRUCTURE_CHAIN], args.out, verbose=False)[0]
    print("Running the sequence path (5 dihydrofolate reductases)...")
    q_dir = run_sequences(model, read_fasta(FASTA_FILE), args.out, FASTA_JOB,
                          f"file {os.path.basename(FASTA_FILE)}", verbose=False)

    if args.update_expected:
        for d in (s_dir, q_dir):
            dest = os.path.join(EXPECTED, os.path.basename(d))
            shutil.rmtree(dest, ignore_errors=True)
            shutil.copytree(d, dest)
        print(f"\nExpected outputs updated in {EXPECTED}")
        return 0

    s_exp = os.path.join(EXPECTED, os.path.basename(s_dir))
    q_exp = os.path.join(EXPECTED, os.path.basename(q_dir))
    s_base, q_base = f"4HHB_{STRUCTURE_CHAIN}_ConSurf-AI", f"{FASTA_JOB}_ConSurf-AI"
    results = []
    print("\n1. Structure path (4HHB, chain A)")
    for ok, lines in (
        compare_tables(os.path.join(s_exp, f"{s_base}_grades.tsv"), os.path.join(s_dir, f"{s_base}_grades.tsv"),
                       ["RESIDUE_NUMBER"], args.tolerance, "grades table"),
        compare_bfactors(os.path.join(s_exp, f"{s_base}_structure.cif"), os.path.join(s_dir, f"{s_base}_structure.cif"),
                         "structure file"),
        compare_texts(s_exp, s_dir, [f"{s_base}_pymol.pml", f"{s_base}_chimerax.cxc"], "PyMOL / ChimeraX scripts"),
    ):
        results.append(ok)
        print("\n".join(lines))
    print("\n2. Sequence path (dihydrofolate_reductases.fasta, 5 sequences)")
    for ok, lines in (
        compare_tables(os.path.join(q_exp, f"{q_base}_all_residues.tsv"), os.path.join(q_dir, f"{q_base}_all_residues.tsv"),
                       ["SEQ_ID", "POS"], args.tolerance, "all-residues table"),
        compare_texts(q_exp, q_dir, [f"{q_base}_grade_strings.fasta"], "grade strings"),
        compare_texts(q_exp, q_dir, [f"{q_base}_summary.tsv"], "summary table"),   # derived from grades only
    ):
        results.append(ok)
        print("\n".join(lines))

    print(f"\nProduced files: {args.out}")
    if all(results):
        print("✅ PASSED: your installation reproduces the reference predictions.")
        return 0
    print("❌ DIFFERENCES FOUND (see ✗ above). Small differences can come from different hardware or library "
          "versions (e.g. a GPU running ESM2 in float16); different grades for many residues point to a problem "
          "with the installation.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
