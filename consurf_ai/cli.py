"""Command-line interface:

    consurf-ai 4HHB --chain A                 # PDB ID (RCSB)
    consurf-ai P69905                         # UniProt accession (AlphaFold DB model)
    consurf-ai my_model.pdb --chain all       # local PDB / mmCIF file (optionally .gz)
    consurf-ai proteins.fasta                 # sequences only -> one CSV per sequence
"""
import argparse
import os
import re
import sys

from . import __version__

_FASTA_EXT = (".fa", ".fasta", ".faa", ".fas", ".fa.gz", ".fasta.gz", ".faa.gz")


def _safe(name):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_") or "query"


def read_fasta(path):
    import gzip
    opener = gzip.open if path.endswith(".gz") else open
    records, name, seq = [], None, []
    with opener(path, "rt") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq)))
                name, seq = line[1:].split()[0] if line[1:].split() else f"seq{len(records) + 1}", []
            elif line:
                seq.append(line.replace("*", ""))
    if name is not None:
        records.append((name, "".join(seq)))
    return records


def run_fasta(model, path, out_dir):
    from .outputs import grades_table, write_csv
    for name, seq in read_fasta(path):
        seq = seq.upper()
        df = grades_table(model.predict(seq), seq)
        out = os.path.join(out_dir, f"{_safe(name)}_ConSurf-AI_grades.csv")
        write_csv(df, out, name)
        print(f"{name}: {len(seq)} residues -> {out}")


def run_structure(model, query, chains, out_dir):
    from .outputs import grades_table, write_csv, write_structure_outputs
    from .structure import chain_residues, load_structure, protein_chains, resolve_input

    files, label, source = resolve_input(query, os.path.join(out_dir, "input"))
    st = load_structure(files.get("cif") or files["pdb"])
    available = protein_chains(st)
    if not available:
        raise SystemExit(f"No protein chain found in {source}.")
    if chains == ["all"]:
        chains = available
    elif not chains:
        chains = available[:1]
    missing = [c for c in chains if c not in available]
    if missing:
        raise SystemExit(f"Chain(s) {', '.join(missing)} not found. Protein chains: {', '.join(available)}")

    for chain in chains:
        res = chain_residues(st, chain)
        seq = "".join(r["aa"] for r in res)
        df = grades_table(model.predict(seq), seq, chain, [r["key"] for r in res], [r["resname"] for r in res])
        base = f"{_safe(label)}_{chain}_ConSurf-AI"
        chain_dir = os.path.join(out_dir, f"{_safe(label)}_{chain}")
        os.makedirs(chain_dir, exist_ok=True)
        csv = os.path.join(chain_dir, f"{base}_grades.csv")
        write_csv(df, csv, f"{source}, chain {chain}")
        written = [csv] + write_structure_outputs(df, files, st, chain, chain_dir, base)
        print(f"{source} chain {chain}: {len(seq)} residues, {(~df.CONFIDENT).sum()} uncertain")
        for p in written:
            print("   " + p)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="consurf-ai", description="ConSurf-AI: per-residue conservation grades (1-9) "
                                 "from a single sequence or structure, without a multiple sequence alignment.")
    ap.add_argument("input", help="PDB ID, UniProt accession, AlphaFold model ID, structure file (PDB/mmCIF) "
                                  "or FASTA file")
    ap.add_argument("-c", "--chain", action="append", default=[],
                    help="chain to predict (repeatable; 'all' for every protein chain; default: first chain)")
    ap.add_argument("-o", "--out", default="consurf_ai_results", help="output directory")
    ap.add_argument("--fasta", action="store_true", help="treat the input as FASTA regardless of its extension")
    ap.add_argument("--device", help="cuda, cuda:1, mps or cpu (default: best available)")
    ap.add_argument("--weights", help="path to consurf_ai_3B_weights.pth (default: repository copy or download)")
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = ap.parse_args(argv)

    from .predictor import ConSurfAI
    os.makedirs(args.out, exist_ok=True)
    model = ConSurfAI(weights=args.weights, device=args.device)
    if model.device.type == "cpu":
        print("Note: running on CPU. This works but is slow (ESM2-3B; ~12 GB of RAM needed).", file=sys.stderr)
    if args.fasta or args.input.lower().endswith(_FASTA_EXT):
        run_fasta(model, args.input, args.out)
    else:
        run_structure(model, args.input, args.chain, args.out)


if __name__ == "__main__":
    main()
