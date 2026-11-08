"""Command-line interface:

    consurf-ai proteins.fasta                       # one or many sequences -> tables (like the sequence notebook)
    consurf-ai --sequence MVLSPADKTNVKAAW --name my_protein
    consurf-ai 4HHB --chain A                       # PDB ID (RCSB)            -> grades mapped on the structure
    consurf-ai P69905                               # UniProt accession (AlphaFold DB model)
    consurf-ai my_model.cif --chain all             # local PDB / mmCIF file (optionally .gz)
"""
import argparse
import os
import re
import sys
import time

from . import __version__

_FASTA_EXT = (".fa", ".fasta", ".faa", ".fas", ".fa.gz", ".fasta.gz", ".faa.gz")


def run_sequences(model, raw_records, out_root, job, source, verbose=True):
    """Predicts every sequence; writes <out_root>/<job>/ (summary, all residues, grade strings, per_sequence/,
    README.txt). Returns the output folder."""
    from .sequences import prepare_sequences, sequence_grades_table, write_sequence_outputs
    records = prepare_sequences(raw_records)
    if not records:
        raise SystemExit("No sequence found in the input.")
    todo = [r for r in records if r["status"] == "ok"]
    tables, t0 = {}, time.time()
    for i, r in enumerate(todo, 1):
        try:
            tables[r["id"]] = sequence_grades_table(r, model.predict(r["sequence"]))
        except Exception as e:                        # keep going (e.g. out of memory on a long sequence)
            r["status"] = f"error: {type(e).__name__}: {str(e)[:150]}"
        if verbose:
            print(f"   {i}/{len(todo)} {r['id']} ({len(r['sequence'])} residues): {r['status']} "
                  f"[{time.time() - t0:.0f} s]")
    out_dir = os.path.join(out_root, job)
    summary, written = write_sequence_outputs(out_dir, job, source, records, tables)
    if verbose:
        failed = [r for r in records if r["status"] != "ok"]
        print(f"Predicted {len(tables)} sequence(s)" + (f"; {len(failed)} skipped or failed (see STATUS)" if failed else "")
              + f". Files written to {out_dir}:")
        for p in written:
            print("   " + os.path.basename(p))
        if tables:
            print(f"   per_sequence/ ({len(tables)} files)")
    return out_dir


def run_structure(model, query, chains, out_root, verbose=True):
    """Predicts the chain(s) of a structure; writes <out_root>/<label>_<chain>/ per chain. Returns the folders."""
    from .outputs import grades_table, write_structure_outputs
    from .sequences import safe_name
    from .structure import chain_residues, load_structure, protein_chains, resolve_input

    files, label, source = resolve_input(query, os.path.join(out_root, "input"))
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
    out_dirs = []
    for chain in chains:
        res = chain_residues(st, chain)
        seq = "".join(r["aa"] for r in res)
        df = grades_table(model.predict(seq), res, chain)
        base = f"{safe_name(label)}_{chain}_ConSurf-AI"
        out_dir = os.path.join(out_root, f"{safe_name(label)}_{chain}")
        written = write_structure_outputs(df, files, st, chain, out_dir, base, source)
        out_dirs.append(out_dir)
        if verbose:
            print(f"{source} chain {chain}: {len(seq)} residues, {(~df.CONFIDENT).sum()} uncertain. "
                  f"Files written to {out_dir}:")
            for p in written:
                print("   " + os.path.basename(p))
    return out_dirs


def main(argv=None):
    ap = argparse.ArgumentParser(prog="consurf-ai", description="ConSurf-AI: per-residue conservation grades (1-9) "
                                 "from a single sequence or structure, without a multiple sequence alignment.")
    ap.add_argument("input", nargs="?", help="FASTA file, PDB ID, UniProt accession, AlphaFold model ID or "
                                             "structure file (PDB/mmCIF)")
    ap.add_argument("-s", "--sequence", help="a single protein sequence (one-letter code) instead of an input file")
    ap.add_argument("-n", "--name", help="name of the results (default: the input file name, or 'query')")
    ap.add_argument("-c", "--chain", action="append", default=[],
                    help="structure chain to predict (repeatable; 'all' for every protein chain; default: first chain)")
    ap.add_argument("-o", "--out", default="consurf_ai_results", help="output directory")
    ap.add_argument("--fasta", action="store_true", help="treat the input as FASTA regardless of its extension")
    ap.add_argument("--device", help="cuda, cuda:1, mps or cpu (default: best available)")
    ap.add_argument("--weights", help="path to consurf_ai_3B_weights.pth (default: repository copy or download)")
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = ap.parse_args(argv)
    if bool(args.input) == bool(args.sequence):
        ap.error("give either an input (FASTA file, PDB ID, UniProt accession or structure file) or --sequence")

    from .predictor import ConSurfAI
    from .sequences import read_fasta, safe_name
    os.makedirs(args.out, exist_ok=True)
    model = ConSurfAI(weights=args.weights, device=args.device)
    if model.device.type == "cpu":
        print("Note: running on CPU. This works but is slow (ESM2-3B; ~12 GB of RAM needed).", file=sys.stderr)

    if args.sequence:
        name = safe_name(args.name or "query")          # whole name as sequence ID (spaces -> _)
        run_sequences(model, [(name, args.sequence)], args.out, name, "sequence given on the command line")
    elif args.fasta or args.input.lower().endswith(_FASTA_EXT):
        job = safe_name(args.name or re.sub(r"\.(fa|fasta|faa|fas)(\.gz)?$", "", os.path.basename(args.input), flags=re.I))
        run_sequences(model, read_fasta(args.input), args.out, job, f"file {os.path.basename(args.input)}")
    else:
        run_structure(model, args.input, args.chain, args.out)


if __name__ == "__main__":
    main()
