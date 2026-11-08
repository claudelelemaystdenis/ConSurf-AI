"""Output files shared by both modes: tab-separated tables + README.txt, and for structures the structure file
with grades as B-factors and PyMOL / ChimeraX colouring scripts (the same files as the Colab notebooks)."""
import datetime
import os

import numpy as np
import pandas as pd

from .structure import cif_with_bfactors, pdb_with_bfactors, query_subchains

CONSURF_COLORS = {1: "#10C8D1", 2: "#8CFFFF", 3: "#D7FFFF", 4: "#EAFFFF", 5: "#FFFFFF",
                  6: "#FCEDF4", 7: "#FAC9DE", 8: "#F07DAB", 9: "#A02560"}
LOWCONF_COLOR = "#FFFF96"

# Description of every column, written to README.txt next to the tables (the tables have no comment lines,
# so they open cleanly in Excel, R and pandas).
COLUMN_HELP = {
    "SEQ_ID": "sequence identifier (first word of the FASTA header, made unique)",
    "POS": "position in the sequence (1 = first residue)",
    "CHAIN": "chain ID in the structure file",
    "RESIDUE_NUMBER": "residue number in the structure file (with insertion code, if any)",
    "RESNAME": "residue name (three-letter code)",
    "AA": "amino acid (one-letter code)",
    "GRADE": "ConSurf-AI conservation grade: 1 = most variable ... 9 = most conserved",
    "UNCERTAINTY_RAW": "predicted uncertainty of the grade, 0-8 (lower = better)",
    "UNCERTAINTY": "UNCERTAINTY_RAW rounded to the nearest integer",
    "CONFIDENT": "True if UNCERTAINTY < 2; False = uncertain prediction (same role as ConSurf's "
                 "'insufficient data' flag)",
    "P_GRADE_1 ... P_GRADE_9": "predicted probability of each grade (they sum to 1)",
    "DESCRIPTION": "full FASTA header (or the name given with --name)",
    "LENGTH": "number of residues after cleaning",
    "STATUS": "ok, skipped (empty or longer than 2000 residues) or error",
    "NOTE": "changes made to the sequence (removed gaps/stops, unknown letters read as X)",
    "N_UNCERTAIN": "number of residues with an uncertain prediction (CONFIDENT = False)",
    "FRACTION_UNCERTAIN": "N_UNCERTAIN / LENGTH",
    "MEAN_GRADE": "mean grade over the sequence",
    "N_GRADE_1 ... N_GRADE_9": "number of residues with each grade",
    "SEQUENCE": "the sequence as used for the prediction",
}


def write_readme(out_dir, title, input_desc, files, column_sets):
    """files: [(name, description)]; column_sets: [(table name, [columns])]."""
    from . import __version__
    w = max(len(n) for n, _ in files) + 2
    out = [title, "=" * len(title), "",
           f"Input:  {input_desc}",
           f"Date:   {datetime.datetime.now(datetime.timezone.utc):%Y-%m-%d %H:%M} UTC",
           f"Model:  ConSurf-AI {__version__} (attention network on ESM2-3B embeddings)", "",
           "Files", "-----"]
    out += [f"{n:<{w}}{d}" for n, d in files]
    for table, cols in column_sets:
        out += ["", f"Columns of {table}", "-" * (11 + len(table))]
        shown = [c for c in cols if not c.startswith("P_GRADE_")]
        if any(c.startswith("P_GRADE_") for c in cols):
            shown.append("P_GRADE_1 ... P_GRADE_9")
        cw = max(len(c) for c in shown) + 2
        out += [f"{c:<{cw}}{COLUMN_HELP.get(c, '')}" for c in shown]
    out += ["", "Tables are tab-separated (.tsv). Excel: File > Open (or drag the file onto Excel). "
                "R: read.delim(file). Python: pandas.read_csv(file, sep='\\t')."]
    path = os.path.join(out_dir, "README.txt")
    with open(path, "w") as f:
        f.write("\n".join(out) + "\n")
    return path


# ── Structure mode ───────────────────────────────────────────────────────────
def grades_table(pred, residues, chain):
    """Per-residue table for a structure chain (residues: list of dicts from structure.chain_residues)."""
    df = pd.DataFrame({
        "POS": np.arange(1, len(residues) + 1),
        "CHAIN": chain,
        "RESIDUE_NUMBER": [r["key"] for r in residues],
        "RESNAME": [r["resname"] for r in residues],
        "AA": [r["aa"] for r in residues],
        "GRADE": pred["grades"].astype(int),
        "UNCERTAINTY_RAW": np.round(pred["uncertainty"], 3),
        "UNCERTAINTY": pred["uncertainty_rounded"].astype(int),
        "CONFIDENT": pred["confident"].astype(bool),
    })
    for g in range(9):
        df[f"P_GRADE_{g + 1}"] = np.round(pred["probs"][:, g], 4)
    return df


def _rgb01(hex_color):
    h = hex_color.lstrip("#")
    return [round(int(h[i:i + 2], 16) / 255, 3) for i in (0, 2, 4)]


def pymol_script(df, chain, structure_file):
    obj, q = "consurf_ai", f"consurf_ai and polymer and chain {chain}"
    lines = [f"# ConSurf-AI colouring. Run from the folder containing {structure_file}:",
             "#   pymol <this_script>.pml", f"load {structure_file}, {obj}",
             f"hide everything, {obj}", f"show cartoon, {obj}", f"color grey80, {obj}",
             f"show sticks, {obj} and organic", f"color green, {obj} and organic and elem C",
             f"color atomic, {obj} and organic and not elem C",
             f"show spheres, {obj} and inorganic", f"set sphere_scale, 0.5, {obj} and inorganic"]
    for g, c in CONSURF_COLORS.items():
        lines.append(f"set_color consurf_{g}, {_rgb01(c)}")
    lines.append(f"set_color consurf_lowconf, {_rgb01(LOWCONF_COLOR)}")
    for g in range(1, 10):
        keys = df.loc[(df.GRADE == g) & df.CONFIDENT, "RESIDUE_NUMBER"].tolist()
        if keys:
            lines.append(f"color consurf_{g}, {q} and resi {'+'.join(map(str, keys))}")
    low = df.loc[~df.CONFIDENT, "RESIDUE_NUMBER"].tolist()
    if low:
        lines.append(f"color consurf_lowconf, {q} and resi {'+'.join(map(str, low))}")
    lines += ["bg_color white", f"orient {q}"]
    return "\n".join(lines) + "\n"


def chimerax_script(df, chain, structure_file):
    lines = [f"# ConSurf-AI colouring for ChimeraX (keep next to {structure_file}); open this .cxc file in ChimeraX",
             f"open {structure_file}", "cartoon", "color #1 lightgray",
             "show ligand", "style ligand stick", "color ligand & C forestgreen", "color ligand byhetero"]
    for g, c in CONSURF_COLORS.items():
        keys = df.loc[(df.GRADE == g) & df.CONFIDENT, "RESIDUE_NUMBER"].tolist()
        if keys:
            lines.append(f"color #1/{chain}:{','.join(map(str, keys))} & protein {c}")
    low = df.loc[~df.CONFIDENT, "RESIDUE_NUMBER"].tolist()
    if low:
        lines.append(f"color #1/{chain}:{','.join(map(str, low))} & protein {LOWCONF_COLOR}")
    lines.append(f"view #1/{chain}")
    return "\n".join(lines) + "\n"


def write_structure_outputs(df, inputs, st, chain, out_dir, base, source):
    """Writes <base>_grades.tsv, the structure file(s) with grades as B-factors, PyMOL / ChimeraX scripts and
    README.txt. `inputs` = {"cif": path, "pdb": path} (original files). Returns the list of written paths."""
    os.makedirs(out_dir, exist_ok=True)
    tsv = os.path.join(out_dir, f"{base}_grades.tsv")
    df.to_csv(tsv, sep="\t", index=False)
    written = [tsv]
    values = dict(zip(df.RESIDUE_NUMBER, df.GRADE.astype(int)))
    struct_name = None
    if "cif" in inputs:
        p = os.path.join(out_dir, f"{base}_structure.cif")
        cif_with_bfactors(inputs["cif"], p, chain, values, query_subchains(st, chain))
        written.append(p)
        struct_name = os.path.basename(p)
    if "pdb" in inputs:
        p = os.path.join(out_dir, f"{base}_structure.pdb")
        pdb_with_bfactors(inputs["pdb"], p, chain, values)
        written.append(p)
        struct_name = os.path.basename(p)      # PDB preferred for the scripts (CONECT is read everywhere)
    for ext, fn, tag in (("pml", pymol_script, "pymol"), ("cxc", chimerax_script, "chimerax")):
        p = os.path.join(out_dir, f"{base}_{tag}.{ext}")
        with open(p, "w") as f:
            f.write(fn(df, chain, struct_name))
        written.append(p)
    desc = {".tsv": "per-residue grades and uncertainties (one row per residue)",
            ".cif": f"input structure file; B-factor column = grade on chain {chain} (0 elsewhere)",
            ".pdb": f"input structure file (PDB format); B-factor column = grade on chain {chain} (0 elsewhere)",
            ".pml": "PyMOL script: ConSurf colouring (run next to the structure file)",
            ".cxc": "ChimeraX script: ConSurf colouring (open next to the structure file)"}
    written.append(write_readme(out_dir, "ConSurf-AI results", f"{source}, chain {chain}",
                                [(os.path.basename(p), desc[os.path.splitext(p)[1]]) for p in written],
                                [(os.path.basename(tsv), list(df.columns))]))
    return written
