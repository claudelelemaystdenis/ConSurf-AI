"""Output files: per-residue CSV, structure with grades as B-factors, PyMOL / ChimeraX colouring scripts."""
import os

import numpy as np
import pandas as pd

from .structure import cif_with_bfactors, pdb_with_bfactors, query_subchains

CONSURF_COLORS = {1: "#10C8D1", 2: "#8CFFFF", 3: "#D7FFFF", 4: "#EAFFFF", 5: "#FFFFFF",
                  6: "#FCEDF4", 7: "#FAC9DE", 8: "#F07DAB", 9: "#A02560"}
LOWCONF_COLOR = "#FFFF96"

CSV_HEADER = ("# GRADE: 1 = most variable ... 9 = most conserved. UNCERTAINTY_RAW: predicted uncertainty "
              "(0-8, lower = better); UNCERTAINTY: rounded value; CONFIDENT = UNCERTAINTY < 2. "
              "P_GRADE_k: predicted probability of grade k.\n")


def grades_table(pred, aa, chain="", residue_numbers=None, resnames=None):
    n = len(aa)
    df = pd.DataFrame({
        "POS": np.arange(1, n + 1),
        "CHAIN": chain,
        "RESIDUE_NUMBER": residue_numbers if residue_numbers is not None else [str(i) for i in range(1, n + 1)],
        "RESNAME": resnames if resnames is not None else [""] * n,
        "AA": list(aa),
        "GRADE": pred["grades"].astype(int),
        "UNCERTAINTY_RAW": np.round(pred["uncertainty"], 3),
        "UNCERTAINTY": pred["uncertainty_rounded"].astype(int),
        "CONFIDENT": pred["confident"].astype(bool),
    })
    for g in range(9):
        df[f"P_GRADE_{g + 1}"] = np.round(pred["probs"][:, g], 4)
    return df


def write_csv(df, path, title):
    with open(path, "w") as f:
        f.write(f"# ConSurf-AI predictions for {title}\n" + CSV_HEADER)
        df.to_csv(f, index=False)


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


def write_structure_outputs(df, inputs, st, chain, out_dir, base):
    """Writes the structure file(s) with grades as B-factors and the PyMOL / ChimeraX scripts.
    `inputs` = {"cif": path, "pdb": path} (original files). Returns the list of written paths."""
    values = dict(zip(df.RESIDUE_NUMBER, df.GRADE.astype(int)))
    written, struct_name = [], None
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
    return written
