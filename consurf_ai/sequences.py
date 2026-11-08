"""Sequence input (one sequence or a FASTA file) and the result files of the sequence mode
(the same files as the ConSurf-AI sequence notebook on Colab)."""
import gzip
import os
import re
import shutil

import numpy as np
import pandas as pd

from .model import MAX_LEN
from .outputs import COLUMN_HELP, write_readme

ALLOWED_AA = set("ACDEFGHIKLMNPQRSTVWYXUOBZ")


def safe_name(name):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_") or "query"


def read_fasta(path):
    """[(header, raw sequence)] from a FASTA file (optionally .gz)."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", errors="replace") as f:
        return parse_fasta_text(f.read())


def parse_fasta_text(text):
    """[(header, raw sequence)]; text without any '>' line is read as one sequence."""
    records, header, seq = [], None, []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            if header is not None or seq:
                records.append((header, "".join(seq)))
            header, seq = line[1:].strip(), []
        else:
            seq.append(line)
    if header is not None or seq:
        records.append((header, "".join(seq)))
    return records


def clean_sequence(raw):
    """Upper-case, drop spaces/digits/gaps/stops; unknown letters -> X. Returns (sequence, notes)."""
    s = re.sub(r"[\s\d]", "", raw).upper()
    notes = []
    n_gap = len(re.findall(r"[-.*]", s))
    if n_gap:
        s = re.sub(r"[-.*]", "", s)
        notes.append(f"{n_gap} gap/stop character(s) removed")
    bad = sorted(set(s) - ALLOWED_AA)
    if bad:
        s = "".join(c if c in ALLOWED_AA else "X" for c in s)
        notes.append(f"unknown character(s) {''.join(bad)} read as X")
    return s, notes


def prepare_sequences(raw_records):
    """List of dicts: id (first word of the header, made unique), description, sequence, status, note."""
    out, used = [], set()
    for i, (header, raw) in enumerate(raw_records, 1):
        header = " ".join((header or "").split()) or f"seq{i}"      # also removes tabs
        base = safe_name(header.split()[0])[:60]
        sid, k = base, 2
        while sid in used:
            sid, k = f"{base}_{k}", k + 1
        used.add(sid)
        seq, notes = clean_sequence(raw)
        status = "ok"
        if not seq:
            status = "skipped: empty sequence"
        elif len(seq) > MAX_LEN:
            status = f"skipped: longer than {MAX_LEN} residues"
        out.append({"id": sid, "description": header, "sequence": seq, "status": status, "note": "; ".join(notes)})
    return out


def sequence_grades_table(rec, pred):
    seq = rec["sequence"]
    df = pd.DataFrame({
        "SEQ_ID": rec["id"],
        "POS": np.arange(1, len(seq) + 1),
        "AA": list(seq),
        "GRADE": pred["grades"].astype(int),
        "UNCERTAINTY_RAW": np.round(pred["uncertainty"], 3),
        "UNCERTAINTY": pred["uncertainty_rounded"].astype(int),
        "CONFIDENT": pred["confident"].astype(bool),
    })
    for g in range(9):
        df[f"P_GRADE_{g + 1}"] = np.round(pred["probs"][:, g], 4)
    return df


def _summary_row(rec, df=None):
    row = {"SEQ_ID": rec["id"], "DESCRIPTION": rec["description"], "LENGTH": len(rec["sequence"]),
           "STATUS": rec["status"], "NOTE": rec["note"]}
    if df is not None:
        row.update({"N_UNCERTAIN": int((~df.CONFIDENT).sum()),
                    "FRACTION_UNCERTAIN": round(float((~df.CONFIDENT).mean()), 3),
                    "MEAN_GRADE": round(float(df.GRADE.mean()), 2)})
        counts = df.GRADE.value_counts()
        row.update({f"N_GRADE_{g}": int(counts.get(g, 0)) for g in range(1, 10)})
    row["SEQUENCE"] = rec["sequence"]
    return row


def summary_table(records, tables):
    df = pd.DataFrame([_summary_row(r, tables.get(r["id"])) for r in records])
    for c in [c for c in df if c.startswith("N_")]:      # whole numbers (empty for skipped sequences)
        df[c] = df[c].astype("Int64")
    return df


def write_sequence_outputs(out_dir, job, source, records, tables):
    """Writes, in out_dir: <job>_ConSurf-AI_summary.tsv, <job>_ConSurf-AI_all_residues.tsv,
    <job>_ConSurf-AI_grade_strings.fasta, per_sequence/<SEQ_ID>_ConSurf-AI_grades.tsv and README.txt."""
    shutil.rmtree(out_dir, ignore_errors=True)
    os.makedirs(os.path.join(out_dir, "per_sequence"))
    summary = summary_table(records, tables)
    files = []
    p = os.path.join(out_dir, f"{job}_ConSurf-AI_summary.tsv")
    summary.to_csv(p, sep="\t", index=False)
    files.append((p, "one row per input sequence: status, length, uncertain residues, grade counts, sequence"))
    all_cols = []
    if tables:
        all_res = pd.concat(tables.values(), ignore_index=True)
        all_cols = list(all_res.columns)
        p = os.path.join(out_dir, f"{job}_ConSurf-AI_all_residues.tsv")
        all_res.to_csv(p, sep="\t", index=False)
        files.append((p, "one row per residue of every sequence (grade, uncertainty, probabilities)"))
        p = os.path.join(out_dir, f"{job}_ConSurf-AI_grade_strings.fasta")
        with open(p, "w") as f:
            for sid, df in tables.items():
                f.write(f">{sid} sequence\n{''.join(df.AA)}\n>{sid} grades\n{''.join(map(str, df.GRADE))}\n"
                        f">{sid} uncertainty\n{''.join(map(str, df.UNCERTAINTY))}\n")
        files.append((p, "per sequence: the sequence, its grades and its rounded uncertainties as aligned strings "
                         "(one character per residue; uncertainty >= 2 = uncertain)"))
        for sid, df in tables.items():
            df.to_csv(os.path.join(out_dir, "per_sequence", f"{sid}_ConSurf-AI_grades.tsv"), sep="\t", index=False)
        files.append(("per_sequence/", "the per-residue table, one file per sequence (<SEQ_ID>_ConSurf-AI_grades.tsv)"))
    sum_cols = [c for c in summary.columns if not c.startswith("N_GRADE_")]
    if any(c.startswith("N_GRADE_") for c in summary.columns):
        sum_cols.insert(sum_cols.index("SEQUENCE"), "N_GRADE_1 ... N_GRADE_9")
    readme = write_readme(out_dir, "ConSurf-AI results (sequences)", f"{source}, {len(records)} sequence(s)",
                          [(os.path.basename(p.rstrip("/")) + ("/" if p.endswith("/") else ""), d) for p, d in files],
                          [("the summary table", sum_cols)] + ([("the residue tables", all_cols)] if all_cols else []))
    return summary, [p for p, _ in files if not p.endswith("/")] + [readme]
