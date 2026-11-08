"""Structure input (PDB ID, UniProt/AlphaFold DB, local file) and B-factor output files.

The structure files are never re-written from a parsed model: outputs are the original files with only the
B-factor column changed, so bonds, CONECT/LINK records, ligand definitions and metadata are preserved."""
import gzip
import json
import os
import re
import shutil
import urllib.error
import urllib.request

import gemmi

_PDB_RE = re.compile(r"^[0-9][A-Z0-9]{3}$")
_UNIPROT_RE = re.compile(r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")
_ALLOWED_AA = set("ACDEFGHIKLMNPQRSTVWYXUOBZ")


# ── Retrieval ────────────────────────────────────────────────────────────────
def _download(url, path, timeout=300):
    req = urllib.request.Request(url, headers={"User-Agent": "ConSurf-AI"})
    tmp = path + ".part"
    with urllib.request.urlopen(req, timeout=timeout) as r, open(tmp, "wb") as f:
        shutil.copyfileobj(r, f)
    os.replace(tmp, path)


def structure_format(path):
    """'cif' for mmCIF files, 'pdb' otherwise (judged from the content, not the extension)."""
    with open(path, errors="replace") as f:
        head = f.read(4096)
    return "cif" if re.search(r"^\s*(data_|loop_|_atom_site\.)", head, re.M) else "pdb"


def fetch_structure(query, out_dir):
    """PDB ID -> RCSB mmCIF (+ PDB-format file when the entry has one); UniProt accession / AF model ID ->
    AlphaFold DB model. Returns (files, label, description), files = {"cif": path, "pdb": path}."""
    os.makedirs(out_dir, exist_ok=True)
    q = query.strip().upper()
    af_id = None
    m = re.match(r"^AF-([A-Z0-9]+)-F(\d+)", q)
    if m:
        q, af_id = m.group(1), f"AF-{m.group(1)}-F{m.group(2)}"
    if _PDB_RE.match(q):
        files = {"cif": os.path.join(out_dir, f"{q}.cif")}
        try:
            _download(f"https://files.rcsb.org/download/{q}.cif", files["cif"])
        except urllib.error.HTTPError as e:
            raise ValueError(f"PDB entry '{q}' could not be downloaded from RCSB (HTTP {e.code}).") from None
        try:                                     # large entries have no PDB-format file
            pdb_path = os.path.join(out_dir, f"{q}.pdb")
            _download(f"https://files.rcsb.org/download/{q}.pdb", pdb_path)
            files["pdb"] = pdb_path
        except urllib.error.URLError:
            pass
        return files, q, f"PDB entry {q} (asymmetric unit, model 1)"
    if _UNIPROT_RE.match(q):
        req = urllib.request.Request(f"https://alphafold.ebi.ac.uk/api/prediction/{q}",
                                     headers={"User-Agent": "ConSurf-AI"})
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                entries = json.loads(r.read())
        except urllib.error.HTTPError as e:
            raise ValueError(f"No AlphaFold DB model found for UniProt accession '{q}' (HTTP {e.code}).") from None
        if not entries:
            raise ValueError(f"No AlphaFold DB model found for UniProt accession '{q}'.")
        if af_id:
            entries = [e for e in entries if (e.get("entryId") or e.get("modelEntityId") or "").upper() == af_id] or entries
        own = [e for e in entries if (e.get("uniprotAccession") or "").upper() == q]
        entry = (own or entries)[0]
        model_id = entry.get("entryId") or entry.get("modelEntityId") or f"AF-{q}-F1"
        files = {}
        for fmt, key in (("cif", "cifUrl"), ("pdb", "pdbUrl")):
            if entry.get(key):
                files[fmt] = os.path.join(out_dir, f"{model_id}.{fmt}")
                _download(entry[key], files[fmt])
        desc = entry.get("uniprotDescription") or ""
        org = entry.get("organismScientificName") or ""
        return files, model_id, f"AlphaFold DB {model_id} – {desc} ({org})".replace(" ()", "")
    raise ValueError(f"'{query}' is not a file, a PDB ID (e.g. 4HHB) or a UniProt accession (e.g. P69905).")


def resolve_input(query, work_dir):
    """Local structure file (optionally .gz) or PDB ID / UniProt accession. Returns (files, label, description)."""
    if os.path.isfile(query):
        name = os.path.basename(query)
        if name.lower().endswith(".gz"):
            name = name[:-3]
            path = os.path.join(work_dir, name)
            os.makedirs(work_dir, exist_ok=True)
            with gzip.open(query, "rb") as fin, open(path, "wb") as fout:
                shutil.copyfileobj(fin, fout)
        else:
            path = os.path.abspath(query)
        label = re.sub(r"\.(pdb|ent|cif|mmcif)$", "", name, flags=re.I)
        return {structure_format(path): path}, label, f"file {name}"
    return fetch_structure(query, work_dir)


# ── Parsing (sequence and numbering only) ────────────────────────────────────
def load_structure(path):
    st = gemmi.read_structure(path)
    st.setup_entities()
    while len(st) > 1:                       # first model only (NMR ensembles)
        del st[len(st) - 1]
    st.remove_alternative_conformations()
    return st


def _residue_letter(res):
    info = gemmi.find_tabulated_residue(res.name)
    if info is None or not info.is_amino_acid():
        return "X"
    letter = info.one_letter_code.upper()     # lowercase = modified residue, e.g. MSE -> 'm' -> 'M'
    return letter if letter in _ALLOWED_AA else "X"


def _protein_polymer(chain):
    poly = chain.get_polymer()
    if len(poly) == 0 or poly.check_polymer_type() not in (gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD):
        return None
    return poly


def res_key(r):
    return f"{r.seqid.num}{r.seqid.icode.strip()}"


def protein_chains(st):
    names = []
    for ch in st[0]:
        if ch.name not in names and _protein_polymer(ch) is not None:
            names.append(ch.name)
    return names


def chain_residues(st, chain_name):
    """Observed residues of a protein chain: list of dicts (key, number, icode, resname, aa)."""
    for ch in st[0]:
        poly = _protein_polymer(ch) if ch.name == chain_name else None
        if poly is not None:
            return [{"key": res_key(r), "number": r.seqid.num, "icode": r.seqid.icode.strip(),
                     "resname": r.name, "aa": _residue_letter(r)} for r in poly]
    return None


def query_subchains(st, chain):
    for ch in st[0]:
        poly = _protein_polymer(ch) if ch.name == chain else None
        if poly is not None:
            return {r.subchain for r in poly}
    return set()


# ── B-factor files ───────────────────────────────────────────────────────────
def cif_with_bfactors(path, out_path, chain, values_by_key, subchains, other=0.0):
    """Copy of an mmCIF file in which only _atom_site.B_iso_or_equiv is changed."""
    doc = gemmi.cif.read(path)
    tab = doc[0].find("_atom_site.", ["auth_asym_id", "auth_seq_id", "?pdbx_PDB_ins_code", "?label_asym_id",
                                      "B_iso_or_equiv"])
    if not tab:
        raise ValueError(f"{path}: no _atom_site table with auth_asym_id / auth_seq_id / B_iso_or_equiv.")
    s = gemmi.cif.as_string
    for row in tab:
        v = None
        if s(row[0]) == chain and (not row.has(3) or s(row[3]) in subchains):
            icode = s(row[2]).strip() if row.has(2) and row[2] not in ("?", ".") else ""
            v = values_by_key.get(s(row[1]) + icode)
        row[4] = f"{float(v if v is not None else other):.2f}"
    doc.write_file(out_path)


def pdb_with_bfactors(path, out_path, chain, values_by_key, other=0.0):
    """Copy of a PDB file in which only columns 61-66 (B-factor) of ATOM/HETATM records are changed."""
    with open(path, errors="replace") as fin, open(out_path, "w") as fout:
        for line in fin:
            stripped = line.rstrip("\r\n")
            if line.startswith(("ATOM  ", "HETATM")) and len(stripped) >= 54:
                body = stripped.ljust(66)
                v = values_by_key.get(body[22:26].strip() + body[26].strip()) if body[21] == chain else None
                line = body[:60] + f"{float(v if v is not None else other):6.2f}" + body[66:] + line[len(stripped):]
            fout.write(line)
