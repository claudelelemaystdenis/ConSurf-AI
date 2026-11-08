# ConSurf-AI

ConSurf-AI predicts **per-residue evolutionary conservation grades** on the ConSurf scale
(**1** = variable … **9** = conserved) from a single protein sequence or structure, **without building a
multiple sequence alignment**. It runs a small attention network on ESM2-3B embeddings, and it also predicts
an uncertainty for each residue (the equivalent of ConSurf's "insufficient data" flag).

## Ways to run it

**Colab: structure.** One structure (PDB ID, UniProt accession → AlphaFold DB model, or your own file) →
grades mapped on the 3D structure, PyMOL / ChimeraX scripts.

<a href="https://colab.research.google.com/github/claudelelemaystdenis/ConSurf-AI/blob/main/colab/ConSurf_AI.ipynb"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open the structure notebook in Colab" height="48"></a>

**Colab: sequences.** One sequence or a FASTA file with many sequences → tables of grades.

<a href="https://colab.research.google.com/github/claudelelemaystdenis/ConSurf-AI/blob/main/colab/ConSurf_AI_sequences.ipynb"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open the sequence notebook in Colab" height="48"></a>

The Colab notebooks need no installation; a free T4 GPU is enough.

**Locally** ([below](#local-installation)). Both, from the command line or Python: many proteins, whole proteomes,
pipelines, or data that must stay on your machine.

## Local installation

```bash
git clone https://github.com/claudelelemaystdenis/ConSurf-AI.git
cd ConSurf-AI
conda env create -f environment.yml      # or, without conda: pip install -e .   (Python ≥ 3.9)
conda activate consurf-ai
```

The ConSurf-AI weights (21 MB) are included in [`weights/`](weights/). The first run also downloads
ESM2-3B (~11 GB) from Hugging Face into `~/.cache/huggingface`.

**Hardware.** The GPU is used automatically when there is one:
- **Linux with an NVIDIA GPU:** works out of the box (pip installs PyTorch with CUDA).
- **Mac with Apple silicon (M1–M4):** uses the Mac's GPU (MPS). At least 16 GB of memory is recommended.
- **Windows with an NVIDIA GPU:** after creating the environment, install the CUDA build of PyTorch with the
  pip command from [pytorch.org](https://pytorch.org/get-started/locally/).

Without a GPU, ConSurf-AI runs on the CPU. It needs ~12 GB of RAM and is much slower (from several seconds to
about a minute per protein, depending on the machine). On NVIDIA GPUs with less than 14 GB of memory, ESM2 runs in fp16.

## Usage

```bash
# Sequences -> tables (the same files as the sequence notebook)
consurf-ai proteins.fasta                                   # one or many sequences
consurf-ai --sequence MVLSPADKTNVKAAWGKVGAH --name my_protein  # a single sequence

# Structure -> grades mapped on the structure (the same files as the structure notebook)
consurf-ai 4HHB --chain A                # PDB ID, fetched from RCSB
consurf-ai P69905                        # UniProt accession, AlphaFold DB model
consurf-ai my_model.cif --chain all      # local PDB / mmCIF file (.gz accepted)

consurf-ai --help
```

Results go to a folder in `consurf_ai_results/` (change it with `-o`): `<name>/` for sequences (named after
the FASTA file or `--name`), `<structure>_<chain>/` for each structure chain.

From Python:

```python
from consurf_ai import ConSurfAI

model = ConSurfAI()                      # loads ESM2-3B + ConSurf-AI once
pred = model.predict("VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF")
pred["grades"], pred["uncertainty"], pred["confident"]
```

## Output files

Tables are tab-separated (`.tsv`): open them in Excel with File → Open, in R with `read.delim()`, in Python with
`pandas.read_csv(file, sep="\t")`. Every results folder has a `README.txt` describing the run and every column.

**Sequences** (`consurf_ai_results/<name>/`)

| File | Content |
|---|---|
| `*_summary.tsv` | One row per input sequence: ID, FASTA header, length, status (ok / skipped / error), changes made to the sequence, number and fraction of uncertain residues, mean grade, number of residues in each grade, sequence |
| `*_all_residues.tsv` | **All predictions together**, one row per residue of every sequence: sequence ID, position, amino acid, **grade (1–9)**, **uncertainty** (raw 0–8 and rounded; lower = better), confidence flag (rounded uncertainty < 2), probability of each grade |
| `*_grade_strings.fasta` | For each sequence: the sequence, its grades and its rounded uncertainties as three aligned FASTA records |
| `per_sequence/*_grades.tsv` | The per-residue table, one file per sequence |

**Structures** (`consurf_ai_results/<structure>_<chain>/`)

| File | Content |
|---|---|
| `*_grades.tsv` | One row per residue: chain numbering, amino acid, grade, uncertainty, confidence flag, probability of each grade |
| `*_structure.pdb` / `.cif` | The input structure file, unchanged except for the **B-factor** column, which holds the grade on the query chain (0 elsewhere) |
| `*_pymol.pml` / `*_chimerax.cxc` | Scripts that colour the structure with the ConSurf colours (open them next to the structure file) |

Predictions use the residues present in the structure (first model, first alternate conformation); missing
residues are not filled in from SEQRES. Sequences and chains of up to 2000 residues are supported.

## Checking your installation

```bash
python tests/run_tests.py            # add --device cpu to force the CPU
```

This runs ConSurf-AI on two inputs in [`tests/input/`](tests/input/) and compares the results with the expected
outputs in [`tests/expected/`](tests/expected/), which were produced on an NVIDIA A100 GPU (ESM2 in float32):

1. **structure path:** `4HHB.cif` (human haemoglobin, chain A);
2. **sequence path:** `dihydrofolate_reductases.fasta` (5 dihydrofolate reductases).

It reports, for each path, whether the grades are identical (and lists the residues whose grade differs) and how
much the uncertainties and probabilities differ. Small numerical differences between machines are normal (CPU vs
GPU: at most 0.001 in our tests) and are accepted up to 0.01.

## Repository layout

```
consurf_ai/        Python package (model, ESM2 embedding, structure and sequence I/O, output files, command line)
weights/           ConSurf-AI network weights + SHA256SUMS (checksum: sha256sum -c SHA256SUMS)
colab/             Google Colab notebooks (structure, sequences)
tests/             Reference test: input/, expected/ and run_tests.py
```

## Citation

If you use ConSurf-AI, please cite:

TO CHANGE
> <Authors>. <Paper title>. *<Journal>* (2026). https://doi.org/<DOI>

```bibtex
@article{consurf_ai_2026,
  author  = {<Family name>, <Given name> and <Family name>, <Given name>},
  title   = {<Paper title>},
  journal = {<Journal>},
  year    = {2026},
  doi     = {<DOI>}
}
```

## License

The code and the model weights are released under the [Apache License 2.0](LICENSE) (see also [`NOTICE`](NOTICE)).
ESM2 (used to compute the embeddings) is © Meta Platforms, Inc., also under the MIT License.
