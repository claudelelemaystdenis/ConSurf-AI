# ConSurf-AI

ConSurf-AI predicts **per-residue evolutionary conservation grades** on the ConSurf scale
(**1** = variable … **9** = conserved) from a single protein sequence or structure, **without building a
multiple sequence alignment**. It runs a small attention network on ESM2-3B embeddings, and it also predicts
an uncertainty for each residue (the equivalent of ConSurf's "insufficient data" flag).

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/claudelelemaystdenis/ConSurf-AI/blob/main/colab/ConSurf_AI.ipynb)

## Two ways to run it

| | When to use it |
|---|---|
| **Google Colab** (badge above) | No installation. One structure at a time, with an interactive 3D view. A free T4 GPU is enough. |
| **Locally** (below) | Many proteins, whole proteomes, pipelines, or data that must stay on your machine. |

## Local installation

```bash
git clone https://github.com/claudelelemaystdenis/ConSurf-AI.git
cd ConSurf-AI
conda env create -f environment.yml      # or: pip install -e .   (into an environment with PyTorch)
conda activate consurf-ai
```

The ConSurf-AI weights (21 MB) are included in [`weights/`](weights/). The first run also downloads
ESM2-3B (~11 GB) from Hugging Face into `~/.cache/huggingface`.

**Hardware.** A GPU with ≥ 8 GB of memory is recommended (ESM2 runs in fp16 on GPUs with < 14 GB). CPU-only
also works: it needs ~12 GB of RAM and takes about a minute for a 150-residue chain.

## Usage

```bash
consurf-ai 4HHB --chain A                # PDB ID, fetched from RCSB
consurf-ai P69905                        # UniProt accession, AlphaFold DB model
consurf-ai my_model.cif --chain all      # local PDB / mmCIF file (.gz accepted)
consurf-ai proteins.fasta                # sequences only: one CSV per sequence
consurf-ai --help
```

Results go to `consurf_ai_results/` (change it with `-o`).

From Python:

```python
from consurf_ai import ConSurfAI

model = ConSurfAI()                      # loads ESM2-3B + ConSurf-AI once
pred = model.predict("VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF")
pred["grades"], pred["uncertainty"], pred["confident"]
```

## Output files

| File | Content |
|---|---|
| `*_grades.csv` | One row per residue: numbering, amino acid, **grade (1–9)**, **uncertainty** (raw 0–8 and rounded; lower = better), confidence flag (rounded uncertainty < 2), probability of each grade |
| `*_structure.pdb` / `.cif` | The input structure file, unchanged except for the **B-factor** column, which holds the grade on the query chain (0 elsewhere) |
| `*_pymol.pml` / `*_chimerax.cxc` | Scripts that colour the structure with the ConSurf colours (open them next to the structure file) |

Predictions use the residues present in the structure (first model, first alternate conformation); missing
residues are not filled in from SEQRES. Chains of up to 2000 residues are supported.

## Repository layout

```
consurf_ai/        Python package (model, ESM2 embedding, structure I/O, output files, command line)
weights/           ConSurf-AI network weights + SHA256 checksum
colab/             Google Colab notebook
tests/             Quick CPU tests (no ESM2 download):  pytest tests
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
