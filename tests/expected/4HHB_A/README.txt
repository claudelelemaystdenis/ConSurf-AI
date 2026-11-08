ConSurf-AI results
==================

Input:  file 4HHB.cif, chain A
Date:   2026-11-08 12:43 UTC
Model:  ConSurf-AI 1.0.0 (attention network on ESM2-3B embeddings)

Files
-----
4HHB_A_ConSurf-AI_grades.tsv     per-residue grades and uncertainties (one row per residue)
4HHB_A_ConSurf-AI_structure.cif  input structure file; B-factor column = grade on chain A (0 elsewhere)
4HHB_A_ConSurf-AI_pymol.pml      PyMOL script: ConSurf colouring (run next to the structure file)
4HHB_A_ConSurf-AI_chimerax.cxc   ChimeraX script: ConSurf colouring (open next to the structure file)

Columns of 4HHB_A_ConSurf-AI_grades.tsv
---------------------------------------
POS                      position in the sequence (1 = first residue)
CHAIN                    chain ID in the structure file
RESIDUE_NUMBER           residue number in the structure file (with insertion code, if any)
RESNAME                  residue name (three-letter code)
AA                       amino acid (one-letter code)
GRADE                    ConSurf-AI conservation grade: 1 = most variable ... 9 = most conserved
UNCERTAINTY_RAW          predicted uncertainty of the grade, 0-8 (lower = better)
UNCERTAINTY              UNCERTAINTY_RAW rounded to the nearest integer
CONFIDENT                True if UNCERTAINTY < 2; False = uncertain prediction (same role as ConSurf's 'insufficient data' flag)
P_GRADE_1 ... P_GRADE_9  predicted probability of each grade (they sum to 1)

Tables are tab-separated (.tsv). Excel: File > Open (or drag the file onto Excel). R: read.delim(file). Python: pandas.read_csv(file, sep='\t').
