ConSurf-AI results (sequences)
==============================

Input:  file dihydrofolate_reductases.fasta, 5 sequence(s)
Date:   2026-11-08 12:43 UTC
Model:  ConSurf-AI 1.0.0 (attention network on ESM2-3B embeddings)

Files
-----
dihydrofolate_reductases_ConSurf-AI_summary.tsv          one row per input sequence: status, length, uncertain residues, grade counts, sequence
dihydrofolate_reductases_ConSurf-AI_all_residues.tsv     one row per residue of every sequence (grade, uncertainty, probabilities)
dihydrofolate_reductases_ConSurf-AI_grade_strings.fasta  per sequence: the sequence, its grades and its rounded uncertainties as aligned strings (one character per residue; uncertainty >= 2 = uncertain)
per_sequence/                                            the per-residue table, one file per sequence (<SEQ_ID>_ConSurf-AI_grades.tsv)

Columns of the summary table
----------------------------
SEQ_ID                   sequence identifier (first word of the FASTA header, made unique)
DESCRIPTION              full FASTA header (or the name given with --name)
LENGTH                   number of residues after cleaning
STATUS                   ok, skipped (empty or longer than 2000 residues) or error
NOTE                     changes made to the sequence (removed gaps/stops, unknown letters read as X)
N_UNCERTAIN              number of residues with an uncertain prediction (CONFIDENT = False)
FRACTION_UNCERTAIN       N_UNCERTAIN / LENGTH
MEAN_GRADE               mean grade over the sequence
N_GRADE_1 ... N_GRADE_9  number of residues with each grade
SEQUENCE                 the sequence as used for the prediction

Columns of the residue tables
-----------------------------
SEQ_ID                   sequence identifier (first word of the FASTA header, made unique)
POS                      position in the sequence (1 = first residue)
AA                       amino acid (one-letter code)
GRADE                    ConSurf-AI conservation grade: 1 = most variable ... 9 = most conserved
UNCERTAINTY_RAW          predicted uncertainty of the grade, 0-8 (lower = better)
UNCERTAINTY              UNCERTAINTY_RAW rounded to the nearest integer
CONFIDENT                True if UNCERTAINTY < 2; False = uncertain prediction (same role as ConSurf's 'insufficient data' flag)
P_GRADE_1 ... P_GRADE_9  predicted probability of each grade (they sum to 1)

Tables are tab-separated (.tsv). Excel: File > Open (or drag the file onto Excel). R: read.delim(file). Python: pandas.read_csv(file, sep='\t').
