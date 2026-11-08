# CRISPR HTS Validator

Analyse and validate **pooled CRISPR screens** (from raw FASTQ reads or a count table to gene hit lists) and
**arrayed screens** (plate-reader data to hit wells and genes), with a plain-language quality check at every step.
A desktop app for the Mac: it runs a small local server and opens in your browser; nothing is uploaded anywhere.

Project page: https://paul-h-kim.com/software/crispr-hts-validator

## Pooled screens

- **Input**: MAGeCK or BAGEL count tables (detected automatically), or FASTQ/FASTQ.gz files plus a guide library.
  Read counting reproduces `mageck count` (exact match, automatic guide position), with optional one-mismatch
  matching, reverse-complement reads and scaffold search.
- **Samples**: sample sheet with reference (plasmid/T0) samples, conditions and replicates; comparisons; BAGEL
  screen-info files can be imported.
- **Quality checks**, each pass / warn / fail with one sentence why: reads per guide, mapped fraction, zero-count
  guides, Gini index, guide coverage and library skew, replicate agreement, sample PCA, and essential vs
  non-essential gene separation (ROC AUC and precision–recall with the Hart lab CEGv2 / NEGv1 lists), including a
  warning when a timepoint is too early.
- **Normalisation**: median ratio (default), total count, or control guides.
- **Scoring**: guide log fold changes and a MAGeCK-style test — mean–variance model, guide p-values, and α-RRA gene
  ranking with FDR — ported step by step from MAGeCK 0.5.9.5, including its random-number generator, so results
  can be checked against `mageck test`.
- **Results**: volcano and rank plots, searchable gene table, per-gene view with every guide, export to CSV,
  MAGeCK-format text, or Excel.

## More scoring, copy number, pathways, PDF

- **BAGEL2 Bayes factors** for essentiality, ported from the Hart lab's BAGEL.py (10-fold cross-validation, same
  random draws for a given seed), with precision–recall against the reference essential / non-essential genes
  (human CEGv2/NEGv1, mouse CEG/NEG).
- **drugZ** for treated vs untreated screens (paired or unpaired), ported from the Hart lab's drugz.py.
- **Copy-number correction** in the style of CRISPRcleanR, from guide positions alone (a CRISPRcleanR library file
  works); off by default.
- **Pathway enrichment** offline: Gene Ontology and Reactome (CC BY 4.0, bundled in `crisprhts/data/genesets`,
  built by `scripts/build_genesets.py`) or your own .gmt; over-representation of hits or a ranked test.
- **PDF report** for pooled and arrayed screens: data, every check with its reason, results, methods.

## Arrayed screens

- **Input**: cellHTS-style sets (plate list, one file per plate and replicate, plate configuration with wildcards,
  annotation, screen log, description), per-well tables (plate, well, one column per replicate or channel), and
  plate-reader grids such as EnVision exports. Messy files are tolerated and what was skipped is reported.
- **Plate map**: well types from the plate configuration, roles (sample / negative / positive control / other /
  empty), extra wells marked by hand, assay direction (lower, higher or either signal), channel.
- **Quality checks**: plate-to-plate medians, Z′ and SSMD of the controls (per plate with ≥ 2 wells of each,
  otherwise pooled, and the app says so), edge effect, row and column effects, replicate agreement; raw plate
  heatmaps and row/column profiles.
- **Normalisation per plate**: robust z (default), B-score (median polish), plate median, percent of control,
  NPI, and the classic plate z-score used for the published KcViab analysis.
- **Hits**: combined score over replicates with a threshold every replicate must also support; gene calls when
  several wells target a gene. Plate maps of scores, replicate scatter, score distribution, hits per plate,
  exports.

## Run it

**Download** the Mac app from the project page (free, not registered with Apple: the first time, open System
Settings → Privacy & Security and click *Open Anyway*; see `macos/HOW-TO-OPEN.txt`). Or, from this repository:

```bash
macos/build_app.sh          # builds "CRISPR HTS Validator.app" here (code and example data inside it)
./launch.sh                 # or run it straight from the folder
```

Needs macOS 12 or newer and Python 3.9 or newer. The first start makes a private Python environment in
`~/Library/Application Support/CRISPR HTS Validator/venv` (numpy, pandas, scipy, statsmodels, openpyxl,
matplotlib, reportlab; needs internet once). The app opens at http://localhost:8788 and quits completely when its
window closes. Eight example screens (pooled and arrayed) open with one click.

## Tests and validation

```bash
"$HOME/Library/Application Support/CRISPR HTS Validator/venv/bin/python" -m pip install pytest
"$HOME/Library/Application Support/CRISPR HTS Validator/venv/bin/python" -m pytest tests
```

- `tests/test_units.py` — every QC metric, normalisation and the RRA pieces on small synthetic data with known answers.
- `tests/test_validation.py` — agreement with the official MAGeCK 0.5.9.5 outputs in `test-data/reference/` (made by
  `scripts/make_reference.sh`): identical counts, statistics equal to MAGeCK's printed precision, identical hit lists.
- `tests/test_arrayed.py` — plate readers, normalisations, B-score, Z′/SSMD, edge and row/column checks and hit
  calls on synthetic plates, and agreement with the published KcViab scores (cellHTS2 `Analysis2003.txt`).
- `tests/test_methods03.py` — BAGEL2, drugZ, copy-number correction, enrichment and the PDF report on synthetic data
  with known answers, and against the official BAGEL2 and drugZ outputs (made by `scripts/make_reference_hart.sh`; identical to their printed precision).
- `tests/test_app.py` — the app end to end, and quitting when the window closes.

## Test data and credits

`test-data/` holds public data only, each with its licence: the Hart lab BAGEL repository (HAP1 TKO and TKOv3
screens, CEGv2/NEGv1 gene lists; MIT), the MAGeCK demos (BSD), and the cellHTS2 example screens (Artistic-2.0).
Methods: Li W et al., *Genome Biology* 15:554 (2014) — MAGeCK; Hart T et al., *G3* 7:2719 (2017) — CEGv2/NEGv1;
Kim E & Hart T, *Genome Medicine* 13:2 (2021) — BAGEL2; Colic M et al., *Genome Medicine* 11:52 (2019) — drugZ;
Iorio F et al., *BMC Genomics* 19:604 (2018) — CRISPRcleanR. Gene sets: Gene Ontology Consortium and Reactome
(CC BY 4.0); see `crisprhts/data/genesets/SOURCES.txt`.

## Licence

MIT — Copyright (c) 2026 Paul H. Kim
