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

## Screen design checks (version 1.1)

The checks that no amount of sequencing can replace, from the numbers you enter on the Samples page:

- **Library representation**: cells carrying each guide after selection, at the tightest passage and at harvest
  (at least 500× at every step; below that, guides are lost by sampling, which looks exactly like a real dropout).
- **MOI**: with the Poisson fraction of transduced cells that carry a single guide (about 86% at MOI 0.3).
- **Guides per gene**: the distribution, and the genes with fewer than four independent guides.
- **Non-targeting controls**: how many are in the library, whether they sit at a log₂ fold change of zero (a shift
  means the normalisation moved, not the cells), and the noise floor their spread sets — the smallest gene effect
  that can be told apart from no effect.
- **Spearman** rank correlation beside Pearson for replicate agreement, so a few huge fold changes cannot carry it.

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

## How it looks

Five styles, chosen from **Style** in the top bar and remembered: **Sequencer** (the original deep teal, with
reads drifting behind the page), **Paper** (a journal page — warm stock, serif headings, good for figures bound
for a manuscript), **Clean** (Apple-like greys and one blue), **Plate** (after the arrayed heatmaps — cool grid,
square wells, cyan and magenta) and **Contrast** (heavy rules, no translucency, no background animation, for a
projector or tired eyes). Each works in light and dark, so the ☀ Light | ☾ Dark pill keeps working, and every
text colour in all ten combinations clears the 4.5:1 contrast ratio. Charts and the background animation take
their colours from whichever style is on.

## Comparing two screens (version 1.5)

Put two of the project's comparisons side by side — a cell line against its isogenic control, or treated
against untreated — and see which dependencies they share and which belong to one of them.

- Each screen keeps **its own reference** sample, so a two-arm experiment is compared the way it was run.
- Both are put on **one scale** before anything is subtracted: the median core-essential gene (CEGv2) at −1 and
  the median non-essential gene (NEGv1) at 0, which removes the difference in dynamic range that timepoint,
  cell line and depth introduce. A screen without enough reference genes is said to be unscaled, not scaled
  silently.
- The **differential test** works guide by guide, so each guide's own efficiency cancels: the difference is
  z-scored against an empirical null estimated in a sliding window of guides ordered by abundance (drugZ's
  method), combined per gene as sumZ / √n, two-sided, with Benjamini–Hochberg FDR.
- **Views**: screen against screen, a differential volcano, a concordance summary (Pearson, Spearman, shared
  dependencies by Jaccard), a sortable gene table with search, and every guide of a chosen gene in both screens.
- **Quality gate**: each screen's essential-gene separation is shown beside the result, with a plain warning
  when one of them is too weak to trust — an early timepoint looks like this.
- **Built-in example**: the HAP1 TKOv3 screen at T18 against T18 starved. It finds 30 starvation-specific
  dependencies led by TSC2, TSC1 and KEAP1, and 36 shared ones that deepen under starvation.
- Exports to CSV, a sheet in the Excel workbook, and a section in the PDF report.
- **Confirmed defaults** (22 November 2026): FDR 0.05, smallest reported change 0.5, at least 3 guides, a
  dependency at a scaled effect of −0.5 or below. Guide agreement is shown in the table, not used to filter.

## Comparing screens in two projects, and two libraries (version 1.6)

A screen in a **separate count table** — another experiment, another lab, another library — can be compared
with this one. Nothing is merged: each keeps its own guides, samples and normalisation, and the two are joined
only by gene name and by the essential-gene scale.

This route is blunter than the one above, unavoidably. Two libraries target a gene with different sequences, so
the guide-by-guide cancellation is gone. Split TKOv3 in half and the cost is plain: **EEF2's guides in one half
read −0.6 and −0.2 while its guides in the other read −5.2 and −7.3** — a measured effect of −0.4 against −5.3,
from one cell line at one timepoint. Everything below follows from taking that seriously.

- **Gene names** are matched exactly after case, spacing and Ensembl version suffixes are normalised, with an
  optional two-column mapping file for the rest. No alias table is bundled and nothing is downloaded, so names
  that do not line up are reported and diagnosed (it names the identifier type on each side), never guessed at.
- **The uncertainty of a gene's effect** is built from the spread of its guides, shrunk toward the spread a gene
  of that strength usually shows — a mean–variance trend, as MAGeCK and limma-voom fit, because every guide of a
  harmless gene reads near zero whether it cuts or not while against an essential gene a dead guide reads 0 and a
  good one −6 — floored by the spread its read abundance implies, and weighted by degrees of freedom with the
  prior strength **estimated from the data** by Smyth's (2004) moment method rather than chosen by hand.
- **A t distribution, not a normal one.** With four guides a side the standard error is itself barely known, and
  treating it as exact manufactures significance.
- **The scale's own error is propagated** into every gene in proportion to how deep that gene is, so a thin
  library's deepest genes are its least certain. Below **30 reference genes** in either set the comparison is
  refused: for two different libraries the scale is the whole of the bridge, not a cosmetic step.
- **The resolution is computed and printed.** A gene with the typical guide spread needs a difference of *x* to
  clear the FDR bar across every gene tested; read a null result as "no difference this large", not "no
  difference". Guides per gene is the lever, and the figure falls roughly as its square root.
- **It tells you when you do not need it.** Guides are matched on the protospacer inside the ID, so the same
  library under two ID conventions is recognised and the app says to use the paired route instead.
- **Built-in example**: the same HAP1 screen held in two tables (TKOv3 `GENE_20mer` IDs against BAGEL bare
  20-mers — no ID in common). Pearson 0.996 and **nothing called different** out of 17,620 genes, with the
  shared-library warning raised.
- Exports to CSV, a sheet in the Excel workbook, and a section in the PDF report.
- **Confirmed defaults** (24 November 2026): FDR 0.05, smallest reported change **1.0** (twice the
  within-project figure, because a smaller cross-library difference is inside the systematic error between
  libraries, and below what four guides a side can resolve), at least 3 guides in each library, a dependency at
  a scaled effect of −0.5 or below. A screen with fewer than 30 reference genes of either kind is **refused**
  rather than warned about, and the empirical null may only ever be **widened, never narrowed** — on two tables
  of one HAP1 screen the statistic scatters at 0.09 where 1 is expected, and the app says so in words rather
  than tightening the test.

Not built: merging three or more projects, and comparing across species.

## Thresholds

Every pass / warn / fail line is a scientific choice. The defaults follow published guidance and were confirmed on
7 November 2026: 300/100 reads per guide; 80%/65% mapped; Gini 0.10/0.20 at baseline; 1%/5% zero-count guides in
the reference; 90/10 skew 10/20; replicate Pearson r 0.60/0.40 (Spearman shown beside it); essential-gene AUC
0.90/0.80; 500×/200× cells per guide; MOI 0.4/0.6; 4/2 guides per gene; 100/20 non-targeting controls. Arrayed:
plate median ratio 2/5, Z′ 0.5/0, |SSMD| 3/2, edge 10%/20%, row-and-column share 10%/25%, replicate r 0.70/0.40,
default normalisation **robust z**, hit when the combined score passes 2 and every replicate passes 1.

Change any of them under *Thresholds* on the Quality page (pooled and arrayed each have their own panel). A changed threshold is remembered in
`~/Library/Application Support/CRISPR HTS Validator/preferences.json` (that file holds nothing else, and
*Restore defaults* removes it). Set `CH_PREFS` to use a different file.

### Confirmed method defaults

- **BAGEL2**: a gene is essential at **BF > 5** (BAGEL2's own tutorial value). On HAP1 TKOv3 T18 that is 1,796 genes,
  94.0% recall, no false positives; on the too-early T3 timepoint it still holds 98% precision, where BF > 1 would
  admit 55 non-essential genes. The results page accepts any cut-off for the view; exports use 5.
- **Copy-number correction**: **off unless asked for**, and segments must cover at least **3 genes and 10 guides**
  (CRISPRcleanR's defaults, so results stay comparable). At those settings it removes planted amplicons of 6-81
  genes to within about 12%, but on the near-diploid HAP1 line it also flattens clusters of neighbouring essential
  genes and the essential-gene AUC falls from 0.998 to 0.994 - hence off by default, for cancer lines.

## Speed

A genome-wide MAGeCK-style test on the 70,877-guide TKOv3 screen takes about **6 seconds** on an idle
recent Mac with several cores, against 25 before version 1.4 (any timing depends on the machine and on what
else it is doing). Nearly all of that time is the permutation null, and nearly all of *that* goes to
the few genes with many guides: in TKOv3 a single gene with 94 guides costs more than all 17,247 genes with four
put together. Those permutation blocks now run on several cores, each starting from the generator state it would
have reached by drawing in turn (jump-ahead), so the random numbers, and therefore the p-values, are exactly
those of a single-threaded run — `CH_RRA_WORKERS=1` forces that path and gives bit-identical output.

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
- `tests/test_design_qc.py` — library representation, MOI (against the Poisson distribution), guides per gene,
  non-targeting controls, and Spearman against scipy.
- `tests/test_prefs.py` — thresholds are remembered, survive a damaged file and a renamed threshold, and can be
  restored; the confirmed defaults are checked by value.
- `tests/test_units.py` also checks that jump-ahead reproduces the generator's stream and that the parallel and
  single-threaded gene tests agree exactly.
- `tests/test_compare.py` — the within-project comparison: planted screen-specific genes recovered, identical
  screens calling nothing, the scaling landing on −1 and 0, and agreement with drugZ when both screens share a
  reference (Spearman 0.96, same top genes).
- `tests/test_cross.py` — the cross-project comparison on two libraries of known truth (planted genes recovered,
  resolution improving with guides per gene, a 1.7× difference in dynamic range absorbed by the scale and nothing
  called), gene-name matching and mapping files, and two real-data checks: **splitting TKOv3's four guides per
  gene into two disjoint halves**, where TSC2 and TSC1 still rank in the top 15 but two guides a side cannot clear
  a genome-wide bar and the app says so instead of pretending; and **one HAP1 screen read from two count tables**
  with no guide ID in common, where nothing at all is called different.
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
