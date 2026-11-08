# Reference outputs from the official tools

Made once with `scripts/make_reference.sh`, which records the exact versions and command lines in
`mageck/COMMANDS.txt`. The validation tests in `tests/test_validation.py` compare this tool against these files.

| Folder | Tool | Version | Licence |
|---|---|---|---|
| `mageck/` | MAGeCK (`mageck count`, `mageck test`) | 0.5.9.5 | BSD |
| `bagel/` (from version 0.3) | BAGEL2 (`BAGEL.py fc`, `bf`, `pr`) | GitHub hart-lab/bagel | MIT |

## Result (2026-11-07, MAGeCK 0.5.9.5)

All 13 validation tests pass. Agreement per comparison:

| Data | Genes | Guides | Counts | Largest relative difference (guide p, RRA score) | Gene p-values equal | Top-100 overlap (depleted / enriched) | Hits at FDR < 0.05, ours = MAGeCK (depl. / enr.) |
|---|---|---|---|---|---|---|---|
| demo 2 (FASTQ) | 2,356 | 2,550 | identical (reads, mapped, zero-count too) | 1.5e-5, 3.6e-5 | 100% | 100 / 100 | 0 / 0 |
| demo 1 | 100 | 999 | — | 4.7e-5, 3.9e-5 | 100% | 100 / 100 | 2 / 0 |
| demo 3 | 1,000 | 9,999 | — | 4.8e-5, 4.7e-5 | 100% | 100 / 100 | 36 / 0 |
| TKOv3 T3 vs T0 | 18,056 | 70,877 | — | 4.9e-5, 4.9e-5 | 100% | 100 / 100 | 107 / 0 |
| TKOv3 T18 vs T0 | 18,056 | 70,877 | — | 5.0e-5, 4.8e-5 | 100% | 100 / 100 | 844 / 4 |
| TKOv3 T18 starved vs T0 | 18,055 | 70,877 | — | 4.9e-5, 4.9e-5 | 100% | 100 / 100 | 737 / 15 |

Differences of up to 5e-5 are MAGeCK printing 5 significant digits. Gene p-values match because the permutation
uses MAGeCK's own random-number generator and seed. Deliberate differences: none.

## Arrayed: KcViab against the published 2003 analysis

`test-data/cellhts2/KcViab/Analysis2003.txt` holds the published per-well scores (Boutros et al., *Science* 2004;
shipped with cellHTS2). The 2003 score turned out to be a **classic plate z-score**: (value − plate mean) ÷ plate SD
over all 384 wells, sign reversed (lower viability = higher score), mean of the two replicates, with the wells
flagged in the screen log kept. With those settings this tool reproduces 99.98% of the 21,060 published scores to
two decimals (r = 0.99995). Hits below are wells with score ≥ 2 (417 published).

| Method in this tool | r with published scores | Hit wells | Published hits recovered | Hit genes |
|---|---|---|---|---|
| Classic plate z, as in 2003 (no replicate rule, flagged wells kept) | 0.9999 | 415 | 415 (99.5%) | 307 |
| Classic plate z + replicate agreement | 0.9998 | 420 | 414 (99.3%) | 305 |
| Robust z (default) | 0.981 | 668 | 415 (99.5%) | 500 |
| Plate median → robust z | 0.982 | 716 | 414 (99.3%) | 506 |
| B-score | 0.705 | 1,411 | 408 (97.8%) | 936 |
| Percent of control | 0.844 | 692 | 386 (92.6%) | 422 |
| NPI | 0.325 | 1,603 | 372 (89.2%) | 1,248 |

Why they differ (deliberate):
- *Robust z* uses the median and MAD of the sample wells. Real hits inflate the plate SD, so the MAD is smaller and
  the same threshold of 2 is less strict than in 2003: it keeps nearly all published hits and adds about 250 wells.
- *B-score* also removes row and column effects. KcViab has them: outer columns read 12–13% above the plate centre
  and rows show a U-shape (row/column effects explain a median 41% of within-plate variance). The 2003 analysis did
  not correct for them, so B-scores rank wells differently.
- *Percent of control* and *NPI* scale each plate by its control wells. KcViab has one positive and one negative
  control well per plate, so this scaling is noisy; these methods suit screens with several control wells per plate.

## BAGEL2 and drugZ (version 0.3)

`scripts/make_reference_hart.sh` runs the official Hart lab programs (BAGEL.py and drugz.py from GitHub, in their
own environment with numpy < 2 and pandas 2.1.4) on HAP1 TKOv3: BAGEL2 `fc`, `bf` (seed 1) and `pr` for T3, T18 and
T18 starved against T0; drugZ for T18 starved (treated) vs T18 (untreated), paired and unpaired. The tests in
`tests/test_methods03.py` then require Bayes factors equal to BAGEL's three printed decimals (same seed, same
cross-validation draws), the same number of genes with BF > 5, drugZ normZ and sumZ equal to two decimals, and the
same FDR < 0.05 calls. Until the script has been run these tests are skipped.

**Result (run 2026-11-07; BAGEL commit 53388ad, drugZ commit eb15d34; numpy 1.26.4, pandas 2.1.4):**

| Output | Genes | Agreement with the official program | Calls |
|---|---|---|---|
| BAGEL2 fold changes | 70,877 guides | identical (max difference 0) | – |
| BAGEL2 BF, T3 | 18,056 | identical at 3 decimals for 100% of genes | BF > 5: 901 = 901, same genes |
| BAGEL2 BF, T18 | 18,056 | identical at 3 decimals for 100% of genes | BF > 5: 1,796 = 1,796, same genes |
| BAGEL2 BF, T18 starved | 18,056 | identical at 3 decimals for 100% of genes | BF > 5: 1,706 = 1,706, same genes |
| drugZ paired | 18,055 | sumZ and normZ identical at 2 decimals; ranks identical; p-values within printing precision (0.5%) | FDR < 0.05: 54 synthetic, 10 suppressor, same genes |
| drugZ unpaired | 18,055 | as above | FDR < 0.05: 49 synthetic, 7 suppressor, same genes |


Data added for 0.3: `hart-bagel/tkov3-example/TKOv3-guide-positions.txt` (TKOv3 guide coordinates from the BAGEL
repository's CRISPRcleanR example, MIT) for copy-number correction.
