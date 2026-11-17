# iCARE — reproducibility code

Analysis code for the three computational parts of the manuscript:

1. **`01_leakage/`** — does the iCARE pipeline actually keep the test set out of
   feature screening and preprocessing? A controlled experiment in which the
   correct and the leaky workflow differ in exactly one line.
2. **`02_head2head/`** — how does iCARE compare with **caret**, **tidymodels** and
   **mlr3** on the same splits: discrimination, calibration, selected features,
   runtime and amount of code.
3. **`03_clinical_scenarios/`** — three NHANES screening scenarios (undiagnosed
   dysglycemia, albuminuria, undiagnosed hypertension), each developed on
   2007–2008 and validated on the independent 2009–2010 cycle.

`R/` holds the shared setup and helpers; nothing else is shared, so any module can
be read and run on its own.

```
.
├── run_all.R                     run everything (or one module)
├── data/                         optional local overrides (see data/README.md)
├── R/
│   ├── setup.R                   root discovery, output folders, AUC / Brier / CI helpers
│   └── stat_replay.R             rebuild and replay an icare Module-1 Stat object
├── 01_leakage/                   [README](01_leakage/README.md)
│   ├── leakage_validation.R      experiments A–E, block by block
│   ├── make_composite_figure.R   the 2 × 2 composite figure
│   └── run_leakage.R             batch entry point
├── 02_head2head/                 [README](02_head2head/README.md)
│   ├── framework_head2head.R     harness: tasks × datasets × frameworks
│   └── blocks/                   one file per task × framework (these are the timed,
│                                 executed, line-counted implementations)
├── 03_clinical_scenarios/        [README](03_clinical_scenarios/README.md)
│   ├── 01_undiagnosed_dysglycemia/  ┐
│   ├── 02_albuminuria/              ├ 01_download.R → 02_cohort.R → 03_model.R
│   ├── 03_undiagnosed_hypertension/ ┘ each folder ships its NHANES cache in data/
│   ├── run_scenarios.R           all three, in order
│   └── derived/                  calibration, NRI/IDI, threshold audit, DeLong, Fig. 2/3
└── outputs/                      everything generated (git-ignored)
```

## Quick start

```r
# from the repository root
source("run_all.R")                 # all three modules, full settings
```

```bash
# or from a shell, one module at a time
Rscript 01_leakage/run_leakage.R
Rscript 02_head2head/framework_head2head.R
Rscript 03_clinical_scenarios/run_scenarios.R
Rscript 03_clinical_scenarios/derived/run_derived.R
```

**Smoke test (minutes, not hours):**

```bash
QUICK=1 Rscript run_all.R
QUICK=1 H2H_DATASETS=Noise Rscript 02_head2head/framework_head2head.R   # single dataset
```

`QUICK=1` keeps the same code path and simply reduces repetitions
(20/20/30 → 5/5/5 for the leakage module, 20 → 3 splits for the head-to-head).
Figures produced in QUICK mode are narrower, because they describe fewer
repetitions — the manuscript numbers come from the default settings.

## Requirements

R ≥ 4.4 and:

```r
install.packages(c(
  "here",                                  # repository-relative paths
  "caret", "pROC", "ggplot2", "readxl", "dplyr",
  "tidymodels", "mlr3", "mlr3pipelines", "mlr3learners", "mlr3filters",
  "glmnet", "ranger", "rpart",
  "patchwork", "scales",                   # composite figures
  "nhanesA"                                # clinical scenarios: NHANES download
))
remotes::install_github("YuLab-SMU/icare")  # icare >= 1.1.0

install.packages("colino")   # optional: the only supervised filter for tidymodels
```

`icare` must be **1.1.0 or newer** (the group label moved into `info.data`; older
versions fail in `01_leakage/leakage_validation.R` part C). The script checks this
and stops with the upgrade command.

`colino` is optional. Without it, two cells are skipped and reported as skipped:
tidymodels in task T2 and in the scalability grid. T1 and T3 are unaffected.

Versions used for the manuscript (written to `outputs/*/sessionInfo.txt` by every
module): R 4.4.0, icare 1.1.0, caret 7.0.1, tidymodels 1.4.1, mlr3 1.6.0,
pROC 1.19.0.1.

## Data

| Input | Where it comes from |
|---|---|
| UTUC plasma proteome and routine laboratory panel (n = 601, PMID 37633276) | shipped inside the `icare` package (`inst/extdata`); a local copy in `data/` takes precedence |
| NHANES 2007–2008 and 2009–2010 tables | `03_clinical_scenarios/*/data/*.rds` (cached, ≈ 3 MB per scenario); delete a folder to force a fresh download with `nhanesA` |

The proteome workbook is used as published. Its preprocessing in these scripts is
label-free (proteins quantified in ≥ 80 % of samples, `log2(x + 1)`, undetected set
to 0), so it cannot carry outcome information; everything that *is* fitted
(scaling, screening, model) is fitted inside each framework, on the training
partition only.

## Outputs

All generated material is written under `outputs/` and git-ignored:

| Path | Content |
|---|---|
| `outputs/01_leakage/tables/exp{A,B,D}_*.csv` | per-repetition results of the leakage experiments |
| `outputs/01_leakage/figures/Figure_leakage_validation.{pdf,png}` | 2 × 2 composite figure |
| `outputs/01_leakage/tables/figure_numbers.csv` | every number printed in that figure |
| `outputs/02_head2head/tables/expE_*.csv` | long results, equivalence vs iCARE, code volume, scaling, selection agreement |
| `outputs/02_head2head/figures/FigE1–E5.*` | AUC, equivalence, runtime, code volume, scaling |
| `outputs/03_clinical_scenarios/tables/*.csv` | calibration, NRI/IDI, confusion matrices, DeLong gain, threshold audit |
| `outputs/03_clinical_scenarios/figures/*` | Fig. 2 / Fig. 3 composites, threshold positions |
| `outputs/<module>/sessionInfo.txt` | exact package versions of the run |

Per-scenario main results stay next to the scenario that produced them
(`03_clinical_scenarios/<scenario>/output/`), also git-ignored.

## How the numbers are pinned down

- Every experiment sets its own seed and re-draws the split from it, so a rerun
  reproduces the same splits from a clean session.
- ROC direction is always passed explicitly (`direction = "<"`). pROC's
  `direction = "auto"` flips any AUC below 0.5 upward, which would turn a true 0.5
  on null data into ≈ 0.55 and hide exactly the bias being measured.
- The equivalence margin (AUC ±0.01, TOST on paired splits) and the risk
  categories for the categorical NRI (2 / 10 / 50 / 95 %) are set in one place at
  the top of the script that uses them.
- External validation never refits anything: thresholds, calibration and
  preprocessing learned on 2007–2008 are applied to 2009–2010 as they are.
- Failures are recorded as rows with an `err` column rather than dropped
  (see `expE_runs_long.csv`).

## Caveats worth knowing before quoting a number

- **Runtime is not a like-for-like comparison in T1/T2.** iCARE's
  `ModelTrainAnalysis()` always cross-validates inside the training partition,
  while the caret / tidymodels / mlr3 blocks fit a single model. T3 is the only
  task where all four frameworks do the same amount of work.
- **iCARE is built on caret**, so identical numbers between those two are
  expected by construction. The informative comparisons are tidymodels and mlr3.
- **The T2 feature sets are not identical across frameworks** by design:
  tidymodels has no built-in supervised filter, so it uses
  `colino::step_select_roc` (univariate ROC ranking) while the other three use an
  ANOVA F-test. The overlap is reported (Jaccard) instead of being hidden.
- The leakage experiment on the real cohort has 23 features and a strong signal,
  so leakage inflates very little there — an honest null result. The noise
  experiment (p = 500, no signal) is where the effect is visible.

## Citation and license

If you use this code, please cite the manuscript (reference to be added on
publication). Code is released under the MIT license; the NHANES data are public
domain and the proteome data follow the terms of the original publication.
