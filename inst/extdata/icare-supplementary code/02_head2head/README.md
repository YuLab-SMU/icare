# 02 — Head-to-head: iCARE vs caret / tidymodels / mlr3

**Question.** Reviewer: "compare iCARE with caret / mlr3 / tidymodels pipelines in
performance and in the amount of code."

**Design.** Four frameworks, four datasets, five axes, the *same* training/test
split for every framework in every repetition.

| Task | Pipeline |
|---|---|
| T1 | z-score fitted on the training partition → logistic regression |
| T2 | z-score → univariate top-K screening → logistic regression, all fitted on the training partition (the leakage-prone case) |
| T3 | glmnet / ranger / rpart with identical fixed hyper-parameters, 5-fold CV on the training partition, best algorithm by CV AUC, evaluated on the held-out partition |
| T4 | scalability: runtime of the T2 pipeline over a grid of n and p |

| Dataset | n | p | Note |
|---|---|---|---|
| `UTUC_lab` | 601 | 23 | routine laboratory, strong signal |
| `UTUC_prot` | 601 | 730 | plasma proteome, p > n_train |
| `UTUC_prot_n120` | 120 | 730 | small-n / high-dimensional |
| `Noise` | 120 | 500 | random labels: true AUC 0.5, any AUC above it means leakage |

**Axes.** discrimination (AUC, direction fixed) and calibration (Brier) · paired
equivalence against iCARE (TOST, ±0.01 AUC, 90 % CI, same splits) · agreement of
selected features (Jaccard) and of the selected algorithm · runtime including the
scaling grid · code volume (executable lines, function calls, distinct functions,
packages needed).

## Why `blocks/` exists

Each implementation is a separate file, `blocks/<task>_<framework>.R`. The harness
*executes* those files and times them, so the lines counted are the lines that ran
and the seconds measured are the seconds those lines took. Nothing is quoted from
documentation and nothing is re-typed for the figure.

A block reads `d_tr`, `d_va`, `s` (seed) and `K` from its environment and leaves
behind `prob` (T1/T2) or `cv_auc` + `probs` (T3).

## Run

```bash
Rscript 02_head2head/framework_head2head.R                                   # ~30-45 min
QUICK=1 Rscript 02_head2head/framework_head2head.R                           # 3 splits, small grid
QUICK=1 H2H_DATASETS=Noise Rscript 02_head2head/framework_head2head.R        # smoke test, ~3 min
```

Outputs: `outputs/02_head2head/tables/expE_*.csv`, `figures/FigE1–E5.*`,
`expE_sessionInfo.txt`.

## Reading the results honestly

- **T1/T2 runtime is not a like-for-like comparison.** iCARE's
  `ModelTrainAnalysis()` always cross-validates inside the training partition; the
  caret / tidymodels / mlr3 blocks fit once. T3 is the only task where all four do
  the same amount of work, and it is the only one we quote for "relative cost".
- **iCARE is built on caret**, so identical numbers between those two are expected
  by construction. The informative comparisons are tidymodels and mlr3.
- **T2 feature sets differ by design.** tidymodels has no built-in supervised
  filter; the closest standard option is `colino::step_select_roc` (univariate ROC
  ranking), while caret / mlr3 / iCARE use an ANOVA F-test. The overlap is
  reported as Jaccard and the consequence is stated, not hidden.
- **The `Noise` dataset is the sanity check.** If any framework returns an AUC
  clearly above 0.5 there, something in that block leaks.
- Failures are rows with an `err` column in `expE_runs_long.csv`; they are never
  silently dropped.
