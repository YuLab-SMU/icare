# 01 — Data-leakage validation

**Question.** Leakage prevention is presented as a core property of iCARE. Does the
pipeline actually keep the test set out of everything that gets fitted?

**Design.** Both arms run the same iCARE functions in the same order:

```
CreateModelObject -> FeatureSelectSBF -> ApplyFeatureSelection ->
FitApplyPreprocess -> ModelTrainAnalysis
```

and differ in **one** thing: the object handed to `FeatureSelectSBF()` /
`FitApplyPreprocess()`.

| Arm | Object passed in | Samples the step can see |
|---|---|---|
| correct | `CreateModelObject(data = train)` | training only |
| leaky | `CreateModelObject(data = all)` | training + the future test set |

Every other parameter — folds, seed, methods, `tuneLength` — is identical, so the
AUC difference is the inflation caused by leakage, and nothing else.

## Experiments

| Part | Data | What it isolates |
|---|---|---|
| **A** | pure noise, n = 120, p = 500, random labels (true AUC 0.5) | screening leak; a sample-level audit counts how many test rows each arm's object actually contains |
| **B** | PMID37633276 cohort, n = 601, p = 23 | 2 × 2: preprocessing fitted on train vs all × screening on train vs all |
| **C** | same cohort, Module-1 `Stat` object | the most realistic path in the tutorials: clean/impute/scale the whole cohort first, split afterwards (optional, needs `R/stat_replay.R`) |
| **D** | noise, 4 algorithms | optimism of picking the algorithm by **test-set** AUC instead of training CV |
| **E** | — | summary numbers and the 2 × 2 composite figure |

## Files

| File | Role |
|---|---|
| `leakage_validation.R` | the experiments; written to be read and run block by block in RStudio |
| `make_composite_figure.R` | draws the composite figure from the tables the first script wrote |
| `run_leakage.R` | batch entry point: runs both, prints the elapsed time |

## Run

```bash
Rscript 01_leakage/run_leakage.R           # ~20 min
QUICK=1 Rscript 01_leakage/run_leakage.R   # 5 repeats per part, ~2 min
```

Outputs: `outputs/01_leakage/{tables,figures}/` plus `sessionInfo.txt`.

## Two details that matter for the numbers

- **ROC direction is fixed to `"<"`.** iCARE's own `@all.results` stores AUC
  computed with pROC's `direction = "auto"`, which flips any value below 0.5
  upward. On null data that turns 0.5 into ≈ 0.56 and would make the *unbiased*
  selection look biased too. Part D therefore recomputes every test AUC with the
  direction fixed and keeps the iCARE values alongside, as a separate column.
- **The real cohort is an honest null result.** With 23 features and a strong
  signal, leakage barely moves the AUC (part B differences are within noise). The
  harm of leakage grows with many features, a weak signal and a small sample —
  which is exactly what part A shows.

Part A's audit columns (`seen_test_correct`, `seen_test_leaky`) are the direct
evidence: 0 test samples seen by screening in the correct arm, all of them in the
leaky arm, in every repetition.
