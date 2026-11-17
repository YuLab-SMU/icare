# 03 — Three NHANES clinical screening scenarios

Three scenarios, one identical pipeline each. Every scenario is developed on the
NHANES **2007–2008** cycle and validated on the **2009–2010** cycle, which is never
touched during development.

| Folder | Label | Clinical action it mimics |
|---|---|---|
| `01_undiagnosed_dysglycemia` | HbA1c ≥ 6.5 % in adults reporting no diabetes (DIQ010 = 2) | who should be sent for an HbA1c test |
| `02_albuminuria` | urine albumin-to-creatinine ratio ≥ 30 mg/g | who should be sent for a confirmatory urine test |
| `03_undiagnosed_hypertension` | mean SBP ≥ 140 or DBP ≥ 90 mmHg, or on treatment | who should have blood pressure measured properly |

In each one the label is a **gold-standard test result a clinician would otherwise
have to order**; the features are only routine, low-cost measurements (CBC,
biochemistry, anthropometry, one blood-pressure reading). The gold-standard
variables never enter the feature set.

## Per-scenario pipeline

```
01_download.R   download the NHANES tables this scenario needs with nhanesA
                (skipped when data/ is already cached; delete data/ to re-download)
02_cohort.R     merge, apply eligibility rules, write output/cohort.rds + cohort_flow.csv
03_model.R      train with iCARE, validate on the other cycle
                output/results.csv, cv_rankings.csv, clinical_utility.csv,
                       roc_external.png, model.rds
```

`03_model.R` departs from the official Quick Start in exactly two places, both
stated in its header:

1. the algorithm is chosen by **training-CV ROC**, not by the maximum test-set AUC
   (the latter is optimistically biased);
2. **no feature selection is done across the split** — if screening is needed it
   happens inside the training partition only.

The fitted model and its preprocessing are packaged with `ModelDeployment()` and
replayed on the external cycle, so nothing is refitted there: the Youden threshold,
the calibration and the preprocessing are the 2007–2008 ones.

## Run

```bash
Rscript 03_clinical_scenarios/run_scenarios.R                    # all three
Rscript 03_clinical_scenarios/run_scenarios.R 02_albuminuria     # one
```

Requires `nhanesA` only if `data/` has been removed.

## Derived analyses (`derived/`)

These read the saved `cohort.rds` / `model.rds` and add the secondary tables. They
never refit a model, so their numbers are consistent with the main results by
construction.

| Script | Produces |
|---|---|
| `01_calibration.R` | calibration-in-the-large and calibration slope per scenario |
| `02_nri_idi_external.R` | continuous / categorical NRI, IDI (bootstrap CI), external confusion matrices |
| `03_threshold_positions.R` | where the Youden threshold sits inside the predicted-risk distribution |
| `04_auc_gain_delong.R` | paired DeLong test for the AUC gain of the full panel over age + sex |
| `05_composite_figures.R` | the Fig. 2 / Fig. 3 composite panels |

```bash
Rscript 03_clinical_scenarios/derived/run_derived.R     # all five, in order
```

Outputs: `outputs/03_clinical_scenarios/{tables,figures,cache}/`.

## Design notes

- **Development / external cycles are disjoint**: the 2007–2008 cycle is split
  70/30 for training and internal validation; the whole 2009–2010 cycle is the
  external test set, used once.
- **Categorical NRI** uses clinical risk strata at 2 / 10 / 50 / 95 %, defined once
  in `derived/02_nri_idi_external.R`.
- **The Youden threshold looks small** (0.04–0.15) because the predicted risks are
  themselves small at 3–10 % prevalence: the threshold is a cut point on the
  ranking, not a probability of 0.5. `derived/03_threshold_positions.R` shows where
  it falls in each distribution, and the operating characteristics at that point are
  in `output/results.csv` (`external_sensitivity`, `external_specificity`,
  `external_ppv`, `external_npv`).
- **Rule-out framing**: at these prevalences the models are used to *exclude*, so
  the relevant threshold is the one that reaches 90 % / 95 % sensitivity, reported
  in `output/clinical_utility.csv`; the PPV stays low by arithmetic, not by model
  failure.
