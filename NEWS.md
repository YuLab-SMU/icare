# Icare 1.1.0

### Breaking change: the group label lives in `info.data` only

`CreateStatObject()` no longer accepts the group/outcome column inside `raw.data`
or `clean.data`. It must be a column of `info.data`, and a label found in a
feature frame is now a hard error (no silent migration) whose message contains a
copy-pasteable fix:

```r
# before
inf$group <- ifelse(index, "0", "1")
mat$group <- inf$group
stat_obj  <- CreateStatObject(raw.data = mat, info.data = inf, group_col = "group")

# after
inf$group <- factor(ifelse(index, "Healthy Control", "Case"),
                    levels = c("Healthy Control", "Case"))
stat_obj  <- CreateStatObject(raw.data = mat, info.data = inf, group_col = "group")
```

**Why.** A label sitting among the features is indistinguishable from a feature
to any code that selects columns by type. A numeric `0`/`1` label was silently
consumed by `PlotPCA()` and `PlotCorrelationHeatmap()`, and because the label
dominated the first component the plotted group separation looked roughly five
times stronger than the features actually justified. The new layout makes that
class of bug unrepresentable rather than merely discouraged.

**Consequences.**

*   The label's R type no longer affects the analysis. `factor`/character labels
    are supported throughout Module 1, and readable labels become the Table 1
    column headers (`Case (N=30)` instead of `1 (N=30)`).
*   `clean.data`, `raw.data` and `scale.data` are guaranteed to contain features
    only; `Stat@group_col` resolves against `info.data`.
*   `info.data` is append-only. Operations that drop samples from `clean.data`
    leave it untouched, and labels are aligned by row name at read time.
*   New internal access layer (`.stat_labels()`, `.stat_features()`,
    `.stat_with_group()`, `.stat_adopt_legacy()`) replaces the per-call-site
    "exclude the label by name" convention, which had already been forgotten in
    three places.
*   `Train_Model` gains an `info.data` slot and `PrognosiX` gains a `group_col`
    slot, so round trips between classes no longer lose the label.
*   Objects saved by icare <= 1.0.1 still load: the label is lifted out of the
    feature frames on first read, with a message.

### Bug Fixes

*   **`stat_var_feature()` produced NaN logFC on data containing negative values.**
    The logFC strategy was inferred from the recorded normalisation method and
    hardcoded to `"log2ratio"` on the clean-data path, but `log2()` of a
    non-positive group mean is undefined. Every feature whose group mean was not
    positive came back with `logFC = NaN` and only a bare `NaNs produced`
    warning -- and because whether a mean is negative depends on the data, the
    same input could succeed or fail depending on the random seed. The guard in
    `batch_Wilcoxon()` that should have caught this was dead code on the main
    path, since `stat_var_feature()` always supplies `logfc_mat`. The strategy is
    now validated against the actual data content and falls back to `"diff"`
    (the mean difference, which is the correct effect size whenever negatives
    are present) with an actionable warning. The applied strategy and its data
    source are recorded on the result as `attr("logfc_type")` and
    `attr("logfc_data_src")`. `.compute_log2_means()` now names the offending
    feature/group instead of propagating NaN.
*   **`ModelTrainAnalysis()` never recorded `preprocess_columns`.** The field
    existed only on the `PreprocessModelData()` path, so `deploy_clinlab_app()`
    fell through to a hardcoded column-name list that could leave the outcome
    (and stray metadata) in the UI inputs and then in `predict(preproc, ...)`.
    Both training branches now record the fitted preprocess columns.
*   **`one_hot_encode()` could one-hot-encode a Train_Model outcome.** With the
    default `group_col = "group"` and an outcome named anything else, the label
    was treated as an ordinary categorical and encoded away, leaving an invalid
    object. The S4 branch now resolves the outcome from `@group_col`.
*   **`.extract_xy()` crashed with a base-R error on an outcome-less
    Train_Model.** `group_col = NULL` is the class prototype, so this was
    reachable; it now raises a clear error (the Stat branch keeps its legal
    `y = NULL`).
*   **`MergeObjects()` dropped metadata.** The Train_Model branch rebuilt via
    `CreateModelObject(data = ...)`, losing `info.data` entirely; the PrognosiX
    branch lost `group_col`. Both are preserved now.
*   **`RenameFeatures()` could orphan the outcome name.** Renaming the label
    column left `group_col` pointing at a nonexistent column; the rename now
    propagates to `group_col` and to info.data, and a collision with an
    existing feature is rejected.
*   **`InspectObject()`** now reports Train_Model's `info.data` and PrognosiX's
    `group_col`.
*   **`one_hot_encode()` on PrognosiX** filters `sub.data` along with
    clean.data / info.data / survival.data (it used to keep all rows).
*   **`CreateStatObject(na.action = "omit")` silently kept samples whose LABEL
    is NA.** Their rows survive by design (info.data is append-only), but a
    warning now says they will drop out of every group comparison.
*   **Roxygen examples still put the group column inside the feature matrix**,
    which the new contract rejects. The first pass caught 14 examples, but R
    CMD check on `main` still failed: 18 *executable* examples (the stat
    module) and 3 `\dontrun{}` examples (the prognosis module) had been missed,
    because the interactive-gated and non-run blocks are not exercised by the
    test suite. All 21 were rewritten to pass the label through info.data, and
    the `@param group_col` missing from `CreateSubtypingObject()` was
    documented, clearing the matching `Rd \usage` WARNING.
*   **`ConvertObject(Stat -> Subtyping)` left `group_col` pointing at a
    column that no longer exists.** Subtyping renames info.data with
    `janitor::make_clean_names()`, a DIFFERENT rule from the
    `modify_column_names()` used by `CreateStatObject()` -- "Diag Group"
    becomes "Diag_Group" in Stat but "diag_group" in Subtyping. The slot kept
    the Stat-side name, which matched nothing: labels were unreadable and
    `stratify_by` could not find the column. The name is now re-resolved
    against the actual info.data columns after conversion, with a message.
*   **`batch_Wilcoxon_MultiClass()` / `_OrderedMultiClass()` crashed when a
    comparison produced nothing significant.** A 0-row frame aborted on scalar
    column assignment ("replacement has 1 row, data has 0"), and an all-empty
    result list made `do.call(rbind, list())` return NULL. Both now skip the
    empty comparison and warn; an entirely empty result is returned as a
    well-formed 0-row frame. (Commit 0c0efdd guarded only `batch_Wilcoxon()`.)
*   **`ConvertObject(Subtyping -> PrognosiX)` lost `group_col`**, so the label
    was not excluded from the survival candidate-feature scan and could be
    tested as a Cox feature.
*   **`stat_convert_variables()` dropped `group_col`** when calling its worker, so
    the label was treated as an ordinary variable and coerced with
    `as.numeric()`. Character labels such as `"HC"`/`"Case"` became `NA` for
    every sample, with only a `NAs introduced by coercion` warning. `one_hot_encode`
    had the same class of omission.
*   **`stat_compute_descriptive()` dropped `group_col`**, silently producing
    ungrouped descriptives whenever the label was not literally named `"group"`;
    `compute_descriptive_stats()` also aborted when `group_col` was `NULL`
    (`NULL %in% colnames(d)` is `logical(0)`).
*   **`batch_Wilcoxon_MultiClass()` / `_OrderedMultiClass()` required the label in
    `logfc_mat`** and raised `subscript out of bounds` on a feature-only frame --
    which is now the normal case. The label is borrowed from `mat` instead.
*   **`one_hot_encode()` desynchronised three frames** on a `PrognosiX` object:
    rows with missing time/status were dropped from `survival.data` only, leaving
    `clean.data` and `info.data` at the original length.
*   **`CreatePrognosiXObject()` crashed on any `time_col` that needed renaming.**
    Column names were normalised and then looked up under the original spelling,
    so `time_col = "Time To Event"` raised `replacement has 0 rows`. Names are
    now tracked by position across the rename.
*   **Outcome leakage in `surv_filter_features_clinical()`.** `info.data` was
    merged wholesale and then scanned for numeric columns, so a `0`/`1`-coded
    group column was tested as a candidate prognostic marker. The label is now
    excluded from both the merge and the candidate set.
*   **`DownsampleObject(per_group = TRUE)` matched groups positionally.**
    `all_samples` came from `clean.data` while the grouping vector could come
    from `info.data`; with differing row orders the wrong samples were sampled,
    silently. Matching is now by row name, and a missing row is an error.
*   **`deploy_clinlab_app()` hardcoded the outcome column names**, so extra label
    columns became UI inputs and were then passed to
    `predict(preproc, newdata)`. The feature set is now taken from
    `process.info$preprocess_columns`.

### Tests

*   Added a `tests/testthat/` suite (18 cases) covering the label contract, label
    type freedom, and regressions for each bug above. The package previously had
    no tests directory and never called `validObject()`, so its validity
    assertions were unreachable.

# Icare 1.0.1

### Dependency Refactoring
*   **Streamlined Dependencies**: Moved heavy/optional packages from `Imports` to `Suggests`, including the `mlr3` ecosystem (`mlr3`, `mlr3proba`, `mlr3tuning`, `mlr3pipelines`, `mlr3tuningspaces`, `paradox`), `shiny`, `shinydashboard`, `DT`, `NMF`, `Rtsne`, `umap`, `mclust`, `NbClust`, `factoextra`, `mice`, `rms`, `timeROC`, `nricens`, `survminer`, `maxstat`, `survcomp`, `risksetROC`, `shapviz`, `survex`, `rBayesianOptimization`, `autoReg`, `officer`, `flextable`, `randomForestSRC`, `Hmisc`, and others. Only core data-handling and plotting packages remain hard dependencies.
*   **Lazy Dependency Checks**: Added internal `.require_pkgs()` helper that validates optional packages at function entry and raises actionable `install.packages()` hints; applied across the model, prognosis, subtyping, stat, and viz modules.
*   **Namespaced Calls**: Removed blanket `@import` / `@importFrom` directives for optional packages and switched to explicit `pkg::fun()` calls (e.g. `shiny::`, `shinydashboard::`, `mlr3::`, `mlr3proba::`, `mlr3pipelines::`, `mlr3tuning::`, `paradox::`) in deployment apps and mlr3-based functions.
*   **Roxygen & Docs**: Upgraded `roxygen2` to 8.1.0 and regenerated `NAMESPACE` and man pages; wrapped runnable examples in `\dontrun{}`; declared `utils::globalVariables(c("BIC", "time"))`; dropped `mdbrown/rmda` from `Remotes`.
*   **`DMwR` via GitHub Remotes**: `DMwR` is archived on CRAN, so it is installed from the official CRAN read-only mirror `cran/DMwR` (added to `Remotes`, kept in `Suggests`). This fixes the GitHub Actions dependency step, which failed because `setup-r-dependencies` could not resolve `DMwR` from CRAN. `R/model_05_train.R`'s SMOTE path now prints `remotes::install_github("cran/DMwR")` when the package is missing.
*   **Housekeeping**: Added `.reasonix` / `reasonix.toml` to both `.Rbuildignore` and `.gitignore`.

# Icare 1.0.0

### New Features & Improvements
*   **Version 1.0 Release**: Completed version 1.0 of the Icare R package — functions and tests refined and officially updated.
*   **Non-binary Downstream Analysis**: Added `R/model_05b_nonbinary_downstream.R` for multiclass and regression downstream tasks (multi-class ROC / confusion matrices, regression fit diagnostics, preprocessing helpers).
*   **Prognosis Extensions**: Added `R/prognosis_04_Extensions.R` covering time-dependent ROC comparison, repeated cross-validation evaluation, DeLong-style model comparison, and Bayesian hyperparameter tuning.
*   **Bundled Reference Data**: Added built-in `allmodel.rda` and `inst/extdata/all_model.csv` with model metadata.
*   **`sub_dir()` Utility**: New helper creating subdirectories under the global output root, with a session-temp fallback when no root is set.

### Bug Fixes
*   **Encoding Cleanup**: Replaced Unicode dashes and non-ASCII characters in code, docs, and messages with ASCII equivalents across deploy and prognosis modules.
*   **Ensemble Prediction**: Hardened error handling in ensemble prediction dispatch and non-binary probability conversion in `ModelDeployment`.
*   **Dependencies**: Added `maxstat` and `survcomp` to `Imports`, and `bbotk` / `mlr3mbo` to `Suggests`.

### Packaging & Housekeeping 
*   **GitHub Packaging**: Added `mlr-org/mlr3proba` to `Remotes` so a current version is installed from GitHub.
*   **Ignore Temporary Files**: Removed stray `.Rhistory` and added an ignore rule.
*   **Merge PR #7**: Resolved merge conflicts, keeping both sides' complementary changes.

# Icare 0.1.2

### Refactoring & Improvements
*   **Code Consolidation**: Merged `run_univariate_cox_analysis` from `module_3` and `module_4` into a single, enhanced implementation in `module_4_pn5_0_cox_univariate.R`.
    *   Added support for **covariate adjustment** and custom **P-value thresholds** (migrated from module 3).
    *   Retained **robust status column handling** and data validation (from module 4).
    *   Deprecated the redundant function in `module_3` to prevent namespace conflicts and installation warnings.

# Icare 0.1.1

### Bug Fixes
*   **Scope Issues**: Fixed variable scope errors in `survfit` and `survdiff` calls across multiple models (Lasso, PLS, RSF, CoxPH, SuperPC) by implementing formula injection.
*   **CoxBoost**: Resolved row mismatch errors in hazard ratio calculation caused by missing feature names.
*   **RSF**: Fixed execution errors in Random Survival Forest model when variable importance is empty; added robust handling for empty result data frames.
*   **SuperPC**: Fixed "logical subscript too long" error by ensuring proper numeric conversion of transposed feature matrices.
*   **Theme Conflicts**: Resolved conflicts between `ggprism` and `ggplot2` themes by standardizing on `theme_classic` and removing unstable dependencies.

### New Features & Improvements
*   **Dependencies**: Added `randomForestSRC`, `plsRcox`, `superpc`, and `rmda` to package dependencies to support advanced modeling and DCA.
*   **DCA Support**: Enhanced Decision Curve Analysis (DCA) integration, ensuring robust extraction of best model results.


# Icare 0.1.0

Icare is a comprehensive R package designed for survival analysis, clinical prediction modeling, and bioinformatics data analysis. It provides a unified framework for:

*   **Data Preprocessing**: Missing value handling, outlier detection, normalization, and descriptive statistics.
*   **Survival Modeling**: Implementation of various survival models including Lasso-Cox, CoxPH, Random Survival Forests (RSF), PLS-Cox, CoxBoost, and SuperPC.
*   **Model Evaluation**: Time-dependent ROC analysis, Kaplan-Meier survival curves, Decision Curve Analysis (DCA), and calibration assessment.
*   **Subtyping Analysis**: Molecular subtyping using K-means, NMF, and consensus clustering with t-SNE/UMAP visualization.
*   **Visualization**: High-quality, publication-ready plots for all analysis steps.


