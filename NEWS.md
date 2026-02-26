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


