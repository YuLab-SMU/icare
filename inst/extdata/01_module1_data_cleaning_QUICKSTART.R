# =============================================================================
# icare Package Vignette — Module 1: Data Cleaning & Exploratory Analysis
# QUICK-START EDITION
# =============================================================================
# Audience : First-time users, or anyone who wants a working, minimal
#            analysis in a few minutes without touching every option.
# This script covers the core workflow only:
#   load data -> create a StatObject -> clean it -> look at one or two
#   summary plots -> run a group comparison -> get a baseline table.
#
# Looking for every option, publication-ready figures, or more diagnostics?
# See "01_module1_data_cleaning_ADVANCED.R" instead.
options(warn = -1)
# =============================================================================
rm(list = ls())
library(icare)
# ---------------------------------------------------------------------------
# 1. Load the data
# ---------------------------------------------------------------------------
library(readxl)
raw_dis=as.data.frame(read_excel("./PMID37633276_DIA_plasmaproomic.xlsx",sheet = 3))
raw_dis=raw_dis[!grepl(x = raw_dis$`Case ID`,pattern = "^KC"),]
dput(colnames(raw_dis))
mat=raw_dis[,c( "White blood cell count", 
               "Lymphocytes%", "Monocytes%", "Neutrophilic granulocyte%", "Eosinophil%", 
               "Basophil%", "Lymphocytes count", "Monocytes count", "Neutrophilic granulocyte count", 
               "Eosinophil count", "Basophil count", "Red blood cell count", 
               "Hemoglobin", "Hematocrit", "Mean red blood cell volume", "Mean cell hemoglobin", 
               "Mean corpuscular hemoglobin concentration", "Red blood cell distribution width CV", 
               "Red blood cell distribution width SD", "Platelet count", "Mean platelet volume", 
               "PCT", "Platelet volume distribution width")]
index=grepl(x = raw_dis$`Case ID`,pattern = "^HC")
inf=raw_dis[,c("Case ID", "Sample type", "Plasma proteome ID", "Tissue proteome ID", 
      "Gender", "Age at diagnosis", "Tobacco smoking history", "Drinking history", 
      "Vital status", "Days to last followup", "Days until death", 
      "Combined days to last followup or death", "Progressive-free survival_Status", 
      "Progressive-free survival_Months", "Tumor location", "Tumor size (cm)", 
      "T-category", "Lymph node involvement", "Distant metastasis", 
      "TNM Stage", "Grade", "Morphology")]
inf$group=ifelse(index,'0','1')
mat$group=inf$group
mat=as.matrix(mat)
storage.mode(mat) <- "numeric"

rownames(mat)=inf$`Case ID`
rownames(inf)=inf$`Case ID`
# ---------------------------------------------------------------------------
# 2. Create the Stat object
# ---------------------------------------------------------------------------
# This is the container icare uses to track your data through cleaning.
# group_col tells icare which column is the outcome you want to compare groups on.
stat_obj <- CreateStatObject(
  raw.data  = mat,info.data = inf,
  group_col = "group",
  na.action = "allow"
)

# ---------------------------------------------------------------------------
# 3. Clean the data (one line per step, sensible defaults)
# ---------------------------------------------------------------------------
stat_obj <- stat_diagnose_variable_type(stat_obj)          # detect numeric vs. categorical
stat_obj <- stat_convert_variables(stat_obj)               # apply those types
stat_obj <- stat_miss_processed(stat_obj, impute_method = "median_mode")  # fill missing values
stat_obj <- stat_detect_and_mark_outliers(stat_obj, method = "iqr")       # flag outliers
stat_obj <- stat_handle_outliers(stat_obj, method = "impute", impute_value = "median") # fix them
stat_obj <- stat_onehot_encode(stat_obj)                   # encode categorical variables
stat_obj <- stat_normalize_process(stat_obj, method = "auto")            # normalize numeric variables

cat("Clean data ready:", nrow(stat_obj@clean.data), "rows x",
    ncol(stat_obj@clean.data), "columns\n")

# ---------------------------------------------------------------------------
# 4. Two quick summary plots
# ---------------------------------------------------------------------------
# A) How do the groups differ across a few key numeric variables?
numeric_vars <- names(stat_obj@clean.data)[sapply(stat_obj@clean.data, is.numeric)]
PlotGroupedDistribution(
  object    = stat_obj,
  features  = head(numeric_vars, 6),
  group_col = "group",
  test      = "wilcox.test",
  save_plot = TRUE
)

# B) Do the two groups separate overall (PCA)?
PlotPCA(
  object   = stat_obj,
  color_by = "group",
  ellipse  = TRUE,
  save_plot = TRUE
)

# ---------------------------------------------------------------------------
# 5. Which features differ most between groups?
# ---------------------------------------------------------------------------
stat_obj <- stat_var_feature(stat_obj, p_threshold = 0.05)
deg_result <- ExtractLastTestSig(stat_obj)
# ---------------------------------------------------------------------------
# 6. Baseline characteristics table ("Table 1")
# ---------------------------------------------------------------------------
stat_obj <- stat_gaze_analysis(stat_obj, save_word = TRUE)
# ---------------------------------------------------------------------------
# 7. Save your work
# ---------------------------------------------------------------------------
# Module 2 (modeling) and Module 3 (subtyping) both can start by loading this file.
saveRDS(stat_obj, file = "stat_obj.rds")

cat("\nDone! stat_obj.rds saved — you're ready for Module 2 (modeling)",
    "or Module 3 (subtyping).\n")

