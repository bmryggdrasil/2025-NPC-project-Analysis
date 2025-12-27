# ==============================================================================
# Script Name: Survival Analysis of ANXA2 Expression in NPC Patients
# Method: Kaplan-Meier Analysis with Optimized Cut-off (maxstat)
# Journal: Nature-style Code Submission
# ==============================================================================

# 1. Load Dependencies
# Ensure packages are installed: install.packages(c("survival", "survminer", "readxl"))
library(survival)
library(survminer)
library(readxl)

# 2. Data Loading
# NOTE: Replace the path below with your local path for testing.
# For GitHub, users should place their data in the 'data/' directory.
data_path <- "data/GSE102234_final_data.xlsx" 
if (!file.exists(data_path)) {
  stop("Data file not found. Please ensure the data file is in the 'data/' folder.")
}

raw_data <- read_excel(data_path)

# 3. Data Pre-processing
# Extracting core variables: Survival Time, Event Status, and ANXA2 Expression
# Columns mapping: 'time to event' -> Time, 'event' -> Status, 'ANXA2' -> ANXA2
analysis_df <- raw_data[, c("time to event", "event", "ANXA2")]
colnames(analysis_df) <- c("Time", "Status", "ANXA2")

# Ensure numeric format
analysis_df$Time <- as.numeric(analysis_df$Time)
analysis_df$Status <- as.numeric(analysis_df$Status)
analysis_df$ANXA2 <- as.numeric(analysis_df$ANXA2)

# ==============================================================================
# 4. Optimal Cut-point Determination
# ==============================================================================
# Calculating the best cut-off value using maximally selected rank statistics
res_cut <- surv_cutpoint(
  analysis_df,
  time = "Time",
  event = "Status",
  variables = "ANXA2"
)

# Output summary of the cut-point calculation
print(summary(res_cut))

# Categorize ANXA2 expression into "high" and "low" groups based on optimal cut-point
surv_cat_df <- surv_categorize(res_cut)

# ==============================================================================
# 5. Survival Curve Fitting (Kaplan-Meier)
# ==============================================================================
fit <- survfit(Surv(Time, Status) ~ ANXA2, data = surv_cat_df)

# Generate the visualization
ggsurvplot(
  fit,
  data = surv_cat_df,
  pval = TRUE,             # Display Log-rank test p-value
  conf.int = TRUE,         # Show confidence intervals
  risk.table = TRUE,       # Add risk table at the bottom
  legend.labs = c("High", "Low"),
  palette = c("#E41A1C", "#377EB8"), # Professional red and blue palette
  ggtheme = theme_minimal(),
  title = "Overall Survival: ANXA2 High vs Low"
)

# ==============================================================================
# End of Script
# ==============================================================================