#############################
## 0. Install and load packages
#############################
# packages <- c(
#   "readxl", "dplyr", "tibble", "purrr", "stringr",
#   "survival", "timeROC", "openxlsx",
#   "survminer", "ggplot2"
# )
# 
# to_install <- packages[!packages %in% installed.packages()[, "Package"]]
# if (length(to_install) > 0) install.packages(to_install)
library(readxl)
library(dplyr)
library(tibble)
library(purrr)
library(stringr)
library(survival)
library(timeROC)
library(openxlsx)
library(survminer)
library(ggplot2)
#############################
## 1. Parameter settings
#############################
# -------- File paths --------
#abundance_file <- "~/Desktop/20260715.total_bracken/MAG_COUNT_ge30_filtered_expression/bracken_relative_Tongue_coating_filtered.MAG_COUNT_ge5.filtered.csv"
abundance_file <- "~/Desktop/20260715.total_bracken/MAG_COUNT_ge30_filtered_expression/bracken_relative_Stool_filtered.MAG_COUNT_ge5.filtered.csv"
#survival_file  <- "~/Desktop/Oral_OS_0729.xlsx"
survival_file  <- "~/Desktop/Stool_OS_0729/Stool_OS_0729.xlsx"
# -------- Sheet name --------
# abundance_file is csv now, no need for abundance_sheet
survival_sheet  <- 1
# -------- Column name settings --------
sample_col_abund <- "SampleID"
sample_col_surv  <- "SampleID"
time_col   <- "OS_time"
status_col <- "OS_status"
# -------- Covariate column names --------
gender_col <- "Gender"
age_col    <- "Age"
# -------- Filter thresholds --------
min_prevalence <- 0.10
min_nonzero_n  <- 10
# -------- Significance threshold --------
fdr_cutoff <- 0.05
# -------- Time points for timeROC evaluation (month) --------
roc_times <- c(12, 24, 36)
setwd("~/Desktop")
# -------- Output files --------
output_rank_file    <- "species_survival_ranking_adjusted_Gender_Age.xlsx"
output_score_file   <- "presence_based_highrisk_score_patients_adjusted_Gender_Age.xlsx"
output_risk_file    <- "presence_based_risk_species_adjusted_Gender_Age.xlsx"
output_protect_file <- "presence_based_protective_species_adjusted_Gender_Age.xlsx"
output_summary_file <- "presence_based_score_summary_adjusted_Gender_Age.xlsx"
output_km_pdf       <- "KM_curve_presence_based_score_adjusted_Gender_Age.pdf"
#############################
## 2. Read data
#############################
# Read stool relative abundance csv
abund_raw0 <- data.table::fread(
  abundance_file,
  sep = ",",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)
surv_raw <- read_excel(
  survival_file,
  sheet = survival_sheet
)
surv_raw <- as.data.frame(surv_raw)
cat("Raw dimension of abundance table:\n")
print(dim(abund_raw0))
cat("First 10 column names of abundance table:\n")
print(colnames(abund_raw0)[1:min(10, ncol(abund_raw0))])
#############################
## 2.1 Clean abundance table format
## Target format:
## Row = sample
## First column = SampleID
## Subsequent columns = species
#############################
# ID cleaning function to avoid unmatched IDs caused by ., -, space, W prefix
clean_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\u00A0", "", x)
  x <- gsub("[[:space:]]+", "", x)
  x <- sub("\\.0$", "", x)
  x <- gsub("\\.", "-", x)
  x <- sub("^W(?=[0-9])", "", x, perl = TRUE)
  x
}
# Check whether abundance_raw0 already contains SampleID column
if (sample_col_abund %in% colnames(abund_raw0)) {
  
  cat("Detected abundance table format: rows = samples, columns = species.\n")
  
  abund_raw <- abund_raw0
  
  colnames(abund_raw)[colnames(abund_raw) == sample_col_abund] <- "SampleID"
  abund_raw$SampleID <- clean_id(abund_raw$SampleID)
  
} else {
  
  cat("SampleID column not detected; default abundance format: 1st column = species, rest columns = samples. Start transposition.\n")
  
  species_col <- colnames(abund_raw0)[1]
  
  # Use first column as species names
  species_names <- as.character(abund_raw0[[species_col]])
  
  # Remove rows with empty species names
  keep_row <- !is.na(species_names) & species_names != ""
  abund_raw0 <- abund_raw0[keep_row, , drop = FALSE]
  species_names <- species_names[keep_row]
  
  # Deduplicate species names to prevent column name conflict after transposition
  species_names <- make.unique(species_names, sep = "__dup")
  
  # Extract abundance matrix: row = species, column = sample
  abund_mat <- abund_raw0[, -1, drop = FALSE]
  
  # Clean sample column names
  sample_ids <- clean_id(colnames(abund_mat))
  
  # Convert to numeric matrix
  abund_mat <- as.matrix(abund_mat)
  mode(abund_mat) <- "numeric"
  
  # Transpose to: row = sample, column = species
  abund_t <- as.data.frame(t(abund_mat), check.names = FALSE)
  colnames(abund_t) <- species_names
  
  abund_t$SampleID <- sample_ids
  
  # Move SampleID to first column
  abund_raw <- abund_t[, c("SampleID", setdiff(colnames(abund_t), "SampleID")), drop = FALSE]
}
cat("Dimension of cleaned abundance table:\n")
print(dim(abund_raw))
cat("First 10 column names of cleaned abundance table:\n")
print(colnames(abund_raw)[1:min(10, ncol(abund_raw))])
# Convert survival table to data.frame
surv_raw <- as.data.frame(surv_raw)
#############################
## 2.2 Column name validation
#############################
if (!"SampleID" %in% colnames(abund_raw)) {
  stop("SampleID column still missing after cleaning abundance table; please check bracken_relative_Stool_filtered.csv format.")
}
required_surv_cols <- c(sample_col_surv, time_col, status_col, gender_col, age_col)
if (!all(required_surv_cols %in% colnames(surv_raw))) {
  missing_cols <- setdiff(required_surv_cols, colnames(surv_raw))
  stop(paste0("Survival table missing columns: ", paste(missing_cols, collapse = ", ")))
}
#############################
## 3. Organize survival data, covariates and merge datasets
#############################
surv_df <- surv_raw %>%
  dplyr::select(
    SampleID = all_of(sample_col_surv),
    time = all_of(time_col),
    status = all_of(status_col),
    Gender = all_of(gender_col),
    Age = all_of(age_col)
  ) %>%
  mutate(
    SampleID = as.character(SampleID),
    time = as.numeric(time),
    status = as.numeric(status),
    Gender = as.character(Gender),
    Age = as.numeric(Age)
  ) %>%
  filter(!is.na(SampleID), !is.na(time), !is.na(status)) %>%
  filter(time > 0, status %in% c(0, 1))
# Treat Unknown / blank as NA
na_like_values <- c(
  "Unknown", "unknown", "UNKNOWN",
  "UnKnown", "unk", "UNK",
  "", " ", "NA", "N/A", "na", "n/a",
  "NULL", "null", "None", "none"
)
surv_df$Gender[surv_df$Gender %in% na_like_values] <- NA
# Impute Age with median
age_median <- median(surv_df$Age, na.rm = TRUE)
if (is.na(age_median)) {
  stop("All Age values are NA; median imputation unavailable.")
}
surv_df$Age[is.na(surv_df$Age)] <- age_median
# Impute Gender with most frequent category
gender_tab <- table(surv_df$Gender, useNA = "no")
if (length(gender_tab) == 0) {
  stop("All Gender values are NA; mode imputation unavailable.")
}
gender_mode <- names(gender_tab)[which.max(gender_tab)]
surv_df$Gender[is.na(surv_df$Gender)] <- gender_mode
surv_df$Gender <- factor(surv_df$Gender)
# Clean abundance table
abund_df <- abund_raw
colnames(abund_df)[colnames(abund_df) == sample_col_abund] <- "SampleID"
abund_df$SampleID <- as.character(abund_df$SampleID)
# Remove duplicate samples
surv_df  <- surv_df %>% distinct(SampleID, .keep_all = TRUE)
abund_df <- abund_df %>% distinct(SampleID, .keep_all = TRUE)
# Merge tables
dat0 <- inner_join(surv_df, abund_df, by = "SampleID")
if (nrow(dat0) < 20) {
  warning("Merged sample size <20; statistical stability may be compromised.")
}
# Identify species columns: exclude SampleID, time, status, Gender, Age
species_cols <- setdiff(
  colnames(dat0),
  c("SampleID", "time", "status", "Gender", "Age")
)
if (length(species_cols) == 0) {
  stop("No species columns identified; check abundance table format.")
}
# Convert species abundance to numeric
dat0[species_cols] <- lapply(dat0[species_cols], function(x) {
  as.numeric(as.character(x))
})
cat("===== Covariate distribution =====\n")
cat("Gender:\n")
print(table(dat0$Gender, useNA = "ifany"))
cat("Age summary:\n")
print(summary(dat0$Age))
#############################
## 4. Pre-filter species
#############################
species_stats <- lapply(species_cols, function(sp) {
  x <- dat0[[sp]]
  data.frame(
    species = sp,
    non_na_n = sum(!is.na(x)),
    nonzero_n = sum(x > 0, na.rm = TRUE),
    prevalence = mean(x > 0, na.rm = TRUE),
    sd_raw = sd(x, na.rm = TRUE)
  )
}) %>% bind_rows()
keep_species <- species_stats %>%
  filter(
    nonzero_n >= min_nonzero_n,
    prevalence >= min_prevalence
  ) %>%
  pull(species)
if (length(keep_species) == 0) {
  stop("No species retained after filtering; relax min_prevalence / min_nonzero_n threshold.")
}
dat <- dat0 %>%
  dplyr::select(
    SampleID, time, status, Gender, Age,
    all_of(keep_species)
  )
#############################
## 5. Data transformation: log10(x + 1e-10) + Z-score normalization
#############################
## HR represents risk change per 1 SD increment of log10-transformed abundance
#############################
transform_one_species <- function(x) {
  x <- as.numeric(x)
  
  # Set negative values to NA
  x[x < 0] <- NA_real_
  
  # log10 transformation with tiny pseudo-value to avoid log10(0)
  x2 <- log10(x + 1e-10)
  
  # Normalize to per-1-SD increment
  sdx <- sd(x2, na.rm = TRUE)
  if (is.na(sdx) || sdx == 0) return(rep(NA_real_, length(x2)))
  
  as.numeric(scale(x2))
}
for (sp in keep_species) {
  dat[[sp]] <- transform_one_species(dat[[sp]])
}
#############################
## 6. Function for single-species adjusted Cox survival analysis
## Adjusted for Gender + Age
#############################
run_adjusted_survival <- function(df, sp, roc_times = c(12, 24, 36)) {
  
  tmp <- df %>%
    dplyr::select(time, status, Gender, Age, all_of(sp)) %>%
    rename(marker = all_of(sp)) %>%
    filter(
      !is.na(time),
      !is.na(status),
      !is.na(marker),
      !is.na(Gender),
      !is.na(Age)
    )
  
  # Skip if insufficient samples
  if (nrow(tmp) < 20) {
    return(data.frame(
      species = sp,
      n = nrow(tmp),
      events = sum(tmp$status == 1),
      beta = NA_real_,
      HR = NA_real_,
      HR_low95 = NA_real_,
      HR_high95 = NA_real_,
      z = NA_real_,
      p = NA_real_,
      c_index = NA_real_,
      auc_12 = NA_real_,
      auc_24 = NA_real_,
      auc_36 = NA_real_
    ))
  }
  
  # Skip if event count too low for stable estimation
  if (sum(tmp$status == 1) < 5) {
    return(data.frame(
      species = sp,
      n = nrow(tmp),
      events = sum(tmp$status == 1),
      beta = NA_real_,
      HR = NA_real_,
      HR_low95 = NA_real_,
      HR_high95 = NA_real_,
      z = NA_real_,
      p = NA_real_,
      c_index = NA_real_,
      auc_12 = NA_real_,
      auc_24 = NA_real_,
      auc_36 = NA_real_
    ))
  }
  
  # Adjust only for Age if Gender has single level in subset
  if (length(unique(tmp$Gender)) >= 2) {
    cox_formula <- Surv(time, status) ~ marker + Gender + Age
  } else {
    cox_formula <- Surv(time, status) ~ marker + Age
  }
  
  # Fit Cox model
  fit <- tryCatch(
    coxph(cox_formula, data = tmp),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(data.frame(
      species = sp,
      n = nrow(tmp),
      events = sum(tmp$status == 1),
      beta = NA_real_,
      HR = NA_real_,
      HR_low95 = NA_real_,
      HR_high95 = NA_real_,
      z = NA_real_,
      p = NA_real_,
      c_index = NA_real_,
      auc_12 = NA_real_,
      auc_24 = NA_real_,
      auc_36 = NA_real_
    ))
  }
  
  s <- summary(fit)
  
  if (!"marker" %in% rownames(s$coefficients)) {
    return(data.frame(
      species = sp,
      n = nrow(tmp),
      events = sum(tmp$status == 1),
      beta = NA_real_,
      HR = NA_real_,
      HR_low95 = NA_real_,
      HR_high95 = NA_real_,
      z = NA_real_,
      p = NA_real_,
      c_index = NA_real_,
      auc_12 = NA_real_,
      auc_24 = NA_real_,
      auc_36 = NA_real_
    ))
  }
  
  beta <- unname(s$coefficients["marker", "coef"])
  HR   <- unname(s$coefficients["marker", "exp(coef)"])
  zval <- unname(s$coefficients["marker", "z"])
  pval <- unname(s$coefficients["marker", "Pr(>|z|)"])
  
  ci <- tryCatch(
    suppressWarnings(confint(fit)),
    error = function(e) NULL
  )
  
  if (!is.null(ci) && "marker" %in% rownames(ci)) {
    HR_l <- exp(ci["marker", 1])
    HR_h <- exp(ci["marker", 2])
  } else {
    HR_l <- NA_real_
    HR_h <- NA_real_
  }
  
  # Linear predictor: full model risk score including marker + covariates
  lp <- tryCatch(
    predict(fit, type = "lp"),
    error = function(e) rep(NA_real_, nrow(tmp))
  )
  
  # C-index
  c_index <- tryCatch({
    survival::concordance(Surv(time, status) ~ lp, data = tmp, reverse = TRUE)$concordance
  }, error = function(e) NA_real_)
  
  # Time-dependent AUC
  auc_vec <- setNames(rep(NA_real_, length(roc_times)), paste0("auc_", roc_times))
  
  tr <- tryCatch({
    timeROC(
      T = tmp$time,
      delta = tmp$status,
      marker = lp,
      cause = 1,
      weighting = "marginal",
      times = roc_times,
      iid = FALSE
    )
  }, error = function(e) NULL)
  
  if (!is.null(tr)) {
    auc_vals <- tr$AUC
    names(auc_vals) <- paste0("auc_", roc_times)
    auc_vec[names(auc_vals)] <- auc_vals
  }
  
  out <- data.frame(
    species = sp,
    n = nrow(tmp),
    events = sum(tmp$status == 1),
    beta = beta,
    HR = HR,
    HR_low95 = HR_l,
    HR_high95 = HR_h,
    z = zval,
    p = pval,
    c_index = c_index,
    auc_12 = ifelse("auc_12" %in% names(auc_vec), auc_vec["auc_12"], NA_real_),
    auc_24 = ifelse("auc_24" %in% names(auc_vec), auc_vec["auc_24"], NA_real_),
    auc_36 = ifelse("auc_36" %in% names(auc_vec), auc_vec["auc_36"], NA_real_)
  )
  
  return(out)
}
#############################
## 7. Iterate Cox analysis over all filtered species
#############################
res <- purrr::map_dfr(
  keep_species,
  ~run_adjusted_survival(dat, .x, roc_times = roc_times)
)
#############################
## 8. Multiple testing correction, annotation and ranking
#############################
res2 <- res %>%
  mutate(
    FDR = p.adjust(p, method = "BH"),
    direction = case_when(
      is.na(HR) ~ NA_character_,
      HR > 1 ~ "Risk",
      HR < 1 ~ "Protective",
      TRUE ~ "Neutral"
    ),
    abs_z = abs(z),
    rank_score = -log10(pmax(FDR, 1e-300)) * abs(beta)
  )
# Primary ranking: FDR + |z| + C-index
rank_main <- res2 %>%
  arrange(FDR, desc(abs_z), desc(c_index))
# Risk-associated species
rank_risk <- res2 %>%
  filter(!is.na(HR), HR > 1) %>%
  arrange(FDR, desc(abs_z), desc(c_index))
# Protective species
rank_protective <- res2 %>%
  filter(!is.na(HR), HR < 1) %>%
  arrange(FDR, desc(abs_z), desc(c_index))
# Statistically significant species
sig_res <- rank_main %>%
  filter(!is.na(FDR), FDR < fdr_cutoff)
#############################
## 9. Export species Cox ranking results
#############################
wb <- createWorkbook()
addWorksheet(wb, "All_results_main_rank")
writeData(wb, "All_results_main_rank", rank_main)
addWorksheet(wb, "Risk_species")
writeData(wb, "Risk_species", rank_risk)
addWorksheet(wb, "Protective_species")
writeData(wb, "Protective_species", rank_protective)
addWorksheet(wb, "Significant_FDR")
writeData(wb, "Significant_FDR", sig_res)
addWorksheet(wb, "Filter_stats")
writeData(wb, "Filter_stats", species_stats)
saveWorkbook(wb, output_rank_file, overwrite = TRUE)
#############################
## 10. Calculate presence-based score using significant risk / protective species
## Higher score = higher mortality risk
##
## Score = (Proportion of risk species present - Proportion of protective species present + 1) / 2
## Presence defined as: relative abundance > 0
## Denominator: total count of species with FDR < fdr_cutoff in each category
#############################
sig_risk_species <- rank_main %>%
  filter(!is.na(FDR), FDR < fdr_cutoff, HR > 1) %>%
  pull(species)
sig_protect_species <- rank_main %>%
  filter(!is.na(FDR), FDR < fdr_cutoff, HR < 1) %>%
  pull(species)
cat("===== Number of significant risk species =====\n")
print(length(sig_risk_species))
cat("===== Number of significant protective species =====\n")
print(length(sig_protect_species))
if (length(sig_risk_species) == 0 & length(sig_protect_species) == 0) {
  stop("No eligible risk or protective species with FDR < 0.05; relax significance threshold.")
}
# Calculate score from raw abundance dat0 (not log10+Z-score transformed data)
score_df <- dat0 %>%
  dplyr::select(
    SampleID, time, status, Gender, Age,
    all_of(keep_species)
  )
# Safety filter: retain only valid significant species columns
sig_risk_species <- intersect(sig_risk_species, colnames(score_df))
sig_protect_species <- intersect(sig_protect_species, colnames(score_df))
# Binarization: presence = abundance > 0
presence_df <- score_df
presence_species <- setdiff(
  colnames(presence_df),
  c("SampleID", "time", "status", "Gender", "Age")
)
presence_df[presence_species] <- lapply(
  presence_df[presence_species],
  function(x) as.numeric(as.numeric(x) > 0)
)
# Proportion of risk species present
if (length(sig_risk_species) > 0) {
  risk_present_n <- rowSums(
    presence_df[, sig_risk_species, drop = FALSE],
    na.rm = TRUE
  )
  risk_present_prop <- risk_present_n / length(sig_risk_species)
} else {
  risk_present_n <- rep(0, nrow(presence_df))
  risk_present_prop <- rep(0, nrow(presence_df))
}
# Proportion of protective species present
if (length(sig_protect_species) > 0) {
  protect_present_n <- rowSums(
    presence_df[, sig_protect_species, drop = FALSE],
    na.rm = TRUE
  )
  protect_present_prop <- protect_present_n / length(sig_protect_species)
} else {
  protect_present_n <- rep(0, nrow(presence_df))
  protect_present_prop <- rep(0, nrow(presence_df))
}
# Core risk score: higher value = higher risk
score_df$Score <- (risk_present_prop - protect_present_prop + 1) / 2
# Save auxiliary variables for interpretation
score_df$risk_present_n <- risk_present_n
score_df$risk_total_n <- length(sig_risk_species)
score_df$risk_present_prop <- risk_present_prop
score_df$protect_present_n <- protect_present_n
score_df$protect_total_n <- length(sig_protect_species)
score_df$protect_present_prop <- protect_present_prop
#############################
## 11. Stratify patients by Score median
#############################
score_cutoff <- median(score_df$Score, na.rm = TRUE)
score_df$Group <- ifelse(score_df$Score >= score_cutoff, "High", "Low")
score_df$Group <- factor(score_df$Group, levels = c("Low", "High"))
cat("===== Score stratification count =====\n")
print(table(score_df$Group))
#############################
## 12. Survival analysis for Score subgroups
## KM curve unadjusted; Cox HR adjusted for Gender + Age
#############################
fit_km <- survfit(Surv(time, status) ~ Group, data = score_df)
# Unadjusted Cox model
fit_cox_unadjusted <- coxph(Surv(time, status) ~ Group, data = score_df)
cox_sum_unadjusted <- summary(fit_cox_unadjusted)
HR_unadj <- cox_sum_unadjusted$coefficients[1, "exp(coef)"]
HR_low_unadj <- cox_sum_unadjusted$conf.int[1, "lower .95"]
HR_high_unadj <- cox_sum_unadjusted$conf.int[1, "upper .95"]
p_cox_unadj <- cox_sum_unadjusted$coefficients[1, "Pr(>|z|)"]
# Adjusted Cox model: Group + Gender + Age
if (length(unique(score_df$Gender)) >= 2) {
  fit_cox_adjusted <- coxph(Surv(time, status) ~ Group + Gender + Age, data = score_df)
} else {
  fit_cox_adjusted <- coxph(Surv(time, status) ~ Group + Age, data = score_df)
}
cox_sum_adjusted <- summary(fit_cox_adjusted)
HR_adj <- cox_sum_adjusted$coefficients["GroupHigh", "exp(coef)"]
HR_low_adj <- cox_sum_adjusted$conf.int["GroupHigh", "lower .95"]
HR_high_adj <- cox_sum_adjusted$conf.int["GroupHigh", "upper .95"]
p_cox_adj <- cox_sum_adjusted$coefficients["GroupHigh", "Pr(>|z|)"]
# Log-rank test
logrank_test <- survdiff(Surv(time, status) ~ Group, data = score_df)
p_logrank <- 1 - pchisq(logrank_test$chisq, df = 1)
cat("===== Cox results for Score subgroups: unadjusted =====\n")
cat(sprintf("HR = %.4f\n", HR_unadj))
cat(sprintf("95%% CI = %.4f ~ %.4f\n", HR_low_unadj, HR_high_unadj))
cat(sprintf("Cox p = %.4g\n", p_cox_unadj))
cat(sprintf("Log-rank p = %.4g\n", p_logrank))
cat("\n===== Cox results for Score subgroups: adjusted for Gender + Age =====\n")
cat(sprintf("Adjusted HR = %.4f\n", HR_adj))
cat(sprintf("95%% CI = %.4f ~ %.4f\n", HR_low_adj, HR_high_adj))
cat(sprintf("Adjusted Cox p = %.4g\n", p_cox_adj))
#############################
## 13. Plot KM survival curve
#############################
label_txt <- sprintf(
  "Adjusted HR = %.2f (95%% CI %.2f-%.2f)\nAdjusted Cox p = %.2e\nLog-rank p = %.2e",
  HR_adj, HR_low_adj, HR_high_adj, p_cox_adj, p_logrank
)
p_score <- ggsurvplot(
  fit_km,
  data = score_df,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.col = "strata",
  surv.median.line = "hv",
  xlab = "Time (months)",
  ylab = "Overall survival probability",
  title = "KM Curve by Presence-based Score",
  legend.title = "Group",
  legend.labs = c("Low", "High"),
  palette = c("#1B9E77", "#D95F02"),
  xlim = c(0, 40),
  break.time.by = 6
)
p_score$plot <- p_score$plot +
  scale_x_continuous(
    limits = c(0, 40),
    breaks = c(0, 6, 12, 18, 24, 30, 36),
    expand = c(0, 0)
  ) +
  annotate(
    "text",
    x = 5,
    y = 0.20,
    label = label_txt,
    hjust = 0,
    size = 4
  ) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
p_score$table <- p_score$table +
  scale_x_continuous(
    limits = c(0, 40),
    breaks = c(0, 6, 12, 18, 24, 30, 36),
    expand = c(0, 0)
  )
print(p_score)
# Save KM curve PDF
pdf(output_km_pdf, width = 7, height = 7)
print(p_score)
dev.off()
#############################
## 14. Export score table and species lists
#############################
write.xlsx(score_df, output_score_file, overwrite = TRUE)
write.xlsx(data.frame(species = sig_risk_species), output_risk_file, overwrite = TRUE)
write.xlsx(data.frame(species = sig_protect_species), output_protect_file, overwrite = TRUE)
wb_score <- createWorkbook()
addWorksheet(wb_score, "patient_score")
writeData(wb_score, "patient_score", score_df)
addWorksheet(wb_score, "risk_species")
writeData(wb_score, "risk_species", data.frame(species = sig_risk_species))
addWorksheet(wb_score, "protective_species")
writeData(wb_score, "protective_species", data.frame(species = sig_protect_species))
addWorksheet(wb_score, "score_summary")
writeData(
  wb_score,
  "score_summary",
  data.frame(
    fdr_cutoff = fdr_cutoff,
    n_risk_species = length(sig_risk_species),
    n_protective_species = length(sig_protect_species),
    score_cutoff_median = score_cutoff,
    
    HR_high_vs_low_unadjusted = HR_unadj,
    HR_low95_unadjusted = HR_low_unadj,
    HR_high95_unadjusted = HR_high_unadj,
    Cox_p_unadjusted = p_cox_unadj,
    
    HR_high_vs_low_adjusted_Gender_Age = HR_adj,
    HR_low95_adjusted_Gender_Age = HR_low_adj,
    HR_high95_adjusted_Gender_Age = HR_high_adj,
    Cox_p_adjusted_Gender_Age = p_cox_adj,
    
    Logrank_p = p_logrank,
    output_rank_file = output_rank_file,
    output_score_file = output_score_file,
    output_km_pdf = output_km_pdf
  )
)
saveWorkbook(wb_score, output_summary_file, overwrite = TRUE)
#############################
## 15. Console summary output
#############################
cat("\nAnalysis completed.\n")
cat("Species abundance transformation: log10(x + 1e-10) + Z-score normalization\n")
cat("Single-species Cox model: Surv(time, status) ~ marker + Gender + Age\n")
cat("Score formula: (risk_present_prop - protect_present_prop + 1) / 2\n")
cat("Score subgroup Cox model: Surv(time, status) ~ Group + Gender + Age\n")
cat("Total raw species count:", length(species_cols), "\n")
cat("Species retained after filtering:", length(keep_species), "\n")
cat("Significant species count (FDR <", fdr_cutoff, "):", nrow(sig_res), "\n")
cat("Significant risk species count:", length(sig_risk_species), "\n")
cat("Significant protective species count:", length(sig_protect_species), "\n")
cat("Score median cutoff:", score_cutoff, "\n")
cat("Species Cox ranking output:", output_rank_file, "\n")
cat("Patient score table output:", output_score_file, "\n")
cat("Score summary output:", output_summary_file, "\n")
cat("KM curve PDF output:", output_km_pdf, "\n\n")
cat("Top 20 species (primary ranking adjusted for Gender + Age):\n")
print(
  rank_main %>%
    dplyr::select(
      species, n, events,
      HR, HR_low95, HR_high95,
      z, p, FDR,
      c_index, auc_12, auc_24, auc_36,
      direction
    ) %>%
    head(20)
)
#############################
## 10. Top10 risk + Top10 protective microbial score
#############################
top_n <- 10
#############################
# 10.1 Extract Top10 risk species
#############################
top_risk_species <- rank_main %>%
  filter(
    !is.na(FDR),
    FDR < fdr_cutoff,
    HR > 1
  ) %>%
  arrange(FDR, desc(abs_z), desc(c_index)) %>%
  head(top_n) %>%
  pull(species)
#############################
# 10.2 Extract Top10 protective species
#############################
top_protect_species <- rank_main %>%
  filter(
    !is.na(FDR),
    FDR < fdr_cutoff,
    HR < 1
  ) %>%
  arrange(FDR, desc(abs_z), desc(c_index)) %>%
  head(top_n) %>%
  pull(species)
cat("Top risk species list:\n")
print(top_risk_species)
cat("\nTop protective species list:\n")
print(top_protect_species)
#############################
# 10.3 Export combined 20 species table
#############################
top20_species <- data.frame(
  
  species=c(
    top_risk_species,
    top_protect_species
  ),
  
  type=c(
    rep("Risk",length(top_risk_species)),
    rep("Protective",length(top_protect_species))
  )
)
write.xlsx(
  top20_species,
  "Top10_Risk_Top10_Protective_species.xlsx",
  overwrite=TRUE
)
#############################
# 11. Calculate presence-based Score using top 20 species
#############################
score_df <- dat0 %>%
  select(
    SampleID,
    time,
    status,
    Gender,
    Age,
    all_of(
      intersect(
        top20_species$species,
        colnames(dat0)
      )
    )
  )
risk_species <- intersect(
  top_risk_species,
  colnames(score_df)
)
protect_species <- intersect(
  top_protect_species,
  colnames(score_df)
)
# Binarize presence/absence
score_binary <- score_df
score_cols <- setdiff(
  colnames(score_binary),
  c(
    "SampleID",
    "time",
    "status",
    "Gender",
    "Age"
  )
)
score_binary[score_cols] <- lapply(
  score_binary[score_cols],
  function(x){
    as.numeric(x>0)
  }
)
#############################
# Risk species presence proportion
#############################
risk_prop <-
  rowSums(
    score_binary[,risk_species,drop=FALSE],
    na.rm=TRUE
  ) /
  length(risk_species)
#############################
# Protective species presence proportion
#############################
protect_prop <-
  rowSums(
    score_binary[,protect_species,drop=FALSE],
    na.rm=TRUE
  ) /
  length(protect_species)
#############################
# Final risk score
#############################
score_df$Score <-
  (risk_prop - protect_prop + 1)/2
#############################
# 12. High / Low risk stratification
#############################
score_cutoff <-
  median(
    score_df$Score,
    na.rm=TRUE
  )
score_df$Group <-
  ifelse(
    score_df$Score >= score_cutoff,
    "High",
    "Low"
  )
score_df$Group <-
  factor(
    score_df$Group,
    levels=c(
      "Low",
      "High"
    )
  )
cat("\nScore subgroup count:\n")
print(table(score_df$Group))
#############################
# 13. Cox survival regression
#############################
cox_fit <- coxph(
  
  Surv(time,status)~
    Group+
    Gender+
    Age,
  
  data=score_df
)
cox_summary <- summary(cox_fit)
HR <- 
  cox_summary$coefficients[
    "GroupHigh",
    "exp(coef)"
  ]
HR_low <-
  cox_summary$conf.int[
    "GroupHigh",
    "lower .95"
  ]
HR_high <-
  cox_summary$conf.int[
    "GroupHigh",
    "upper .95"
  ]
cox_p <-
  cox_summary$coefficients[
    "GroupHigh",
    "Pr(>|z|)"
  ]
#############################
# Log-rank test
#############################
logrank <- survdiff(
  Surv(time,status)~
    Group,
  data=score_df
)
logrank_p <-
  1-pchisq(
    logrank$chisq,
    df=1
  )
#############################
# Export HR summary table
#############################
score_summary <- data.frame(
  
  top_species=10,
  
  risk_species_number=
    length(risk_species),
  
  protective_species_number=
    length(protect_species),
  
  HR=HR,
  
  HR_low95=HR_low,
  
  HR_high95=HR_high,
  
  Cox_p=cox_p,
  
  Logrank_p=logrank_p
  
)
write.xlsx(
  score_summary,
  "Top10_score_HR_summary.xlsx",
  overwrite=TRUE
)
print(score_summary)
#############################
# 14. Draw KM survival curve
#############################
fit_km <- survfit(
  Surv(time,status)~
    Group,
  data=score_df
)
km <- ggsurvplot(
  
  fit_km,
  
  data=score_df,
  
  pval=TRUE,
  
  risk.table=TRUE,
  
  palette=c(
    "#4C78A8",
    "#D95F59"
  ),
  
  legend.title="Score",
  
  legend.labs=c(
    "Low",
    "High"
  ),
  
  xlab="Time (months)",
  
  ylab="Overall survival probability",
  
  title="Top10 microbial score"
  
)
print(km)
pdf(
  "Top10_microbiome_score_KM.pdf",
  width=7,
  height=7
)
print(km)
dev.off()
#############################
# 15. Export patient-level score table
#############################
write.xlsx(
  score_df,
  "Top10_microbiome_score_patient.xlsx",
  overwrite=TRUE
)
cat("\nTop10 score analysis finished\n")
#############################
# Tertile stratification + Cox HR comparison
#############################
library(survival)
# ==========================
# Tertile split
# ==========================
score_cutoff <- quantile(
  score_df$Score,
  probs = c(0, 1/3, 2/3, 1),
  na.rm = TRUE
)
score_df$Group <- cut(
  score_df$Score,
  breaks = score_cutoff,
  include.lowest = TRUE,
  labels = c(
    "Low",
    "Middle",
    "High"
  )
)
score_df$Group <- factor(
  score_df$Group,
  levels = c(
    "Low",
    "Middle",
    "High"
  )
)
cat("Patient count per tertile group:\n")
print(table(score_df$Group))
# ==========================
# Cox regression model
# ==========================
cox_fit <- coxph(
  Surv(time, status) ~ 
    Group +
    Gender +
    Age,
  data = score_df
)
cox_summary <- summary(cox_fit)
# ==========================
# Extract HR estimates
# ==========================
hr_result <- data.frame(
  
  Comparison = c(
    "Middle_vs_Low",
    "High_vs_Low"
  ),
  
  HR = c(
    cox_summary$conf.int["GroupMiddle", "exp(coef)"],
    cox_summary$conf.int["GroupHigh", "exp(coef)"]
  ),
  
  CI_low95 = c(
    cox_summary$conf.int["GroupMiddle", "lower .95"],
    cox_summary$conf.int["GroupHigh", "lower .95"]
  ),
  
  CI_high95 = c(
    cox_summary$conf.int["GroupMiddle", "upper .95"],
    cox_summary$conf.int["GroupHigh", "upper .95"]
  ),
  
  P_value = c(
    cox_summary$coefficients["GroupMiddle", "Pr(>|z|)"],
    cox_summary$coefficients["GroupHigh", "Pr(>|z|)"]
  )
  
)
# ==========================
# Linear trend test
# ==========================
score_df$Group_num <- as.numeric(score_df$Group)
cox_trend <- coxph(
  Surv(time,status) ~
    Group_num +
    Gender +
    Age,
  data = score_df
)
trend_p <- summary(cox_trend)$coefficients[
  "Group_num",
  "Pr(>|z|)"
]
hr_result$Trend_P <- trend_p
# ==========================
# Print HR results
# ==========================
print(hr_result)
#############################
# Export tertile grouping metadata (preserve original sample order)
#############################
# Sample count per tertile
cat("Sample size per tertile:\n")
print(table(score_df$Group))
# Score distribution statistics per group
group_info <- score_df %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    min_score = min(Score, na.rm = TRUE),
    Q1_score = quantile(Score, 0.25, na.rm = TRUE),
    median_score = median(Score, na.rm = TRUE),
    Q3_score = quantile(Score, 0.75, na.rm = TRUE),
    max_score = max(Score, na.rm = TRUE)
  )
cat("\nScore distribution statistics by tertile:\n")
print(group_info)
# Retain original sample row order
group_samples <- score_df %>%
  select(
    SampleID,
    Score,
    Group
  )
cat("\nPatient list with tertile grouping (original order):\n")
print(group_samples)
# Export to Excel
write.xlsx(
  group_samples,
  "microbial_score_tertile_group_information.xlsx",
  overwrite = TRUE
)
