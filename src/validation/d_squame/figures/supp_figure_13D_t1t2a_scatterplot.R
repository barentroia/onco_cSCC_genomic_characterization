#===============================================================================
# Aim:      Generate Supplementary Figure 13D: scatter plot of weighted PPV
#           versus weighted sensitivity in the BWH T1-T2a subset of the
#           D-SQUAME validation dataset (val102), for SCCore-GEP at the
#           operating point chosen on discovery and for single adverse
#           features (BWH/AJCC stage, NCCN 1.2026 very-high-risk criteria and
#           the additional European Guideline 2026 criteria).
# Input:    - Merged imputed validation clinical data and mice object
#             (results_dir/intermediate/clinical)
#           - NCC weights (e7_weights)
#           - Recalibrated SCCore-GEP predictions (output of recalibration.R)
#           - Discovery threshold metrics (output of Discovery_operating_point.R)
# Output:   Written to results_dir/publication (exactly 3 files):
#           - supp_figure_13D_t1t2a_scatterplot_PPV_vs_SE_scatter.pdf
#           - supp_figure_13D_t1t2a_scatterplot_PPV_vs_SE_scatter_sourcedata.csv
#           - SeSp_supp_figure_13D_t1t2a_scatterplot.csv
# Requires: results_dir, project_dir, dir_results_intermediate_sccore_gep (main.R)
# Usage:    Sourced from main.R, before figure_5d_t1t2a_precision_recall_plot.R
#===============================================================================

# Background
#-------------------------------------------------------------------------------
# One point per rule: x = weighted sensitivity, y = weighted PPV. Weighted
# PPV = TP_w / (TP_w + FP_w), confusion matrix pooled across imputations.
# Because the weights restore the population case/control ratio, this is the
# PPV at the cohort's own weighted prevalence (stated in the y-axis label).
#
# Minimum flagged patients: an adverse feature is dropped from the figure
# unless it flags more than MIN_FLAGGED_FEATURE patients (mean across
# imputations). The operating point is never dropped by this rule.
#
# Operating point: the threshold with the maximal finite weighted LR+ in the
# discovery T1-T2a threshold analysis (0.0698, i.e. 6.98% 5-year absolute risk).
# High LR+ == 10% PPV for this specific threshold matrix
#
# Dataset: val102 is a mice `mids` object with m completed datasets. Cohort
# membership is itself imputed (a patient with missing staging is in T1-T2a
# in some imputations and not in others), so every metric is pooled across
# all m imputations rather than taken from one dataset.
#-------------------------------------------------------------------------------

# 1. Packages
#-------------------------------------------------------------------------------
library(conflicted)
library(tidyverse)
library(survival)
library(mice)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)
conflicts_prefer(base::setdiff)
conflicts_prefer(base::intersect)
conflicts_prefer(base::rbind)
conflicts_prefer(base::cbind)

## ggrepel is optional: used for non-overlapping point labels when installed
HAS_GGREPEL <- requireNamespace("ggrepel", quietly = TRUE)
## Prefix for every output file name
FIGURE_PREFIX <- "supp_figure_13D_t1t2a_scatterplot"

#-------------------------------------------------------------------------------

# 2. Paths and settings
#-------------------------------------------------------------------------------
## val102: the 102-patient D-SQUAME validation set
p102 <- list(
  clinical = file.path(results_dir, "intermediate", "clinical", "merged_clinical_data_imputed_d_squame_validation.csv"),
  mimp = file.path(results_dir, "intermediate", "clinical",
                   "clinical_data_imputed_d_squame_validation.Rdata"),
  weights = e7_weights,
  preds = file.path(results_dir, "intermediate", "validation", "d_squame",
                    "recalibrated_SCCoreGEP_predictions.csv")
)

## Output directory
output_dir <- file.path(results_dir, "publication")
if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = TRUE)}

pred_var <- "sccore_gep"
outcome_var <- "Metastasis_num"
weights_var <- "Weight_rescaled"
thickness.var <- "Breslow_thickness"

### Get operating point from discovery
thr_file <- file.path(dir_results_intermediate_sccore_gep,
                      "discovery_T1_T2a_threshold_metrics_discoveryCal_OOB.csv")
thr_tab <- readr::read_csv(thr_file)

# Criterion: maximum finite weighted LR+; PPV around 10%
op_row <- thr_tab[is.finite(thr_tab$LR_pos_w), ]
op_row <- op_row[which.max(op_row$LR_pos_w), ]

OPERATING_POINTS <- op_row$threshold
#OPERATING_POINTS <- 0.069823237
stopifnot(abs(OPERATING_POINTS - 0.069823237) < 1e-8)  # guard against silent changes

BWH_LEVELS <- c("T1", "T2a")
ID_CANDIDATES <- c("Patient_ID_SKY", "Skyline_ID", "Sample_id", "Skyline_Sample_ID", "Set_id")
SUPPLEMENT_GROUPS <- c("Set_id", "Set_info", "Type_of_material_bin", "Metastasis",
                       "FU_metastasis_years", "OTR_at_cSCC", "HM_at_cSCC", "Tumor_location",
                       "Resection_margin_excision", "Bone_invasion", "Differntiation_broder",
                       "Number_of_cSCC_before_culprit", "BWH_staging", "AJCC_staging")
MIN_FLAGGED_FEATURE <- 5
## Upper limit of the x-axis (weighted sensitivity)
X_MAX <- 0.6

THICKNESS_FEATURE <- "thick_tumor_gt6mm"

## NCCN 1.2026, VERY-HIGH-RISK criteria (+ BWH T2a as staging comparator)
##   Adenosquamous / sarcomatoid histology not buildable (no Morphology_subtype column).
NCCN_FEATURES <- c("very_high_tumor_diameter",
                   THICKNESS_FEATURE,
                   "beyond_subcutaneous_fat",
                   "lymphovascular_invasion",
                   "poor_differentiation",
                   "perineural_invasion",
                   "bwh_T2a",
                   "ajcc_T2_T3")

## European Guideline 2026 - criteria added on top of NCCN
##   Lip/ear/temple location and desmoplasia not buildable in these files.
EUROPEAN_EXTRA_FEATURES <- c("high_tumor_diameter",
                             "bone_invasion",
                             "immunosuppressed",
                             "incomplete_excision")

NCCN_EU_FEATURES <- c(NCCN_FEATURES, EUROPEAN_EXTRA_FEATURES)

## Every adverse feature that can be built (only NCCN_EU_FEATURES are plotted)
ADVERSE_FEATURES <- c("bwh_T2a",
                      "ajcc_T2_T3",
                      "prior_cSCC",
                      "high_tumor_diameter",
                      "very_high_tumor_diameter",
                      "poor_differentiation",
                      "poor_differentiation_broder",
                      "moderate_or_poor_differentiation_broder",
                      "head_neck_location",
                      "perineural_invasion",
                      "lymphovascular_invasion",
                      "immunosuppressed",
                      "incomplete_excision",
                      "thick_tumor_gt6mm",
                      "moderate_thickness_2to6mm",
                      "beyond_subcutaneous_fat",
                      "deep_invasion_gt6mm",
                      "deep_invasion_ge4mm",
                      "bone_invasion")

## Readable labels for the figure
RULE_LABELS <- c(OP1 = "SCCore-GEP",
                 bwh_T2a = "BWH T2a",
                 ajcc_T2_T3 = "AJCC T2-T3",
                 very_high_tumor_diameter = "Diameter \u2265 4 cm",
                 thick_tumor_gt6mm = "Thickness > 6 mm",
                 beyond_subcutaneous_fat = "Beyond subcutaneous fat",
                 lymphovascular_invasion = "Lymphovascular invasion",
                 poor_differentiation = "Poor differentiation",
                 perineural_invasion = "Perineural invasion",
                 high_tumor_diameter = "Diameter \u2265 2 cm",
                 bone_invasion = "Bone invasion",
                 immunosuppressed = "Immunosuppressed",
                 incomplete_excision = "Incomplete excision")

GROUP_LEVELS <- c("SCCore-GEP",
                  "Staging (vs T1)",
                  "NCCN 1.2026 very-high-risk",
                  "European Guideline 2026 (additional)")
GROUP_COLOURS <- setNames(c("#1f77b4", "#d62728", "#2ca02c", "#ff7f0e"), GROUP_LEVELS)
GROUP_SHAPES <- setNames(c(18, 15, 16, 17), GROUP_LEVELS)


#-------------------------------------------------------------------------------

# 3. Check every input path before loading anything
#-------------------------------------------------------------------------------

#' Check that input files exist
#'
#' @param paths named list/character. File paths.
#' @param label character. Label used in the messages.
#'
#' @return character vector of missing paths (empty if all exist).
checkPaths <- function(paths, label) {
  missing <- character(0)
  for (nm in names(paths)) {
    if (!file.exists(paths[[nm]])) {
      missing <- c(missing, paste0("  ", label, "$", nm, ": ", paths[[nm]]))
    }
  }
  return(missing)
}

missing_all <- checkPaths(p102, "p102")

if (length(missing_all) > 0) {
  stop("These input files were not found:\n", paste(missing_all, collapse = "\n"),
       "\n\nMost of the val102 paths hang off results_dir, currently:\n  ", results_dir,
       "\nIf that directory is wrong, fixing the results_dir line is enough.")
}

print("All input files found:")
for (nm in names(p102)) print(paste0("  p102$", nm, ": ", p102[[nm]]))


#-------------------------------------------------------------------------------

# 4. Metric functions
#-------------------------------------------------------------------------------

#' (Weighted) threshold-based metrics for a binary classification
#'
#' Samples with prediction > cutoff are classified as high risk.
#' From L. Pozza, unchanged.
#'
#' @param cutoff numeric. Threshold on the predictions.
#' @param outcome_num numeric vector. Outcome (1 = event, 0 = no event).
#' @param predictions numeric vector. Predicted risks.
#' @param weights numeric vector or NULL. Sample weights (NULL = unweighted).
#'
#' @return Named numeric vector of weighted 2x2 cells, SE, SP, PPV, NPV,
#'   likelihood ratios and proportions at low/high risk.
computeDiscMetricsBinary <- function(cutoff, outcome_num, predictions, weights = NULL){
  
  labels_cutoff <- ifelse(predictions > cutoff, 1, 0)
  
  if (is.null(weights)){
    weights <- rep(1, length(predictions))
  }
  
  correctly_assigned <- which(outcome_num == labels_cutoff)
  wrongly_assigned <-  which(outcome_num != labels_cutoff)
  positives <- which(outcome_num == 1)
  negatives <- which(outcome_num == 0)
  TN <- sum(weights[intersect(correctly_assigned, negatives)])
  TP <- sum(weights[intersect(correctly_assigned, positives)])
  FN <- sum(weights[intersect(wrongly_assigned, positives)])
  FP <- sum(weights[intersect(wrongly_assigned, negatives)])
  Pos <- sum(weights[positives])
  Neg <- sum(weights[negatives])
  
  SE <- if ((TP + FN) > 0) TP / (TP + FN) else NA
  SP <- if ((TN + FP) > 0) TN / (TN + FP) else NA
  PPV <- if ((TP + FP) > 0) TP / (TP + FP) else NA
  NPV <- if ((TN + FN) > 0) TN / (TN + FN) else NA
  
  LR_pos <- NA
  LR_neg <- NA
  if (!is.na(SE) & !is.na(SP)){
    LR_pos <- SE /(1 - SP)
    LR_neg <- (1 - SE) / SP
  }
  
  low_risk_n <- FN + TN
  high_risk_n <- FP + TP
  low_risk_prop <- low_risk_n / (low_risk_n + high_risk_n)
  high_risk_prop <- 1 - low_risk_prop
  
  performance_metrics <- c("FP_w" = FP, "FN_w" = FN, "TP_w" = TP, "TN_w" = TN,
                           "SE_w" = SE, "SP_w" = SP, "PPV_w" = PPV, "NPV_w" = NPV,
                           "LR_pos_w" = LR_pos, "LR_neg_w" = LR_neg,
                           "low_risk_w" = low_risk_prop, "high_risk_w" = high_risk_prop,
                           "positives_w" = Pos, "negatives_w" = Neg,
                           "n_low_risk_w" = low_risk_n, "n_high_risk_w" = high_risk_n)
  
  return(performance_metrics)
}

#' PPV at an assumed prevalence (Bayes' rule)
#'
#' Despite the name, this returns Se*pi / (Se*pi + (1 - Sp)*(1 - pi)).
#'
#' @param sensitivity numeric. Sensitivity.
#' @param specificity numeric. Specificity.
#' @param prevalence numeric. Assumed prevalence (pi).
#'
#' @return numeric. Positive predictive value.
computeEER <- function(sensitivity, specificity, prevalence) {
  return(sensitivity * prevalence /
           (sensitivity * prevalence + (1 - specificity) * (1 - prevalence)))
}

#' Division that returns NA when the denominator is 0 or NA
#'
#' @param a numeric. Numerator.
#' @param b numeric. Denominator.
#'
#' @return numeric.
safeDiv <- function(a, b) {
  if (is.na(a) || is.na(b) || b <= 0) return(NA_real_)
  return(a / b)
}

#' Minimum ignoring NA; NA if all values are NA
#'
#' @param x numeric vector.
#'
#' @return numeric.
safeMin <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
#' Maximum ignoring NA; NA if all values are NA
#'
#' @param x numeric vector.
#'
#' @return numeric.
safeMax <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)

#' Membership test that keeps NA as NA
#'
#' @param x vector. Values to test.
#' @param set vector. Allowed values.
#'
#' @return logical vector (NA where x is NA).
inSet <- function(x, set) {
  out <- as.integer(as.character(x) %in% set)
  out[is.na(x)] <- NA
  return(out)
}

#' First candidate column present in a data frame
#'
#' @param df data.frame.
#' @param candidates character vector. Candidate column names.
#'
#' @return character (column name) or NULL if none is present.
pickCol <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(NULL)
  return(hit[1])
}

#' Numeric version of a column, or all-NA of the right length if absent
#'
#' @param df data.frame.
#' @param candidates character vector. Candidate column names.
#'
#' @return numeric vector of length nrow(df).
numCol <- function(df, candidates) {
  cn <- pickCol(df, candidates)
  if (is.null(cn)) return(rep(NA_real_, nrow(df)))
  return(suppressWarnings(as.numeric(as.character(df[[cn]]))))
}

#' Weighted and unweighted Se/Sp per rule for one completed dataset
#'
#' @param df data.frame. One completed dataset with binary rule columns.
#' @param rule_labels character vector. Rule (column) names to evaluate.
#' @param outcome_var character. Outcome column.
#' @param weights_var character. Weights column.
#'
#' @return data.frame with one row per rule: 2x2 cells, SE, SP, PPV, NPV, ...
getSeSpTable <- function(df, rule_labels, outcome_var, weights_var) {
  rule_labels <- rule_labels[rule_labels %in% colnames(df)]
  
  se_sp_table <- do.call(rbind, lapply(rule_labels, function(rule_label) {
    rule <- df[[rule_label]]
    outcome <- df[[outcome_var]]
    weights <- df[[weights_var]]
    
    i.ok <- which(!is.na(rule) & !is.na(outcome) & !is.na(weights))
    
    m <- computeDiscMetricsBinary(0.5, outcome[i.ok], rule[i.ok], weights[i.ok])
    mu <- computeDiscMetricsBinary(0.5, outcome[i.ok], rule[i.ok])
    
    data.frame(
      Rule = rule_label,
      SE_w = m[["SE_w"]], SP_w = m[["SP_w"]], PPV_w = m[["PPV_w"]], NPV_w = m[["NPV_w"]],
      LR_pos_w = m[["LR_pos_w"]], LR_neg_w = m[["LR_neg_w"]],
      TP_w = m[["TP_w"]], FP_w = m[["FP_w"]], FN_w = m[["FN_w"]], TN_w = m[["TN_w"]],
      SE_unw = mu[["SE_w"]], SP_unw = mu[["SP_w"]], PPV_unw = mu[["PPV_w"]],
      NPV_unw = mu[["NPV_w"]],
      TP_unw = mu[["TP_w"]], FP_unw = mu[["FP_w"]], FN_unw = mu[["FN_w"]], TN_unw = mu[["TN_w"]],
      n_flagged = sum(rule == 1, na.rm = TRUE),
      n_flagged_cases = sum(rule == 1 & outcome == 1, na.rm = TRUE),
      n_flagged_controls = sum(rule == 1 & outcome == 0, na.rm = TRUE),
      n_missing = sum(is.na(rule)),
      stringsAsFactors = FALSE
    )
  }))
  
  return(se_sp_table)
}

#' Pool the Se/Sp tables across the completed datasets
#'
#' The four confusion matrices are averaged and Se, Sp, PPV, NPV and LRs are
#' recomputed from those averages. The per-imputation min/max of SE_w and
#' PPV_w are kept in the source data.
#'
#' @param datasets list of data.frames. Completed (imputed) datasets.
#' @param rule_labels character vector. Rules to evaluate.
#' @param outcome_var character. Outcome column.
#' @param weights_var character. Weights column.
#'
#' @return data.frame with one row per rule (pooled metrics, ranges and
#'   mean number of flagged patients).
getSeSpTablePooled <- function(datasets, rule_labels, outcome_var, weights_var) {
  tabs <- lapply(datasets, getSeSpTable, rule_labels = rule_labels,
                 outcome_var = outcome_var, weights_var = weights_var)
  tabs <- tabs[!vapply(tabs, is.null, logical(1))]
  
  if (length(tabs) == 0) return(NULL)
  
  rules <- tabs[[1]]$Rule
  
  pooled <- do.call(rbind, lapply(rules, function(r) {
    rows <- do.call(rbind, lapply(tabs, function(t) t[t$Rule == r, , drop = FALSE]))
    
    TP <- mean(rows$TP_w, na.rm = TRUE); FP <- mean(rows$FP_w, na.rm = TRUE)
    FN <- mean(rows$FN_w, na.rm = TRUE); TN <- mean(rows$TN_w, na.rm = TRUE)
    
    TPu <- mean(rows$TP_unw, na.rm = TRUE); FPu <- mean(rows$FP_unw, na.rm = TRUE)
    FNu <- mean(rows$FN_unw, na.rm = TRUE); TNu <- mean(rows$TN_unw, na.rm = TRUE)
    
    SE <- safeDiv(TP, TP + FN); SP <- safeDiv(TN, TN + FP)
    SEu <- safeDiv(TPu, TPu + FNu); SPu <- safeDiv(TNu, TNu + FPu)
    
    data.frame(
      Rule = r,
      SE_w = SE, SP_w = SP,
      PPV_w = safeDiv(TP, TP + FP), NPV_w = safeDiv(TN, TN + FN),
      LR_pos_w = if (!is.na(SE) && !is.na(SP) && SP < 1) SE / (1 - SP) else NA_real_,
      LR_neg_w = if (!is.na(SE) && !is.na(SP) && SP > 0) (1 - SE) / SP else NA_real_,
      TP_w = TP, FP_w = FP, FN_w = FN, TN_w = TN,
      SE_unw = SEu, SP_unw = SPu,
      PPV_unw = safeDiv(TPu, TPu + FPu), NPV_unw = safeDiv(TNu, TNu + FNu),
      TP_unw = TPu, FP_unw = FPu, FN_unw = FNu, TN_unw = TNu,
      n_flagged = mean(rows$n_flagged, na.rm = TRUE),
      n_flagged_cases = mean(rows$n_flagged_cases, na.rm = TRUE),
      n_flagged_controls = mean(rows$n_flagged_controls, na.rm = TRUE),
      n_missing = mean(rows$n_missing, na.rm = TRUE),
      SE_w_min = safeMin(rows$SE_w), SE_w_max = safeMax(rows$SE_w),
      SP_w_min = safeMin(rows$SP_w), SP_w_max = safeMax(rows$SP_w),
      PPV_w_min = safeMin(rows$PPV_w), PPV_w_max = safeMax(rows$PPV_w),
      n_datasets = nrow(rows),
      stringsAsFactors = FALSE
    )
  }))
  
  return(pooled)
}

#-------------------------------------------------------------------------------

# 5. Plotting: weighted PPV vs weighted sensitivity
#-------------------------------------------------------------------------------

#' Figure group of each rule (SCCore-GEP, staging, NCCN, European Guideline)
#'
#' @param rule character vector. Rule names.
#'
#' @return factor with levels GROUP_LEVELS.
ruleGroup <- function(rule) {
  out <- ifelse(rule %in% c("OP1", "OP-low", "OP-high"), GROUP_LEVELS[1],
                ifelse(rule %in% c("bwh_T2a", "ajcc_T2_T3"), GROUP_LEVELS[2],
                       ifelse(rule %in% NCCN_FEATURES, GROUP_LEVELS[3],
                              ifelse(rule %in% EUROPEAN_EXTRA_FEATURES, GROUP_LEVELS[4], NA_character_))))
  return(factor(out, levels = GROUP_LEVELS))
}

#' Readable figure label of each rule (falls back to the rule name)
#'
#' @param rule character vector. Rule names.
#'
#' @return character vector.
ruleLabel <- function(rule) {
  out <- unname(RULE_LABELS[rule])
  out[is.na(out)] <- rule[is.na(out)]
  return(out)
}

#' Scatter plot of weighted PPV versus weighted sensitivity
#'
#' @param se_sp_table data.frame. Pooled metrics per rule.
#' @param dir.output character. Output folder.
#' @param file.ID character. File name prefix.
#' @param title_plot character. Plot title.
#' @param prevalence_line numeric. Weighted prevalence, shown in the y label.
#'
#' @return Invisibly, the source data; saves the pdf and source data csv.
plot_ppv_vs_se <- function(se_sp_table, dir.output, file.ID, title_plot,
                           prevalence_line = NA) {
  
  df <- se_sp_table
  df$Group <- ruleGroup(df$Rule)
  df$Label <- ruleLabel(df$Rule)
  
  i.drop <- is.na(df$SE_w) | is.na(df$PPV_w)
  if (any(i.drop)) {
    print(paste0("  not plotted (SE or PPV undefined): ",
                 paste(df$Rule[i.drop], collapse = ", ")))
  }
  df <- df[!i.drop, , drop = FALSE]
  
  if (nrow(df) == 0) {
    print("  no usable rules; no figure produced.")
    return(invisible(NULL))
  }
  
  ## Flag rules whose weighted 2x2 has an empty cell: the point is still drawn, but SE or PPV
  ## sits on a boundary (0 or 1) and should be read with care.
  df$empty_2x2_cell <- with(df, TP_w == 0 | FP_w == 0 | FN_w == 0 | TN_w == 0)
  
  df <- df[order(df$Group, -df$PPV_w), ]
  
  fn <- file.path(dir.output, paste0(file.ID, "_PPV_vs_SE_scatter.pdf"))
  
  ## ---- source data: exactly the values drawn ----
  ##   SE_w            x position of the point
  ##   PPV_w           y position of the point
  ##   *_min / *_max   range across imputations (not drawn; for reference only)
  ##   cohort_prevalence_weighted   the prevalence at which PPV_w holds (stated in the y label)
  ##   ref_line_vertical_SE         dashed vertical line (SCCore-GEP sensitivity)
  ##   ref_line_horizontal_PPV      dashed horizontal line (SCCore-GEP PPV)
  src <- data.frame(
    Rule = df$Rule,
    Label = df$Label,
    Group = as.character(df$Group),
    SE_w = df$SE_w,
    PPV_w = df$PPV_w,
    NPV_w = df$NPV_w,
    SE_w_min = df$SE_w_min,
    SE_w_max = df$SE_w_max,
    PPV_w_min = df$PPV_w_min,
    PPV_w_max = df$PPV_w_max,
    PPV_at_cohort_prevalence = computeEER(df$SE_w, df$SP_w, prevalence_line),
    SP_w = df$SP_w,
    LR_pos_w = df$LR_pos_w,
    LR_neg_w = df$LR_neg_w,
    TP_w = df$TP_w, FP_w = df$FP_w, FN_w = df$FN_w, TN_w = df$TN_w,
    n_flagged = df$n_flagged,
    n_flagged_cases = df$n_flagged_cases,
    n_flagged_controls = df$n_flagged_controls,
    empty_2x2_cell = df$empty_2x2_cell,
    n_datasets = df$n_datasets,
    cohort_prevalence_weighted = prevalence_line,
    ref_line_vertical_SE = if (any(df$Rule == "OP1")) df$SE_w[df$Rule == "OP1"][1] else NA_real_,
    ref_line_horizontal_PPV = if (any(df$Rule == "OP1")) df$PPV_w[df$Rule == "OP1"][1] else NA_real_,
    stringsAsFactors = FALSE
  )
  fn_src <- sub("\\.pdf$", "_sourcedata.csv", fn)
  write.csv(src, fn_src, row.names = FALSE)
  print(paste0("  wrote ", basename(fn_src), " (", nrow(src), " rows)"))
  
  ## ---- figure ----
  ## All points plain black. Dashed reference lines through the SCCore-GEP point:
  ## vertical at its sensitivity, horizontal at its PPV.
  y_top <- max(df$PPV_w, na.rm = TRUE) * 1.15
  prev_txt <- if (is.na(prevalence_line)) "" else
    paste0(" (at cohort prevalence ", signif(100 * prevalence_line, 2), "%)")
  
  i.gep <- which(df$Rule == "OP1")
  
  outside <- df$Rule[df$SE_w > X_MAX]
  if (length(outside) > 0) {
    warning("Outside the x-axis (sensitivity > ", 100 * X_MAX, "%), not shown: ",
            paste(outside, collapse = ", "))
  }
  
  p <- ggplot(df, aes(x = SE_w, y = PPV_w))
  
  if (length(i.gep) == 1) {
    p <- p +
      geom_vline(xintercept = df$SE_w[i.gep], linetype = 2, colour = "grey40") +
      geom_hline(yintercept = df$PPV_w[i.gep], linetype = 2, colour = "grey40")
  }
  
  p <- p +
    geom_point(colour = "black", shape = 16, size = 3) +
    scale_x_continuous(limits = c(0, X_MAX), breaks = seq(0, X_MAX, 0.1),
                       labels = scales::label_percent(accuracy = 1),
                       expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(limits = c(0, y_top),
                       labels = scales::label_percent(accuracy = 0.1),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Weighted sensitivity",
         y = paste0("Weighted PPV"),
         title = title_plot,
         caption = paste0("Dashed lines: sensitivity and PPV of SCCore-GEP. ",
                          "Features flagging \u2264 ", MIN_FLAGGED_FEATURE,
                          " patients (mean) are omitted.")) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 11),
          plot.caption = element_text(size = 8, colour = "grey30"))
  
  if (HAS_GGREPEL) {
    p <- p + ggrepel::geom_text_repel(aes(label = Label), colour = "black", size = 3,
                                      min.segment.length = 0, box.padding = 0.4,
                                      max.overlaps = Inf, seed = 1)
  } else {
    p <- p + geom_text(aes(label = Label), colour = "black", size = 3, vjust = -1)
  }
  
  ggsave(fn, p, width = 5, height = 6, device = cairo_pdf)
  print(paste0("  wrote ", basename(fn)))
  
  return(invisible(src))
}

#-------------------------------------------------------------------------------

# 6. Adverse features
#-------------------------------------------------------------------------------

#' Add binary adverse-feature columns (size, depth, differentiation, ...)
#'
#' @param df data.frame. One completed dataset.
#' @param thickness.var character. Column with the tumor thickness (mm).
#' @param verbose logical. Print a summary of the built features.
#'
#' @return df with one 0/1 column per feature in ADVERSE_FEATURES.
buildAdverseFeatures <- function(df, thickness.var = "Breslow_thickness", verbose = TRUE) {
  
  df <- as.data.frame(df)
  
  ## Size
  dia <- numCol(df, "Tumor_diameter")
  df[["high_tumor_diameter"]] <- as.integer(dia >= 2)
  df[["very_high_tumor_diameter"]] <- as.integer(dia >= 4)
  
  df[["head_neck_location"]] <- inSet(df$Tumor_location, c("Face", "Scalp/neck"))
  
  ## Histology
  df[["poor_differentiation"]] <- inSet(df$Differentiation, "Poor/undifferentiated")
  
  cn_broder <- pickCol(df, c("Differntiation_broder", "Differentiation_broder",
                             "Differntiation_Broder", "Differentiation_Broder"))
  if (!is.null(cn_broder)) {
    lv <- sort(unique(as.character(df[[cn_broder]])))
    lv_poor <- lv[grepl("poor", lv, ignore.case = TRUE)]
    lv_mod <- lv[grepl("moderate", lv, ignore.case = TRUE)]
    
    if (verbose) {
      print(paste0("  Broder column: ", cn_broder, " | levels: ", paste(lv, collapse = " | ")))
    }
    if (length(lv_poor) > 0) {
      df[["poor_differentiation_broder"]] <- inSet(df[[cn_broder]], lv_poor)
    }
    if (length(c(lv_poor, lv_mod)) > 0) {
      df[["moderate_or_poor_differentiation_broder"]] <- inSet(df[[cn_broder]],
                                                               c(lv_poor, lv_mod))
    }
  } else if (verbose) {
    print("  no Broder differentiation column found")
  }
  
  ## Invasion
  df[["perineural_invasion"]] <- inSet(df$PNI_bin, "Yes")
  
  if (!("LVI_bin" %in% colnames(df))) {
    df$LVI_bin <- ifelse(df$PNI_or_LVI == "Yes" & df$PNI_bin == "No", "Yes",
                         ifelse(df$PNI_or_LVI == "No", "No", NA))
  }
  df[["lymphovascular_invasion"]] <- inSet(df$LVI_bin, "Yes")
  
  df[["beyond_subcutaneous_fat"]] <- inSet(df$Tissue_involvement, "Beyond subcutaneous fat")
  
  cn_bone <- pickCol(df, c("Bone_invasion", "Invasion_of_bones"))
  df[["bone_invasion"]] <- if (is.null(cn_bone)) NA_integer_ else inSet(df[[cn_bone]], "Yes")
  
  ## Thickness and depth
  thick <- numCol(df, c(thickness.var, "Breslow_thickness"))
  df[["thick_tumor_gt6mm"]] <- as.integer(thick > 6)
  df[["moderate_thickness_2to6mm"]] <- as.integer(thick >= 2 & thick <= 6)
  
  doi <- numCol(df, c("Depth_of_invasion", "Depth_of_Invasion"))
  df[["deep_invasion_gt6mm"]] <- as.integer(doi > 6)
  df[["deep_invasion_ge4mm"]] <- as.integer(doi >= 4)
  
  ## BWH stage within the cohort: 1 for T2a, 0 for T1 (staging comparator)
  cn_bwh <- pickCol(df, c("BWH_staging", "BWH"))
  df[["bwh_T2a"]] <- if (is.null(cn_bwh)) NA_integer_ else inSet(df[[cn_bwh]], "T2a")
  ## AJCC stage within the cohort: 1 for T2 or T3, 0 for T1 (staging comparator)
  cn_ajcc <- pickCol(df, c("AJCC_staging", "AJCC"))
  df[["ajcc_T2_T3"]] <- if (is.null(cn_ajcc)) NA_integer_ else inSet(df[[cn_ajcc]], c("T2", "T3"))
  
  ## Other
  df[["prior_cSCC"]] <- as.integer(numCol(df, "Number_of_cSCC_before_culprit") > 0)
  df[["immunosuppressed"]] <- inSet(df$Immunosuppressed, "Yes")
  df[["incomplete_excision"]] <- inSet(df$Resection_margin_excision, "Incomplete")
  
  return(df)
}


#-------------------------------------------------------------------------------

# 7. Loaders
#-------------------------------------------------------------------------------

#' First ID column shared by two data frames
#'
#' @param a data.frame.
#' @param b data.frame.
#' @param candidates character vector. Candidate ID columns.
#' @param what character. Description used in the error message.
#'
#' @return character. Column name (error if none is shared).
shared_id <- function(a, b, candidates, what) {
  hit <- candidates[candidates %in% colnames(a) & candidates %in% colnames(b)]
  if (length(hit) == 0) {
    stop(sprintf("No shared %s. Left has [%s], right has [%s]", what,
                 paste(intersect(candidates, colnames(a)), collapse = ", "),
                 paste(intersect(candidates, colnames(b)), collapse = ", ")))
  }
  return(hit[1])
}

#' Read SCCore-GEP predictions on the percent scale
#'
#' The prediction column and the scale (proportion or percent) are detected
#' from the data.
#'
#' @param path character. Path to the predictions csv.
#' @param id_candidates character vector. Candidate ID columns.
#'
#' @return data.frame with the ID column and sccore_gep (%), one row per ID.
readPreds <- function(path, id_candidates) {
  raw <- read.csv(path)
  
  pred_col <- pickCol(raw, c("sccore_gep_cox_calibrated", "cox", "sccore_gep"))
  if (is.null(pred_col)) {
    stop("No prediction column in ", path, ". Columns present: ",
         paste(colnames(raw), collapse = ", "))
  }
  
  vals <- suppressWarnings(as.numeric(raw[[pred_col]]))
  scale_factor <- if (max(vals, na.rm = TRUE) <= 1) 100 else 1
  
  print(paste0("  prediction column: ", pred_col,
               " | range ", signif(min(vals, na.rm = TRUE), 3), " - ",
               signif(max(vals, na.rm = TRUE), 3),
               " -> x", scale_factor, " to reach the percent scale"))
  
  id_col <- id_candidates[id_candidates %in% colnames(raw)][1]
  out <- data.frame(raw[[id_col]], vals * scale_factor, stringsAsFactors = FALSE)
  names(out) <- c(id_col, "sccore_gep")
  
  return(out[!duplicated(out[[id_col]]), ])
}

#' Add the NCC weights and the predictions to a completed dataset
#'
#' @param d data.frame. Completed dataset.
#' @param weights data.frame. NCC weights.
#' @param preds data.frame. Output of readPreds().
#'
#' @return data.frame with Weight_rescaled and sccore_gep added.
attachWeightsAndPreds <- function(d, weights, preds) {
  id_w <- shared_id(d, weights, ID_CANDIDATES, "ID for weights")
  id_p <- shared_id(d, preds, ID_CANDIDATES, "ID for predictions")
  
  w_small <- weights[, c(id_w, "Weight_rescaled"), drop = FALSE]
  w_small <- w_small[!is.na(w_small$Weight_rescaled), , drop = FALSE]
  w_small <- w_small[!duplicated(w_small[[id_w]]), , drop = FALSE]
  
  d <- merge(d, w_small, by = id_w, all.x = TRUE, sort = FALSE)
  d <- merge(d, preds, by = id_p, all.x = TRUE, sort = FALSE)
  
  return(as.data.frame(d))
}

#' Add Metastasis_num and Immunosuppressed if they are missing
#'
#' @param d data.frame. Completed dataset.
#'
#' @return d with the derived columns.
addDerived <- function(d) {
  if (!("Metastasis_num" %in% colnames(d)) && "Metastasis" %in% colnames(d)) {
    d$Metastasis_num <- ifelse(as.character(d$Metastasis) == "Case", 1, 0)
  }
  if (!("Immunosuppressed" %in% colnames(d)) &&
      all(c("OTR_at_cSCC", "HM_at_cSCC") %in% colnames(d))) {
    d$Immunosuppressed <- ifelse(d$OTR_at_cSCC == "Yes" | d$HM_at_cSCC == "Yes", "Yes", "No")
  }
  return(d)
}

#' Load the first mids object from an .Rdata file
#'
#' @param path character. Path to the .Rdata file.
#'
#' @return mids object.
load_mids <- function(path){
  e <- new.env()
  load(path, envir = e)
  objs <- ls(e)
  is_mids <- vapply(objs, function(o) inherits(get(o, envir = e), "mids"), logical(1))
  if (!any(is_mids)){
    stop(sprintf("No mids object in %s (objects: %s)", path, paste(objs, collapse = ", ")))
  }
  return(get(objs[is_mids][1], envir = e))
}

#' Load the D-SQUAME validation cohort (val102) as completed datasets
#'
#' Completes every imputation, supplements columns dropped by mice from the
#' raw clinical file, and adds derived variables, weights, predictions and
#' adverse features. Uses the paths in p102.
#'
#' @return list with label ("val102") and datasets (list of data.frames).
load_val102 <- function() {
  message("Loading cohort: 102 D-SQUAME validation")
  
  mimp <- load_mids(p102$mimp)
  print(mimp)
  
  clinical_raw <- read.csv(p102$clinical)
  weights <- read.csv(p102$weights)
  preds <- readPreds(p102$preds, ID_CANDIDATES)
  
  datasets <- lapply(seq_len(mimp$m), function(i) as.data.frame(complete(mimp, i)))
  print(paste0("  m = ", length(datasets), " completed datasets, n = ", nrow(datasets[[1]])))
  
  ## Put back columns mice dropped. Nothing already present is overwritten.
  join_id <- shared_id(datasets[[1]], clinical_raw, ID_CANDIDATES, "ID for the raw clinical file")
  add <- setdiff(intersect(SUPPLEMENT_GROUPS, colnames(clinical_raw)), colnames(datasets[[1]]))
  if (length(add) > 0) {
    print(paste0("  supplementing from the raw clinical file: ", paste(add, collapse = ", ")))
    raw_small <- clinical_raw[, c(join_id, add), drop = FALSE]
    raw_small <- raw_small[!duplicated(raw_small[[join_id]]), , drop = FALSE]
    datasets <- lapply(datasets, function(d) merge(d, raw_small, by = join_id,
                                                   all.x = TRUE, sort = FALSE))
  }
  
  datasets <- lapply(seq_along(datasets), function(i) {
    d <- addDerived(datasets[[i]])
    d <- attachWeightsAndPreds(d, weights, preds)
    d <- buildAdverseFeatures(d, thickness.var = thickness.var, verbose = (i == 1))
    return(d)
  })
  
  d1 <- datasets[[1]]
  print(paste0("  n = ", nrow(d1),
               " | without a prediction: ", sum(is.na(d1[[pred_var]])),
               " | without a weight: ", sum(is.na(d1[[weights_var]]))))
  print("  BWH_staging, imputation 1:")
  print(table(d1$BWH_staging, useNA = "ifany"))
  
  n_t1t2a <- vapply(datasets, function(d)
    sum(d$BWH_staging %in% BWH_LEVELS, na.rm = TRUE), integer(1))
  print(paste0("  T1-T2a size per imputation: mean ", round(mean(n_t1t2a), 1),
               " (range ", min(n_t1t2a), "-", max(n_t1t2a), ")"))
  
  return(list(label = "val102", datasets = datasets))
}


#-------------------------------------------------------------------------------

# 8. Feature-set and minimum-flagged filters
#-------------------------------------------------------------------------------

#' Keep the requested features that flag enough patients
#'
#' Operating point rules are always kept.
#'
#' @param se_sp_table data.frame. Pooled metrics per rule.
#' @param keep_features character vector. Features to show.
#' @param min_flagged numeric. Features flagging <= min_flagged patients
#'   (mean across imputations) are dropped.
#' @param label character. Label used in the messages.
#'
#' @return data.frame with the operating point rows and the kept features.
selectFeatures <- function(se_sp_table, keep_features, min_flagged, label) {
  
  op_rule <- se_sp_table$Rule %in% c("OP1", "OP-low", "OP-high")
  op_rows <- se_sp_table[op_rule, , drop = FALSE]
  feat_rows <- se_sp_table[!op_rule, , drop = FALSE]
  
  missing <- setdiff(keep_features, feat_rows$Rule)
  if (length(missing) > 0) {
    print(paste0("  ", label, ": not available in these data - ",
                 paste(missing, collapse = ", ")))
  }
  
  feat_rows <- feat_rows[feat_rows$Rule %in% keep_features, , drop = FALSE]
  
  too_few <- feat_rows$n_flagged <= min_flagged
  if (any(too_few)) {
    print(paste0("  ", label, ": dropped for flagging <= ", min_flagged, " patients - ",
                 paste(sprintf("%s (%.1f)", feat_rows$Rule[too_few],
                               feat_rows$n_flagged[too_few]), collapse = ", ")))
  }
  
  feat_rows <- feat_rows[!too_few, , drop = FALSE]
  
  out <- rbind(op_rows, feat_rows)
  print(paste0("  ", label, ": ", nrow(feat_rows), " features kept plus ",
               nrow(op_rows), " operating point(s)"))
  
  return(out)
}

#-------------------------------------------------------------------------------

# 9. One cohort, one operating point
#-------------------------------------------------------------------------------

#' Run the analysis for one cohort at one operating point
#'
#' Restricts each imputation to BWH T1-T2a, flags patients above the
#' operating point, pools Se/Sp for SCCore-GEP and the adverse features,
#' writes the SeSp table and draws the figure.
#'
#' @param cohort list. Output of load_val102().
#' @param op_proportion numeric. Operating point as a proportion (5-year risk).
#'
#' @return list with file.ID, se_sp (pooled table), sourcedata and more.
runOnePoint <- function(cohort, op_proportion) {
  
  ## sccore_gep is a PERCENTAGE; the operating point is quoted as a proportion
  op_percent <- 100 * op_proportion
  file.ID <- FIGURE_PREFIX
  
  print(paste0("=================== ", file.ID, " ==================="))
  
  ## Restrict every completed dataset to T1-T2a, then apply the cut
  subsets <- lapply(cohort$datasets, function(d) {
    i <- which(!is.na(d$BWH_staging) & d$BWH_staging %in% BWH_LEVELS)
    s <- d[i, , drop = FALSE]
    s[["OP1"]] <- as.integer(suppressWarnings(as.numeric(s[[pred_var]])) > op_percent)
    return(s)
  })
  
  n_i <- vapply(subsets, nrow, integer(1))
  pos_i <- vapply(subsets, function(d) sum(d[[outcome_var]] == 1, na.rm = TRUE), integer(1))
  neg_i <- vapply(subsets, function(d) sum(d[[outcome_var]] == 0, na.rm = TRUE), integer(1))
  flag_i <- vapply(subsets, function(d) sum(d[["OP1"]] == 1, na.rm = TRUE), integer(1))
  
  print(paste0("  operating point ", op_proportion, " as a proportion = ", op_percent,
               "% on the sccore_gep scale"))
  print(paste0("  n: mean ", round(mean(n_i), 1), " (", min(n_i), "-", max(n_i), ")",
               " | cases ", round(mean(pos_i), 1),
               " | controls ", round(mean(neg_i), 1),
               " | flagged ", round(mean(flag_i), 1),
               " (", min(flag_i), "-", max(flag_i), ")"))
  
  if (mean(flag_i) < 5) {
    warning(file.ID, ": the cut flags a mean of only ", round(mean(flag_i), 1), " patients.")
  }
  
  ## Prevalence, weighted and unweighted, averaged across the completed datasets
  prev_w_i <- vapply(subsets, function(d)
    sum(d[[weights_var]] * d[[outcome_var]], na.rm = TRUE) /
      sum(d[[weights_var]], na.rm = TRUE), numeric(1))
  prev_unw_i <- vapply(subsets, function(d)
    mean(d[[outcome_var]] == 1, na.rm = TRUE), numeric(1))
  
  prev_w <- mean(prev_w_i)
  prev_unw <- mean(prev_unw_i)
  
  print(paste0("  weighted prevalence ", round(prev_w, 5),
               " (range ", round(min(prev_w_i), 5), "-", round(max(prev_w_i), 5), ")",
               " | unweighted ", round(prev_unw, 4)))
  
  ## Pooled weighted Se/Sp/PPV for OP1 and every adverse feature
  se_sp_table <- getSeSpTablePooled(subsets, c("OP1", ADVERSE_FEATURES),
                                    outcome_var, weights_var)
  se_sp_table <- se_sp_table[order(-se_sp_table$LR_pos_w), ]
  
  ## OP1 + BWH T2a + NCCN 1.2026 + European Guideline 2026 features
  tab_fs <- selectFeatures(se_sp_table, NCCN_EU_FEATURES, MIN_FLAGGED_FEATURE,
                           "NCCN + European")
  
  print(tab_fs[, c("Rule", "SE_w", "PPV_w", "SP_w", "LR_pos_w", "n_flagged",
                   "SE_w_min", "SE_w_max", "PPV_w_min", "PPV_w_max")],
        digits = 4, row.names = FALSE)
  
  write.csv(tab_fs, file.path(output_dir, paste0("SeSp_", file.ID, ".csv")), row.names = FALSE)
  
  ttl <- sprintf("%s  BWH T1-T2a  (SCCore-GEP cut-off %.4f%%)\nSCCore-GEP, BWH T2a, AJCC T2-T3, NCCN 1.2026 + European Guideline 2026",
                 cohort$label, op_percent)
  
  src <- plot_ppv_vs_se(tab_fs, output_dir, file.ID, ttl, prevalence_line = prev_w)
  
  return(list(file.ID = file.ID, se_sp = tab_fs, sourcedata = src,
              prevalence_weighted = prev_w, prevalence_unweighted = prev_unw))
}



#-------------------------------------------------------------------------------

# 10. Run VAL102 at the selected operating point
#-------------------------------------------------------------------------------
cohorts <- list(load_val102())

results <- list()
for (co in cohorts) {
  for (op in OPERATING_POINTS) {
    results[[length(results) + 1]] <- runOnePoint(co, op)
  }
}
#-------------------------------------------------------------------------------