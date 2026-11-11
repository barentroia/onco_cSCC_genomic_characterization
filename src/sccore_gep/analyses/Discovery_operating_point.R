#===============================================================================
# Aim:      Compute the threshold-based metrics used to choose the SCCore-GEP
#           operating point on the BWH T1-T2a subset of the D-SQUAME discovery
#           dataset. Thresholds are expressed as 5-year absolute metastasis
#           risk, the same units as in validation. Only discovery data are
#           used; no validation data are read in this script.
# Input:    - Final SCCore-GEP model with OOB predictions
#             (Discovery_models/coxnet_Late_feature_integration_robj.rds)
#           - Merged imputed discovery clinical data
#             (results_dir/intermediate/clinical)
#           - NCC weights (p10_weights)
# Output:   Written to dir_results_intermediate_sccore_gep:
#           - discovery_T1_T2a_threshold_metrics_discoveryCal_<source>.csv:
#             weighted SE, SP, PPV, NPV, ... at every candidate threshold
#           - precision_recall_T1T2a_and_entire_cohort_discoveryCal_<source>.pdf:
#             precision-recall curves with the operating point (section 5)
# Requires: dir_scripts_sccore_gep, dir_results_intermediate_sccore_gep,
#           results_dir (main.R) and config.R paths
# Usage:    Sourced from main.R
# Notes:    The operating point itself  is read
#           from the threshold metrics by supp_figure_13D_t1t2a_scatterplot.R.
#           The precision-recall curves in section 5 support that choice.
#===============================================================================

# Background
#-------------------------------------------------------------------------------
# What "the same space as validation" means here:
# It means the same UNITS: 5-year absolute metastasis risk. A threshold of 6%
# means "6% predicted five-year risk" in either cohort, because both sides
# express their predictions as that same quantity. It does not mean the two
# cohorts share one fitted mapping.
#
# The recalibration is fitted on DISCOVERY, so nothing about validation enters
# the choice of threshold. Two separately fitted calibrations are
# not guaranteed to agree, so the threshold transfers only to the extent that
# both calibrations are good.
#
# Which discovery predictions:
# OUT-OF-BAG predictions by default (PREDICTION_SOURCE). Each patient's
# held-out scores across the bootstrap fits are averaged into one score. This
# removes optimism from (1) the threshold and (2) the calibration slope beta,
# which in-sample scores would push toward zero.
#
# Caveat for the methods: the averaged OOB score is a bagged estimate, a
# slightly different object from the single final model.
#-------------------------------------------------------------------------------

# 1. Packages
#-------------------------------------------------------------------------------
library(dplyr)
library(tibble)
library(survival)
library(readr)
library(ggplot2)
#-------------------------------------------------------------------------------

# 2. Paths and settings
#-------------------------------------------------------------------------------
path.model.Data <- file.path(dir_results_intermediate_sccore_gep,  "Discovery_models", "coxnet_Late_feature_integration_robj.rds")

path.clinical.Data <- file.path(results_dir,  "intermediate",
                                "clinical", "merged_clinical_data_imputed_d_squame_discovery.csv")

## Sampling weights (366 samples). Columns: Patient_ID_SKY, Samp_prob_patient, Weight,
## Weight_rescaled, Metastasis, SkylineDx.ID
path.weights.Data <- p10_weights
## Column in the weights file that matches clinical$SkylineDx.ID (and the OOB SkylineDx.ID)
WEIGHTS_ID_COL <- "SkylineDx.ID"

## "OOB"       averaged out-of-bag scores
## "Model_fit" the final model's own in-sample predictions
PREDICTION_SOURCE <- "OOB"
slot.final.model  <- "Model_fit"

output_dir <- file.path(dir_results_intermediate_sccore_gep)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

outcome_var <- "Metastasis_num"
time_var    <- "FU_metastasis_years"
weights_var <- "Weight_rescaled"
tp          <- 5

BWH_LEVELS <- c("T1", "T2a")

## "full"   recalibrate on the whole discovery cohort, then restrict to T1-T2a (mirrors validation)
## "cohort" recalibrate on T1-T2a only
CALIBRATE_ON <- "full"

#-------------------------------------------------------------------------------

# 3. Functions
#-------------------------------------------------------------------------------

#' Recalibrate predictions with the Cox method
#'
#' Fits a weighted Cox model of the outcome on log(-log(predictions)) and
#' returns the recalibrated event probabilities at tp.
#' Unchanged from L. Pozza / B. Rentroia Pacheco.
#'
#' @param predictions numeric vector. Predicted event probabilities.
#' @param y data.frame. Outcome, follow-up and weights.
#' @param outcome_n character. Outcome column in y.
#' @param fup_n character. Follow-up time column in y.
#' @param tp numeric. Time point (years).
#' @param weights_n character. Weights column in y.
#'
#' @return list with uncalibrated_predictions, recalibrated_predictions,
#'   recalibration_cox_lp, recalibration_model and baseline_surv.
recalib_preds <- function(predictions, y, outcome_n, fup_n, tp, weights_n){
  
  predictions <- ifelse(1 - predictions < 1e-16, 1 - 1e-16, predictions)
  predictions <- ifelse(predictions < 1e-16, 1e-16, predictions)
  
  loglog_predictions <- log(-log(predictions))
  
  recalib_model <- coxph(Surv(y[, fup_n], y[, outcome_n]) ~ loglog_predictions,
                         weights = y[, weights_n])
  
  lp <- predict(recalib_model, type = "lp")
  baseline_surv <- summary(survfit(recalib_model), times = tp)$surv
  
  recalib_predictions <- as.vector(1 - baseline_surv^exp(lp))
  
  list(uncalibrated_predictions = predictions,
       recalibrated_predictions = recalib_predictions,
       recalibration_cox_lp = lp,
       recalibration_model = recalib_model,
       baseline_surv = baseline_surv)
}

#' Apply a fitted calibration model to raw scores
#'
#' predict() centres on fit$means, which matches the baseline survival from
#' survfit(). Do not replace with a hand-written beta * x.
#'
#' @param scores numeric vector. Raw predicted event probabilities.
#' @param calibration_model coxph. Calibration model from recalib_preds().
#' @param baseline_surv numeric. Baseline survival at the time point.
#'
#' @return list with loglog (transformed scores), lp and risk (calibrated
#'   event probabilities).
applyCalibration <- function(scores, calibration_model, baseline_surv) {
  
  scores_clipped <- pmin(pmax(scores, 1e-16), 1 - 1e-16)
  newdat <- data.frame(loglog_predictions = log(-log(scores_clipped)))
  lp <- as.vector(predict(calibration_model, newdata = newdat, type = "lp"))
  
  list(loglog = newdat$loglog_predictions,
       lp = lp,
       risk = as.vector(1 - baseline_surv^exp(lp)))
}

#' (Weighted) threshold-based metrics for a binary classification
#'
#' Samples with prediction > cutoff are classified as high risk.
#' From L. Pozza, unchanged.
#'
#' @param cutoff numeric. Threshold on the predictions.
#' @param outcome_num numeric vector. Outcome (1 = case, 0 = control).
#' @param predictions numeric vector. Predicted risks.
#' @param weights numeric vector or NULL. Sample weights (NULL = unweighted).
#'
#' @return Named numeric vector: weighted FP, FN, TP, TN, SE, SP, PPV, NPV,
#'   likelihood ratios, proportions/numbers at low and high risk, positives
#'   and negatives.
computeDiscMetricsBinary <- function(cutoff, outcome_num, predictions, weights = NULL){
  
  labels_cutoff <- ifelse(predictions > cutoff, 1, 0)
  if (is.null(weights)) weights <- rep(1, length(predictions))
  
  correctly_assigned <- which(outcome_num == labels_cutoff)
  wrongly_assigned   <- which(outcome_num != labels_cutoff)
  positives <- which(outcome_num == 1)
  negatives <- which(outcome_num == 0)
  TN <- sum(weights[intersect(correctly_assigned, negatives)])
  TP <- sum(weights[intersect(correctly_assigned, positives)])
  FN <- sum(weights[intersect(wrongly_assigned, positives)])
  FP <- sum(weights[intersect(wrongly_assigned, negatives)])
  Pos <- sum(weights[positives])
  Neg <- sum(weights[negatives])
  
  SE  <- if ((TP + FN) > 0) TP / (TP + FN) else NA
  SP  <- if ((TN + FP) > 0) TN / (TN + FP) else NA
  PPV <- if ((TP + FP) > 0) TP / (TP + FP) else NA
  NPV <- if ((TN + FN) > 0) TN / (TN + FN) else NA
  
  LR_pos <- NA
  LR_neg <- NA
  if (!is.na(SE) & !is.na(SP)){
    LR_pos <- SE / (1 - SP)
    LR_neg <- (1 - SE) / SP
  }
  
  low_risk_n  <- FN + TN
  high_risk_n <- FP + TP
  low_risk_prop  <- low_risk_n / (low_risk_n + high_risk_n)
  high_risk_prop <- 1 - low_risk_prop
  
  c("FP_w" = FP, "FN_w" = FN, "TP_w" = TP, "TN_w" = TN,
    "SE_w" = SE, "SP_w" = SP, "PPV_w" = PPV, "NPV_w" = NPV,
    "LR_pos_w" = LR_pos, "LR_neg_w" = LR_neg,
    "low_risk_w" = low_risk_prop, "high_risk_w" = high_risk_prop,
    "positives_w" = Pos, "negatives_w" = Neg,
    "n_low_risk_w" = low_risk_n, "n_high_risk_w" = high_risk_n)
}

#' Named numeric vector of predictions
#'
#' @param x Named numeric vector, or data.frame with an ID and a prediction column.
#' @param id_col character. ID column (data.frame input).
#' @param pred_col character or NULL. Prediction column (data.frame input).
#'
#' @return Numeric vector of predictions named by sample ID.
as_named_preds <- function(x, id_col = "SkylineDx.ID", pred_col = NULL) {
  if (is.data.frame(x)) {
    if (is.null(pred_col)) pred_col <- setdiff(colnames(x), id_col)[1]
    out <- setNames(as.numeric(x[[pred_col]]), as.character(x[[id_col]]))
  } else {
    out <- setNames(as.numeric(x), names(x))
  }
  out[!is.na(names(out))]
}

#' Recode the outcome (Case/Control, 1/0, Yes/No) to 1/0
#'
#' @param x vector. Outcome values.
#'
#' @return Numeric vector (1 = metastasis, 0 = no metastasis).
norm_metastasis <- function(x) {
  x <- trimws(as.character(x))
  ifelse(x %in% c("Case", "1", "TRUE", "Yes"), 1L,
         ifelse(x %in% c("Control", "0", "FALSE", "No"), 0L, NA_integer_))
}

#' Read the NCC sampling weights
#'
#' @param path character. Path to the weights file.
#' @param id_col character. Sample ID column.
#' @param weight_col character. Weight column.
#'
#' @return data.frame with SkylineDx.ID, the weight and (if available) the
#'   outcome in the file, used as a cross-check.
read_weights <- function(path, id_col, weight_col) {
  w <- read.csv(path, stringsAsFactors = FALSE)
  
  missing_cols <- setdiff(c(id_col, weight_col), names(w))
  if (length(missing_cols) > 0) {
    stop("Weights file lacks column(s): ", paste(missing_cols, collapse = ", "))
  }
  
  w[[id_col]] <- trimws(as.character(w[[id_col]]))
  dup <- unique(w[[id_col]][duplicated(w[[id_col]])])
  if (length(dup) > 0) {
    stop(length(dup), " duplicated IDs in the weights file, e.g. ",
         paste(head(dup), collapse = ", "))
  }
  
  out <- data.frame(SkylineDx.ID = w[[id_col]], stringsAsFactors = FALSE)
  out[[weight_col]] <- as.numeric(w[[weight_col]])
  if ("Metastasis" %in% names(w)) out$Metastasis_weightsfile <- norm_metastasis(w$Metastasis)
  
  cat("Weights file:", nrow(out), "samples |",
      weight_col, "range", signif(range(out[[weight_col]], na.rm = TRUE), 4), "\n")
  out
}

#' One SCCore-GEP score per discovery patient
#'
#' @param modelGEP list. Saved SCCore-GEP model object.
#' @param source character. "OOB" (mean out-of-bag prediction per patient)
#'   or "Model_fit" (in-sample predictions of the final model).
#' @param slot_final character. Slot of the final model in modelGEP.
#'
#' @return Numeric vector of scores named by SkylineDx.ID.
getDiscoveryScores <- function(modelGEP, source, slot_final) {
  
  if (source == "Model_fit") {
    out <- as_named_preds(modelGEP[[slot_final]]$predictions)
    cat("  source: ", slot_final, " (in-sample), ", length(out), " patients\n", sep = "")
    return(out)
  }
  
  if (!("OOB_predicted_probs" %in% names(modelGEP))) {
    stop("PREDICTION_SOURCE is \"OOB\" but the model object has no OOB_predicted_probs.")
  }
  
  oob <- modelGEP$OOB_predicted_probs
  oob$SkylineDx.ID <- trimws(as.character(oob$SkylineDx.ID))
  
  avg <- oob %>%
    group_by(SkylineDx.ID) %>%
    summarise(score = mean(predicted_probs), n_times_oob = n(), .groups = "drop")
  
  cat("  source: out-of-bag averages, ", nrow(avg), " patients\n", sep = "")
  cat("  times held out per patient: median ", median(avg$n_times_oob),
      " (range ", min(avg$n_times_oob), "-", max(avg$n_times_oob), ")\n", sep = "")
  
  n_model_fit <- length(as_named_preds(modelGEP[[slot_final]]$predictions))
  if (nrow(avg) < n_model_fit) {
    warning(n_model_fit - nrow(avg), " patients have no out-of-bag score (never held out) and ",
            "drop out of this analysis.")
  }
  
  setNames(avg$score, avg$SkylineDx.ID)
}

#' Combine clinical data, scores and weights into the analysis data
#'
#' @param clinical data.frame. Discovery clinical data.
#' @param preds_named numeric vector. Scores named by SkylineDx.ID.
#' @param weights_df data.frame. Output of read_weights().
#'
#' @return data.frame of complete cases with outcome, follow-up, weights,
#'   BWH stage and score.
build_discovery_data <- function(clinical, preds_named, weights_df) {
  
  clinical$SkylineDx.ID <- trimws(as.character(clinical$SkylineDx.ID))
  
  if (weights_var %in% names(clinical)) {
    message("Clinical data already has '", weights_var, "'; it is replaced by the weights file.")
    clinical[[weights_var]] <- NULL
  }
  
  clinical[[outcome_var]] <- norm_metastasis(clinical$Metastasis)
  if (all(is.na(clinical[[outcome_var]]))) {
    stop("clinical$Metastasis could not be coded to 0/1 (values: ",
         paste(head(unique(clinical$Metastasis)), collapse = ", "), ").")
  }
  
  df <- clinical %>%
    filter(SkylineDx.ID %in% names(preds_named)) %>%
    left_join(tibble(SkylineDx.ID = names(preds_named),
                     sccore_gep = as.numeric(preds_named)),
              by = "SkylineDx.ID") %>%
    left_join(weights_df, by = "SkylineDx.ID") %>%
    as.data.frame()
  
  cat("  unique scored IDs:", length(unique(names(preds_named))), "\n")
  cat("  matched clinical samples:", nrow(df), "\n")
  cat("  with a weight:", sum(!is.na(df[[weights_var]])), "\n")
  
  if (sum(!is.na(df[[weights_var]])) == 0) {
    stop("No clinical sample matched the weights file. Check WEIGHTS_ID_COL ('",
         WEIGHTS_ID_COL, "') against clinical$SkylineDx.ID.")
  }
  
  ## Outcome in the weights file must agree with the clinical outcome
  if ("Metastasis_weightsfile" %in% names(df)) {
    both <- !is.na(df[[outcome_var]]) & !is.na(df$Metastasis_weightsfile)
    n_mismatch <- sum(df[[outcome_var]][both] != df$Metastasis_weightsfile[both])
    cat("  outcome mismatches clinical vs weights file:", n_mismatch, "\n")
    if (n_mismatch > 0) {
      warning(n_mismatch, " patients have a different metastasis status in the weights file ",
              "than in the clinical data. The weights may belong to a different sampling.")
    }
    df$Metastasis_weightsfile <- NULL
  }
  
  cat("  BWH distribution:\n")
  print(table(df$BWH, useNA = "ifany"))
  
  n_before <- nrow(df)
  df <- df %>%
    filter(!is.na(sccore_gep),
           !is.na(.data[[outcome_var]]),
           !is.na(.data[[time_var]]),
           !is.na(.data[[weights_var]]))
  
  cat("\nComplete rows for the calibration:", nrow(df),
      "(dropped", n_before - nrow(df), ")\n")
  cat("  raw score range:", signif(range(df$sccore_gep), 4), "\n")
  cat("  follow-up: median", round(median(df[[time_var]]), 2),
      "max", round(max(df[[time_var]]), 2), "\n")
  
  if (max(df[[time_var]]) < tp) {
    stop("Follow-up does not reach ", tp, " years; the baseline survival cannot be taken there.")
  }
  
  ## Orientation: a score running the wrong way would silently invert the risks
  cor_raw <- suppressWarnings(cor(df$sccore_gep, df[[outcome_var]]))
  cat("  correlation of the raw score with metastasis:", signif(cor_raw, 3), "\n")
  if (cor_raw < 0) {
    warning("The raw score is NEGATIVELY correlated with metastasis. If sccore_gep is a survival ",
            "probability rather than a risk, the calibrated risks will be inverted.")
  }
  
  df
}

#' Weighted metrics at one cutoff of the calibrated risk (strict >)
#'
#' @param cutoff numeric. Threshold on the calibrated 5-year risk.
#' @param df data.frame. Analysis data with sccore_gep_cox_calibrated.
#'
#' @return One-row data.frame with the threshold, numbers flagged and the
#'   metrics of computeDiscMetricsBinary().
metrics_row <- function(cutoff, df) {
  r <- df$sccore_gep_cox_calibrated
  y <- df[[outcome_var]]
  w <- df[[weights_var]]
  flagged <- r > cutoff
  
  m <- computeDiscMetricsBinary(cutoff, y, r, w)
  
  data.frame(threshold = cutoff,
             threshold_pct = 100 * cutoff,
             n_flagged = sum(flagged),
             n_flagged_cases = sum(flagged & y == 1),
             n_flagged_controls = sum(flagged & y == 0),
             as.list(m),
             check.names = FALSE)
}

#' Run sections 5-9 for one prediction source
#'
#' Builds the discovery data, fits the calibration, computes calibrated
#' 5-year risks, restricts to BWH T1-T2a and computes metrics at every
#' threshold.
#'
#' @param source character. Prediction source ("OOB" or "Model_fit").
#' @param modelGEP list. Saved SCCore-GEP model object.
#' @param clinical data.frame. Discovery clinical data.
#' @param weights_df data.frame. Output of read_weights().
#'
#' @return list with source, df_all, df_metrics, thresh_based_metrics and
#'   the calibration (model, baseline survival, beta, x_bar).
run_pipeline <- function(source, modelGEP, clinical, weights_df) {
  

  cat("PREDICTION SOURCE:", source, "\n")
  # 3.1. Build the discovery dataset
  #-----------------------------------------------------------------------------
  preds_named <- getDiscoveryScores(modelGEP, source, slot.final.model)
  df_all <- build_discovery_data(clinical, preds_named, weights_df)
  
  # 3.2. Fit the calibration on discovery
  #-----------------------------------------------------------------------------
  df_fit <- if (CALIBRATE_ON == "cohort") df_all %>% filter(BWH %in% BWH_LEVELS) else df_all
  
  cat("\nFitting the calibration on", CALIBRATE_ON, "-", nrow(df_fit), "patients\n")
  
  recalib_res <- recalib_preds(df_fit$sccore_gep, df_fit,
                               outcome_var, time_var, tp, weights_var)
  
  calibration_model <- recalib_res$recalibration_model
  baseline_surv_5yr <- recalib_res$baseline_surv
  beta_cal <- coef(calibration_model)[["loglog_predictions"]]
  x_bar    <- calibration_model$means[["loglog_predictions"]]
  
  cat("  beta (calibration slope) :", signif(beta_cal, 6), "\n")
  cat("  baseline survival at", tp, "yr:", signif(baseline_surv_5yr, 6), "\n")
  cat("  centring value fit$means :", signif(x_bar, 6), "\n")
  cat("  |beta| =", signif(abs(beta_cal), 4),
      if (abs(beta_cal) > 1.5) " - the raw predictions were substantially too compressed"
      else if (abs(beta_cal) < 0.7) " - the raw predictions were substantially too spread out"
      else " - close to 1, the raw spread was about right", "\n")
  
  # 3.3. Calibrated 5-year risk for every discovery patient
  #-----------------------------------------------------------------------------
  cal_all <- applyCalibration(df_all$sccore_gep, calibration_model, baseline_surv_5yr)
  df_all$sccore_gep_cox_calibrated <- cal_all$risk
  df_all$recalibration_lp <- cal_all$lp
  
  ## Self-check: the mapping must reproduce itself on the patients it was fitted on
  i.fit <- match(df_fit$SkylineDx.ID, df_all$SkylineDx.ID)
  max_diff <- max(abs(df_all$sccore_gep_cox_calibrated[i.fit] -
                        recalib_res$recalibrated_predictions), na.rm = TRUE)
  
  cat("\nSelf-check, calibration reproduced through applyCalibration()\n")
  cat("  max absolute difference:", signif(max_diff, 3), "\n")
  if (max_diff > 1e-10) {
    stop("applyCalibration() does not reproduce recalib_preds() (max difference ",
         signif(max_diff, 3), "). Do not continue.")
  }
  cat("  PASSED\n")
  
  cat("\nCalibrated 5-year risk, all discovery patients:\n")
  print(signif(quantile(df_all$sccore_gep_cox_calibrated, c(0, .25, .5, .75, 1)), 4))
  
  w_all <- df_all[[weights_var]]
  cat("Weighted observed metastasis rate:",
      signif(sum(w_all * df_all[[outcome_var]]) / sum(w_all), 5),
      "| weighted mean predicted risk:",
      signif(sum(w_all * df_all$sccore_gep_cox_calibrated) / sum(w_all), 5), "\n")
  
  # 3.4. Restrict to T1-T2a
  #-----------------------------------------------------------------------------
  df_metrics <- df_all %>% filter(BWH %in% BWH_LEVELS)
  
  cat("\nDiscovery T1-T2a\n")
  cat("  n:", nrow(df_metrics),
      "| cases:", sum(df_metrics[[outcome_var]] == 1),
      "| controls:", sum(df_metrics[[outcome_var]] == 0), "\n")
  cat("  weighted prevalence:",
      signif(sum(df_metrics[[weights_var]] * df_metrics[[outcome_var]]) /
               sum(df_metrics[[weights_var]]), 5), "\n")
  cat("  calibrated risk range:", signif(range(df_metrics$sccore_gep_cox_calibrated), 4), "\n")
  
  if (nrow(df_metrics) == 0) stop("No T1-T2a patients left after filtering.")
  
  # 3.5. Metrics at every threshold
  #-----------------------------------------------------------------------------
  # Every observed calibrated risk is a candidate cut. Strict >, so a threshold equal to a
  # patient's risk puts that patient BELOW the cut.
  thresholds <- sort(unique(df_metrics$sccore_gep_cox_calibrated))
  cat("\nCandidate thresholds:", length(thresholds), "\n")
  
  thresh_based_metrics <- bind_rows(lapply(thresholds, metrics_row, df = df_metrics))

  list(source = source,
       df_all = df_all,
       df_metrics = df_metrics,
       thresh_based_metrics = thresh_based_metrics,
       
       calibration_model = calibration_model,
       baseline_surv = baseline_surv_5yr,
       beta = beta_cal,
       x_bar = x_bar)
}

#' Write the outputs of one run
#'
#' @param res list. Output of run_pipeline().
#'
#' @return Called for its side effect: writes the threshold metrics (csv).
write_outputs <- function(res) {
  sfx <- paste0("discoveryCal_", res$source)
  
  write.csv(res$thresh_based_metrics,
            file.path(output_dir, paste0("discovery_T1_T2a_threshold_metrics_", sfx, ".csv")),
            row.names = FALSE)
  
  
}


#-------------------------------------------------------------------------------

# 4. Read the discovery model, clinical data and weights
#-------------------------------------------------------------------------------

modelGEP <- readRDS(path.model.Data)
stopifnot(slot.final.model %in% names(modelGEP))

#clinical <- attr(get(load(path.clinical.Data)), "clinical")
clinical <- read_csv(path.clinical.Data)
weights_df <- read_weights(path.weights.Data, WEIGHTS_ID_COL, weights_var)


res_main <- run_pipeline(PREDICTION_SOURCE, modelGEP, clinical, weights_df)
write_outputs(res_main)


#-------------------------------------------------------------------------------
# 5. PRECISION-RECALL CURVES: T1-T2a AND ENTIRE DEVELOPMENT COHORT
#-------------------------------------------------------------------------------
OP_THRESHOLD <- 0.0698   # selected operating point (5-year absolute risk)

#' Precision-recall data for one patient group
#'
#' @param df data.frame. Patients with sccore_gep_cox_calibrated.
#' @param label character. Name of the group for the plot.
#'
#' @return list with the curve, the operating point and the weighted prevalence.
pr_group <- function(df, label) {
  thr <- sort(unique(df$sccore_gep_cox_calibrated))
  curve <- bind_rows(lapply(thr, metrics_row, df = df)) %>%
    filter(!is.na(PPV_w)) %>%
    arrange(desc(threshold)) %>%
    mutate(group = label)
  op <- metrics_row(OP_THRESHOLD, df) %>% mutate(group = label)
  prev <- sum(df[[weights_var]] * df[[outcome_var]]) / sum(df[[weights_var]])
  list(curve = curve, op = op,
       prev = data.frame(group = label, prev = prev))
}

pr_t1t2a <- pr_group(res_main$df_metrics, "Discovery BWH T1-T2a")
pr_all   <- pr_group(res_main$df_all,     "Entire development cohort")

curves <- bind_rows(pr_t1t2a$curve, pr_all$curve)
ops    <- bind_rows(pr_t1t2a$op,    pr_all$op)
prevs  <- bind_rows(pr_t1t2a$prev,  pr_all$prev)

p_pr <- ggplot(curves, aes(x = SE_w, y = PPV_w, colour = group)) +
  geom_path() +
  geom_hline(data = prevs, aes(yintercept = prev, colour = group), linetype = "dashed") +
  geom_point(data = ops, size = 3, shape = 18) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_colour_manual(values = c("Discovery BWH T1-T2a" = "steelblue",
                                 "Entire development cohort" = "black")) +
  labs(x = "Weighted sensitivity (recall)", y = "Weighted PPV (precision)", colour = NULL,
       title = "Precision-recall curves: discovery T1-T2a and entire development cohort",
       subtitle = paste0("Calibrated averaged OOB scores (", res_main$source,
                         "); diamond = operating point (", 100 * OP_THRESHOLD,
                         "% 5-year risk); dashed lines = weighted prevalence")) +
  theme_bw() + theme(legend.position = "bottom")

ggsave(file.path(output_dir, paste0("precision_recall_T1T2a_and_entire_cohort_discoveryCal_",
                                    res_main$source, ".pdf")),
       p_pr, width = 7, height = 6)
p_pr

cat("\nOperating point (", OP_THRESHOLD, ") per group:\n")
print(ops %>% select(group, PPV_w, SE_w, n_flagged, n_flagged_cases, n_flagged_controls) %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))))
print(prevs %>% mutate(prev = signif(prev, 3)))
#-------------------------------------------------------------------------------