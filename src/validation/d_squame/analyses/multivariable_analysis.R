#-------------------------------------------------------------------------------
# Aim: Peform multivariable analyses (BWH, AJCC8, EMC model) on D-SQUAME validation dataset
# Author: L.Pozza
# Input: Recalibrated SCCore-GEP predictions + clinical data + weights of D-SQUAME validation
# Output: Multivariable analyses results
#-------------------------------------------------------------------------------

# Functions
#-------------------------------------------------------------------------------
pool_logHR_rubin <- function(est_metric, est_se, n = NULL, k = NULL){
  #' Pool log HR and its standard error, and compute pooled p-value with Rubin's rule
  #'
  #' @param est_metric numeric vector. Vector with log HR estimates for each imputed datatset
  #' @param est_se numeric vector. Vector with SE of log HR for each imputed datatset
  #' @param n integer or NULL. Number of events
  #' @param k integer or NULL. Number of parameters to fit
  
  m <- length(est_metric)
  # Pooled estimate
  mean_metric <- mean(est_metric)
  # Within-imputation variance
  w_metric <- mean(est_se^2)
  # Betweeen-imputation variance
  b_metric <- var(est_metric)
  # Total variance
  tv_metric <- w_metric + (1 + (1/m)) * b_metric
  # Pooled standard error
  se_total <- sqrt(tv_metric)
  # Fraction of missing information (FMI)
  lambda <- (b_metric + (b_metric / m)) / tv_metric
  # Rubin's degrees of freedom
  v_rubin <- (m - 1) / (lambda ^ 2)
  # Bernard-Rubin adjusted degrees of freedom
  if (!(is.null(n) & is.null(k))){
    v_com <- n - k
    v_obs <- (v_com + 1) / (v_com + 3) * v_com * (1 - lambda)
    v_bernard_rubin <- (v_rubin * v_obs) / (v_rubin + v_obs)
    v <- v_bernard_rubin
  } else {
    v <- v_rubin
  }
  # t-quantile
  t <- qt(0.975, v)
  # 95% CI on log-HR scale
  metrics_l <- mean_metric - t * se_total
  metrics_u <- mean_metric + t * se_total
  # p-value
  pvalue <- pt(q = abs(mean_metric / se_total), df = v, lower.tail = FALSE) * 2
  metrics_res <- c(metrics_l, mean_metric, metrics_u, pvalue)
  names(metrics_res) <- c("95% Low", "Metric", "95% Up", "p-value")
  metrics_res
}
#-------------------------------------------------------------------------------

# Load libraries
#-------------------------------------------------------------------------------
library(conflicted)
library(openxlsx)
library(tidyverse)
library(survival)
library(mice)
library(patchwork)
conflicted::conflicts_prefer(dplyr::filter)
#-------------------------------------------------------------------------------

# Filenames and settings
#-------------------------------------------------------------------------------
source(file.path(code_dir, "functions", "boot_metrics_se_ci.R"))
# Clinical data
clinical_fn <- e3_clinical
clinical_mimp_fn <- file.path(results_dir, "intermediate", "clinical", "clinical_data_imputed_d_squame_validation.Rdata")
# Weights
lr_weights_fn <- e7_weights
# Recalibrated SCCore-GEP predictions
sccoregep_fn <- file.path(results_dir, "intermediate", "validation", "d_squame", "recalibrated_SCCoreGEP_predictions.csv")
# Time point
tp <- 5
# Number of bootstrap repetitions
n_boot <- 100
# Subset of cohort to consider
subset_oi <- ifelse(exists("subset_oi"), subset_oi, "entiredataset")
if (subset_oi == "t1t2a"){
  bwh_subset_to_keep <- c("T1", "T2a")
  ajcc8_subset_to_keep <- NULL
  pair_var <- NULL
} else if (subset_oi == "t1t2") {
  bwh_subset_to_keep <- NULL
  ajcc8_subset_to_keep <- c("T1", "T2")
  pair_var <- NULL
} else {
  bwh_subset_to_keep <- c("T1", "T2a", "T2b", "T3")
  ajcc8_subset_to_keep <- NULL
  pair_var <- "Set_id"
}
# Output folder
output_dir <- file.path(results_dir, "intermediate", "validation", "d_squame")
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}
#-------------------------------------------------------------------------------

# Source functions
#-------------------------------------------------------------------------------
source(file.path(code_dir, "functions", "colors.R"))
#-------------------------------------------------------------------------------

# Read data
#-------------------------------------------------------------------------------
# Clinical data
clinical <- read.csv(clinical_fn)
load(clinical_mimp_fn)
# Weights
lr_weights <- read.csv(lr_weights_fn)
# Recalibrated SCCore-GEP predictions
sccoregep <- read.csv(sccoregep_fn)
#-------------------------------------------------------------------------------

# Preprocess data
#-------------------------------------------------------------------------------
# We use calibrated, logit-transformed SCCore-GEP probabilities
sccoregep <- sccoregep %>%
  left_join(clinical %>%
              mutate(Metastasis_num = ifelse(Metastasis == "Case", 1, 0)) %>%
              select(Skyline_ID, Metastasis_num, Set_id),
            by = "Skyline_ID") %>%
  left_join(lr_weights %>% select(Patient_ID_SKY, Weight_rescaled) %>% unique(), by = "Patient_ID_SKY") %>%
  mutate(sccore_gep_cox_calibrated = qlogis(sccore_gep_cox_calibrated)) %>%
  select(Patient_ID_SKY, sccore_gep_cox_calibrated, Metastasis_num, Set_id, Weight_rescaled)
#-------------------------------------------------------------------------------

# Multivariable analysis: SCCore-GEP vs BWH
#-------------------------------------------------------------------------------
# Compute point estimate and bootstrap SE of log HR
res_cox_sccoregep_bwh <- lapply(1 : mimp$m,
                                function(i, preds_df, status_var, time_var, strat_vars, ids_var, weight_var) {
                                  df <- complete(mimp, i) %>%
                                    mutate(BWH_staging = factor(BWH_staging, levels = c("T1", "T2a", "T2b", "T3"))) %>%
                                    left_join(preds_df, by = ids_var)
                                  if (!is.null(bwh_subset_to_keep)){
                                    df <- df %>%
                                      filter(BWH_staging %in% bwh_subset_to_keep)
                                  } else {
                                    df <- df %>%
                                      filter(AJCC_staging %in% ajcc8_subset_to_keep)              
                                  }
                                  if ("T2b" %in% unique(df$BWH_staging)){
                                    df <- df %>%
                                      mutate(BWH_staging = ifelse(BWH_staging %in% c("T2a", "T2b"), "T2a-T2b", BWH_staging))
                                  }
                                  df <- df %>% droplevels()
                                  boot_res <- boot_hr_se_ci(df, status_var, time_var, weight_var, strat_vars, pair_var, n_boot, "whr")
                                  data.frame(log_hr = boot_res$est,
                                             se = boot_res$se,
                                             n_imp = i,
                                             n_events = df %>% filter(Metastasis == "Case") %>% nrow()) %>%
                                    rownames_to_column(var = "term")},
                                preds_df = sccoregep,
                                status_var = "Metastasis_num",
                                time_var = "FU_metastasis_years",
                                strat_vars = c("sccore_gep_cox_calibrated", "BWH_staging"),
                                ids_var = "Patient_ID_SKY",
                                weight_var = "Weight_rescaled")
boot_params_bwh <- res_cox_sccoregep_bwh %>%
  bind_rows()
# Pool results using Rubin's rules
pooled_params_bwh <- boot_params_bwh %>%
  mutate(n_preds = length(unique(term)),
         avg_events = round(mean(n_events))) %>%
  group_by(term) %>%
  dplyr::summarize(res = list(pool_logHR_rubin(log_hr, se, unique(avg_events), unique(n_preds)))) %>%
  unnest_wider(res) %>%
  as.data.frame() %>%
  mutate(term = recode(term,
                       "BWH_stagingT2a" = "BWH:T2a",
                       "BWH_stagingT2b" = "BWH:T2b",
                       "BWH_stagingT2a-T2b" = "BWH:T2a-T2b",
                       "sccore_gep_cox_calibrated" = "SCCore-GEP"),
         hr_ci95_low = round(exp(`95% Low`), 4),
         hr = round(exp(Metric), 4),
         hr_ci95_up = round(exp(`95% Up`), 4),
         pvalue = `p-value`) %>%
  select(term, hr, hr_ci95_low, hr_ci95_up, pvalue)
#-------------------------------------------------------------------------------

# Multivariable analysis: SCCore-GEP vs AJCC8
#-------------------------------------------------------------------------------
# Compute point estimate and bootstrap SE of log HR
res_cox_sccoregep_ajcc8 <- lapply(1 : mimp$m,
                                  function(i, preds_df, status_var, time_var, strat_vars, ids_var, weight_var) {
                                    df <- complete(mimp, i) %>%
                                      mutate(BWH_staging = factor(BWH_staging, levels = c("T1", "T2a", "T2b", "T3"))) %>%
                                      left_join(preds_df, by = ids_var)
                                    if (!is.null(bwh_subset_to_keep)){
                                      df <- df %>%
                                        filter(BWH_staging %in% bwh_subset_to_keep)
                                    } else {
                                      df <- df %>%
                                        filter(AJCC_staging %in% ajcc8_subset_to_keep)              
                                    }
                                    df <- df %>%
                                      mutate(AJCC_staging = ifelse(AJCC_staging %in% c("T2", "T3"), "T2-T3", AJCC_staging)) %>%
                                      droplevels()
                                    boot_res <- boot_hr_se_ci(df, status_var, time_var, weight_var, strat_vars, pair_var, n_boot, "whr")
                                    data.frame(log_hr = boot_res$est,
                                               se = boot_res$se,
                                               n_imp = i,
                                               n_events = df %>% filter(Metastasis == "Case") %>% nrow()) %>%
                                      rownames_to_column(var = "term")},
                                  preds_df = sccoregep,
                                  status_var = "Metastasis_num",
                                  time_var = "FU_metastasis_years",
                                  strat_vars = c("sccore_gep_cox_calibrated", "AJCC_staging"),
                                  ids_var = "Patient_ID_SKY",
                                  weight_var = "Weight_rescaled")
boot_params_ajcc8 <- res_cox_sccoregep_ajcc8 %>%
  bind_rows()
# Pool results using Rubin's rules
pooled_params_ajcc8 <- boot_params_ajcc8 %>%
  mutate(n_preds = length(unique(term)),
         avg_events = round(mean(n_events))) %>%
  group_by(term) %>%
  dplyr::summarize(res = list(pool_logHR_rubin(log_hr, se, unique(avg_events), unique(n_preds)))) %>%
  unnest_wider(res) %>%
  as.data.frame() %>%
  mutate(term = recode(term,
                       "AJCC_stagingT2" = "AJCC8:T2",
                       "AJCC_stagingT3" = "AJCC8:T3",
                       "AJCC_stagingT2-T3" = "AJCC8:T2-T3",
                       "sccore_gep_cox_calibrated" = "SCCore-GEP"),
         hr_ci95_low = round(exp(`95% Low`), 4),
         hr = round(exp(Metric), 4),
         hr_ci95_up = round(exp(`95% Up`), 4),
         pvalue = `p-value`) %>%
  select(term, hr, hr_ci95_low, hr_ci95_up, pvalue)
#-------------------------------------------------------------------------------

# Multivariable analysis: EMC model
#-------------------------------------------------------------------------------
# Compute point estimate and bootstrap SE of log HR
res_cox_sccoregep_emcmodel <- lapply(1 : mimp$m,
                                     function(i, preds_df, status_var, time_var, strat_vars, ids_var, weight_var) {
                                       df <- complete(mimp, i) %>%
                                         mutate(BWH_staging = factor(BWH_staging, levels = c("T1", "T2a", "T2b", "T3")),
                                                CP_risk = qlogis(CP_risk/100)) %>%
                                         left_join(preds_df, by = ids_var)
                                       if (!is.null(bwh_subset_to_keep)){
                                         df <- df %>%
                                           filter(BWH_staging %in% bwh_subset_to_keep)
                                       } else {
                                         df <- df %>%
                                           filter(AJCC_staging %in% ajcc8_subset_to_keep)              
                                       }
                                       df <- df %>% droplevels()
                                       boot_res <- boot_hr_se_ci(df, status_var, time_var, weight_var, strat_vars, pair_var, n_boot, "whr")
                                       data.frame(log_hr = boot_res$est,
                                                  se = boot_res$se,
                                                  n_imp = i,
                                                  n_events = df %>% filter(Metastasis == "Case") %>% nrow()) %>%
                                         rownames_to_column(var = "term")},
                                     preds_df = sccoregep,
                                     status_var = "Metastasis_num",
                                     time_var = "FU_metastasis_years",
                                     strat_vars = c("sccore_gep_cox_calibrated", "CP_risk"),
                                     ids_var = "Patient_ID_SKY",
                                     weight_var = "Weight_rescaled")
boot_params_emcmodel <- res_cox_sccoregep_emcmodel %>%
  bind_rows()
# Pool results using Rubin's rules
pooled_params_emcmodel <- boot_params_emcmodel %>%
  mutate(n_preds = length(unique(term)),
         avg_events = round(mean(n_events))) %>%
  group_by(term) %>%
  dplyr::summarize(res = list(pool_logHR_rubin(log_hr, se, unique(avg_events), unique(n_preds)))) %>%
  unnest_wider(res) %>%
  as.data.frame() %>%
  mutate(term = recode(term,
                       "CP_risk" = "EMC model",
                       "sccore_gep_cox_calibrated" = "SCCore-GEP"),
         hr_ci95_low = round(exp(`95% Low`), 4),
         hr = round(exp(Metric), 4),
         hr_ci95_up = round(exp(`95% Up`), 4),
         pvalue = `p-value`) %>%
  select(term, hr, hr_ci95_low, hr_ci95_up, pvalue)
#-------------------------------------------------------------------------------

# Save results
#-------------------------------------------------------------------------------
write.csv(pooled_params_bwh,
          file.path(output_dir, paste0("multivariable_analysis_", ifelse(is.null(subset_oi), "", paste0(subset_oi, "_")), "SCCoreGEP_vs_BWH.csv")),
          row.names = F)
write.csv(pooled_params_ajcc8,
          file.path(output_dir, paste0("multivariable_analysis_", ifelse(is.null(subset_oi), "", paste0(subset_oi, "_")), "SCCoreGEP_vs_AJCC8.csv")),
          row.names = F)
write.csv(pooled_params_emcmodel,
          file.path(output_dir, paste0("multivariable_analysis_", ifelse(is.null(subset_oi), "", paste0(subset_oi, "_")), "SCCoreGEP_vs_EMCmodel.csv")),
          row.names = F)
#-------------------------------------------------------------
