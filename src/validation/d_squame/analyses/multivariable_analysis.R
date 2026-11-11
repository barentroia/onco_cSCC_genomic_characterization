#===============================================================================
# Aim:      Multivariable Cox analyses in the D-SQUAME validation dataset:
#           SCCore-GEP (recalibrated, logit-transformed) combined with BWH,
#           AJCC8, the EMC model or the DvP score. Log HRs are estimated per
#           imputed dataset with bootstrap SEs and pooled with Rubin's rules.
# Author:   Adapted from code of L. Pozza
# Input:    - D-SQUAME validation clinical data (e3_clinical), multiple
#             imputation object and NCC weights (e7_weights)
#           - Recalibrated SCCore-GEP predictions (output of recalibration.R)
#           - D_SQUAME_validation_data_set.rds and DvP signature (p7_bailey_xsls)
# Output:   Written to results_dir/intermediate/validation/d_squame:
#           - multivariable_analysis_<subset>_SCCoreGEP_vs_<BWH|AJCC8|EMCmodel|DvPscore>.csv
#           - DvP scores per sample (for reference)
# Requires: code_dir, results_dir, output_data_dir (main.R); subset_oi
#           ("t1t2a" = BWH T1-T2a, "t1t2" = AJCC8 T1-T2, otherwise entire
#           dataset, analysed with paired bootstrap); config.R paths
# Usage:    Sourced from main.R (once per subset)
#===============================================================================

# 1. Load libraries
#-------------------------------------------------------------------------------
library(conflicted)
library(openxlsx)
library(tidyverse)
library(survival)
library(mice)
library(patchwork)
# genefu (Bioconductor) is called via genefu:: to avoid masking dplyr functions
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::rename)
conflicts_prefer(base::cbind)
conflicts_prefer(base::rbind)

#-------------------------------------------------------------------------------

# 2. Functions
#-------------------------------------------------------------------------------
#' Pool log HRs over imputed datasets with Rubin's rules
#'
#' Uses Barnard-Rubin adjusted degrees of freedom when n and k are given.
#'
#' @param est_metric numeric vector. Log HR estimate per imputed dataset.
#' @param est_se numeric vector. SE of the log HR per imputed dataset.
#' @param n integer or NULL. Number of events.
#' @param k integer or NULL. Number of parameters in the model.
#'
#' @return Named numeric vector: "95% Low", "Metric" (pooled log HR),
#'   "95% Up" and "p-value".
pool_logHR_rubin <- function(est_metric, est_se, n = NULL, k = NULL){
  m <- length(est_metric)
  # Pooled estimate
  mean_metric <- mean(est_metric)
  # Within-imputation variance
  w_metric <- mean(est_se^2)
  # Betweeen-imputation variance (0 for a single dataset)
  b_metric <- if (m > 1) var(est_metric) else 0
  # Total variance
  tv_metric <- w_metric + (1 + (1/m)) * b_metric
  # Pooled standard error
  se_total <- sqrt(tv_metric)
  # Fraction of missing information (FMI)
  lambda <- (b_metric + (b_metric / m)) / tv_metric
  # Rubin's degrees of freedom (infinite when there is no between-imputation variance)
  v_rubin <- if (lambda > 0) (m - 1) / (lambda ^ 2) else Inf
  # Bernard-Rubin adjusted degrees of freedom
  if (!(is.null(n) & is.null(k))){
    v_com <- n - k
    v_obs <- (v_com + 1) / (v_com + 3) * v_com * (1 - lambda)
    # for v_rubin -> Inf the Barnard-Rubin df tends to v_obs
    v_bernard_rubin <- if (is.finite(v_rubin)) (v_rubin * v_obs) / (v_rubin + v_obs) else v_obs
    v <- v_bernard_rubin
  } else {
    v <- v_rubin
  }
  if (b_metric == 0 && m > 1){
    message("  pool_logHR_rubin: estimates identical across the ", m,
            " imputed datasets (no imputed variable in the model) -> single-dataset result")
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

# 3. Filenames and settings
#-------------------------------------------------------------------------------
source(file.path(code_dir, "functions", "boot_metrics_se_ci.R"))
source(file.path(code_dir, "functions", "sig.score.R"))
# Clinical data
clinical_fn <- e3_clinical
clinical_mimp_fn <- file.path(results_dir, "intermediate", "clinical", "clinical_data_imputed_d_squame_validation.Rdata")
# Weights
lr_weights_fn <- e7_weights
# Recalibrated SCCore-GEP predictions
sccoregep_fn <- file.path(results_dir, "intermediate", "validation", "d_squame", "recalibrated_SCCoreGEP_predictions.csv")
# Expression data validation cohort
output_data_dir<- file.path(results_dir, "data_rnaseq")
validation_data_fn <- file.path(output_data_dir, "D_SQUAME_validation_data_set.rds")  
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

# 4. Source functions
#-------------------------------------------------------------------------------
source(file.path(code_dir, "functions", "colors.R"))
#-------------------------------------------------------------------------------

# 5. Read data
#-------------------------------------------------------------------------------
# Clinical data
clinical <- read.csv(clinical_fn)
load(clinical_mimp_fn)
# Weights
lr_weights <- read.csv(lr_weights_fn)
# Recalibrated SCCore-GEP predictions
sccoregep <- read.csv(sccoregep_fn)
# Expression data
validation_data <- readRDS(validation_data_fn)
# DvP signature genes
dvp_genes <- readxl::read_xlsx(p7_bailey_xsls, sheet = "DvP_Signature")
#-------------------------------------------------------------------------------

# 6. DvP score
#-------------------------------------------------------------------------------
gtf <- validation_data@out
# Signature gene symbols that differ from the GTF gene names
replacement_gene_ids <- list(LOR = "ENSG00000203782.6",
                             GRAMD3 = "ENSG00000250602.6",
                             SEPP1 = "ENSG00000286271.2",
                             NOV = "ENSG00000136999.5",
                             FAM84A = "ENSG00000162981.14",
                             HIST1H2BC = "ENSG00000180596.7",
                             HIST1H2BD = "ENSG00000158373.8",
                             WDR66 = "ENSG00000158023.10",
                             KIAA0101 = "ENSG00000166803.14",
                             GSG2 = "ENSG00000177602.5",
                             DFNA5 = "ENSG00000105928.16",
                             FAM64A = "ENSG00000129195.16",
                             KIAA1524 = "ENSG00000163507.14")
dvp_genes$Gene_ID <- gtf[match(dvp_genes$Gene_symbol, gtf$gene_name), ]$gene_id
dvp_genes[match(names(replacement_gene_ids), dvp_genes$Gene_symbol), ]$Gene_ID <- unlist(replacement_gene_ids)
stopifnot(!anyNA(dvp_genes$Gene_ID))

log2tpms <- log2(validation_data@sequencing_data$tpms + 1)
# all signature genes must be present (also catches Ensembl version mismatches)
stopifnot(all(dvp_genes$Gene_ID %in% rownames(log2tpms)))

# Same scoring as in discovery: z-scored log2(TPM+1), continuous coefficients
zlog2tpms <- as.data.frame(scale(log2tpms, center = TRUE, scale = TRUE))
zlog2tpms_sel <- t(zlog2tpms[match(dvp_genes$Gene_ID, rownames(zlog2tpms)), ])

dvp_genes_prepared <- dvp_genes %>%
  rename(EntrezGene.ID = Gene_ID,
         coefficient = Coefficient) %>%
  mutate(probe = EntrezGene.ID) %>%
  select(probe, EntrezGene.ID, coefficient)

annot_prepared <- dvp_genes_prepared %>%
  select(EntrezGene.ID) %>%
  mutate(row.names = EntrezGene.ID)

sig_score_con <- sig.score(x = dvp_genes_prepared,
                                   data = zlog2tpms_sel,
                                   annot = annot_prepared,
                                   do.mapping = FALSE,
                                   signed = FALSE,
                                   verbose = TRUE)

# Sample IDs of the expression data are assumed to be Skyline IDs
dvp_scores <- data.frame(Skyline_ID = rownames(zlog2tpms_sel),
                         DvP_score = as.numeric(sig_score_con$score)*-1)

#-------------------------------------------------------------------------------

# 7. Preprocess data
#-------------------------------------------------------------------------------
# We use calibrated, logit-transformed SCCore-GEP probabilities
sccoregep <- sccoregep %>%
  left_join(clinical %>%
              mutate(Metastasis_num = ifelse(Metastasis == "Case", 1, 0)) %>%
              select(Skyline_ID, Metastasis_num, Set_id),
            by = "Skyline_ID") %>%
  left_join(lr_weights %>% select(Patient_ID_SKY, Weight_rescaled) %>% unique(), by = "Patient_ID_SKY") %>%
  left_join(dvp_scores, by = "Skyline_ID") %>%
  mutate(sccore_gep_cox_calibrated = qlogis(sccore_gep_cox_calibrated)) %>%
  select(Patient_ID_SKY, sccore_gep_cox_calibrated, DvP_score, Metastasis_num, Set_id, Weight_rescaled)

# every sample with a SCCore-GEP prediction should have a DvP score
stopifnot(!anyNA(sccoregep$DvP_score))

#-------------------------------------------------------------------------------

# 8. Multivariable analysis: SCCore-GEP vs BWH
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

# 9. Multivariable analysis: SCCore-GEP vs AJCC8
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

# 10. Multivariable analysis: EMC model
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

# 11. Multivariable analysis: SCCore-GEP vs DvP score
#-------------------------------------------------------------------------------
# Compute point estimate and bootstrap SE of log HR
res_cox_sccoregep_dvp <- lapply(1 : mimp$m,
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
                                  # Standardise DvP score (mean 0, SD 1) within this imputed dataset -> HR per 1 SD
                                  df <- df %>%
                                    mutate(DvP_score = as.numeric(scale(DvP_score))) %>%
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
                                strat_vars = c("sccore_gep_cox_calibrated", "DvP_score"),
                                ids_var = "Patient_ID_SKY",
                                weight_var = "Weight_rescaled")
boot_params_dvp <- res_cox_sccoregep_dvp %>%
  bind_rows()
# Pool results using Rubin's rules
pooled_params_dvp <- boot_params_dvp %>%
  mutate(n_preds = length(unique(term)),
         avg_events = round(mean(n_events))) %>%
  group_by(term) %>%
  dplyr::summarize(res = list(pool_logHR_rubin(log_hr, se, unique(avg_events), unique(n_preds)))) %>%
  unnest_wider(res) %>%
  as.data.frame() %>%
  mutate(term = recode(term,
                       "DvP_score" = "DvP score",
                       "sccore_gep_cox_calibrated" = "SCCore-GEP"),
         hr_ci95_low = round(exp(`95% Low`), 4),
         hr = round(exp(Metric), 4),
         hr_ci95_up = round(exp(`95% Up`), 4),
         pvalue = `p-value`) %>%
  select(term, hr, hr_ci95_low, hr_ci95_up, pvalue)
#-------------------------------------------------------------------------------

# 12. Save results
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
write.csv(pooled_params_dvp,
          file.path(output_dir, paste0("multivariable_analysis_", ifelse(is.null(subset_oi), "", paste0(subset_oi, "_")), "SCCoreGEP_vs_DvPscore.csv")),
          row.names = F)
# DvP scores, for reference
write.csv(dvp_scores,
          file.path(output_dir, "DvP_scores_d_squame_validation.csv"),
          row.names = F)
#-------------------------------------------------------------------------------


