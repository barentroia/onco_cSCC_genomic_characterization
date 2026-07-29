#-------------------------------------------------------------------------------
# Aim: Compute SCCore-GEP, BWH, AJCC8, EMC model performances in D-SQUAME validation dataset
# Author: L.Pozza, adapted from J. Traets
# Input: D-SQUAME validation TPM & imputed clinical data + SCCore-GEP model
# Output: Performances of SCCore-GEP, BWH, AJCC8, EMC model in D-SQUAME validation dataset
#-------------------------------------------------------------------------------

# Functions
#-------------------------------------------------------------------------------
compute_perfs <- function(df, preds_df, ids_var, setids_var, status_var, time_var, weight_var, strat_vars, nboot, bwh_subset_to_keep, ajcc8_subset_to_keep){
  #' Compute AUC and C-index stratified by (imputed) variables
  #'
  #' @param df dataframe. Dataframe with imputed values
  #' @param preds_df dataframe. Dataframe with predictions and weights
  #' @param ids_var string. String indicating column with samples IDs
  #' @param setids_var string. String indicating column with set IDs
  #' @param status_var string. String indicating column with status
  #' @param time_var string. String indicating column with time to status
  #' @param weight_var string. String indicating column with weights
  #' @param strat_var string vector. Strings indicating stratification variables
  #' @param nboot integer. Number of bootstrap repetitions
  #' @param bwh_subset_to_keep vector. String indicating BWH stages to keep
  #' @param ajcc8_subset_to_keep vector. String indicating AJCC8 stages to keep
  #' 
  #' @return weighted/unweighted AUC and C-index and corresponding standard errors
  
  if (!is.null(bwh_subset_to_keep)){
    df_mod <- df %>%
       filter(BWH_staging %in% bwh_subset_to_keep)
  } else {
    df_mod <- df %>%
      filter(AJCC_staging %in% ajcc8_subset_to_keep)              
  }
  # Preprocess data
  df_gathered <- df_mod %>%
     mutate(BWH_staging = factor(BWH_staging, levels = c("T1", "T2a", "T2b", "T3")),
            AJCC_staging = factor(AJCC_staging, levels = c("T1", "T2", "T3", "T4"))) %>%
     droplevels() %>%
     mutate(BWH_model = as.numeric(BWH_staging),
            AJCC8_model = as.numeric(AJCC_staging)) %>%
     select(all_of(c(ids_var, setdiff(strat_vars, c("All", "Type_of_material_bin")), "BWH_model", "AJCC8_model", "CP_risk"))) %>%
     left_join(preds_df, by = ids_var) %>%
     gather(model, score, -all_of(c(ids_var, setids_var, status_var, time_var, weight_var, strat_vars))) %>%
     gather(variable, group, -score, -model, -all_of(c(ids_var, setids_var, status_var, weight_var, time_var))) %>%
     filter(!((model == "AJCC8_model" & variable %in% c("AJCC_staging", "BWH_staging")) |
             (model == "BWH_model" & variable %in% c("AJCC_staging", "BWH_staging"))))
  
  # Compute AUC and C-index for each model
  perf_df <- df_gathered %>%
     group_by(model, variable, group) %>%
     summarise(N = n(),
               auc = ifelse(length(unique(Metastasis_num)) == 1,
                            NA,
                            as.numeric(gsub(".*: ", "", roc(cur_data()[[status_var]], score, quiet = TRUE, direction = "<")$auc))),
               auc_se = ifelse(length(unique(Metastasis_num)) == 1,
                               NA,
                               ifelse(all(variable == "All"),
                                      boot_metrics_se_ci(cur_data(), "score", status_var, NULL, weight_var, unique(setids_var), nboot, "auc")$se,
                                      boot_metrics_se_ci(cur_data(), "score", status_var, NULL, weight_var, NULL, nboot, "auc")$se)),
               wauc = ifelse(length(unique(Metastasis_num)) == 1,
                             NA,
                             WeightedAUC(WeightedROC(score, cur_data()[[status_var]], weight = cur_data()[[weight_var]]))),
               wauc_se = ifelse(length(unique(Metastasis_num)) == 1,
                                NA,
                                ifelse(all(variable == "All"),
                                       boot_metrics_se_ci(cur_data(), "score", status_var, NULL, weight_var, unique(setids_var), nboot, "wauc")$se,
                                       boot_metrics_se_ci(cur_data(), "score", status_var, NULL, weight_var, NULL, nboot, "wauc")$se)), 
               cindex = ifelse(length(unique(Metastasis_num)) == 1,
                               NA,
                               1 - rcorr.cens(score, Surv(cur_data()[[time_var]], cur_data()[[status_var]]))["C Index"] %>% `names<-`(NULL)),
               cindex_se = ifelse(length(unique(Metastasis_num)) == 1,
                                  NA,
                                  ifelse(all(variable == "All"),
                                         boot_metrics_se_ci(cur_data(), "score", status_var, time_var, weight_var, unique(setids_var), nboot, "cindex")$se,
                                         boot_metrics_se_ci(cur_data(), "score", status_var, time_var, weight_var, NULL, nboot, "cindex")$se)),
               wcindex = ifelse(length(unique(Metastasis_num)) == 1,
                                NA,
                                cIndex(cur_data()[[time_var]], event = cur_data()[[status_var]], score, weight = cur_data()[[weight_var]])[[1]]),
               wcindex_se = ifelse(length(unique(Metastasis_num)) == 1,
                                   NA,
                                   ifelse(all(variable == "All"),
                                          boot_metrics_se_ci(cur_data(), "score", status_var, time_var, weight_var, unique(setids_var), nboot, "wcindex")$se,
                                          boot_metrics_se_ci(cur_data(), "score", status_var, time_var, weight_var, NULL, nboot, "wcindex")$se)),
               .groups = "drop")
  perf_df
}
#-------------------------------------------------------------------------------

# Load libraries
#-------------------------------------------------------------------------------
#url <- "https://cran.r-project.org/src/contrib/psfmi_1.4.0.tar.gz"
#install.packages(url, repos = NULL, type = "source")
library(conflicted)
library(openxlsx)
library(tidyverse)
library(survival)
library(glmnet)
library(Hmisc)
library(mice)
library(psfmi)
library(pROC)
library(boot)
library(WeightedROC)
library(intsurv)
conflicted::conflicts_prefer(dplyr::filter)
#-------------------------------------------------------------------------------

# Filenames and settings
#-------------------------------------------------------------------------------
# Clinical data
clinical_fn <- e3_clinical
clinical_mimp_fn <- file.path(results_dir, "intermediate", "clinical", "clinical_data_imputed_d_squame_validation.Rdata")
# Weights
lr_weights_fn <- e7_weights
# SCCore-GEP predictions
sccoregep_preds_fn <- file.path(results_dir, "intermediate", "validation", "d_squame", "SCCoreGEP_predictions.csv")
# Time point at which c-index is computed
tp <- 5
# Number of bootstrap repetition
n_boots <- 100
# Subset of cohort to consider
subset_oi <- ifelse(exists("subset_oi"), subset_oi, "entiredataset")
if (subset_oi == "t1t2a"){
       bwh_subset_to_keep <- c("T1", "T2a")
       ajcc8_subset_to_keep <- NULL
  } else if (subset_oi == "t1t2") {
       bwh_subset_to_keep <- NULL
       ajcc8_subset_to_keep <- c("T1", "T2")
  } else {
       bwh_subset_to_keep <- c("T1", "T2a", "T2b", "T3")
       ajcc8_subset_to_keep <- NULL
}
# Output folder
output_dir <- file.path(results_dir, "intermediate", "validation", "d_squame")
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}
#-------------------------------------------------------------------------------

# Source functions
#-------------------------------------------------------------------------------
source(file.path(code_dir, "functions", "boot_metrics_se_ci.R"))
#-------------------------------------------------------------------------------

# Read data
#-------------------------------------------------------------------------------
# Clinical data
clinical <- read.csv(clinical_fn)
load(clinical_mimp_fn)
## Keep only biopsy for patients with duplicated sample
clinical <- clinical %>%
      filter(!(Set_info == "complete_set_double_sampletype" & Type_of_material_bin != "Biopsy"))
# Weights
lr_weights <- read.csv(lr_weights_fn)
# SCCore-GEP predictions
sccoregep_preds <- read.csv(sccoregep_preds_fn)
#-------------------------------------------------------------------------------

# Preprocess data
#-------------------------------------------------------------------------------
model_predictions_df <- sccoregep_preds %>%
        full_join(clinical %>%
                     select(Skyline_ID, Set_id, Metastasis, FU_metastasis_years, Type_of_material_bin),
                  by = "Skyline_ID") %>%
        full_join(lr_weights %>% select(Patient_ID_SKY, Weight_rescaled) %>% unique(),
                  by = "Patient_ID_SKY") %>%
        mutate(Metastasis_num = ifelse(Metastasis == "Case", 1, 0)) %>%
        select(-Metastasis, -Skyline_ID)
#-------------------------------------------------------------------------------

# Compute discrimination performances in imputed datasets
#-------------------------------------------------------------------------------
# Identify stratification variables (in imputed datasets)
strat_vars_mimp <- c("All", "Type_of_material_bin", "AJCC_staging", "BWH_staging")

# Compute performances (AUC and C-index) in imputed datasets
perf_mimp <- sapply(1:mimp$m, function(i) {
       compute_perfs(complete(mimp, i),
                     model_predictions_df  %>%
                         mutate(All = "All"),
                     "Patient_ID_SKY",
                     "Set_id",
                     "Metastasis_num",
                     "FU_metastasis_years",
                     "Weight_rescaled",
                     strat_vars_mimp,
                     n_boots,
                     bwh_subset_to_keep,
                     ajcc8_subset_to_keep) %>%
             mutate(nimp = i)
       }, simplify = F) %>%
       bind_rows()

# Pool performances
perf_mimp_strata <- perf_mimp %>% select(model, variable, group) %>% unique()
pooled_perf_mimp <- vector("list", nrow(perf_mimp_strata))
for (i in 1:nrow(perf_mimp_strata)){
      # Define subset
      strata <- perf_mimp_strata[i,]
      model_strata <- strata$model
      variable_strata <- strata$variable
      group_strata <- strata$group
      # Subset dataframe
      df_strata <- perf_mimp %>% dplyr::filter(model == model_strata &
                                               variable == variable_strata &
                                               group == group_strata)
      # Pool AUC
      pooled_auc_res <- pool_auc(df_strata$auc, df_strata$auc_se, nimp = mimp$m)
      # Pool C-index
      pooled_cindex_res <- pool_auc(df_strata$cindex, df_strata$cindex_se, nimp = mimp$m)
      # Pool wAUC
      pooled_wauc_res <- pool_auc(df_strata$wauc, df_strata$wauc_se, nimp = mimp$m)
      # Pool wC-index
      pooled_wcindex_res <- pool_auc(df_strata$wcindex, df_strata$wcindex_se, nimp = mimp$m)
      # Save results in dataframe
      pooled_perf_mimp[[i]] <- data.frame(model = model_strata,
                                          variable = variable_strata,
                                          group = group_strata,
                                          N = mean(df_strata$N),
                                          auc = round(ifelse(is.nan(pooled_auc_res[[2]]), mean(df_strata$auc), pooled_auc_res[[2]]), 4),
                                          auc_ci95low = round(pooled_auc_res[[1]], 4),
                                          auc_ci95up = round(pooled_auc_res[[3]], 4),
                                          cindex = round(ifelse(is.nan(pooled_cindex_res[[2]]), mean(df_strata$cindex), pooled_cindex_res[[2]]), 4),
                                          cindex_ci95low = round(pooled_cindex_res[[1]], 4),
                                          cindex_ci95up = round(pooled_cindex_res[[3]], 4),
                                          wauc = round(ifelse(is.nan(pooled_wauc_res[[2]]), mean(df_strata$wauc), pooled_wauc_res[[2]]), 4),
                                          wauc_ci95low = round(pooled_wauc_res[[1]], 4),
                                          wauc_ci95up = round(pooled_wauc_res[[3]], 4),
                                          wcindex = round(ifelse(is.nan(pooled_wcindex_res[[2]]), mean(df_strata$wcindex), pooled_wcindex_res[[2]]), 4),
                                          wcindex_ci95low = round(pooled_wcindex_res[[1]], 4),
                                          wcindex_ci95up = round(pooled_wcindex_res[[3]], 4)) %>%
            mutate(auc = paste0(auc, " (", auc_ci95low, "-", auc_ci95up, ")"),
                   cindex = paste0(cindex, " (", cindex_ci95low, "-", cindex_ci95up, ")"),
                   wauc = paste0(wauc, " (", wauc_ci95low, "-", wauc_ci95up, ")"),
                   wcindex = paste0(wcindex, " (", wcindex_ci95low, "-", wcindex_ci95up, ")")) %>%
            select(-auc_ci95low, -auc_ci95up, -cindex_ci95low, -cindex_ci95up,
                   -wauc_ci95low, -wauc_ci95up, -wcindex_ci95low, -wcindex_ci95up)
}
pooled_perf_mimp <- pooled_perf_mimp %>% bind_rows()
#-------------------------------------------------------------------------------

# Compute discrimination performances in entire dataset for complete stratification variables
# and combine with performances from multiple imputation
#-------------------------------------------------------------------------------
if (is.null(subset_oi)){
       # Identify stratification variables
       strat_vars <- setdiff(strat_vars_mimp, c("AJCC_staging", "BWH_staging"))
       # Compute performances (AUC and C-index)
       perf <- model_predictions_df %>%
              mutate(All = "All",
                     Immunosuppressed = ifelse(HM_at_cSCC == "Yes" | OTR_at_cSCC == "Yes", "Yes", "No"),
                     Prior_cSCCs = ifelse(Number_of_cSCC_before_culprit == 0, "No", "Yes")) %>%
              select(-HM_at_cSCC, -OTR_at_cSCC, -Number_of_cSCC_before_culprit) %>%
              gather(model, score, -all_of(c("Patient_ID_SKY", "Set_id", "Metastasis_num", "FU_metastasis_years", "Weight_rescaled", strat_vars))) %>%
              gather(variable, group, -model, -score, -Patient_ID_SKY, -Set_id, -Metastasis_num, -FU_metastasis_years, -Weight_rescaled) %>%
              group_by(model, variable, group) %>%
              summarise(N = n(),
                        auc = ifelse(length(unique(Metastasis_num)) == 1, NA, round(as.numeric(gsub(".*: ", "", roc(Metastasis_num, as.numeric(score), quiet = TRUE, direction = "<")$auc)), 4)),
                        auc_ci95 = ifelse(length(unique(Metastasis_num)) == 1,
                                          "(NA-NA)",
                                          ifelse(all(variable == "All"),
                                                 boot_metrics_se_ci(cur_data(), "score", "Metastasis_num", NULL, "Weight_rescaled", "Set_id", n_boots, "auc")$ci95,
                                                 boot_metrics_se_ci(cur_data(), "score", "Metastasis_num", NULL, "Weight_rescaled", NULL, n_boots, "auc")$ci95)),,
                        cindex = ifelse(length(unique(Metastasis_num)) == 1, NA, 1 - round(rcorr.cens(score, Surv(FU_metastasis_years, Metastasis_num))["C Index"] %>% `names<-`(NULL), 4)),
                        cindex_ci95 = ifelse(length(unique(Metastasis_num)) == 1,
                                            "(NA-NA)",
                                            ifelse(all(variable == "All"),
                                                   boot_metrics_se_ci(cur_data(), "score", "Metastasis_num", "FU_metastasis_years", "Weight_rescaled", "Set_id", n_boots, "cindex")$ci95,
                                                   boot_metrics_se_ci(cur_data(), "score", "Metastasis_num", "FU_metastasis_years", "Weight_rescaled", NULL, n_boots, "cindex")$ci95)),
                        wauc = ifelse(length(unique(Metastasis_num)) == 1, NA, round(WeightedAUC(WeightedROC(score, Metastasis_num, weight = Weight_rescaled)), 4)),
                        wauc_ci95 = ifelse(length(unique(Metastasis_num)) == 1,
                                          "(NA-NA)",
                                          ifelse(all(variable == "All"),
                                                 boot_metrics_se_ci(cur_data(), "score", "Metastasis_num", NULL, "Weight_rescaled", "Set_id", n_boots, "wauc")$ci95,
                                                 boot_metrics_se_ci(cur_data(), "score", "Metastasis_num", NULL, "Weight_rescaled", NULL, n_boots, "wauc")$ci95)),
                        wcindex = ifelse(length(unique(Metastasis_num)) == 1, NA, round(cIndex(FU_metastasis_years, event = Metastasis_num, score, weight = Weight_rescaled)[[1]], 4)),
                        wcindex_ci95 = ifelse(length(unique(Metastasis_num)) == 1,
                                             "(NA-NA)",
                                             ifelse(all(variable == "All"),
                                                    boot_metrics_se_ci(cur_data(), "score", "Metastasis_num", "FU_metastasis_years", "Weight_rescaled", "Set_id", n_boots, "wcindex")$ci95,
                                                    boot_metrics_se_ci(cur_data(), "score", "Metastasis_num", "FU_metastasis_years", "Weight_rescaled", NULL, n_boots, "wcindex")$ci95)),
                        .groups = "drop") %>%
              mutate(auc = paste0(auc, " ", auc_ci95),
                     cindex = paste0(cindex, " ", cindex_ci95),
                     wauc = paste0(wauc, " ", wauc_ci95),
                     wcindex = paste0(wcindex, " ", wcindex_ci95)) %>%
              select(-auc_ci95, -cindex_ci95, -wauc_ci95, -wcindex_ci95)
       # Combine dataframe with pooled performance dataframe
       pooled_perf_mimp <- perf %>%
              rbind(pooled_perf_mimp) %>%
              arrange(model, variable, group)
}
#-------------------------------------------------------------------------------

# Compute discrimination performances in entire dataset (complete cases and max samples)
#-------------------------------------------------------------------------------
model_predictions_df_mod <- model_predictions_df %>%
        left_join(clinical %>% select(Patient_ID_SKY, AJCC_staging, BWH_staging, CP_risk),
                  by = "Patient_ID_SKY") %>%
        mutate(All = "All",
               BWH_staging = factor(BWH_staging, levels = c("T1", "T2a", "T2b", "T3")),
               AJCC_staging = factor(AJCC_staging, levels = c("T1", "T2", "T3", "T4"))) %>%
        droplevels() %>%
        mutate(BWH_model = as.numeric(BWH_staging),
               AJCC8_model = as.numeric(AJCC_staging))
if (!is.null(bwh_subset_to_keep)){
       model_predictions_df_mod <- model_predictions_df_mod %>%
         filter(BWH_staging %in% bwh_subset_to_keep)
} else {
       model_predictions_df_mod <- model_predictions_df_mod %>%
         filter(AJCC_staging %in% ajcc8_subset_to_keep)
}
# For maximum samples of each model
perf_max_samples <- model_predictions_df_mod %>%
        gather(variable, group, -Patient_ID_SKY, -Set_id, -Metastasis_num, -FU_metastasis_years,
               -Weight_rescaled, -sccore_gep, -AJCC8_model, -BWH_model, -CP_risk) %>%
        gather(model, score, -variable, -group, -Patient_ID_SKY, -Metastasis_num, -FU_metastasis_years, -Weight_rescaled) %>%
        filter(!((model == "AJCC8_model" & variable == "AJCC_staging") | (model == "BWH_model" & variable == "BWH_staging") | is.na(score))) %>%
        group_by(model, variable, group) %>%
        summarise(N = n(),
                  auc = ifelse(length(unique(Metastasis_num)) == 1, NA, round(as.numeric(gsub(".*: ", "", roc(Metastasis_num, as.numeric(score), quiet = TRUE, direction = "<")$auc)), 4)),
                  auc_ci95 = ifelse(length(unique(Metastasis_num)) == 1,
                                    "(NA-NA)",
                                    boot_metrics_se_ci(cur_data(), "score", "Metastasis_num", NULL, "Weight_rescaled", NULL, n_boots, "auc")$ci95),
                  cindex = ifelse(length(unique(Metastasis_num)) == 1, NA, 1 - round(rcorr.cens(score, Surv(FU_metastasis_years, Metastasis_num))["C Index"] %>% `names<-`(NULL), 4)),
                  cindex_ci95 = ifelse(length(unique(Metastasis_num)) == 1,
                                      "(NA-NA)",
                                      boot_metrics_se_ci(cur_data(), "score", "Metastasis_num", "FU_metastasis_years", "Weight_rescaled", NULL, n_boots, "cindex")$ci95),
                  wauc = ifelse(length(unique(Metastasis_num)) == 1, NA, round(WeightedAUC(WeightedROC(score, Metastasis_num, weight = Weight_rescaled)), 4)),
                  wauc_ci95 = ifelse(length(unique(Metastasis_num)) == 1,
                                     "(NA-NA)",
                                     boot_metrics_se_ci(cur_data(), "score", "Metastasis_num", NULL, "Weight_rescaled", NULL, n_boots, "wauc")$ci95),
                  wcindex = ifelse(length(unique(Metastasis_num)) == 1, NA, round(cIndex(FU_metastasis_years, event = Metastasis_num, score, weight = Weight_rescaled)[[1]], 4)),
                  wcindex_ci95 = ifelse(length(unique(Metastasis_num)) == 1,
                                       "(NA-NA)",
                                        boot_metrics_se_ci(cur_data(), "score", "Metastasis_num", "FU_metastasis_years", "Weight_rescaled", NULL, n_boots, "wcindex")$ci95),
               .groups = "drop") %>%
     mutate(auc = paste0(auc, " ", auc_ci95),
            cindex = paste0(cindex, " ", cindex_ci95),
            wauc = paste0(wauc, " ", wauc_ci95),
            wcindex = paste0(wcindex, " ", wcindex_ci95)) %>%
     select(-auc_ci95, -cindex_ci95, -wauc_ci95, -wcindex_ci95)
# For complete cases of all models
perf_complete_cases <- model_predictions_df_mod %>%
        filter(!(is.na(AJCC8_model) | is.na(BWH_model) | is.na(CP_risk))) %>%
        gather(variable, group, -Patient_ID_SKY, -Set_id, -Metastasis_num, -FU_metastasis_years,
               -Weight_rescaled, -sccore_gep, -AJCC8_model, -BWH_model, -CP_risk) %>%
        gather(model, score, -variable, -group, -Patient_ID_SKY, -Metastasis_num, -FU_metastasis_years, -Weight_rescaled) %>%
        filter(!((model == "AJCC8_model" & variable == "AJCC_staging") | (model == "BWH_model" & variable == "BWH_staging") | is.na(score))) %>%
        group_by(model, variable, group) %>%
        summarise(N = n(),
                  auc = ifelse(length(unique(Metastasis_num)) == 1, NA, round(as.numeric(gsub(".*: ", "", roc(Metastasis_num, as.numeric(score), quiet = TRUE, direction = "<")$auc)), 4)),
                  auc_ci95 = ifelse(length(unique(Metastasis_num)) == 1,
                                    "(NA-NA)",
                                    boot_metrics_se_ci(cur_data(), "score", "Metastasis_num", NULL, "Weight_rescaled", NULL, n_boots, "auc")$ci95),
                  cindex = ifelse(length(unique(Metastasis_num)) == 1, NA, 1 - round(rcorr.cens(score, Surv(FU_metastasis_years, Metastasis_num))["C Index"] %>% `names<-`(NULL), 4)),
                  cindex_ci95 = ifelse(length(unique(Metastasis_num)) == 1,
                                      "(NA-NA)",
                                      boot_metrics_se_ci(cur_data(), "score", "Metastasis_num", "FU_metastasis_years", "Weight_rescaled", NULL, n_boots, "cindex")$ci95),
                  wauc = ifelse(length(unique(Metastasis_num)) == 1, NA, round(WeightedAUC(WeightedROC(score, Metastasis_num, weight = Weight_rescaled)), 4)),
                  wauc_ci95 = ifelse(length(unique(Metastasis_num)) == 1,
                                     "(NA-NA)",
                                     boot_metrics_se_ci(cur_data(), "score", "Metastasis_num", NULL, "Weight_rescaled", NULL, n_boots, "wauc")$ci95),
                  wcindex = ifelse(length(unique(Metastasis_num)) == 1, NA, round(cIndex(FU_metastasis_years, event = Metastasis_num, score, weight = Weight_rescaled)[[1]], 4)),
                  wcindex_ci95 = ifelse(length(unique(Metastasis_num)) == 1,
                                       "(NA-NA)",
                                        boot_metrics_se_ci(cur_data(), "score", "Metastasis_num", "FU_metastasis_years", "Weight_rescaled", NULL, n_boots, "wcindex")$ci95),
               .groups = "drop") %>%
     mutate(auc = paste0(auc, " ", auc_ci95),
            cindex = paste0(cindex, " ", cindex_ci95),
            wauc = paste0(wauc, " ", wauc_ci95),
            wcindex = paste0(wcindex, " ", wcindex_ci95)) %>%
     select(-auc_ci95, -cindex_ci95, -wauc_ci95, -wcindex_ci95)
#-------------------------------------------------------------------------------

# Save results
#-------------------------------------------------------------------------------
write.csv(model_predictions_df  %>%
              full_join(clinical %>% select(Patient_ID_SKY, Skyline_ID),
                        by = "Patient_ID_SKY") %>%
              select(Skyline_ID, Patient_ID_SKY, sccore_gep),
          file.path(output_dir, "SCCoreGEP_predictions.csv"),
          row.names = F)
write.csv(pooled_perf_mimp ,
          file.path(output_dir, paste0("SCCoreGEP_BWH_AJCC8_EMCmodel_performances_", ifelse(is.null(subset_oi), "", paste0(subset_oi, "_")), "with_imputation.csv")),
          row.names = F)
write.csv(perf_max_samples,
          file.path(output_dir, paste0("SCCoreGEP_BWH_AJCC8_EMCmodel_performances_", ifelse(is.null(subset_oi), "", paste0(subset_oi, "_")), "max_samples.csv")),
          row.names = F)
write.csv(perf_complete_cases,
          file.path(output_dir, paste0("SCCoreGEP_BWH_AJCC8_EMCmodel_performances_", ifelse(is.null(subset_oi), "", paste0(subset_oi, "_")), "complete_cases.csv")),
          row.names = F)
#-------------------------------------------------------------------------------
