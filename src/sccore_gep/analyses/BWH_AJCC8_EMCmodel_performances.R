#-------------------------------------------------------------------------------
# Aim: Compute performances of BWH, AJCC8, EMC model in D-SQUAME discovery cohort
# Author: L.Pozza
# Input: Clinical data of D-SQUAME discovery cohort
# Output: Performances of BWH, AJCC8, EMC model in D-SQUAME discovery cohort
#-------------------------------------------------------------------------------

# Functions
#-------------------------------------------------------------------------------
compute_perfs <- function(df, weights_df, ids_var, setids_var, status_var, time_var, weight_var, strat_vars, nboot){
  #' Compute AUC and C-index stratified by (imputed) variables
  #'
  #' @param df dataframe. Dataframe with imputed values
  #' @param weights_df dataframe. Dataframe with weights
  #' @param ids_var string. String indicating column with samples IDs
  #' @param setids_var string. String indicating column with set IDs
  #' @param status_var string. String indicating column with status
  #' @param time_var string. String indicating column with time to status
  #' @param weight_var string. String indicating column with weights
  #' @param strat_var string vector. Strings indicating stratification variables
  #' @param nboot integer. Number of bootstrap repetitions
  #' 
  #' @return weighted/unweighted AUC and C-index and corresponding standard errors
  
  # Preprocess data
  df_gathered <- df %>%
    inner_join(weights_df, by = ids_var) %>%
    mutate(Metastasis_num = ifelse(Metastasis == "Case", 1, 0),
           Immunosuppressed = ifelse(OTR_at_cSCC == "Yes" | HM_at_cSCC == "Yes", "Yes", "No"),
           Sample_type_bin = ifelse(Biopsy_excision == "Biopsy", "Biopsy", "Excision"),
           All = "All",
           BWH_model = as.numeric(factor(BWH, levels = c("T1", "T2a", "T2b", "T3"))),
           AJCC8_model = as.numeric(factor(AJCC_8, levels = c("T1", "T2", "T3", "T4")))) %>%
    select(all_of(c(ids_var, setids_var, status_var, time_var, weight_var, strat_vars, "BWH_model", "AJCC8_model", "CP_score"))) %>%
    gather(model, score, -all_of(c(ids_var, setids_var, status_var, time_var, weight_var, strat_vars))) %>%
    gather(variable, group, -score, -model, -all_of(c(ids_var, setids_var, status_var, weight_var, time_var))) %>%
    filter(!(model %in% c("BWH_model", "AJCC8_model") & variable %in% c("BWH_imputed_merged", "AJCC8_imputed_merged")))
  
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
              wauc =  ifelse(length(unique(Metastasis_num)) == 1,
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
library(Hmisc)
library(mice)
library(psfmi)
library(pROC)
library(boot)
library(WeightedROC)
library(intsurv)
conflicted::conflicts_prefer(dplyr::filter)
conflicts_prefer(pROC::roc)
conflicts_prefer(pROC::var)
#-------------------------------------------------------------------------------

# Filenames and settings
#-------------------------------------------------------------------------------
# Clinical data
clinical_fn <- p3_clinical
clinical_mimp_fn <- file.path(results_dir, "intermediate", "clinical", "clinical_data_imputed_d_squame_discovery.Rdata")
clinical_mimp_merged_fn <- file.path(results_dir, "intermediate", "clinical", "merged_clinical_data_imputed_d_squame_discovery.csv")
# Weights
lr_weights_fn <- p10_weights
# Number of bootstrap repetition
n_boots <- 100
# Output folder
output_dir <- file.path(results_dir, "intermediate", "sccore_gep")
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}
#-------------------------------------------------------------------------------

# Source functions
#-------------------------------------------------------------------------------
source(file.path(code_dir, "functions", "boot_metrics_se_ci.R"))
#-------------------------------------------------------------------------------

# Read data
#-------------------------------------------------------------------------------
# Clinical
clinical_samples <- openxlsx::read.xlsx(clinical_fn, sheet = "Request Lara - Samples (values)")
# Clinical imputed
load(clinical_mimp_fn)
# Clinical imputed merged
clinical_mimp_merged <- read.csv(clinical_mimp_merged_fn)
# Weights
lr_weights <- read.csv(lr_weights_fn)
#-------------------------------------------------------------------------------

# Preprocess data
#-------------------------------------------------------------------------------
lr_weights <- lr_weights %>%
  rename(Skyline_ID = SkylineDx.ID) %>%
  left_join(clinical_samples %>% select(Skyline_ID, Set_id),
            by = "Skyline_ID") %>%
  left_join(clinical_mimp_merged %>% select(SkylineDx.ID, BWH, AJCC_8), by = c("Skyline_ID" = "SkylineDx.ID")) %>%
  select(Skyline_ID, Set_id, Weight_rescaled, BWH, AJCC_8) %>%
  mutate(BWH_bin_imputed_merged = ifelse(BWH %in% c("T1", "T2a"), "T1-T2a", "T2b-T3"),
         AJCC8_bin_imputed_merged = ifelse(AJCC_8 %in% c("T1", "T2"), "T1-T2", "T3-T4")) %>%
  rename(BWH_imputed_merged = BWH,
         AJCC8_imputed_merged = AJCC_8)
#-------------------------------------------------------------------------------

# Compute discrimination performances in imputed datasets
#-------------------------------------------------------------------------------
# Identify stratification variables (in imputed datasets)
strat_vars <- c("All", "Immunosuppressed", "Sample_type_bin", "BWH_imputed_merged", "AJCC8_imputed_merged", "BWH_bin_imputed_merged", "AJCC8_bin_imputed_merged")
# Compute performances (AUC and C-index) in imputed datasets
perf_mimp <- sapply(1:mimp$m, function(i) {
  compute_perfs(complete(mimp, i),
                lr_weights,
                "Skyline_ID",
                "Set_id",
                "Metastasis_num",
                "FU_metastasis_years",
                "Weight_rescaled",
                strat_vars,
                n_boots) %>%
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
                                      auc_95ci_low = round(pooled_auc_res[[1]], 4),
                                      auc_95ci_up = round(pooled_auc_res[[3]], 4),
                                      cindex = round(ifelse(is.nan(pooled_cindex_res[[2]]), mean(df_strata$cindex), pooled_cindex_res[[2]]), 4),
                                      cindex_95ci_low = round(pooled_cindex_res[[1]], 4),
                                      cindex_95ci_up = round(pooled_cindex_res[[3]], 4),
                                      wauc = round(ifelse(is.nan(pooled_wauc_res[[2]]), mean(df_strata$wauc), pooled_wauc_res[[2]]), 4),
                                      wauc_95ci_low = round(pooled_wauc_res[[1]], 4),
                                      wauc_95ci_up = round(pooled_wauc_res[[3]], 4),
                                      wcindex = round(ifelse(is.nan(pooled_wcindex_res[[2]]), mean(df_strata$wcindex), pooled_wcindex_res[[2]]), 4),
                                      wcindex_95ci_low = round(pooled_wcindex_res[[1]], 4),
                                      wcindex_95ci_up = round(pooled_wcindex_res[[3]], 4)) %>%
    mutate(auc = paste0(auc, " (", auc_95ci_low, "-", auc_95ci_up, ")"),
           cindex = paste0(cindex, " (", cindex_95ci_low, "-", cindex_95ci_up, ")"),
           wauc = paste0(wauc, " (", wauc_95ci_low, "-", wauc_95ci_up, ")"),
           wcindex = paste0(wcindex, " (", wcindex_95ci_low, "-", wcindex_95ci_up, ")")) %>%
    select(-auc_95ci_low, -auc_95ci_up, -cindex_95ci_low, -cindex_95ci_up,
           -wauc_95ci_low, -wauc_95ci_up, -wcindex_95ci_low, -wcindex_95ci_up)
  
}
pooled_perf_mimp <- pooled_perf_mimp %>% bind_rows()
#-------------------------------------------------------------------------------

# Compute discrimination performances in entire dataset (complete cases and max samples)
#-------------------------------------------------------------------------------
clinical_df_mod <- complete(mimp, 0) %>%
  mutate(Metastasis_num = ifelse(Metastasis == "Case", 1, 0),
         Immunosuppressed = ifelse(OTR_at_cSCC == "Yes" | HM_at_cSCC == "Yes", "Yes", "No"),
         Sample_type_bin = ifelse(Biopsy_excision == "Biopsy", "Biopsy", "Excision"),
         All = "All",
         BWH_model = as.numeric(factor(BWH, levels = c("T1", "T2a", "T2b", "T3"))),
         AJCC8_model = as.numeric(factor(AJCC_8, levels = c("T1", "T2", "T3", "T4")))) %>%
  inner_join(lr_weights, by = "Skyline_ID") %>%
  select(Skyline_ID, Metastasis_num, FU_metastasis_years, Weight_rescaled, All, Immunosuppressed,
         Sample_type_bin, BWH_imputed_merged, AJCC8_imputed_merged, BWH_bin_imputed_merged,
         AJCC8_bin_imputed_merged, BWH_model, AJCC8_model, CP_score)
# For maximum samples of each model
perf_max_samples <- clinical_df_mod %>%
  gather(variable, group, -Skyline_ID, -Metastasis_num, -FU_metastasis_years, -Weight_rescaled, -AJCC8_model, -BWH_model, -CP_score) %>%
  gather(model, score, -variable, -group, -Skyline_ID, -Metastasis_num, -FU_metastasis_years, -Weight_rescaled) %>%
  filter(!(is.na(score) | (model %in% c("BWH_model", "AJCC8_model") & variable %in% c("BWH_imputed_merged", "AJCC8_imputed_merged")))) %>%
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
perf_complete_cases <- clinical_df_mod %>%
  filter(!(is.na(AJCC8_model) | is.na(BWH_model) | is.na(CP_score))) %>%
  gather(variable, group, -Skyline_ID, -Metastasis_num, -FU_metastasis_years, -Weight_rescaled, -AJCC8_model, -BWH_model, -CP_score) %>%
  gather(model, score, -variable, -group, -Skyline_ID, -Metastasis_num, -FU_metastasis_years, -Weight_rescaled) %>%
  filter(!(is.na(score) | (model %in% c("BWH_model", "AJCC8_model") & variable %in% c("BWH_imputed_merged", "AJCC8_imputed_merged")))) %>%
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
write.csv(pooled_perf_mimp,
          file.path(output_dir, "BWH_AJCC8_EMCmodel_performances_with_imputation.csv"),
          row.names = F)
write.csv(perf_max_samples,
          file.path(output_dir, "BWH_AJCC8_EMCmodel_performances_max_samples.csv"),
          row.names = F)
write.csv(perf_complete_cases,
          file.path(output_dir, "BWH_AJCC8_EMCmodel_performances_complete_cases.csv"),
          row.names = F)
#--------
