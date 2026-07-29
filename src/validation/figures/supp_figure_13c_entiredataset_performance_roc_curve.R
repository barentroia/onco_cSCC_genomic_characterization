# Load libraries
#-------------------------------------------------------------------------------
library(conflicted)
library(tidyverse)
library(patchwork)
library(WeightedROC)
library(mice)
conflicted::conflicts_prefer(dplyr::filter)
#-------------------------------------------------------------------------------

# Filenames and settings
#-------------------------------------------------------------------------------
# D-SQUAME validation dataset
## Clinical data
clinical_d_squame_fn <- e3_clinical
clinical_mimp_d_squame_fn <- file.path(results_dir, "intermediate", "clinical", "clinical_data_imputed_d_squame_validation.Rdata")
# Weights
lr_weights_d_squame_fn <- e7_weights
# SCCore-GEP predictions
predictions_d_squame_fn <- file.path(results_dir, "intermediate", "validation", "d_squame", "SCCoreGEP_predictions.csv")
# SCCore-GEP performances in entire dataset
performances_entiredataset_d_squame_fn <- file.path(results_dir, "intermediate", "validation", "d_squame", "SCCoreGEP_BWH_AJCC8_EMCmodel_performances_entiredataset_with_imputation.csv")
# Nassir et al. validation dataset
## Clinical data
clinical_nassir_fn <- n3_clinical
## SCCore-GEP predictions
predictions_nassir_fn <- file.path(results_dir, "intermediate", "validation", "nassir", "SCCoreGEP_predictions.csv")
# SCCore-GEP performances
performances_nassir_fn <- file.path(results_dir, "intermediate", "validation", "nassir", "SCCoreGEP_BWH_performances_max_samples.csv")

# Output folder
output_dir <- file.path(results_dir, "publication")
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}
#-------------------------------------------------------------------------------

# Source functions
#-------------------------------------------------------------------------------
source(file.path(code_dir, "functions", "colors.R"))
#-------------------------------------------------------------------------------

# Read data
#-------------------------------------------------------------------------------
# D-SQUAME validation data
## Clinical data
clinical_d_squame <- read.csv(clinical_d_squame_fn)
load(clinical_mimp_d_squame_fn)
## Weights
lr_weights_d_squame <- read.csv(lr_weights_d_squame_fn)
## Predictions
predictions_d_squame <- read.csv(predictions_d_squame_fn)
## Performances in entire dataset
performances_entiredataset_d_squame <- read.csv(performances_entiredataset_d_squame_fn)

# Nassir data
## Clinical data
clinical_nassir <- read.csv(clinical_nassir_fn)
## Predictions
predictions_nassir <- read.csv(predictions_nassir_fn)
## Performances
performances_nassir <- read.csv(performances_nassir_fn)
#-------------------------------------------------------------------------------

# Preprocess clinical and predictions data
#-------------------------------------------------------------------------------
# D-SQUAME validation data: keep only biopsy for patients with double sample type, add column with weights
d_squame_df <- clinical_d_squame %>%
      filter(!(Set_info == "complete_set_double_sampletype" & Type_of_material_bin != "Biopsy")) %>%
      mutate(Metastasis_num = ifelse(Metastasis == "Case", 1, 0)) %>%
      inner_join(predictions_d_squame, by = c("Skyline_ID", "Patient_ID_SKY")) %>%
      left_join(lr_weights_d_squame %>% select(Patient_ID_SKY, Weight_rescaled) %>% unique(), by = "Patient_ID_SKY") %>%
      select(Skyline_ID, Patient_ID_SKY, Metastasis_num, Weight_rescaled, sccore_gep)


# Nassir: combine predictions with clinical data
nassir_df <- predictions_nassir %>%
    inner_join(clinical_nassir, by = c("Sample_id" = "Run")) %>%
    mutate(BWH = str_extract(SampleID, "T[0-9]+[a-z]?"),
           Metastasis_num = ifelse(MetNoMet == "Met", 1, 0),
           BWH = ifelse(is.na(BWH), "T2b", BWH),
           BWH_num = as.numeric(factor(BWH, levels = c("T1", "T2a", "T2b", "T3")))) %>%
    select(Sample_id, Metastasis_num, BWH, sccore_gep, BWH_num)
#-------------------------------------------------------------------------------

# ROC plots
#-------------------------------------------------------------------------------
# Entire dataset
## ROC
d_squame_roc_entiredataset_sccore_gep_df <- WeightedROC(d_squame_df$sccore_gep, d_squame_df$Metastasis_num, d_squame_df$Weight_rescaled)
d_squame_roc_entiredataset_bwh_df <- sapply(1:mimp$m, function(i, df, id, outcome_var, pred_var, weight_var) {
        mimp_df <- complete(mimp, i) %>%
            select(all_of(c(id, "BWH_staging"))) %>%
            full_join(df, by = id) %>%
            mutate(BWH_num = as.numeric(factor(BWH_staging, levels = c("T1", "T2a", "T2b", "T3"))))
        mimp_roc_df <- WeightedROC(mimp_df[[pred_var]], mimp_df[[outcome_var]], mimp_df[[weight_var]])
        },
        df = d_squame_df,
        id = "Patient_ID_SKY",
        outcome_var = "Metastasis_num",
        pred_var = "BWH_num",
        weight_var = "Weight_rescaled",
        simplify = F) %>%
    bind_rows() %>%
    gather(metric, value, -threshold) %>%
    group_by(metric, threshold) %>%
    summarise(value = mean(value), .groups = "keep") %>%
    ungroup() %>%
    spread(metric, value)
d_squame_roc_entiredataset_ajcc8_df <- sapply(1:mimp$m, function(i, df, id, outcome_var, pred_var, weight_var) {
        mimp_df <- complete(mimp, i) %>%
            select(all_of(c(id, "BWH_staging", "AJCC_staging"))) %>%
            full_join(df, by = id) %>%
            mutate(AJCC_num = as.numeric(factor(AJCC_staging, levels = c("T1", "T2", "T3", "T4"))))
        mimp_roc_df <- WeightedROC(mimp_df[[pred_var]], mimp_df[[outcome_var]], mimp_df[[weight_var]])
        },
        df = d_squame_df,
        id = "Patient_ID_SKY",
        outcome_var = "Metastasis_num",
        pred_var = "AJCC_num",
        weight_var = "Weight_rescaled",
        simplify = F) %>%
    bind_rows() %>%
    gather(metric, value, -threshold) %>%
    group_by(metric, threshold) %>%
    summarise(value = mean(value), .groups = "keep") %>%
    ungroup() %>%
    spread(metric, value)
nassir_roc_entiredataset_sccore_gep_df <- WeightedROC(nassir_df$sccore_gep, nassir_df$Metastasis_num)
nassir_roc_entiredataset_bwh_df <- WeightedROC(nassir_df$BWH_num, nassir_df$Metastasis_num)
## AUC
d_squame_auc_entiredataset_sccore_gep <- performances_entiredataset_d_squame %>%
    filter(model == "sccore_gep" & variable == "All")
d_squame_auc_entiredataset_bwh <- performances_entiredataset_d_squame %>%
    filter(model == "BWH_model" & variable == "All")
d_squame_auc_entiredataset_ajcc8 <- performances_entiredataset_d_squame %>%
    filter(model == "AJCC8_model" & variable == "All")
nassir_auc_entiredataset_sccore_gep <- performances_nassir %>%
    filter(model == "sccore_gep" & variable == "All")
nassir_auc_entiredataset_bwh <- performances_nassir %>%
    filter(model == "bwh" & variable == "All")
## Combine ROC and AUC
d_squame_roc_entiredataset_df <- d_squame_roc_entiredataset_sccore_gep_df %>% mutate(strat_system = paste0("SCCore-GEP, Weighted AUC = ",
                                                             round(as.numeric(sapply(strsplit(d_squame_auc_entiredataset_sccore_gep$wauc, "[(|)|-]"), "[[", 1)), 2),
                                                             " (",
                                                             round(as.numeric(sapply(strsplit(d_squame_auc_entiredataset_sccore_gep$wauc, "[(|)|-]"), "[[", 2)), 2),
                                                             "-",
                                                             round(as.numeric(sapply(strsplit(d_squame_auc_entiredataset_sccore_gep$wauc, "[(|)|-]"), "[[", 3)), 2),
                                                             ")"), cohort = "D-SQUAME validation (N=102)") %>%
    rbind(d_squame_roc_entiredataset_bwh_df %>% mutate(strat_system = paste0("BWH, Weighted AUC = ",
                                                          round(as.numeric(sapply(strsplit(d_squame_auc_entiredataset_bwh$wauc, "[(|)|-]"), "[[", 1)), 2),
                                                          " (",
                                                          round(as.numeric(sapply(strsplit(d_squame_auc_entiredataset_bwh$wauc, "[(|)|-]"), "[[", 2)), 2),
                                                          "-",
                                                          round(as.numeric(sapply(strsplit(d_squame_auc_entiredataset_bwh$wauc, "[(|)|-]"), "[[", 3)), 2),
                                                          ")"), cohort = "D-SQUAME validation (N=102)")) %>%
    rbind(d_squame_roc_entiredataset_ajcc8_df %>% mutate(strat_system = paste0("AJCC8, Weighted AUC = ",
                                                          round(as.numeric(sapply(strsplit(d_squame_auc_entiredataset_ajcc8$wauc, "[(|)|-]"), "[[", 1)), 2),
                                                          " (",
                                                          round(as.numeric(sapply(strsplit(d_squame_auc_entiredataset_ajcc8$wauc, "[(|)|-]"), "[[", 2)), 2),
                                                          "-",
                                                          round(as.numeric(sapply(strsplit(d_squame_auc_entiredataset_ajcc8$wauc, "[(|)|-]"), "[[", 3)), 2),
                                                          ")"), cohort = "D-SQUAME validation (N=102)")) %>%
    mutate(strat_system = factor(strat_system, levels = c(paste0("SCCore-GEP, Weighted AUC = ",
                                                             round(as.numeric(sapply(strsplit(d_squame_auc_entiredataset_sccore_gep$wauc, "[(|)|-]"), "[[", 1)), 2),
                                                             " (",
                                                             round(as.numeric(sapply(strsplit(d_squame_auc_entiredataset_sccore_gep$wauc, "[(|)|-]"), "[[", 2)), 2),
                                                             "-",
                                                             round(as.numeric(sapply(strsplit(d_squame_auc_entiredataset_sccore_gep$wauc, "[(|)|-]"), "[[", 3)), 2),
                                                             ")"),
                                                          paste0("BWH, Weighted AUC = ",
                                                             round(as.numeric(sapply(strsplit(d_squame_auc_entiredataset_bwh$wauc, "[(|)|-]"), "[[", 1)), 2),
                                                             " (",
                                                             round(as.numeric(sapply(strsplit(d_squame_auc_entiredataset_bwh$wauc, "[(|)|-]"), "[[", 2)), 2),
                                                             "-",
                                                             round(as.numeric(sapply(strsplit(d_squame_auc_entiredataset_bwh$wauc, "[(|)|-]"), "[[", 3)), 2),
                                                             ")"),
                                                          paste0("AJCC8, Weighted AUC = ",
                                                          round(as.numeric(sapply(strsplit(d_squame_auc_entiredataset_ajcc8$wauc, "[(|)|-]"), "[[", 1)), 2),
                                                          " (",
                                                          round(as.numeric(sapply(strsplit(d_squame_auc_entiredataset_ajcc8$wauc, "[(|)|-]"), "[[", 2)), 2),
                                                          "-",
                                                          round(as.numeric(sapply(strsplit(d_squame_auc_entiredataset_ajcc8$wauc, "[(|)|-]"), "[[", 3)), 2),
                                                          ")"))))
nassir_roc_entiredataset_df <- nassir_roc_entiredataset_sccore_gep_df %>% mutate(strat_system = paste0("SCCore-GEP, AUC = ",
                                                             round(as.numeric(sapply(strsplit(nassir_auc_entiredataset_sccore_gep$auc, "[(|)|-]"), "[[", 1)), 2),
                                                             " (",
                                                             round(as.numeric(sapply(strsplit(nassir_auc_entiredataset_sccore_gep$auc, "[(|)|-]"), "[[", 2)), 2),
                                                             "-",
                                                             round(as.numeric(sapply(strsplit(nassir_auc_entiredataset_sccore_gep$auc, "[(|)|-]"), "[[", 3)), 2),
                                                             ")"), cohort = "Nassir et al. (N=52)")  %>%
    rbind(nassir_roc_entiredataset_bwh_df %>% mutate(strat_system = paste0("BWH, AUC = ",
                                                          round(as.numeric(sapply(strsplit(nassir_auc_entiredataset_bwh$auc, "[(|)|-]"), "[[", 1)), 2),
                                                          " (",
                                                          round(as.numeric(sapply(strsplit(nassir_auc_entiredataset_bwh$auc, "[(|)|-]"), "[[", 2)), 2),
                                                          "-",
                                                          round(as.numeric(sapply(strsplit(nassir_auc_entiredataset_bwh$auc, "[(|)|-]"), "[[", 3)), 2),
                                                          ")"), cohort = "Nassir et al. (N=52)")) %>%
    mutate(strat_system = factor(strat_system, levels = c(paste0("SCCore-GEP, AUC = ",
                                                             round(as.numeric(sapply(strsplit(nassir_auc_entiredataset_sccore_gep$auc, "[(|)|-]"), "[[", 1)), 2),
                                                             " (",
                                                             round(as.numeric(sapply(strsplit(nassir_auc_entiredataset_sccore_gep$auc, "[(|)|-]"), "[[", 2)), 2),
                                                             "-",
                                                             round(as.numeric(sapply(strsplit(nassir_auc_entiredataset_sccore_gep$auc, "[(|)|-]"), "[[", 3)), 2),
                                                             ")"),
                                                          paste0("BWH, AUC = ",
                                                             round(as.numeric(sapply(strsplit(nassir_auc_entiredataset_bwh$auc, "[(|)|-]"), "[[", 1)), 2),
                                                             " (",
                                                             round(as.numeric(sapply(strsplit(nassir_auc_entiredataset_bwh$auc, "[(|)|-]"), "[[", 2)), 2),
                                                             "-",
                                                             round(as.numeric(sapply(strsplit(nassir_auc_entiredataset_bwh$auc, "[(|)|-]"), "[[", 3)), 2),
                                                             ")"))))

# Plots
## Nassir
nassir_roc_entiredataset_p <- ggplot(data = nassir_roc_entiredataset_df,
                               aes(x = FPR,
                                   y = TPR,
                                   color = strat_system))+
                geom_path(linewidth = 0.8)+
                scale_color_manual(values = c(cols_strat_systems[["SCCore-GEP"]], cols_strat_systems[["BWH"]])) +
                labs(color = "", x = "False positive rate", y = "")+
                theme_bw() +
                theme(axis.text.y = element_blank(),
                      axis.title.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      panel.grid.major = element_line(colour = "grey85", linewidth = 0.1),
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black"),
                      strip.background = element_rect(fill = "white"),
                      legend.position = c(0.99, 0.01),
                      legend.justification = c("right", "bottom"),
                      legend.box.just = "right",
                      legend.margin = margin(6, 6, 6, 6),
                      legend.title = element_blank())+
                facet_wrap("cohort")
## D-SQUAME
d_squame_roc_entiredataset_p <- ggplot(data = d_squame_roc_entiredataset_df,
                               aes(x = FPR,
                                   y = TPR,
                                   color = strat_system))+
                geom_path(linewidth = 0.8)+
                scale_color_manual(values = c(cols_strat_systems[["SCCore-GEP"]], cols_strat_systems[["BWH"]], cols_strat_systems[["AJCC8"]])) +
                labs(color = "", x = "False positive rate", y = "True positive rate")+
                theme_bw() +
                theme(panel.grid.major = element_line(colour = "grey85", linewidth = 0.1),
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black"),
                      strip.background = element_rect(fill = "white"),
                      legend.position = c(0.99, 0.01),
                      legend.justification = c("right", "bottom"),
                      legend.box.just = "right",
                      legend.margin = margin(6, 6, 6, 6),
                      legend.title = element_blank())+
                facet_wrap("cohort")
roc_p <- d_squame_roc_entiredataset_p + nassir_roc_entiredataset_p
#-------------------------------------------------------------------------------

# Source data that needs to be saved
#-------------------------------------------------------------------------------
source_data_df <- d_squame_roc_entiredataset_df %>%
    rbind(nassir_roc_entiredataset_df) %>%
    mutate(strat_system = sub(",.*$", "", strat_system)) %>%
    rename(Model = strat_system,
           Dataset = cohort,
           Threshold = threshold) %>%
    relocate(Dataset, Model, Threshold)
#-------------------------------------------------------------------------------

# Save results
#-------------------------------------------------------------------------------
write.csv(source_data_df, file.path(output_dir, "supp_figure_13c_entiredataset_performance_roc_curve_source_data.csv"), row.names = F)
pdf(file.path(output_dir, "supp_figure_13c_entiredataset_performance_roc_curve.pdf"), height = 4.5, width = 9.5)
print(roc_p)
dev.off()
#-------------------------------------------------------------------------------