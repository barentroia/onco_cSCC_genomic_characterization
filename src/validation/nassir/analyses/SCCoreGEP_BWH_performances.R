#-------------------------------------------------------------------------------
# Aim: Test SCCore-GEP model on Nassir et al. validation dataset
# Author: L.Pozza, adapted from J. Traets
# Input: D-SQUAME discovery and Nassir et al. TPM & clinical data + SCCore-GEP model
# Output: Performance of SCCore-GEP model in Nassir et al. validation dataset
#-------------------------------------------------------------------------------

# Load libraries
#-------------------------------------------------------------------------------
library(conflicted)
library(openxlsx)
library(tidyverse)
library(survival)
library(glmnet)
library(pROC)
library(sva)
conflicted::conflicts_prefer(dplyr::filter)
#-------------------------------------------------------------------------------

# Filenames and settings
#-------------------------------------------------------------------------------
# Nassir et al. data
## Clinical data
clinical_nassir_fn <- n3_clinical
tpms_nassir_fn <- n1_tpms
## List of bad samples
bad_samples <-  c("SRR31748002","SRR31748007") #only 007 is a tumor

# GEP discovery data
## Clinical data
clinical_gep_discovery_fn <- file.path(results_dir, "intermediate", "clinical", "merged_clinical_data_imputed_d_squame_discovery.csv")
## Expression data
tpms_gep_discovery_fn <- p1_tpms
## List of bad samples
bad_samples_gep_discovery_fn <- p7_badQC

# SCCore-GEP
sccoregep_fn <- file.path(results_dir, "intermediate", "sccore_gep","Discovery_models",
                          "coxnet_Late_feature_integration_robj.rds")
# Time point at which c-index is computed
tp <- 5
# Number of bootstrap sampling to compute 95% CI
n_boots <- 100
# Output folder
output_dir <- file.path(results_dir, "intermediate", "validation", "nassir")
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}
#-------------------------------------------------------------------------------

# Source functions
#-------------------------------------------------------------------------------
source(file.path(code_dir, "functions", "boot_metrics_se_ci.R"))
#-------------------------------------------------------------------------------

# Read data
#-------------------------------------------------------------------------------
# Nassir et al.
clinical_nassir <- read.csv(clinical_nassir_fn)
tpms_nassir <- as.data.frame(read_tsv(tpms_nassir_fn))
# GEP discovery
bad_samples_gep_discovery <- unique(read.csv(bad_samples_gep_discovery_fn)$Set_id)
clinical_gep_discovery <- read.csv(clinical_gep_discovery_fn)
tpms_gep_discovery <- read.csv(tpms_gep_discovery_fn)
# SCCore-GEP
sccoregep <- readRDS(sccoregep_fn)
#-------------------------------------------------------------------------------

# Preprocess Nassir et al. data
#-------------------------------------------------------------------------------
clinical_nassir <- clinical_nassir[!(clinical_nassir$Run %in% bad_samples),]
clinical_nassir <- clinical_nassir %>% filter(Tumor != "CTRL")
clinical_nassir$BWH <- str_extract(clinical_nassir$SampleID, "T[0-9]+[a-z]?")
clinical_nassir[is.na(clinical_nassir$BWH),]$BWH <- "T2b"
tpms_nassir <- tpms_nassir[,c("gene_id","gene_name",clinical_nassir$Run)]
#-------------------------------------------------------------------------------

# Preprocess GEP discovery data
#-------------------------------------------------------------------------------
clinical_gep_discovery <- clinical_gep_discovery %>%
        filter(!Set_id %in% bad_samples_gep_discovery) %>%
        mutate(Biopsy_excision = ifelse(Biopsy_excision == "Biopsy", "Biopsy", "Excision"))
tpms_gep_discovery <- tpms_gep_discovery %>% select(all_of(c("ID", clinical_gep_discovery$SkylineDx.ID)))
#-------------------------------------------------------------------------------

# Combine dataset and perform batch correction (combat)
#-------------------------------------------------------------------------------
# Combine clinical and gene expression datasets
tpms_comb <- tpms_gep_discovery %>%
        inner_join(tpms_nassir, by = c("ID" = "gene_id")) %>%
        column_to_rownames(var = "ID") %>%
        select(-gene_name)
logtpms_comb <- log2(tpms_comb + 1)
clinical_comb <- rbind(clinical_gep_discovery %>% 
                         mutate(Sample_id = SkylineDx.ID,
                                Dataset = "Discovery") %>%
                         select(Sample_id, Metastasis, BWH, Dataset),
                       clinical_nassir %>% 
                         mutate(Sample_id = Run,
                                Metastasis = MetNoMet,
                                Dataset = "Nassir") %>%
                         select(Sample_id, Metastasis, BWH, Dataset))
clinical_comb <- clinical_comb[match(colnames(tpms_comb),clinical_comb$Sample_id),]
clinical_comb$Metastasis <- ifelse(clinical_comb$Metastasis %in% c("Case","Met"),"Case","Control")
# Perform batch correction (combat)
modcombat <- model.matrix(~Metastasis+BWH, data = clinical_comb)
stopifnot(clinical_comb$Sample_id == colnames(logtpms_comb))
batch <- clinical_comb$Dataset
logtpms_comb_corrected <- ComBat(dat = logtpms_comb,
                                mod = modcombat,
                                batch = batch,
                                ref.batch = "Discovery")
# Extract only Nassir et al. batch corrected gene expression
logtpms_comb_corrected <- logtpms_comb_corrected[, clinical_nassir$Run]
#-------------------------------------------------------------------------------

# Test SCCore-GEP
#-------------------------------------------------------------------------------
# Prepare data
## Removing .# behind gene ids, not in feats_to_keep/x_fitted_model
newx_data <- t(logtpms_comb_corrected)
colnames(newx_data) <- sub("\\..*", "", colnames(newx_data))
final_fitted_model_sccoregep <- sccoregep$Final_model2$fitted_model
feats_to_keep_sccoregep <- colnames(sccoregep$Final_model2$x_fitted_model)
newx_data_sccoregep <- as.matrix(newx_data[,feats_to_keep_sccoregep])
# Apply SCCore-GEP
sccore_gep_preds <- 1 - as.vector(summary(survfit(final_fitted_model_sccoregep, 
                                     s = "lambda.1se", 
                                     x = sccoregep$Final_model2$x_fitted_model, 
                                     y = sccoregep$Final_model2$y_fitted_model, 
                                     weights = sccoregep$Final_model2$w_fitted_model, 
                                     newx = as.matrix(newx_data_sccoregep)), times = tp)$surv)
names(sccore_gep_preds) <- rownames(newx_data)
# Save predictions in dataframe
model_predictions_df <- as.data.frame(sccore_gep_preds) %>%
        rownames_to_column(var = "Sample_id") %>%
        `colnames<-`(gsub("_preds", "", colnames(.))) %>%
        left_join(clinical_comb,
                  by = "Sample_id") %>%
        mutate(Metastasis_num = ifelse(Metastasis == "Case", 1, 0),
               bwh = as.numeric(factor(ifelse(BWH == "Unknown", NA, BWH), levels = c("T1", "T2a", "T2b", "T3"))),
               All = "All",
               BWH_staging_T1T2a_T2bT3 = ifelse(BWH %in% c("T1", "T2a"), "T1-T2a",
                                                ifelse(BWH %in% c("T2b", "T3"), "T2b-T3", "Unknown"))) %>%
        select(-Metastasis, -Dataset)
#-------------------------------------------------------------------------------

# Compute discrimination performances in entire dataset (complete cases and max samples)
#-------------------------------------------------------------------------------
# For maximum samples of each model
perf_max_samples <- model_predictions_df %>%
        gather(variable, group, -Sample_id, -Metastasis_num, -sccore_gep, -bwh) %>%
        gather(model, score, -variable, -group, -Sample_id, -Metastasis_num) %>%
        filter(!((model == "bwh" & variable %in% c("BWH")) | is.na(score))) %>%
        group_by(model, variable, group) %>%
        summarise(N = n(),
                  auc = ifelse(length(unique(Metastasis_num)) == 1,
                               NA,
                               round(as.numeric(gsub(".*: ", "", pROC::roc(Metastasis_num, as.numeric(score), quiet = TRUE, direction = "<")$auc)), 4)),
                  auc_ci95 = ifelse(length(unique(Metastasis_num)) == 1,
                                     "(NA-NA)",
                                     boot_metrics_se_ci(cur_data(), "score", "Metastasis_num", NULL, "Weight_rescaled", NULL, n_boots, "auc")$ci95)) %>%
        mutate(auc = paste0(auc, " ", auc_ci95)) %>%
        select(-auc_ci95)
# For complete cases of all models
perf_complete_cases <- model_predictions_df %>%
        filter(!is.na(bwh)) %>%
        gather(variable, group, -Sample_id, -Metastasis_num, -sccore_gep, -bwh) %>%
        gather(model, score, -variable, -group, -Sample_id, -Metastasis_num) %>%
        filter(!((model == "bwh" & variable %in% c("BWH")) | is.na(score))) %>%
        group_by(model, variable, group) %>%
        summarise(N = n(),
                  auc = ifelse(length(unique(Metastasis_num)) == 1,
                               NA,
                               round(as.numeric(gsub(".*: ", "", pROC::roc(Metastasis_num, as.numeric(score), quiet = TRUE, direction = "<")$auc)), 4)),
                  auc_ci95 = ifelse(length(unique(Metastasis_num)) == 1,
                                     "(NA-NA)",
                                     boot_metrics_se_ci(cur_data(), "score", "Metastasis_num", NULL, "Weight_rescaled", NULL, n_boots, "auc")$ci95)) %>%
        mutate(auc = paste0(auc, " ", auc_ci95)) %>%
        select(-auc_ci95)
#-------------------------------------------------------------------------------

# Save results
#-------------------------------------------------------------------------------
write.csv(model_predictions_df %>%
              select(Sample_id, sccore_gep),
          file.path(output_dir, paste0("SCCoreGEP_predictions.csv")),
          row.names = F)
write.csv(perf_max_samples,
          file.path(output_dir, "SCCoreGEP_BWH_performances_max_samples.csv"),
          row.names = F)
write.csv(perf_complete_cases,
          file.path(output_dir, "SCCoreGEP_BWH_performances_complete_cases.csv"),
          row.names = F)
#-------------------------------------------------------------------------------