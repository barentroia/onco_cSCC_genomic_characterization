#-------------------------------------------------------------------------------
# Aim: Perform multiple imputation on D-SQUAME validation clinical dataset
# Author: L.Pozza
# Input: D-SQUAME validation clinical dataset
# Output: Multiple imputation object and merged multiple imputed dataset
#-------------------------------------------------------------------------------

# Load libraries
#-------------------------------------------------------------------------------
library(conflicted)
library(tidyverse)
library(mice)
library(VIM)
library(DescTools)
#-------------------------------------------------------------------------------

# Filenames and settings
#-------------------------------------------------------------------------------
# Clinical data
clinical_fn <- e3_clinical
# Output folder
output_dir <- file.path(results_dir, "intermediate", "clinical")
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}
#-------------------------------------------------------------------------------

# Source functions
#-------------------------------------------------------------------------------
source(file.path(code_dir, "functions", "stagings_emcmodel_funs.R"))
source(file.path(code_dir, "functions", "propplot_fun.R"))
#-------------------------------------------------------------------------------

# Read data and prepare output folder
#-------------------------------------------------------------------------------
clinical <- read.csv(clinical_fn)
clinical <- clinical %>%
  mutate(Type_of_material_scored_bin = ifelse(Type_of_material_scored == "Biopsy", "Biopsy", "Excision"))
#-------------------------------------------------------------------------------

# Check pattern of missing data
#-------------------------------------------------------------------------------
pdf(file.path(output_dir, paste0("multiple_imputation_missing_data_d_squame_validation.pdf")), height = 6, width = 20)
aggr(clinical, oma = c(15,5,5,5))
dev.off()
#-------------------------------------------------------------------------------

# Direct imputation of variables and transformation into right class
#-------------------------------------------------------------------------------
# Identify variables to exclude
vars_to_exclude <- c("Skyline_ID", # Not predictive
                     "Sample_id", # Not predictive
                     "Set_id", # Not predictive
                     "random_set", # Not predictive
                     "Type_of_material", # We use Type_of_material_scored_bin
                     "Block_received", # Not predictive
                     "Batchnumber", # Not predictive
                     "Excluded_block", # Not predictive
                     "Excluded_path", # Not predictive
                     "Seq_info", # Not predictive
                     "RNA_seq", # Not predictive
                     "control_type", # Contains some information, but probably we can get them from other variables too
                     "H_E_performed", # Not predictive
                     "Time_between_first_sec_proc_samplessheet", # Not predictive
                     "Path_scoring", # Not predictive
                     "Year", # Not predictive
                     "Tumor_purity", # Not predictive,
                     "First_50_sets", # Not predictive
                     "LAB_ID_SKY", # Not predictive
                     "first_procedure", # We use Type_of_material_scored_bin
                     "Type_of_material_scored", # We use Type_of_material_scored_bin
                     "Time_between_first_sec_proc_tumorssheet", # Not predictive
                     "Scored", # Not predictive
                     "localisation_nm", # Contains some information about location, but location is missing for we samples
                     "metastasis_baseline", # Contains some information, but we've aslo seen that these metastases are not that different form others
                     "Differntiation_broder", # We use differentiation (same as done in discovery)
                     "PNI_or_LVI_in_report", # We use PNI_bin and PNI_or_LVI
                     "PNI_scored", # We use PNI_bin and PNI_or_LVI
                     "PNI_diameter", # Contains limited information
                     "PNI_tissue_involvement", # Contains limited information
                     "LVI_scored", # We use PNI_or_LVI
                     "Morphological_subtype", # Contains information, but it's a variable difficult to score and has too many levels
                     "Differentation_worstpattern", # We use differentiation (same as done in discovery)
                     "Doubt_CSCC", # Contains some information, but limited and these doubt samples look similar to others
                     "Vitfup_metastasis_years_matchvar", # We use FU_metastasis_years
                     "Bad_seq_QC", # Not predictive
                     "Set_info", # Not predictive
                     "Type_of_material_bin" # We use Type_of_material_scored_bin
)
# Subset dataset
clinical_sub <- clinical %>%
  select(-all_of(vars_to_exclude))
clinical_sub <- clinical_sub %>%
  mutate(Type_of_material_scored_bin = factor(Type_of_material_scored_bin, levels = c("Excision", "Biopsy")),
         Metastasis = factor(Metastasis, levels = c("Control", "Case")),
         FU_metastasis_years = as.numeric(FU_metastasis_years),
         Vitstat2023 = factor(Vitstat2023, levels = c("Alive", "Dead")),
         Vitfup_2023 = as.numeric(Vitfup_2023),
         Sex = factor(Sex, levels = c("Female", "Male")),
         Immunosuppressed = factor(ifelse(HM_at_cSCC == "Yes" | OTR_at_cSCC == "Yes", "Yes", "No"), levels = c("No", "Yes")),
         HM_at_cSCC = factor(HM_at_cSCC, levels = c("No", "Yes")),
         OTR_at_cSCC = factor(OTR_at_cSCC, levels = c("No", "Yes")),
         Number_of_cSCC_before_culprit = as.numeric(as.character(Number_of_cSCC_before_culprit)),
         Tumor_location = factor(Tumor_location, levels = c("Scalp/neck", "Trunk/extremities", "Face")),
         AJCC_staging = factor(AJCC_staging, levels = c("T1", "T2", "T3", "T4")),
         BWH_staging = factor(BWH_staging, levels = c("T1", "T2a", "T2b", "T3")),
         CP_risk = as.numeric(CP_risk),
         Tumor_diameter = as.numeric(Tumor_diameter),
         Tissue_involvement = factor(Tissue_involvement, levels = c("Dermis", "Subcutaneous fat", "Beyond subcutaneous fat")),
         Differentiation = factor(Differentiation, levels = c("Good/moderate", "Poor/undifferentiated")),
         PNI_bin = factor(PNI_bin, levels = c("No", "Yes")),
         PNI_or_LVI = factor(PNI_or_LVI, levels = c("No", "Yes")),
         Depth_of_invasion = as.numeric(Depth_of_invasion),
         Breslow_thickness = as.numeric(Breslow_thickness),
         Resection_margin_excision = factor(Resection_margin_excision, levels = c("Complete", "Incomplete")),
         Bone_invasion = factor(Bone_invasion, levels = c("No", "Yes")),
         Solar_elastosis = factor(Solar_elastosis, levels = c("No/nihil", "Moderate", "Extensive")),
         Peritumoral_infiltrate = factor(Peritumoral_infiltrate, levels = c("Absent/mild", "Moderate", "Abundant"))
  ) %>%
  droplevels() %>%
  unique()
#-------------------------------------------------------------------------------

# Check pattern of missing data
#-------------------------------------------------------------------------------
pdf(file.path(output_dir, paste0("multiple_imputation_missing_data_sel_vars_d_squame_validation.pdf")), height = 6, width = 20)
aggr(clinical_sub, oma = c(15,5,5,5))
dev.off()
#-------------------------------------------------------------------------------

# Multiple imputation
#-------------------------------------------------------------------------------
# Initialization of multiple imputation
## As method, use the default options:
## - pmm: for numerical variables
## - logreg: for binary variables
## - polyreg: for categorical variables
init_df <- mice(clinical_sub, maxit = 0, print = F)
print(init_df$loggedEvents %>% as.data.frame())

# Initial prediction matrix
pred <- init_df$pred

# Removal of variables not used in imputation of other variables
pred[, "Patient_ID_SKY"] <- 0
pred[, "Bone_invasion"] <- 0 # Bone invasion recording is unreliable, so it will not be used
pred[, "HM_at_cSCC"] <- 0 # Only 1 sample in the validation set has HM --> we decided to use Immunosuppressed
pred[, "OTR_at_cSCC"] <- 0 # we decided to use Immunosuppressed
pred[c("Age", "Sex", "Number_of_cSCC_before_culprit", "Tumor_location", "Tumor_diameter",
       "Tissue_involvement", "Differentiation", "PNI_or_LVI", "PNI_bin"), "CP_risk"] <- 0
pred[c("Tumor_diameter", "PNI_bin", "Tissue_involvement", "Depth_of_invasion", "Bone_invasion"), "AJCC_staging"] <- 0
pred[c("Tumor_diameter", "Differentiation", "PNI_bin", "Tissue_involvement", "Bone_invasion"), "BWH_staging"] <- 0
pred[c("Breslow_thickness", "Solar_elastosis", "Peritumoral_infiltrate"), "BWH_staging"] <- 0 # BWH is usually missing in case of missing values in BT, SE, PI and gives issues

# Adjust method, for derived variables
## CP score
init_df$meth["CP_risk"] <- "~ I(100 * (1 - 0.973 ^ exp(predict_lp(
  Age, Sex, Number_of_cSCC_before_culprit, Tumor_location, Tumor_diameter,
  Tissue_involvement, Differentiation, PNI_or_LVI) - 1.86)))"
## BWH_staging
init_df$meth["BWH_staging"] <- "~ I(bwh_staging(Tumor_diameter, Differentiation, PNI_bin, Tissue_involvement, Bone_invasion))"
## AJCC
init_df$meth["AJCC_staging"] <- "~ I(ajcc_staging(Tumor_diameter, PNI_bin, Tissue_involvement, Depth_of_invasion, Bone_invasion))"

# Run multiple imputation
n_imp <- 25
mimp <- mice(clinical_sub,
             m = n_imp,
             seed = 567,
             maxit = 40,
             method = init_df$meth,
             pred = pred,
             vis = "monotone",
             print = F)
#-------------------------------------------------------------------------------

# Merge multiple imputed dataset in a single dataset
#-------------------------------------------------------------------------------
merged_mimp <- sjmisc::merge_imputations(complete(mimp, 0), mimp, summary = "none")
# Adjust stagings and CP risk
merged_mimp <- merged_mimp %>%
  cbind(clinical_sub %>% select(-all_of(colnames(merged_mimp)))) %>%
  mutate(CP_risk = 100 * (1 - 0.973 ^ exp(predict_lp(Age, Sex, Number_of_cSCC_before_culprit, Tumor_location, Tumor_diameter,
                                                     Tissue_involvement, Differentiation, PNI_or_LVI) - 1.86)),
         BWH_staging = bwh_staging(Tumor_diameter, Differentiation, PNI_bin, Tissue_involvement, Bone_invasion),
         AJCC_staging = ajcc_staging(Tumor_diameter, PNI_bin, Tissue_involvement, Depth_of_invasion, Bone_invasion)) %>%
  select(all_of(colnames(clinical_sub))) %>%
  relocate(Patient_ID_SKY)
#-------------------------------------------------------------------------------

# Save results
#-------------------------------------------------------------------------------
# Save logged events
write.csv(mimp$loggedEvents %>% as.data.frame(),
          file = file.path(output_dir, paste0("multiple_imputation_logged_events_d_squame_validation.log")),
          row.names = F)
# Plot convergence plots
pdf(file.path(output_dir, paste0("multiple_imputation_convergence_plots_d_squame_validation.pdf")), width = 16, height = 16)
print(plot(mimp, layout = c(4, 10)))
dev.off()
# Plot density plots
pdf(file.path(output_dir, paste0("multiple_imputation_density_plots_d_squame_validation.pdf")))
print(densityplot(mimp))
dev.off()
# Plot proportions plots
pdf(file.path(output_dir, paste0("multiple_imputation_proportion_plots_d_squame_validation.pdf")), width = 19)
print(propplot(mimp))
dev.off()
# Save multiple imputation dataset
save(mimp, file = file.path(output_dir, paste0("clinical_data_imputed_d_squame_validation.Rdata")))
# Save merged multiple imputation dataset
write.csv(merged_mimp,
          file = file.path(output_dir, paste0("merged_clinical_data_imputed_d_squame_validation.csv")),
          row.names = F)
#-------------------------------------------------------------------------------
