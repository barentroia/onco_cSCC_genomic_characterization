#-------------------------------------------------------------------------------
# Aim: Perform multiple imputation on D-SQUAME discovery clinical dataset
# Author: L.Pozza
# Input: D-SQUAME discovery clinical dataset
# Output: Multiple imputation object and merged multiple imputed dataset
#-------------------------------------------------------------------------------
# Load libraries
#-------------------------------------------------------------------------------
library(conflicted)
library(tidyverse)
library(VIM)
library(mice)
library(DescTools)
#-------------------------------------------------------------------------------

# Filenames and settings
#-------------------------------------------------------------------------------
# Clinical data
clinical_fn <- p3_clinical
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
clinical_samples <- openxlsx::read.xlsx(clinical_fn, sheet = "Request Lara - Samples (values)")
clinical_tumors <- openxlsx::read.xlsx(clinical_fn, sheet = "Request Lara - Tumor (values)")
# PNI_bin column is repeated, drop one (column 28)
clinical_tumors <- clinical_tumors[, -28]
stopifnot(!any(duplicated(colnames(clinical_tumors))))
clinical <- clinical_samples %>%
  full_join(clinical_tumors,
            by = c("Skyline_ID",
                   "Patient_ID_SKY",
                   "Set_id",
                   c("Sample_id" = "Sample_ID")),
            suffix = c("_samplessheet", "_tumorsheet"))
#-------------------------------------------------------------------------------
  
# Check pattern of missing data
#-------------------------------------------------------------------------------
pdf(file.path(output_dir,"multiple_imputation_missing_data_d_squame_discovery.pdf"), height = 6, width = 20)
aggr(clinical, oma = c(15,5,5,5))
dev.off()
#-------------------------------------------------------------------------------
  
# Direct imputation of variables and transformation into right class
#-------------------------------------------------------------------------------
clinical <- clinical %>%
              mutate(Invasion_of_bones = ifelse(is.na(Invasion_of_bones), "No", Invasion_of_bones), # Assume that if missing, bones invasion is not present
                     PNI_bin = ifelse(is.na(PNI_bin), "No", PNI_bin), # Assume that if missing, PNI is not present
                     Lymphovascular_invasion_bin = ifelse(is.na(Lymphovascular_invasion_bin), "No", Lymphovascular_invasion_bin), # Assume that if missing, LVI is not present
                     PNI_or_LVI = ifelse(is.na(PNI_or_LVI), "No", PNI_or_LVI), # Assume that if missing, PNI or LVI are not present
                     CP_score = 100 * (1 - 0.973 ^ exp(predict_lp(Age, Sex, Number_of_cSCC_before_culprit, Tumor_location, Tumor_diameter, Tissue_involvement, Differentiation, PNI_or_LVI) - 1.86)), # Recompute CP score after changing PNI_or_LVI
                     AJCC_8 = ajcc_staging(Tumor_diameter, PNI_bin, Tissue_involvement, Depth_of_Invasion, Invasion_of_bones), # Recompute AJCC after changing PNI_bin and Invasion_of_bones
                     BWH = bwh_staging(Tumor_diameter, Differentiation, PNI_bin, Tissue_involvement, Invasion_of_bones) # Recompute BWH after changing PNI_bin and Invasion_of_bones
                     )
# Identify variables to exclude
vars_to_exclude <- c("Set_id",
                     "Sample_id",
                     "Path_scoring",
                     "Sample_id", # I guess it's not a predictive variable
                     "Request", # Can the year be predictive? I guess not
                     "Procedure_number", # Not predictive
                     "Time_between_first_sec_proc_samplessheet",
                     "Time_between_first_sec_proc_tumorsheet",
                     "Year_of_obtained_material", # Not predictive
                     "RNA_seq", # Constant
                     "WES_performed", # Constant
                     "Block_received", # Constant
                     "H_E_performed", # Constant
                     "Tumor_purity", # Not predictive
                     "Patient_ID_SKY", # Not predictive
                     "Batch.summary_1.slide", # Not predictive
                     "Batch.summary_28.slides", # Not predictive
                     "Scoring", # Not predictive
                     "Tumor_location_morecats", # Lots of levels --> giving issues
                     "PNI_diameter", # Associated to PNI_bin
                     "Tumor_budding", # ?
                     "Mitotic_rate", # ?
                     "Morphology_subtype", # Lots of levels --> giving issues
                     "Morphology_less_30", # Lots of levels and would need to be uniformized
                     "Notes_subtype_session_Antien", # Not consistent
                     "Notes_pathologists", # Not consistent
                     "Biopsy_scored", # We already use sample type
                     "Doubt_cscc", # Few samples are doubt CSCC, captured by other variables
                     "No_tumor_on_HE" # Few samples have no tumor on HE, captured by other variables
  )
  # Subset dataset
  clinical_sub <- clinical %>%
    select(-all_of(vars_to_exclude))
  clinical_sub <- clinical_sub %>%
    mutate(Biopsy_excision = factor(ifelse(Biopsy_excision == "Biopsy", "Biopsy", "Excision"), levels = c("Excision", "Biopsy")),
           Metastasis = factor(Metastasis, levels = c("Control", "Case")),
           FU_metastasis_years = as.numeric(FU_metastasis_years),
           vitstat2022 = factor(vitstat2022, levels = c("Alive", "Dead")),
           Vitfup_2022_years = as.numeric(Vitfup_2022_years),
           Sex = factor(Sex, levels = c("Female", "Male")),
           HM_at_cSCC = factor(HM_at_cSCC, levels = c("No", "Yes")),
           OTR_at_cSCC = factor(OTR_at_cSCC, levels = c("No", "Yes")),
           Number_of_cSCC_before_culprit = as.numeric(as.character(Number_of_cSCC_before_culprit)),
           Tumor_location = factor(Tumor_location, levels = c("Scalp/neck", "Trunk/Extremities", "Face")),
           AJCC_8 = factor(AJCC_8, levels = c("T1", "T2", "T3", "T4")),
           BWH = factor(BWH, levels = c("T1", "T2a", "T2b", "T3")),
           CP_score = as.numeric(CP_score),
           Tumor_diameter = as.numeric(Tumor_diameter),
           Tumor_width = as.numeric(Tumor_width),
           Tissue_involvement = factor(Tissue_involvement, levels = c("Dermis", "Subcutaneous fat", "Beyond subcutaneous fat")),
           Differentiation = factor(Differentiation, levels = c("Good/moderate", "Poor/undifferentiated")),
           PNI_bin = factor(PNI_bin, levels = c("No", "Yes")),
           Lymphovascular_invasion_bin = factor(Lymphovascular_invasion_bin, levels = c("No", "Yes")),
           PNI_or_LVI = factor(PNI_or_LVI, levels = c("No", "Yes")),
           Depth_of_Invasion = as.numeric(Depth_of_Invasion),
           Breslow_thickness = as.numeric(Breslow_thickness),
           Resection_margin_cat = factor(Resection_margin_cat, levels = c("R0 (Radical)", "R1/R2 (irradicaal)")),
           Invasion_of_bones = factor(Invasion_of_bones, levels = c("No", "Yes")),
           Solar_elastosis = factor(Solar_elastosis, levels = c("Absent/moderate", "Extensive")),
           Peritumoral_infiltration = factor(Peritumoral_infiltration, levels = c("Absent/moderate", "Abundant")))
#-------------------------------------------------------------------------------

# Check pattern of missing data
#-------------------------------------------------------------------------------
pdf(file.path(output_dir, "multiple_imputation_missing_data_sel_vars_d_squame_discovery.pdf"), height = 6, width = 20)
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
pred[, "Skyline_ID"] <- 0
pred[, "Invasion_of_bones"] <- 0 # Bone invasion recording is unreliable, so it will not be used
pred[c("PNI_bin", "Lymphovascular_invasion_bin"), "PNI_or_LVI"] <- 0
pred[c("Age", "Sex", "Number_of_cSCC_before_culprit", "Tumor_location", "Tumor_diameter",
       "Tissue_involvement", "Differentiation", "PNI_or_LVI", "PNI_bin", "Lymphovascular_invasion_bin"), "CP_score"] <- 0
pred[c("Tumor_diameter", "PNI_bin", "Tissue_involvement", "Depth_of_Invasion", "Invasion_of_bones"), "AJCC_8"] <- 0
pred[c("Tumor_diameter", "Differentiation", "PNI_bin", "Tissue_involvement", "Invasion_of_bones"), "BWH"] <- 0
  
# Adjust method, for derived variables
## PNI or LVI
init_df$meth["PNI_or_LVI"] <- "~ I(ifelse(PNI_bin == 'Yes' | Lymphovascular_invasion_bin =='Yes', 'Yes','No'))"
## CP score
init_df$meth["CP_score"] <- "~ I(100 * (1 - 0.973 ^ exp(predict_lp(
  Age, Sex, Number_of_cSCC_before_culprit, Tumor_location, Tumor_diameter,
  Tissue_involvement, Differentiation, PNI_or_LVI) - 1.86)))"
## BWH
init_df$meth["BWH"] <- "~ I(bwh_staging(Tumor_diameter, Differentiation, PNI_bin, Tissue_involvement, Invasion_of_bones))"
## AJCC
init_df$meth["AJCC_8"] <- "~ I(ajcc_staging(Tumor_diameter, PNI_bin, Tissue_involvement, Depth_of_Invasion, Invasion_of_bones))"
 
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

# Extract imputed rows and for each varaible get median (numeric) or mode (categorical)
imp_vars <- c("Tumor_location", "Tumor_diameter", "Tumor_width", "Tissue_involvement",
              "Differentiation", "Depth_of_Invasion", "Breslow_thickness", "Resection_margin_cat",
              "Solar_elastosis", "Peritumoral_infiltration")
imputed_ds <- sapply(1:n_imp,
                     function(i) complete(mimp, i) %>% mutate(n_imp = i),
                     simplify = F) %>%
  bind_rows() %>%
  select(all_of(c("Skyline_ID", imp_vars))) %>%
  group_by(Skyline_ID) %>%
  reframe(Tumor_location = Mode(Tumor_location),
          Tumor_diameter = median(Tumor_diameter),
          Tumor_width = median(Tumor_width),
          Tissue_involvement = Mode(Tissue_involvement),
          Differentiation = Mode(Differentiation),
          Depth_of_Invasion = median(Depth_of_Invasion),
          Breslow_thickness = median(Breslow_thickness),
          Resection_margin_cat = Mode(Resection_margin_cat),
          Solar_elastosis = Mode(Solar_elastosis),
          Peritumoral_infiltration = Mode(Peritumoral_infiltration)) %>%
   full_join(clinical %>%
              select(all_of(setdiff(colnames(clinical), imp_vars)))%>%
              dplyr::rename(CP_score_orig = CP_score,
                            AJCC_8_orig = AJCC_8,
                            BWH_orig = BWH),
            by = "Skyline_ID") %>%
  mutate(CP_score = 100 * (1 - 0.973 ^ exp(predict_lp(
            Age, Sex, Number_of_cSCC_before_culprit, Tumor_location, Tumor_diameter,
            Tissue_involvement, Differentiation, PNI_or_LVI) - 1.86)),
         AJCC_8 = ajcc_staging(Tumor_diameter, PNI_bin, Tissue_involvement, Depth_of_Invasion, Invasion_of_bones),
         BWH = bwh_staging(Tumor_diameter, Differentiation, PNI_bin, Tissue_involvement, Invasion_of_bones),
         CP_score_check = round(CP_score, 2) == round(CP_score_orig, 2),
         AJCC_8_check = AJCC_8 == AJCC_8_orig,
         BWH_check = BWH == BWH_orig)
stopifnot(all(imputed_ds$CP_score_check[!is.na(imputed_ds$CP_score_check)]))
stopifnot(all(imputed_ds$AJCC_8_check[!is.na(imputed_ds$AJCC_8_check)]))
stopifnot(all(imputed_ds$BWH_check[!is.na(imputed_ds$BWH_check)]))
imputed_ds <- imputed_ds %>%
  select(-CP_score_check, -CP_score_orig, -AJCC_8_check, -AJCC_8_orig, -BWH_check, -BWH_orig) %>%
  arrange(Set_id) %>%
  rename("SkylineDx.ID" = "Skyline_ID")
stopifnot(!any(duplicated(imputed_ds$SkylineDx.ID)))
#-------------------------------------------------------------------------------
  
# Save results
#-------------------------------------------------------------------------------
# Save logged events
write.csv(mimp$loggedEvents %>% as.data.frame(),
          file = file.path(output_dir, "multiple_imputation_logged_events_d_squame_discovery.log"),
          row.names = F)
# Plot convergence plots
pdf(file.path(output_dir, "multiple_imputation_convergence_plots_d_squame_discovery.pdf"), width = 16, height = 16)
print(plot(mimp, layout = c(4, 10)))
dev.off()
# Plot density plots
pdf(file.path(output_dir, "multiple_imputation_density_plots_d_squame_discovery.pdf"))
print(densityplot(mimp))
dev.off()
# Plot proportions plots
pdf(file.path(output_dir, "multiple_imputation_proportion_plots_d_squame_discovery.pdf"), width = 19)
print(propplot(mimp))
dev.off()
# Save multiple imputation dataset
save(mimp, file = file.path(output_dir, "clinical_data_imputed_d_squame_discovery.Rdata"))
# Save merged imputed dataset
write.csv(imputed_ds, file = file.path(output_dir, "merged_clinical_data_imputed_d_squame_discovery.csv"), row.names = F)
#-------------------------------------------------------------------------------
