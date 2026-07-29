#-------------------------------------------------------------------------------
# Aim: Compute weights for patients in D-SQUAME validation dataset
# Author: L. Pozza, functions to compute weights were coded by B. Rentroia-Pacheco
# Input:
# - Dataset with clinical data from the samples in D-SQUAME validation full dataset
# - Samples IDs of samples in D-SQUAME NCC dataset
# Output: Sampling probabilities and weights for samples in D-SQUAME validation dataset
#-------------------------------------------------------------------------------

# Load libraries
#-------------------------------------------------------------------------------
library(conflicted)
library(tidyverse)
#-------------------------------------------------------------------------------

# Functions
#-------------------------------------------------------------------------------
# Sampling probability based on all records
compute_sampling_prob_based_on_all_records <- function(df_all_records_controls, sampled_records_bin, matching_vars = c("Vitfup2023", "Type_of_material1")){
  # Logistic regression to model sampling probability of each record
  fit <- glm(sampled_records_bin~., family = "binomial", data = df_all_records_controls[, matching_vars])
  
  # Get record sampling probabilities
  rec_sampling_probs <- predict(fit, type = "response", newdata = df_all_records_controls)
  df_predictions <- data.frame("administratienummer" = df_all_records_controls$administratienummer,
                               "Samp_prob" = rec_sampling_probs)
  
  # Aggregate sampling probabilities per patient
  df_predictions_per_patient <- df_predictions %>%
    mutate(Samp_prob = 1 - Samp_prob) %>%
    group_by(administratienummer)%>%
    summarise(Samp_prob = 1 - prod(Samp_prob)) %>%
    ungroup() %>%
    as.data.frame()
  
  # Compute weights
  df_predictions_per_patient$Weight <- 1/df_predictions_per_patient$Samp_prob
  weight_ncc_sum <- sum(df_predictions_per_patient$Weight[df_predictions_per_patient$administratienummer %in% df_all_records_controls$administratienummer[df_all_records_controls$Sampled==1]])
  
  # Rescale weights
  df_predictions_per_patient$Weight_rescaled <- df_predictions_per_patient$Weight*nrow(df_predictions_per_patient)/weight_ncc_sum
  return(df_predictions_per_patient)
}
#-------------------------------------------------------------------------------

# Filenames and settings
#-------------------------------------------------------------------------------
# Project directory
project_dir <- file.path("W:", "Data", "ONDERZOEK", "Studies", "2021",
                         "EMCD21054 [StepIdent_CP_G_refined]", "2-Data",
                         "Analyses", "Validation cohort", "Weighted metrics")
# Output directory
output_dir <-  file.path(project_dir, "output_weights_calculation", "random_controls")
# Filename of full D-SQUAME validation dataset
full_ds_fn <- file.path(project_dir, "full_dataset_WP3_v2.csv")
# Filename with NCC samples
ncc_ds_fn <- file.path(project_dir, "fromSkyline", "parsed_clinical_data_v3lp_valset_complete_sets_red.csv")
#-------------------------------------------------------------------------------

# Create output directories
#-------------------------------------------------------------------------------
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}
#-------------------------------------------------------------------------------

# Read data
#-------------------------------------------------------------------------------
full_ds <- read.csv(full_ds_fn) %>% select(-X)
ncc_ds <- read.csv(ncc_ds_fn)
# In NCC dataset, keep only biopsy samples for patients with double sample type
ncc_ds_f <- ncc_ds %>%
  mutate(selection_priority = ifelse(Type_of_material == "Biopsy", 1, 2)) %>%
  group_by(Patient_ID_SKY) %>% 
  arrange(selection_priority, by_group = T) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(-selection_priority)
print(paste0("There are ", nrow(ncc_ds_f), " patients and samples."))
#-------------------------------------------------------------------------------

# First checks and separate cases and controls
#-------------------------------------------------------------------------------
# Check that all samples in the NCC dataset are in the full dataset
samples_ncc_not_in_full <- setdiff(ncc_ds_f$Sample_id, full_ds$Sample_id)
if (length(samples_ncc_not_in_full) > 0){
  stop(paste0("Sample(s) ", paste(samples_ncc_not_in_full, collapse = ", "), " not in full dataset."))
}
# Check that there are not duplicated sample IDs
duplicated_samples <- full_ds$Sample_id[duplicated(full_ds$Sample_id) & !is.na(full_ds$Sample_id) & !full_ds$Sample_id %in% c("Vx", "VX")]
if (length(duplicated_samples) > 0){
  warning(paste0("Sample(s) ", paste(duplicated_samples, collapse = ", "), " is/are duplicated."))
}
# Full cohort
print(paste0("In the full cohort, there are ", length(unique(full_ds$administratienummer)),
             " patients and ", nrow(full_ds), " records."))
# Cases: nr of total cases should be 279 (cases that were requested)
pats_cases <- full_ds %>%
  dplyr::filter(meta_patient == "yes") %>%
  pull(administratienummer) %>%
  unique()
print(paste0("There are ", length(pats_cases), " cases (patients number)."))
cases_ds <- full_ds %>%
  dplyr::filter(administratienummer %in% pats_cases)
print(paste0("There are ", nrow(cases_ds), " case records."))
cases_not_in_ncc <- setdiff(ncc_ds_f %>% dplyr::filter(Metastasis == "Case") %>% pull(Sample_id), cases_ds$Sample_id)
if (length(cases_not_in_ncc) > 0){
  stop(paste0("Sample(s) ", paste(cases_not_in_ncc, collapse = ", "), " not in cases dataset."))
}
# Controls
pats_controls <- full_ds %>%
  dplyr::filter(meta_patient != "yes") %>%
  pull(administratienummer) %>%
  unique()
print(paste0("There are ", length(pats_controls), " controls (patients number)."))
controls_ds <- full_ds %>%
  dplyr::filter(administratienummer %in% pats_controls)
print(paste0("There are ", nrow(controls_ds), " control records."))
controls_not_in_ncc <- setdiff(ncc_ds_f %>% dplyr::filter(Metastasis == "Control") %>% pull(Sample_id), controls_ds$Sample_id)
if (length(controls_not_in_ncc) > 0){
  stop(paste0("Sample(s) ", paste(controls_not_in_ncc, collapse = ", "), " not in controls dataset."))
}
#-------------------------------------------------------------------------------

# Preprocess cases dataset and perform some checks
#-------------------------------------------------------------------------------
# Check that information is complete
pats_cases_wo_fu_or_sid <- cases_ds %>%
  group_by(administratienummer)  %>%
  mutate(at_least_one_sample_id = any(!is.na(Sample_id)),
         at_least_one_met_fu = any(!is.na(followup_untilmeta))) %>%
  ungroup() %>%
  mutate(missing_sampleid_in_sample_with_met_fu = is.na(Sample_id) & !is.na(followup_untilmeta)) %>%
  group_by(administratienummer) %>%
  mutate(all_missing_sampleid_in_sample_with_met_fu = all(missing_sampleid_in_sample_with_met_fu | is.na(followup_untilmeta))) %>%
  dplyr::filter(!at_least_one_sample_id | !at_least_one_met_fu | all_missing_sampleid_in_sample_with_met_fu) %>%
  pull(administratienummer) %>%
  unique()
if (length(pats_cases_wo_fu_or_sid) > 0){
  stop(paste0("Patient(s) ", paste(as.character(pats_cases_wo_fu_or_sid), collapse = ", "),
              " don't have any metastasis FU or any Sample_id
              or have metastasis FU but Sample_id = NA."))
}
# There are multiple tumors per patient
# Select the one that metastasize (culprit) based on metastasis_case column
cases_ds_f <- cases_ds %>%
  dplyr::filter(metastasis_case == 1)
# There are multiple samples per tumor, due to different procedures.
# Sample type is chosen taking into account the following priority: 
# Excision/Reexcision/Resectie > Mohs > Biopsy
# In case of multiple entries with same sample type, the first one will be kept (based on palgaexcerptid).
## Check annotation of sample type
if (!all(unique(cases_ds_f$Type_of_material1) %in%
         c("Excision/Reexcision/Resectie", "Mohs", "Biopsy"))){
  stop("Type of material in dataset is different from expected. Adjust code accordingly.")
}
## Check that there are maximum 2 entries for each tumor
pats_cases_more_2_entries_tumor <- cases_ds_f %>%
  group_by(administratienummer) %>%
  summarise(n = n()) %>%
  dplyr::filter(n > 2) %>%
  pull(administratienummer)
if (length(pats_cases_more_2_entries_tumor) > 0){
  stop(paste0("Patient(s) ", paste(as.character(pats_cases_more_2_entries_tumor), collapse = ", "),
              " have more than 2 entries per culprit CSCC."))
}
## Do selection based on sample type, palgaexcerptid
## Add column with IDs of both excision and biopsy
cases_ds_f <- cases_ds_f %>%
  mutate(selection_priority = ifelse(Type_of_material1 == "Excision/Reexcision/Resectie", 1,
                                            ifelse(Type_of_material1 == "Mohs", 2, 3)),
         palgaexcerptid_order = sapply(strsplit(palgaexcerptid, "[-]"), '[[', 2)) %>%
  group_by(administratienummer) %>%
  mutate(Sample_ids = paste(Sample_id, collapse = ","),
         N = n()) %>% 
  select(-N) %>%
  arrange(selection_priority, palgaexcerptid_order, by_group = T) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(Sample_id1 = sapply(strsplit(Sample_ids, "[,]"), "[[", 1),
         Sample_id1 = ifelse(Sample_id1 == "NA", NA, Sample_id1),
         Sample_id2 = sapply(Sample_ids, function(x){
           s <- unlist(strsplit(x, "[,]"))
           ifelse(length(s) > 1, s[[2]], NA)
         }),
         Sample_id2 = ifelse(Sample_id2 == "NA", NA, Sample_id2),
         Sample_id_add = ifelse(is.na(Sample_id) & !is.na(Sample_id1), Sample_id1,
                                ifelse(is.na(Sample_id) & !is.na(Sample_id2), Sample_id2,
                                       ifelse(Sample_id == Sample_id1, Sample_id2, Sample_id1)))) %>%
  select(-selection_priority, -palgaexcerptid_order, -Sample_id1, -Sample_id2)
# Check that all sample IDs in NCC dataset are also in cases dataset
cases_not_in_ncc <- setdiff(ncc_ds_f %>%
                           dplyr::filter(Metastasis == "Case") %>%
                           pull(Sample_id),
                         c(cases_ds_f$Sample_id, cases_ds_f$Sample_id_add))
if (length(cases_not_in_ncc) > 0){
  stop(paste0("Sample(s) ", paste(cases_not_in_ncc, collapse = ", "), " not in cases dataset."))
}
#-------------------------------------------------------------------------------

# Preprocess controls dataset and perform some checks
#-------------------------------------------------------------------------------
# Remove tumors that have 0 follow-up, they are not part of the study,
# and cannot be considered for matching
controls_ds_f <- controls_ds %>% dplyr::filter(Vitfup2023 > 0)
# Some patients have multiple entries, there can be 2 reasons why this happens:
# - 1. patient has multiple CSCCs (= multiple tumors)
# - 2. biopsy/excision entries for a same tumor
# In the second case, for weights computation purpose we just want to keep one entry
# as they belong to the same tumor.
# Sample types is chose taking into account the following priority: 
# Excision/Reexcision/Resectie > Mohs > Biopsy
# In case of multiple entries with same sample type, the first one will be kept (based on palgaexcerptid).
## Identify patients with multiple tumors
multiple_entries_per_tumor_controls_ds_f <- controls_ds_f %>%
  group_by(volg_id) %>%
  mutate(n_entries_per_tumor = n()) %>%
  ungroup() %>%
  dplyr::filter(n_entries_per_tumor > 1)
## Check that for the same tumor follow-up is the same
pats_multiple_entries_diff_fu_per_tumor_controls_ds_f <- multiple_entries_per_tumor_controls_ds_f %>%
  select(administratienummer, SCC_volgnummer, volg_id, Vitstat2023, n_entries_per_tumor) %>%
  unique() %>%
  group_by(volg_id) %>%
  mutate(n_diff_fu = n()) %>%
  dplyr::filter(n_diff_fu > 1) %>%
  pull(administratienummer)
if (length(pats_multiple_entries_diff_fu_per_tumor_controls_ds_f) > 0){
  stop(paste0("Patient(s) ", paste(pats_multiple_entries_diff_fu_per_tumor_controls_ds_f, collapse = ", "),
              " have multiple entries for same tumor with different follow-up times."))
}
## Select entry based on sample type
if (!all(unique(multiple_entries_per_tumor_controls_ds_f$Type_of_material1) %in%
         c("Excision/Reexcision/Resectie", "Mohs", "Biopsy"))){
  stop("Type of material in dataset is different from expected. Adjust code accordingly.")
}
### Identify samples to keep
samples_to_keep_multiple_entries_per_tumor_controls <- multiple_entries_per_tumor_controls_ds_f %>%
  mutate(selection_priority = ifelse(Type_of_material1 == "Excision/Reexcision/Resectie", 1,
                                     ifelse(Type_of_material1 == "Mohs", 2, 3)),
         palgaexcerptid_order = sapply(strsplit(palgaexcerptid, "[-]"), '[[', 2)) %>%
  group_by(volg_id) %>%
  arrange(selection_priority, palgaexcerptid_order, by_group = T) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  pull(palgaexcerptid)
### Identify samples to remove
samples_to_remove_multiple_entries_per_tumor_controls <- setdiff(
  multiple_entries_per_tumor_controls_ds_f$palgaexcerptid,
  samples_to_keep_multiple_entries_per_tumor_controls)
stopifnot(length(samples_to_keep_multiple_entries_per_tumor_controls) +
            length(samples_to_remove_multiple_entries_per_tumor_controls) ==
            nrow(multiple_entries_per_tumor_controls_ds_f))
### Store all Sample_id in a column, and remove samples identified above
controls_ds_f <- controls_ds_f  %>%
  group_by(volg_id) %>%
  mutate(Sample_ids = paste(Sample_id, collapse = ",")) %>%
  ungroup() %>%
  mutate(Sample_id1 = sapply(strsplit(Sample_ids, "[,]"), "[[", 1),
         Sample_id1 = ifelse(Sample_id1 == "NA", NA, Sample_id1),
         Sample_id2 = sapply(Sample_ids, function(x){
           s <- unlist(strsplit(x, "[,]"))
           ifelse(length(s) > 1, s[[2]], NA)
         }),
         Sample_id2 = ifelse(Sample_id2 == "NA", NA, Sample_id2),
         Sample_id_add = ifelse(is.na(Sample_id) & !is.na(Sample_id1), Sample_id1,
                                ifelse(is.na(Sample_id) & !is.na(Sample_id2), Sample_id2,
                                       ifelse(Sample_id == Sample_id1, Sample_id2, Sample_id1)))) %>%
  select(-Sample_id1, -Sample_id2) %>%
  dplyr::filter(!palgaexcerptid %in% samples_to_remove_multiple_entries_per_tumor_controls)

# For LR weights: Add variable (Sampled) to indicate if samples are in NCC cohort or not
# Make type of procedure a categorical variable
controls_ds_f <- controls_ds_f %>%
  mutate(Sampled = ifelse(Sample_id %in% ncc_ds_f$Sample_id | Sample_id_add %in% ncc_ds_f$Sample_id, 1, 0),
         first_procedure = as.factor(first_procedure))
if (!all(ncc_ds_f %>% dplyr::filter(Metastasis == "Control") %>% pull(Sample_id) %in%
    unique(c(controls_ds_f %>% dplyr::filter(Sampled == 1) %>% pull(Sample_id),
       controls_ds_f %>% dplyr::filter(Sampled == 1) %>% pull(Sample_id_add)))) |
    sum(controls_ds_f$Sampled) != nrow(ncc_ds_f %>% dplyr::filter(Metastasis == "Control"))){
  stop("Not all NCC samples are in the sampled controls in controls dataset
       or number of sampled controls doesn't match.")
}
#-------------------------------------------------------------------------------

# Compute weights
#-------------------------------------------------------------------------------
# Compute weights for cases
## Number of cases in full cohort
tot_cases <- length(pats_cases)
## Compute sampling probabilities and weights
weights_cases_df <- ncc_ds_f %>%
  dplyr::filter(Metastasis == "Case") %>%
  mutate(n_cases = n(),
         Samp_prob = n_cases / tot_cases,
         Weight = tot_cases / n_cases,
         Weight_rescaled = tot_cases / n_cases) %>%
  left_join(cases_ds %>% select(Sample_id, administratienummer) %>% unique(),
            by = "Sample_id") %>%
  select(-n_cases) %>%
  relocate(administratienummer, Sample_id, Patient_ID_SKY, Metastasis)
## Get weights for all samples
weights_cases_ncc_df <- weights_cases_df %>%
  select(Patient_ID_SKY, Metastasis, first_procedure, FU_metastasis_years,
         Vitfup_metastasis_years_matchvar, Samp_prob, Weight, Weight_rescaled) %>%
  right_join(ncc_ds %>%
               dplyr::filter(Metastasis == "Case"),
             by = c("Patient_ID_SKY", "Metastasis", "FU_metastasis_years",
                    "Vitfup_metastasis_years_matchvar", "first_procedure")) %>%
  left_join(full_ds %>%
              select(administratienummer, Sample_id) %>%
              unique(),
            by = "Sample_id") %>%
  relocate(administratienummer, Sample_id, Patient_ID_SKY, Metastasis)

# Compute LR weights for controls
lr_fu_tp <- compute_sampling_prob_based_on_all_records(controls_ds_f,
                                                       controls_ds_f$Sampled,
                                                       c("Vitfup2023", "first_procedure"))
## Get weights for controls in NCC only
lr_fu_tp_ncc_controls_df <- lr_fu_tp %>%
  left_join(full_ds %>%
              select(administratienummer, Sample_id) %>%
              unique(),
            by = "administratienummer") %>%
  right_join(ncc_ds %>%
               dplyr::filter(Metastasis == "Control"),
             by = "Sample_id") %>%
  relocate(administratienummer, Sample_id, Patient_ID_SKY, Metastasis)

# Combine cases and controls weights
lr_weights_ncc <- rbind(weights_cases_ncc_df, lr_fu_tp_ncc_controls_df) %>% arrange(Sample_id)
#-------------------------------------------------------------------------------

# Save results
#-------------------------------------------------------------------------------
write.csv(lr_fu_tp, file.path(output_dir, "WC_LR_allcontrols_valset.csv"), row.names = F)
write.csv(lr_weights_ncc, file.path(output_dir, "WC_LR_NCC_cohort_valset.csv"), row.names = F)
write.csv(lr_weights_ncc %>% select(-administratienummer), file.path(output_dir, "WC_LR_NCC_cohort_valset_noADMIN.csv"), row.names = F)
#-------------------------------------------------------------------------------