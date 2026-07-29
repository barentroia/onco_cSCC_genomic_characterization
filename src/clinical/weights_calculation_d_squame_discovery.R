#-------------------------------------------------------------------------------
# Aim: Compute weights for development and evaluation of gene expression models
# Author: B. Rentoia Pacheco, based also on code of Lara Pozza, Joleen Traets
# Input:
# - Dataset with clinical data from the samples in WP2 cohort
# - Patients and set IDs of samples with bad RNA-seq QCs
# Output: Sampling probabilities for controls in WP2 cohort, in different settings and corresponding weights
#-------------------------------------------------------------------------------

# 0. Setup
#--------------------
library(openxlsx)
library(dplyr)
library(tidyr)
library(naniar)

# Directories: 
# Directory WP2 project
dir_project_WP2 = file.path("W:", "Data", "ONDERZOEK",
                            "Studies", "2021", "EMCD21054 [StepIdent_CP_G_refined]")
# Directory WP1 project
dir_project_WP1 = file.path("W:", "Data", "ONDERZOEK",
                            "Studies", "2021", "EMCD21053 [StepIdent_CP_refined]")
# Results directory
dir_results =  file.path(dir_project_WP2, "2-Data", "Intermediate datafiles",
                         "WP2 clinical data", "Weight computation",
                         "09_final_weights_september_2025")
if (!dir.exists(dir_results)){dir.create(dir_results, recursive = T)}

# 1. Load required datasets
#--------------------
# File with IDs of samples having bad RNA-seq QCs
samples_bad_rnaseq_qcs_fn <- file.path(dir_project_WP2, "2-Data", "Intermediate datafiles",
                                       "WP2 clinical data", "Sent to SkylineDx",
                                       "Skyline sent to EMC", "sample_ids_bad_QC.csv")

samples_bad_rnaseq_qcs <- read.csv(samples_bad_rnaseq_qcs_fn)

# Filename with mapping between `administratienummer` and patient skyline ID
mapping_fn <- file.path(dir_project_WP2, "2-Data", "Intermediate datafiles",
                        "WP2 clinical data","Intermediate datafiles",
                        "2023_11_13 RS_WP2 WP2 clinical data ALL records_v1_1.xlsx")

# Filename with IDs of samples in the 390 NCC cohort
ncc_390_ids_fn <- file.path(dir_project_WP2, "2-Data", "Intermediate datafiles",
                            "WP2 clinical data", "Sent to SkylineDx",
                            "Skyline sent to EMC", "sample_ids_complete_v2.csv")

# File with mapping between `administratienummer` and patient skyline ID
mapping <- read.xlsx(mapping_fn, sheet = "Samples")
# File with IDs of samples in 390 NCC cohort
ncc_390_ids <- read.csv(ncc_390_ids_fn)

# 2. Preprocessing that accounts for the fact that samples with bad QC were not sampled and computes CP risk for all pathology records in the nationwide database:
#----------------------------------------
# This code links databases to obtain the full cohort of 19,120 patients, it uses the entire database and contains IDs from the cancer registry that cannot be made publicly available
# The code produces the dataframes: wp2_dataset_sample_info and ds_potential_controls_0709
source(file.path(dir_project_WP2, "2-Data", "Intermediate datafiles","WP2 clinical data","Script","ncc_weights","WC_00_wp2_preprocessing_weights_computation.R"))

# Pairs with bad RNA-seq QCs are not included in the final NCC cohort,
# so their `Sampled` variable should be 0 (instead of 1):

## Extract Palga excerpt IDs of samples to exclude
palga_excerptids_samples_to_exclude <- wp2_dataset_sample_info$PALGAexcerptid[wp2_dataset_sample_info$Patient_ID_SKY%in%samples_bad_rnaseq_qcs$Patient_ID_SKY]
ds_potential_controls_0709_removed_bad_rnaseq_qcs <- ds_potential_controls_0709 %>%
  mutate(Sampled = ifelse(PALGAexcerptid %in% palga_excerptids_samples_to_exclude, 0, Sampled))

## Check that numbers match: 
### 195 controls in NCC, before removing pairs with bad RNA-seq QCs
### 183 controls in NCC, after removing pairs with bad RNA-seq QCs
stopifnot(sum(ds_potential_controls_0709$Sampled) == sum(ds_potential_controls_0709_removed_bad_rnaseq_qcs$Sampled) + nrow(samples_bad_rnaseq_qcs)/2)
stopifnot(nrow(ds_potential_controls_0709)==nrow(ds_potential_controls_0709_removed_bad_rnaseq_qcs))

# Many records are missing pathological variables, so we fill them in with what is the most likely value:
# Fill in records with missing data in a simple way:
sim_imp_records_0709_removed_bad_rnaseq_qcs = ds_potential_controls_0709_removed_bad_rnaseq_qcs %>%
  select(Sampled,FU_metastasis_years,Lab.name,Sex,Age,Tumor_location_cats,Differentiation,Tumor_diameter,PNI_or_LVI,Number_of_cSCC_before_culprit,Type_of_material,Tissue_involvement,PALGAexcerptid,administratienummer)%>%
  mutate(Differentiation = replace_na(Differentiation,"0"),
         Tumor_diameter = replace_na(Tumor_diameter,median(Tumor_diameter,na.rm=TRUE)),
         Tumor_location_cats = replace_na(Tumor_location_cats,"Face"),
         Type_of_material = replace_na(Type_of_material,"Excision/Reexcision/Resectie")) %>%as.data.frame()
cp_risk = compute_CP_BWH(sim_imp_records_0709_removed_bad_rnaseq_qcs)$Met_risk
sim_imp_records_0709_removed_bad_rnaseq_qcs$CP_risk = cp_risk*100 # Note: this step modifies 12 records/5328, in which the CP risk of the biopsy and the excision were different and the other CP risk was chosen. I decided to keep it this way, because in theory it could also have happened to records in which the CP risk score was missing both in the biopsy and the excision. Therefore, the CP risk score is obtained in the same way for all records.


# 3. Weight computation
#----------------------------------------
# 3.1: Functions to compute different sets of weights:
# Sampling probability based on all records:
compute_sampling_prob_based_on_all_records = function(df_all_records,sampled_records_bin){
  # Logistic regression to model sampling probability of each record:
  fit = glm(sampled_records_bin~.,family="binomial",data=df_all_records[,c("FU_metastasis_years","Sex","Age","Tumor_location_cats","Differentiation","Tumor_diameter","Tissue_involvement","PNI_or_LVI","Number_of_cSCC_before_culprit","Type_of_material")])
  
  # Get record sampling probabilities:
  rec_sampling_probs = predict(fit,type="response",newdata = df_all_records)
  
  df_predictions = data.frame("administratienummer"= df_all_records$administratienummer,"Samp_prob"= rec_sampling_probs)
  
  # Aggregate sampling probabilities per patient:
  df_predictions_per_patient = (1-df_predictions) %>% group_by(administratienummer)%>%
    summarise(Samp_prob_patient = 1-prod(Samp_prob))%>%
    ungroup() %>%
    mutate(administratienummer=1-administratienummer)%>%
    as.data.frame()
  
  df_predictions_per_patient$Weight = 1/df_predictions_per_patient$Samp_prob_patient
  
  weight_ncc_sum = sum(df_predictions_per_patient$Weight[df_predictions_per_patient$administratienummer%in%df_all_records$administratienummer[df_all_records$Sampled==1]])
  # Rescaling:
  df_predictions_per_patient$Weight_rescaled = df_predictions_per_patient$Weight*nrow(df_predictions_per_patient)/weight_ncc_sum
  df_predictions_per_patient$Metastasis = 0
  return(df_predictions_per_patient)
}

# 3.2 Compute weights:
#----------------------------------------
# Define full cohort and sampled cases:
# Full Cohort:
df_all_records_full_control_cohort = sim_imp_records_0709_removed_bad_rnaseq_qcs

# Cases info:
for(n_events in c(183,195)){
  df_cases = wp2_dataset_tumors_info[which(wp2_dataset_tumors_info$Metastasis=="Case"),]
  total_events = 305
  if(n_events == 183){
    df_cases = df_cases%>%filter(Patient_ID_SKY%in%ncc_390_ids$Patient_ID_SKY & !Patient_ID_SKY%in%samples_bad_rnaseq_qcs$Patient_ID_SKY)%>%as.data.frame()
    df_all_records_full_control_cohort$Sampled = sim_imp_records_0709_removed_bad_rnaseq_qcs$Sampled
    admin_ncc_controls = df_all_records_full_control_cohort$administratienummer[which(df_all_records_full_control_cohort$Sampled==1)]
  }else if (n_events==195){
    df_cases = df_cases%>%filter(Patient_ID_SKY%in%ncc_390_ids$Patient_ID_SKY )%>%as.data.frame()
    stopifnot(ds_potential_controls_0709_removed_bad_rnaseq_qcs$PALGAexcerptid==ds_potential_controls_0709$PALGAexcerptid)
    df_all_records_full_control_cohort$Sampled =ds_potential_controls_0709$Sampled
    table(df_all_records_full_control_cohort$Sampled)
    admin_ncc_controls = df_all_records_full_control_cohort$administratienummer[which(df_all_records_full_control_cohort$Sampled==1)]
  }
  
  print(nrow(df_cases))
  print(table(df_all_records_full_control_cohort$Sampled))
  
  # Build dataframe for cases:
  df_weights_cases = data.frame("administratienummer"=df_cases$administratienummer,"Samp_prob_patient"=n_events/total_events,"Weight"= total_events/n_events,"Weight_rescaled"= total_events/n_events,"Metastasis"=1)
  
  # Controls: 
  #Sampling probabilities of Controls:
  
  df_sampling_probs = compute_sampling_prob_based_on_all_records(df_all_records_full_control_cohort,df_all_records_full_control_cohort$Sampled)
  
  
  # Save them (all controls):
  write.csv(df_sampling_probs,file.path(dir_results,paste0("WC_01_",weight_method,"_all_controls_samp_prob_weights",n_events,".csv")), row.names = F)
  
  # Save the probabilities of cSCCs in the NCC cohort only:
  ncc_cohort_weights = rbind(df_sampling_probs[df_sampling_probs$administratienummer%in%admin_ncc_controls,],df_weights_cases)
  ncc_cohort_weights = merge(ncc_cohort_weights ,wp2_dataset_sample_info[,c("administratienummer","Patient_ID_SKY")],by="administratienummer",all.x=TRUE)
  ncc_cohort_weights = ncc_cohort_weights[!duplicated(ncc_cohort_weights$administratienummer),]
  ncc_cohort_weights = merge(ncc_cohort_weights,ncc_390_ids[,c("Patient_ID_SKY","SkylineDx.ID")],by.x = "Patient_ID_SKY")
  write.csv(ncc_cohort_weights,file.path(dir_results,paste0("WC_01_",weight_method,"_samp_prob_weights",n_events,"_NCC_cohort.csv")), row.names = F)
  write.csv(ncc_cohort_weights[,which(colnames(ncc_cohort_weights)!="administratienummer")],file.path(dir_results,paste0("WC_01_",weight_method,"_samp_prob_weights",n_events,"_NCC_cohort_noADMIN.csv")), row.names = F)
    
  
}
