#---------------------------------------------------
# Aim: Prepare WES model development pipeline unmatched with 126 samples
# Author: R.Ruiter
# Input: Clinical dataset
#        WES mutation dataset
#        Bailey DvP gene list
#        QC exclusion list
#        samples that overlap between WES and GEP data
# Output: WES model on 126 samples




#-------------------------------------------------------------------------------

# 0. Load R library
#-------------------------------------------------------------------------------
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("setdiff", "base")
conflicted::conflict_prefer("intersect", "base")


strp <- function(x) substr(x, 1, 15)
#-------------------------------------------------------------------------------
# Helper function to subset gene expression matrix for a given signature
subset_expression_matrix <- function(signature_gene_ids, expression_matrix) {
  genes_in_matrix <- intersect(signature_gene_ids, colnames(expression_matrix))
  return(expression_matrix[, genes_in_matrix, drop = FALSE])
}

#-------------------------------------------------------------------------------

# 1. Input parameters and global functions
#-------------------------------------------------------------------------------
# Removing X and Y chr cnv
remove_sex_chr_cnv <- function(df) {
  df %>% select(-matches("^chr[XY]"))
}

# Mutation feature prefiltering
# pathways: keep if altered in >= 3 samples
# CNAs: keep if altered in >= 10 samples
# CNA features: keep if prevalence >= 20%
prefilter_mutation_data <- function(df,
                                    clinical_data,
                                    samples_bad_qcs,
                                    remove_badQCsamples = TRUE,
                                    pathway_cols,
                                    min_pathway_samples = 3,
                                    min_cna_samples = 10,
                                    min_final_prop = 0.20,
                                    keep_covariates = c("Tumor.cellularity.avg.pct")) {
  
  # Store original full dataset (with bad QC samples)
  df_original <- df
  # Remove bad QC samples temporarily for filtering
  if (remove_badQCsamples) {
    print("remove bad qc samples")
    sample_ids_bad_qcs <- clinical_data %>%
      dplyr::filter(SkylineDx.ID %in% samples_bad_qcs$Skyline_ID) %>%
      pull(SkylineDx.ID)
    df_filtered <- df[
      !(rownames(df) %in% sample_ids_bad_qcs),
      ,
      drop = FALSE
    ]
  } else {
    sample_ids_bad_qcs <- c()
    df_filtered <- df
  }
  # Covariates
  existing_covariates <- intersect(
    keep_covariates,
    colnames(df_filtered)
  )
  # Identify feature columns
  feature_cols <- setdiff(
    colnames(df_filtered),
    existing_covariates
  )
  # Split pathways vs CNA columns
  pathway_features <- intersect(feature_cols, pathway_cols)
  cna_features <- setdiff(feature_cols, pathway_features)
  # 1. Pathway filtering (NO prevalence filter)
  pathway_counts <- colSums(
    df_filtered[, pathway_features, drop = FALSE] > 0,
    na.rm = TRUE
  )
  keep_pathways <- names(pathway_counts)[
    pathway_counts >= min_pathway_samples
  ]
  remove_pathways <- names(pathway_counts)[
    pathway_counts < min_pathway_samples
  ]
  # 2. CNA filtering (sample count first)
  cna_counts <- sapply(cna_features, function(feat) {
    x <- df_filtered[[feat]]
    ux <- sort(unique(na.omit(x)))
    # Binary 0/1
    if (all(ux %in% c(0, 1))) {
      sum(x == 1, na.rm = TRUE)
      # CN state 1/2/3
    } else if (all(ux %in% c(1, 2, 3))) {
      sum(x != 2, na.rm = TRUE)
    } else {
      stop(paste(
        "Unknown CNA encoding for feature:",
        feat,
        "| unique values:",
        paste(ux, collapse = ", ")
      ))
    }
  })
  keep_cnas <- names(cna_counts)[
    cna_counts >= min_cna_samples
  ]
  remove_cnas <- names(cna_counts)[
    cna_counts < min_cna_samples
  ]
  # 3. Final prevalence filter(only CNAs)
  cna_prevalence <- sapply(keep_cnas, function(feat) {
    x <- df_filtered[[feat]]
    ux <- sort(unique(na.omit(x)))
    if (all(ux %in% c(0, 1))) {
      mean(x == 1, na.rm = TRUE)
    } else if (all(ux %in% c(1, 2, 3))) {
      mean(x != 2, na.rm = TRUE)
    } else {
      mean(x > 0, na.rm = TRUE)
    }
  })
  cna_prevalence <- unlist(cna_prevalence)
  keep_cnas_final <- names(cna_prevalence)[
    cna_prevalence >= min_final_prop
  ]
  remove_cnas_final <- names(cna_prevalence)[
    cna_prevalence < min_final_prop
  ]
  # Final feature set
  keep_features_final <- c(
    keep_pathways,
    keep_cnas_final
  )
  # Final dataset (keep covariates + filtered features)
  df_final <- df_original %>%
    dplyr::select(all_of(c(
      existing_covariates,
      keep_features_final
    )))
  # Reporting
  message("Bad QC samples removed temporarily: ",
          length(sample_ids_bad_qcs))
  if (length(remove_pathways) > 0) {
    message("Pathways removed (< ", min_pathway_samples, " altered samples):")
    message(paste(remove_pathways, collapse = ", "))
  }
  if (length(remove_cnas) > 0) {
    message("CNA removed (< ", min_cna_samples, " altered samples):")
    message(paste(remove_cnas, collapse = ", "))
  }
  if (length(remove_cnas_final) > 0) {
    message("CNA removed by prevalence (< ",
            min_final_prop * 100, "%):")
    message(paste(remove_cnas_final, collapse = ", "))
  }
  message("Final feature count: ", length(keep_features_final))
  message("Final matrix dimensions: ",
          nrow(df_final), " x ", ncol(df_final))
  message(paste(colnames(df_final), collapse = ", "))
  message("======================================")
  return(df_final)
}


# Directory where results will be stored (might need to be created)
if (!dir.exists(dir_results_intermediate_integration)){dir.create(dir_results_intermediate_integration, recursive = T)}

#-------------------------------------------------------------------------------

# 1.1 Input parameters
#-------------------------------------------------------------------------------
# Workflow settings
## Remove samples wityh bad QCs
remove_badQCsamples <- T
## Input genes of interest
genes_of_interest <- NULL
## Seed value
seed_val <- 42
## Number of outer bootstrap repetitions
outer_breps <- NULL
## Compute boostrap CI
if (is.null(outer_breps)) {
  compute_cis <- F
} else {
  compute_cis <- T
}
compute_cis<-F
## Number boostrap repetitions
breps <- 200

# Model buildling settings
## gene selection method(s) (multiple options accepted)
# mutation feature selection uses coxnet
fselections <-  c("deseq2")
## Number of input features (multiple options accepted)
nfeats <- NA
# Model(s) only for SCCore-GEP and combination of GEP+mutation
models <- c("coxnet")
## Name of outcome variable
outcome_var <- "Metastasis_num"
## Name of follow-up variable
follow_up_var <- "FU_metastasis_years"
## Name of weights' variable
weights_var <- "Weight_rescaled"
## Time point
time_point <- 5
## Weights'type
weights_type <- "record_based"
## If early or late integration script need to be sourced
source_type <- "late_wes"
## Samples matched or unmatched
samples_matched <- "unmatched"

#-------------------------------------------------------------------------------

# 1.2 Input parameters for model
#-------------------------------------------------------------------------------

## Additional info for model builing for SCCore-GEP (named list or NULL)
### For features selection with DESeq2:
### - deseq2_design: string indicating the design
### - deseq2_padj_thresh: maximum value of padj
### - deseq2_logfc_thresh: minimum absolute value of logFC
### - deseq2_feats_thresh_by_logfc: maximum number of features allowed, they will be filtered by increasing logFC threshold
### - alphas2: list of alpha options for the combined model
###   - c(0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)
###   - c(0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5)
###   - c(0)
### - model_final: which second model needs to be used to stack the combined selected genes
###   for WES+GEP, only model2 has been tested as SCCore-GEP uses coxnet
###   - model2: coxnet
###   - model3: glmnet
###   - RSF: RSF
###   - customRF: customRF
### - combi_type == "Stacked", to combine predictions, or "stackedGEP" to combine selected genes
###   for WES+GEP, only stackedGEP has been tested as SCCore-GEP uses stackedGEP

## Additional info for model builing for WES (named list or NULL)
### - combi_mut_type: how mutation and GEP will be combined
###   - late_feature_mut to combine mutation and GEP features
### - cellularity_check: Cellularity adjustment or not. Empty or "Yes"
### - mut_alphas_numbers: alpha selection for WES models 
###   - 0.5 as paper
### - unmatched: dont use matched set_id samples
## Final use:

models_params_list <-list(

  list(deseq2_design = paste0("~Biopsy_excision_2f+Sex+", outcome_var),
       deseq2_padj_thresh = 0.05,
       deseq2_logfc_thresh = 0.5,
       deseq2_feats_thresh_by_logfc = 50,
       model_final = "model2",
       alphas2 =  0,
       combi_type = "stackedGEP",
       combi_mut_type = "late_feature_mut",
       mnumbers = c(0.5),
       unmatched = "Yes"
))

#-------------------------------------------------------------------------------

# 2. Source file with experiment set-up
#-------------------------------------------------------------------------------
# expression and clinical information is needed to overlap sampels with mutation data
source(file.path(dir_scripts_functions, "Model_input_processing.R"))

#-------------------------------------------------------------------------------

# 3. Make input mutation data
#-------------------------------------------------------------------------------
library(dplyr)
library(stringr)
# remove Chr X and Y
mutation_data_amp_del <- remove_sex_chr_cnv(mutation_data_amp_del)

mutation_data_list <- list(
  amp_del   = mutation_data_amp_del

)

mutation_data_list_no_cellularity <- lapply(
  mutation_data_list,
  function(df) df %>% select(-Tumor.cellularity.avg.pct)
)
#-------------------------------------------------------------------------------

# 6. Get DvP and non-DvP
#-------------------------------------------------------------------------------
x_data_bailey <- subset_expression_matrix(bailey_gene_ids$ensembl_gene_id, x_data)
x_data_counts_bailey <- as.data.frame(t(subset_expression_matrix(bailey_gene_ids$ensembl_gene_id, t(x_data_counts))))

other_genes <- setdiff(colnames(x_data), bailey_gene_ids$ensembl_gene_id)
x_data_other <- subset_expression_matrix(other_genes, x_data)
x_data_counts_other <- as.data.frame(t(subset_expression_matrix(other_genes, t(x_data_counts))))

# Align datasets
stopifnot(all(colnames(x_data_counts) == rownames(x_data)))
stopifnot(all(rownames(x_data_bailey) == clinical_data$SkylineDx.ID))
stopifnot(all(colnames(x_data_counts_bailey) == clinical_data$SkylineDx.ID))
stopifnot(all(colnames(x_data_counts_bailey) == rownames(x_data_bailey)))
print(paste0("gene matrix: ", nrow(x_data), " samples x ", ncol(x_data), " features."))
print(paste0("Bailey gene matrix: ", nrow(x_data_bailey), " samples x ", ncol(x_data_bailey), " features."))
print(paste0("Non-Bailey gene matrix: ", nrow(x_data_other), " samples x ", ncol(x_data_other), " features."))
gene_base <- gsub("\\..*", "", colnames(x_data))
#-------------------------------------------------------------------------------

# 5. Run workflow and save results
#-------------------------------------------------------------------------------
prop_prevalence <- 0.20
for (mut_type in names(mutation_data_list_no_cellularity)) {
  
  message("Running mutation type: ", mut_type)
  
  mutation_data <- mutation_data_list_no_cellularity[[mut_type]]
  
  # Literature-based mutation feature prefiltering
  mutation_data <- prefilter_mutation_data(
    df = mutation_data,
    clinical_data = clinical_data,
    samples_bad_qcs = samples_bad_qcs,
    remove_badQCsamples = remove_badQCsamples,
    pathway_cols = pathway_cols,
    min_pathway_samples = 3,
    min_cna_samples = 10,
    min_final_prop = prop_prevalence,
    keep_covariates = c("Tumor.cellularity.avg.pct")
  )
  
  message("Final mutation matrix dimensions: ",
          nrow(mutation_data), " x ", ncol(mutation_data))

for (md in models){
  for (fselection in fselections){
    for (models_params in models_params_list){
      
      if (md == "coxnet") {
        models_params$model_final <- "model2"
        
      } 

    for (nfeat in nfeats){
      bootstrap_results <- bootstrap_validation(feat_sel_method = fselection,
                                                topfeats = nfeat,
                                                model_method = md,
                                                ds_x_list = list(x_data_bailey, x_data_other) ,  # x_data_other for remaining, x_data for datadriven
                                                ds_additional = clinical_data,
                                                ds_x_counts_list  = list(x_data_counts_bailey, x_data_counts_other), #x_data_counts_other for remaining, x_data_counts for datadriven
                                                ds_mutation = mutation_data,
                                                outcome_name = outcome_var,
                                                fup_name = follow_up_var,
                                                weights_name = weights_var,
                                                tp = time_point,
                                                B_inner = breps,
                                                seed = 123,
                                                cis = compute_cis,
                                                B_outer = outer_breps,
                                                model.params = models_params,
                                                bad_samples = sample_ids_bad_qcs,
                                                remove_bad_samples = remove_badQCsamples,
                                                genes_oi = genes_of_interest,
                                                dir_results = dir_results_intermediate_integration)
      save_model_bootstrap(bootstrap_results, "Combined_model", dir_results_intermediate_integration)
    }
  }
}
}}