#---------------------------------------------------
# Aim: Prepare model development pipeline with
# late prediction integration and late feature integration
# Author: R.Ruiter
# Input: Clinical dataset
#        TPM gene expression matrix
#        Raw count gene expression matrix
#        Bailey DvP gene list
#        QC exclusion list
# Output: Bootstrap-validated prediction models for
#        Late prediction integration
#        Late feature integration (SCCore-GEP)
#---------------------------------------------------

#-------------------------------------------------------------------------------

# 1. Input parameters and global functions
#-------------------------------------------------------------------------------

strp <- function(x) substr(x, 1, 15)
# Helper function to subset gene expression matrix for a given signature
subset_expression_matrix <- function(signature_gene_ids, expression_matrix) {
  genes_in_matrix <- intersect(signature_gene_ids, colnames(expression_matrix))
  return(expression_matrix[, genes_in_matrix, drop = FALSE])
}

## Number of outer bootstrap repetitions if not NULL
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  outer <- NULL
} else {
  outer <- as.numeric(args[1])
  dir_scripts_sccore_gep <- args[2]
  dir_results_intermediate_sccore_gep <- args[3]
  dir_scripts_functions <- args[4]
  p1_tpms <- args[5]
  p2_counts <- args[6]
  p7_badQC <- args[7]
  project_dir <- args[8]
  p11_weights_unmatched <- args[9]
  p10_weights <- args[10]
  p12_genes_to_keep <- args[11]
  p6_bailey <- args[12]
  p13_mutation_data <- args[13]
  p14_all_pathways <- args[14]
  p15_input138 <- args[15]
  experiment_setup <- args[16]
  results_dir <- args[17]
  p9_gtf <- args[18]
}

outer_breps <- outer


# Experiment name
if (is.null(outer_breps)) {
experiment_outer <- experiment_sccore_gep
}else{
  experiment_outer <- paste0("CI/",outer)
}
# Directory where results will be stored (might need to be created)

if (!dir.exists(file.path(dir_results_intermediate_sccore_gep,  experiment_outer))){
  dir.create(file.path(dir_results_intermediate_sccore_gep,  experiment_outer), recursive = T)}


# Workflow settings
## Remove samples wityh bad QCs
remove_badQCsamples <- T
## Input genes of interest
genes_of_interest <- NULL
## Seed value
seed_val <- 42
## Compute bootstrap CI
if (is.null(outer_breps)) {
  compute_cis <- F
} else {
  compute_cis <- T
}
## Number bootstrap repetitions
breps <- 200
## Samples matched or unmatched
samples_matched <- "matched"

# Model buildling settings
## Features selection method(s) (multiple options accepted)
fselection <- "deseq2"

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
source_type <- "late"
## Additional info for model builing
### For features selection with DESeq2:
### - deseq2_design: string indicating the design
### - deseq2_padj_thresh: maximum value of padj
### - deseq2_logfc_thresh: minimum absolute value of logFC
### - deseq2_feats_thresh_by_logfc: maximum number of features allowed, they will be filtered by increasing logFC threshold
### - alphas2: alpha for the combined DvP and non-DvP gene model
### - model_final (included later): which second model needs to be used to stack the combined selected genes
###   - model2: coxnet
###   - RSF: RSF
### - combi_type == "Stacked", to combine predictions, or "stackedGEP" to combine selected genes
## Late feature integration (SCCore-GEP) and late prediction integration, respectively
## For confidence interval, only SCCore-GEP is run
if (is.null(outer_breps)){
  ## Model(s)
  models <- c("coxnet","RSF")
  models_params_list <-list(
    list(deseq2_design = paste0("~Biopsy_excision_2f+Sex+", outcome_var),
         deseq2_padj_thresh = 0.05,
         deseq2_logfc_thresh = 0.5,
         deseq2_feats_thresh_by_logfc = 50,
         alphas2 = c(0),
         combi_type = "stackedGEP"),
    list(deseq2_design = paste0("~Biopsy_excision_2f+Sex+", outcome_var),
         deseq2_padj_thresh = 0.05,
         deseq2_logfc_thresh = 0.5,
         deseq2_feats_thresh_by_logfc = 50,
         alphas2 = c(0),
         combi_type = "Stacked"))
}else{
  ## Model(s)
  models <- c("coxnet")
  models_params_list <-list(
    list(deseq2_design = paste0("~Biopsy_excision_2f+Sex+", outcome_var),
         deseq2_padj_thresh = 0.05,
         deseq2_logfc_thresh = 0.5,
         deseq2_feats_thresh_by_logfc = 50,
         alphas2 = c(0),
         combi_type = "stackedGEP"))
}


#-------------------------------------------------------------------------------

# 2. Source file with experiment set-up
#-------------------------------------------------------------------------------
source(file.path(dir_scripts_functions, "Model_input_processing.R"))

# Subset TPM and counts data to only have relevant genes
# Create expression matrices for all the required signatures
colnames(x_data) <- strp(colnames(x_data))
rownames(x_data_counts) <- strp(rownames(x_data_counts))
# Get DvP and non-DvP
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
if (is.null(outer)) {
print(paste0("gene matrix: ", nrow(x_data), " samples x ", ncol(x_data), " features."))
print(paste0("Bailey gene matrix: ", nrow(x_data_bailey), " samples x ", ncol(x_data_bailey), " features."))
print(paste0("Non-Bailey gene matrix: ", nrow(x_data_other), " samples x ", ncol(x_data_other), " features."))
}
gene_base <- gsub("\\..*", "", colnames(x_data))


#-------------------------------------------------------------------------------

# 5. Run workflow and save results
#-------------------------------------------------------------------------------
for (md in models){
    for (models_params in models_params_list){
    
      # Select final stacking model based on primary model
      models_params$model_final <- ifelse(
        md == "RSF",
        "RSF",
        "model2"
      )
      
      integration_name <- ifelse(
        models_params$combi_type == "stackedGEP",
        "Late_feature_integration",
        "Late_prediction_integration"
      )
      print(paste0("Training: ", md, " - ",integration_name ))
      bootstrap_results <- bootstrap_validation(feat_sel_method = fselection,
                                                topfeats = NA,
                                                model_method = md,
                                                ds_x_list = list(x_data_bailey, x_data_other) ,  
                                                ds_additional = clinical_data,
                                                ds_x_counts_list  = list(x_data_counts_bailey, x_data_counts_other),
                                                outcome_name = outcome_var,
                                                fup_name = follow_up_var,
                                                weights_name = weights_var,
                                                tp = time_point,
                                                B_inner = breps,
                                                seed = seed_val,
                                                cis = compute_cis,
                                                B_outer = outer_breps,
                                                model.params = models_params,
                                                bad_samples = sample_ids_bad_qcs,
                                                remove_bad_samples = remove_badQCsamples,
                                                genes_oi = genes_of_interest,
                                                dir_results = file.path(dir_results_intermediate_sccore_gep,  experiment_outer))
      save_model_bootstrap(
        bootstrap_results,
        paste0(md, "_", integration_name),
        file.path(dir_results_intermediate_sccore_gep,  experiment_outer)
      )
    }
  }
