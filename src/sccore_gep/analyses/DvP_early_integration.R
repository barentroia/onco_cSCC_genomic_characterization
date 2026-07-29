#---------------------------------------------------
# Aim: Prepare model development pipeline with DvP
# non-DvP or early integration
# Author: R.Ruiter
# Input: Clinical dataset
#        TPM gene expression matrix
#        Raw count gene expression matrix
#        Bailey DvP gene list
#        QC exclusion list
# Output: Bootstrap-validated prediction models for
#        DvP
#        Non-DvP
#        Early integration
#---------------------------------------------------


#-------------------------------------------------------------------------------

# 1. Input parameters and global functions
#-------------------------------------------------------------------------------
# To remove versions from genes
strp <- function(x) substr(x, 1, 15)
# To subset expression matrix with genes of interest
subset_expression_matrix <- function(signature_gene_ids, expression_matrix) {
  genes_in_matrix <- intersect(signature_gene_ids, colnames(expression_matrix))
  return(expression_matrix[, genes_in_matrix, drop = FALSE])
}


# Directory where results will be stored (might need to be created)
if (!dir.exists(file.path(dir_results_intermediate_sccore_gep,  experiment_sccore_gep))){
  dir.create(file.path(dir_results_intermediate_sccore_gep,  experiment_sccore_gep), recursive = T)}


# Workflow settings
## Remove samples wityh bad QCs
remove_badQCsamples <- T
## Input genes of interest
genes_of_interest <- NULL
## Seed value
seed_val <- 42
## Number of outer bootstrap repetitions
outer_breps <- NULL
## Compute bootstrap CI, not done for these integration strategies
compute_cis <- F

## Number bootstrap repetitions
breps <- 200
## Samples matched or unmatched
samples_matched <- "matched"

# Model building settings
## Features selection method
fselection <- c("deseq2")
## Number of input features
nfeats <- NA
## Model(s)
models <- c("coxnet","RSF")
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
source_type <- "early"
## Additional info for model builing
### For features selection with DESeq2:
### - deseq2_design: string indicating the design
### - deseq2_padj_thresh: maximum value of padj
### - deseq2_logfc_thresh: minimum absolute value of logFC
### - deseq2_feats_thresh_by_logfc_by_logfc: maximum number of features allowed, they will be filtered by increasing logFC threshold
models_params_list <-list(
                          list(deseq2_design = paste0("~Biopsy_excision_2f+Sex+", outcome_var),
                               deseq2_padj_thresh = 0.05,
                               deseq2_logfc_thresh = 0.5,
                               deseq2_feats_thresh_by_logfc = 50))
                     
#-------------------------------------------------------------------------------

# 2. Source file with experiment set-up
#-------------------------------------------------------------------------------
# Sourcing script for processing the clinical data and prefiltering the gene expression data
source(file.path(dir_scripts_functions, "Model_input_processing.R"))

# Remove versions from genes and only keep DvP genes in gene expression
colnames(x_data) <- strp(colnames(x_data))
x_data_bailey <- subset_expression_matrix(bailey_gene_ids$ensembl_gene_id, x_data)
rownames(x_data_counts) <- strp(rownames(x_data_counts))
x_data_counts_bailey <- as.data.frame(t(subset_expression_matrix(bailey_gene_ids$ensembl_gene_id, t(x_data_counts))))
# Remove versions from genes and only keep non-DvP genes in gene expression
other_genes <- setdiff(colnames(x_data), bailey_gene_ids$ensembl_gene_id)
x_data_other <- subset_expression_matrix(other_genes, x_data)
x_data_counts_other <- as.data.frame(t(subset_expression_matrix(other_genes, t(x_data_counts))))
print(paste0("Early integration gene matrix: ", nrow(x_data), " samples x ", ncol(x_data), " features."))
print(paste0("DvP gene matrix: ", nrow(x_data_bailey), " samples x ", ncol(x_data_bailey), " features."))
print(paste0("Non-DvP gene matrix: ", nrow(x_data_other), " samples x ", ncol(x_data_other), " features."))

#-------------------------------------------------------------------------------
# 3. Select datasets to run
#-------------------------------------------------------------------------------
datasets <- list(
  DvP = list(
    data = x_data_bailey,
    counts = x_data_counts_bailey
  ),
  non_DvP = list(
    data = x_data_other,
    counts = x_data_counts_other
  ),
  Early_integration = list(
    data = x_data,
    counts = x_data_counts
  )
)

#-------------------------------------------------------------------------------
# 4. Run workflow and save results
#-------------------------------------------------------------------------------
for (dataset_name in names(datasets)) {
  
  # Select current dataset
  data <- datasets[[dataset_name]]$data
  counts <- datasets[[dataset_name]]$counts
  
  
  
  for (md in models){
    message("Running integration strategy: ", dataset_name, " on model " ,md)
      for (models_params in models_params_list){ 
          bootstrap_results <- bootstrap_validation(
            feat_sel_method = fselection,
            topfeats = NA,
            model_method = md,
            ds_x = data,
            ds_additional = clinical_data,
            ds_x_counts = counts,
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
            genes_oi = genes_of_interest
          )
          save_model_bootstrap(
            bootstrap_results,
            paste0(md,"_",dataset_name),file.path(dir_results_intermediate_sccore_gep,  experiment_sccore_gep))
          }}}