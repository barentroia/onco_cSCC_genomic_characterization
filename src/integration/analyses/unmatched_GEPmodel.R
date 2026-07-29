#---------------------------------------------------
# Aim: Prepare RNAseq model development pipeline unmatched with 126 samples
# Author: R.Ruiter
# Input: Clinical dataset
#        TPM gene expression matrix
#        Raw count gene expression matrix
#        Bailey DvP gene list
#        QC exclusion list
#        samples that overlap between WES and GEP data
# Output: RNAseq model on 126 samples



#-------------------------------------------------------------------------------

# 0. Load R library
#-------------------------------------------------------------------------------
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("setdiff", "base")
conflicted::conflict_prefer("intersect", "base")

#-------------------------------------------------------------------------------

# 1. Input parameters and global functions
#-------------------------------------------------------------------------------
strp <- function(x) substr(x, 1, 15)
# Helper function to subset gene expression matrix for a given signature
subset_expression_matrix <- function(signature_gene_ids, expression_matrix) {
  genes_in_matrix <- intersect(signature_gene_ids, colnames(expression_matrix))
  return(expression_matrix[, genes_in_matrix, drop = FALSE])
}

# Directory where results will be stored (might need to be created)
if (!dir.exists(dir_results_intermediate_integration)){dir.create(dir_results_intermediate_integration, recursive = T)}

#-------------------------------------------------------------------------------

# 1. Input parameters
#-------------------------------------------------------------------------------
# Workflow settings
## Remove samples with bad QCs
remove_badQCsamples <- T
## Input genes of interest
genes_of_interest <- NULL
## Seed value
seed_val <- 42
## Number of outer bootstrap repetitions
outer_breps <- NULL
## Compute boostrap CI
compute_cis<-F
## Number boostrap repetitions
breps <-200
## Samples matched or unmatched
samples_matched <- "unmatched"
# Model buildling settings
## Features selection method(s) (multiple options accepted)
fselections <- "deseq2"
## Number of input features (multiple options accepted)
nfeats <- NA
## Model(s)
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
### - combi_type ==  "stackedGEP" to combine selected genes
### - unmatched: to account for unmatches samples
## FEATURES
models_params_list <-list(
  list(deseq2_design = paste0("~Biopsy_excision_2f+Sex+", outcome_var),
       deseq2_padj_thresh = 0.05,
       deseq2_logfc_thresh = 0.5,
       deseq2_feats_thresh_by_logfc = 50,
       model_final = "model2",
       alphas2 = c(0),
       combi_type = "stackedGEP",
       unmatched = "Yes"))


#-------------------------------------------------------------------------------

# 2. Source file with experiment set-up
#-------------------------------------------------------------------------------
source(file.path(dir_scripts_functions, "Model_input_processing.R"))

#-------------------------------------------------------------------------------

# 3. Subset TPM and counts data to only have relevant genes
#-------------------------------------------------------------------------------
# Create expression matrices for all the required signatures
colnames(x_data) <- strp(colnames(x_data))
rownames(x_data_counts) <- strp(rownames(x_data_counts))

# 4. Get DvP and non-DvP
#-------------------------------------------------------------------------------

# Align datasets
x_data_bailey <- subset_expression_matrix(bailey_gene_ids$ensembl_gene_id, x_data)
x_data_counts_bailey <- as.data.frame(t(subset_expression_matrix(bailey_gene_ids$ensembl_gene_id, t(x_data_counts))))

other_genes <- setdiff(colnames(x_data), bailey_gene_ids$ensembl_gene_id)
x_data_other <- subset_expression_matrix(other_genes, x_data)
x_data_counts_other <- as.data.frame(t(subset_expression_matrix(other_genes, t(x_data_counts))))
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
max_features <- ncol(x_data_bailey)
# Run workflow in for-loop and save results
for (md in models){
  for (fselection in fselections){
    for (models_params in models_params_list){
      if (md == "coxnet") {
        models_params$model_final <- "model2"
      } else if (md == "RSF") {
        models_params$model_final <- "RSF"
      }
    for (nfeat in nfeats){
      bootstrap_results <- bootstrap_validation(feat_sel_method = fselection,
                                                topfeats = nfeat,
                                                model_method = md,
                                                ds_x_list = list(x_data_bailey, x_data_other) ,  # x_data_other for remaining, x_data for datadriven
                                                ds_additional = clinical_data,
                                                ds_x_counts_list  = list(x_data_counts_bailey, x_data_counts_other), #x_data_counts_other for remaining, x_data_counts for datadriven
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
                                                dir_results = dir_results_intermediate_integration)
      save_model_bootstrap(bootstrap_results, "GEP_model", dir_results_intermediate_integration)
    }
  }
}
}