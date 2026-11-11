#===============================================================================
# Aim:      Develop and bootstrap-validate the WES-only model on the 126
#           unmatched samples: prefiltered WES mutation (pathway) and copy
#           number features in a Coxnet model (alpha = 0.5).
# Author:   R. Ruiter
# Input:    - Clinical, RNA-seq and WES data, loaded via Model_input_processing.R
#           - Bad QC sample list
# Output:   - WES_model_robj.rds in dir_results_intermediate_integration
# Requires: dir_scripts_functions, dir_results_intermediate_integration,
#           experiment_setup = "integration" (main.R) and config.R paths
# Usage:    Sourced from main.R
#===============================================================================

# 1. Settings
#-------------------------------------------------------------------------------
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("setdiff", "base")
conflicted::conflict_prefer("intersect", "base")


#-------------------------------------------------------------------------------

# 2. Helper functions and output directory
#-------------------------------------------------------------------------------
#' Strip the version suffix from Ensembl gene IDs (keep first 15 characters)
#'
#' @param x character vector. Ensembl gene IDs, e.g. "ENSG00000141510.18".
#'
#' @return character vector of unversioned IDs.
strp <- function(x) substr(x, 1, 15)


#' Subset an expression matrix to the genes of a signature
#'
#' @param signature_gene_ids character vector. Gene IDs of the signature.
#' @param expression_matrix matrix/data.frame. Expression data (samples x genes).
#'
#' @return The expression matrix restricted to the signature genes present.
subset_expression_matrix <- function(signature_gene_ids, expression_matrix) {
  genes_in_matrix <- intersect(signature_gene_ids, colnames(expression_matrix))
  return(expression_matrix[, genes_in_matrix, drop = FALSE])
}

#' Remove copy number features of the sex chromosomes (chrX, chrY)
#'
#' @param df data.frame. Mutation/copy number features.
#'
#' @return df without columns starting with "chrX" or "chrY".
remove_sex_chr_cnv <- function(df) {
  df %>% select(-matches("^chr[XY]"))
}

#' Prefilter WES mutation and copy number features
#'
#' Filtering is based on good QC samples only: pathways are kept if altered
#' in at least min_pathway_samples samples; copy number alterations (CNAs)
#' are kept if altered in at least min_cna_samples samples and in at least
#' min_final_prop of the samples. Covariates in keep_covariates are kept.
#'
#' @param df data.frame. Mutation/CNA features (samples x features).
#' @param clinical_data data.frame. Clinical data with SkylineDx.ID.
#' @param samples_bad_qcs data.frame. Bad QC samples (Skyline_ID column).
#' @param remove_badQCsamples logical. Exclude bad QC samples when filtering.
#' @param pathway_cols character vector. Columns that are pathway features.
#' @param min_pathway_samples integer. Minimum number of altered samples per pathway.
#' @param min_cna_samples integer. Minimum number of altered samples per CNA.
#' @param min_final_prop numeric. Minimum proportion of altered samples per CNA.
#' @param keep_covariates character vector. Columns always kept.
#'
#' @return data.frame with all samples, the covariates and the kept features.
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

# 3. Input parameters
#-------------------------------------------------------------------------------
# Workflow settings
## Remove samples with bad QCs
remove_badQCsamples <- T
## Input genes of interest
genes_of_interest <- NULL
## Seed value (the same as WES)
seed_val <- 123
## Number of outer bootstrap repetitions
outer_breps <- NULL
## Compute boostrap CI
compute_cis <- F
## Number boostrap repetitions
breps <- 200
# Model buildling settings
## Features selection method(s) (multiple options accepted)
fselections <- "none"
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
## Samples matched or unmatched
samples_matched <- "unmatched"
## If early or late integration script need to be sourced
source_type <- "early"

#-------------------------------------------------------------------------------

# 4. Source file with experiment set-up
#-------------------------------------------------------------------------------
# Expression and clinical information is needed to overlap samples with mutation data
source(file.path(dir_scripts_functions, "Model_input_processing.R"))

#-------------------------------------------------------------------------------

# 5. Make input mutation data
#-------------------------------------------------------------------------------
# Remove chrX and chrY copy number features
mutation_data_amp_del <- remove_sex_chr_cnv(mutation_data_amp_del)

mutation_data_list <- list(
  amp_del   = mutation_data_amp_del
)

mutation_data_list_no_cellularity <- lapply(
  mutation_data_list,
  function(df) df %>% select(-Tumor.cellularity.avg.pct)
)

# Less stricter regularization for mutation model
# Literature-based alpha was 0.5
models_params <- list(mut_alphas = c(0.5),
                      unmatched = "Yes")

#-------------------------------------------------------------------------------

# 6. Run workflow and save results
#-------------------------------------------------------------------------------
for (mut_type in names(mutation_data_list_no_cellularity)) {
  message("Running mutation type: ", mut_type)
  mutation_data <- mutation_data_list_no_cellularity[[mut_type]]
  prop_prevalence <- 0.20
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

# Run workflow in for-loop and save results
for (md in models){
  for (fselection in fselections){
    for (nfeat in nfeats){
      bootstrap_results <- bootstrap_validation(feat_sel_method = fselection,
                                                topfeats = nfeat,
                                                model_method = md,
                                                ds_x = mutation_data,
                                                ds_additional = clinical_data,
                                                ds_x_counts = NULL,
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
                                                genes_oi = genes_of_interest)
      save_model_bootstrap(
        bootstrap_results,
        "WES_model",
        dir_results_intermediate_integration
      )
    
  }
}}}
#-------------------------------------------------------------------------------
