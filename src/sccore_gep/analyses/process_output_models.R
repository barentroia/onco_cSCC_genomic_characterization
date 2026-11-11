#===============================================================================
# Aim:      Extract and standardize the performance metrics and the number of
#           selected genes of all SCCore-GEP development models (DvP, non-DvP,
#           early integration, late feature and late prediction integration)
#           for the figure scripts.
# Author:   R. Ruiter
# Input:    - Model .rds files in dir_results_intermediate_sccore_gep/<experiment_sccore_gep>
# Output:   - cox23_linear_predictors.csv: linear predictor of SCCore-GEP per sample
#           - perf_OOB_apparent_bootstrap.csv, perf_632plus.csv,
#             perf_apparent.csv and model_gene_counts.csv, used by
#             supp_figure_10A_performance_boxplot.R
#           All written to dir_results_intermediate_sccore_gep.
# Requires: dir_results_intermediate_sccore_gep, experiment_sccore_gep (main.R)
# Usage:    Sourced from main.R
#===============================================================================

#-------------------------------------------------------------------------------

# 1. Find all RDS model files
#-------------------------------------------------------------------------------
folders <- tibble(
  path = list.files(
    file.path(dir_results_intermediate_sccore_gep, experiment_sccore_gep) ,
    pattern = "\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
) %>%
  mutate(
    model_file = tools::file_path_sans_ext(basename(path)),
    model = model_file
  )

#-------------------------------------------------------------------------------

# 2. Lists for storing results
#-------------------------------------------------------------------------------
perf_all_list <- list()
perf_632_list <- list()
perf_app_list <- list()
gene_count_list <- list()

n_models <- nrow(folders)

#-------------------------------------------------------------------------------

# 3. Read models and extract performance
#-------------------------------------------------------------------------------
for(i in seq_len(n_models)){
  
  model <- readRDS(folders$path[i])
  model_name <- folders$model[i]
  
  message(sprintf(
    "[%d/%d] Reading model: %s",
    i,
    n_models,
    model_name
  ))
  

  # 3.1 Linear predictor of SCCore-GEP
  if (grepl("coxnet_Late_feature_integration", model_name)) {
    fit <- model$Final_model2
    
    X    <- fit$x_fitted_model          # 366 x 25
    beta <- fit$coefficients            # length 25, named
    
    # Make sure columns line up by gene ID
    stopifnot(all(colnames(X) %in% names(beta)))
    lp <- drop(X[, names(beta)] %*% beta) #sum of each gene log tpm multiplied by the coefficient
    names(lp) <- rownames(X)
    lp_df <- data.frame(sample = names(lp), lp = lp)
    write.csv(lp_df, file.path(
      dir_results_intermediate_sccore_gep ,
      "cox23_linear_predictors.csv"
    ) , row.names = FALSE)
  }

  # 3.2 OOB + bootstrap Apparent bootstrap
  oob_table_name <- if (!is.null(model$`All out of bag`)) {
    "All out of bag"
  } else {
    "All out of bootstrap"
  }
  
  oob <- model[[oob_table_name]] %>%
    mutate(type = "OOB")
  
  apparent_boot <- model$`All apparent bootstrap` %>%
    mutate(type = "Apparent")
  
  perf_all_list[[i]] <-
    bind_rows(oob, apparent_boot) %>%
    mutate(model = model_name) %>%
    filter(metric == "wcind_model") %>%
    dplyr::select(
      variable,
      group,
      metric,
      value,
      boots_sample,
      type,
      model
    )
  
  # 3.3 0.632+
  if(!is.null(model$`Optimism-Corrected0632plus`)){
    perf_632_list[[i]] <-
      model$`Optimism-Corrected0632plus` %>%
      filter(metric == "wcind_model") %>%
      dplyr::select(
        variable,
        group,
        metric,
        value
      ) %>%
      mutate(model = model_name)
  }
  
  # 3.4 Apparent performance
  if(!is.null(model$Apparent_performance)){
    perf_app_list[[i]] <-
      model$Apparent_performance %>%
      filter(metric == "wcind_model") %>%
      dplyr::select(
        variable,
        group,
        metric,
        value
      ) %>%
      mutate(model = model_name)
  }
  
  # 3.5 Extract final fitted model
  fm <- model$Final_model2
  if(is.null(fm)){
    fm <- model$Model_fit_opt
  }
  if(is.null(fm)){
    fm <- model$Model_fit
  }
  
  # 3.6 Count genes/features
  gene_count <- length(fm$coefficients)
  gene_count_list[[i]] <- tibble(
    model = model_name,
    gene_count = gene_count
  )
  
  rm(model)
  gc()
}

#-------------------------------------------------------------------------------

# 4. Combine results
#-------------------------------------------------------------------------------

perf_all_filtered <- bind_rows(perf_all_list)
perf_632plus <- bind_rows(perf_632_list)
perf_apparent <- bind_rows(perf_app_list)
gene_counts <- bind_rows(gene_count_list)

#-------------------------------------------------------------------------------

# 5. For Late prediction integration models
#-------------------------------------------------------------------------------

late_prediction_idx <- grepl(
  "Late_prediction_integration",
  gene_counts$model,
  ignore.case = TRUE
)

if(any(late_prediction_idx)){
  
  # Create corresponding Late feature integration
  # model names
  feature_model_names <- gene_counts$model[late_prediction_idx] %>%
    str_replace(
      "Late_prediction_integration",
      "Late_feature_integration"
    )
  
  # Find matching gene counts
  matching_gene_counts <- gene_counts$gene_count[
    match(
      feature_model_names,
      gene_counts$model
    )
  ]
  
  # Replace gene counts for Late prediction integration
  gene_counts$gene_count[late_prediction_idx] <-
    matching_gene_counts
}

#-------------------------------------------------------------------------------

# 6. Output directory
#-------------------------------------------------------------------------------
# Create output directory if needed
dir.create(
  dir_results_intermediate_sccore_gep ,
  recursive = TRUE,
  showWarnings = FALSE
)

#-------------------------------------------------------------------------------

# 7. Save processed performance datasets
#-------------------------------------------------------------------------------
write_csv(
  perf_all_filtered,
  file.path(
    dir_results_intermediate_sccore_gep ,
    "perf_OOB_apparent_bootstrap.csv"
  )
)

write_csv(
  perf_632plus,
  file.path(
    dir_results_intermediate_sccore_gep ,
    "perf_632plus.csv"
  )
)

write_csv(
  perf_apparent,
  file.path(
    dir_results_intermediate_sccore_gep ,
    "perf_apparent.csv"
  )
)

write_csv(
  gene_counts,
  file.path(
    dir_results_intermediate_sccore_gep ,
    "model_gene_counts.csv"
  )
)
#-------------------------------------------------------------------------------

