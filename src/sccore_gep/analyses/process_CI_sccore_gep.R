#---------------------------------------------------
# Aim: Prepare confidence interval SCCore-GEP performance
# Author: R.Ruiter
# Input: all batches of confidence performance of SCCore-GEP
# Output: One performance file for confidence performance of SCCore-GEP
#---------------------------------------------------
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflicts_prefer(dplyr::filter)
#---------------------------------------------------
# 1. Load all CI batches
#---------------------------------------------------

batches_list <- Sys.glob(
  file.path(dir_results_intermediate_sccore_gep, "CI/*")
)


# Extract model name from folder path
model_name <- function(path) {
  basename(normalizePath(path))
}


# Load the first .robj.rds file found in each folder
models_list <- map(batches_list, function(folder) {
  
  message("Processing: ", folder)
  
  # Find only existing .robj.rds files
  rds_file <- list.files(
    folder,
    pattern = "robj\\.rds$",
    full.names = TRUE
  )
  
  # If no file is found, return NULL
  if (length(rds_file) == 0) {
    warning("No .robj.rds file found in: ", folder)
    return(NULL)
  }
  
  # Read first matching RDS file
  readRDS(rds_file[1])
})


# Assign model names
names(models_list) <- map_chr(batches_list, model_name)


# Remove models where no RDS file was found
models_list <- models_list[!sapply(models_list, is.null)]


#---------------------------------------------------
# 2. Remove potential duplicates
#---------------------------------------------------
#
# Each bootstrap replicate (n_rep) should only occur once
# across all input models.
#
# The order of models in models_list determines which model
# keeps a duplicated n_rep.
#---------------------------------------------------

seen <- c()


models_list_632plus_clean <- map(models_list, function(model) {
  
  # Extract 0.632+ bootstrap estimates
  df_632plus <- model$`Estimates Optimism-Corrected0632plus`
  
  # Check that the expected object exists
  if (is.null(df_632plus)) {
    warning("Missing 'Estimates Optimism-Corrected0632plus' in model")
    return(NULL)
  }
  
  # Keep only n_rep values not already encountered
  df_632plus_clean <- df_632plus %>%
    filter(!(n_rep %in% seen))
  
  # Update global list of already-used n_rep values
  seen <<- c(seen, df_632plus_clean$n_rep)
  
  return(df_632plus_clean)
})


# Remove NULL entries
models_list_632plus_clean <- models_list_632plus_clean[
  !sapply(models_list_632plus_clean, is.null)
]


#---------------------------------------------------
# 3. combine cleaned 0.632+ outer bootstrap results
#---------------------------------------------------

combined_632plus <- bind_rows(
  models_list_632plus_clean,
  .id = "source_model"
)


#---------------------------------------------------
# 4. Calculate median and 95% CI
#---------------------------------------------------

bootstrap_632plus_summary <- combined_632plus %>%
  group_by(variable, group, metric) %>%
  summarise(
    median = median(value, na.rm = TRUE),
    ci_025 = quantile(
      value,
      probs = 0.025,
      na.rm = TRUE,
      names = FALSE
    ),
    ci_975 = quantile(
      value,
      probs = 0.975,
      na.rm = TRUE,
      names = FALSE
    ),
    .groups = "drop"
  )


#---------------------------------------------------
# 5. Load whole SCCore-GEP
#---------------------------------------------------

GEP23 <- readRDS(file.path(dir_results_intermediate_sccore_gep, experiment_sccore_gep, "coxnet_Late_feature_integration_robj.rds"))


#---------------------------------------------------
# 6. Extract 632+ of SCCor-GEP
#---------------------------------------------------

final_vals_632plus <- GEP23$`Optimism-Corrected0632plus` %>%
  filter(
    metric %in% c("wcind_model", "wauc_model")
  ) %>%
  select(
    variable,
    group,
    metric,
    final_value = value
  )


#---------------------------------------------------
# 7. Merge CI with SCCore-GEP
#---------------------------------------------------

plot_df_632plus <- bootstrap_632plus_summary %>%
  filter(
    metric %in% c("wcind_model", "wauc_model")
  ) %>%
  left_join(
    final_vals_632plus,
    by = c("variable", "group", "metric")
  ) %>%
  select(
    variable,
    group,
    metric,
    median,
    ci_025,
    ci_975,
    final_value
  )


#---------------------------------------------------
# 8. Write final CV
#---------------------------------------------------
write_csv(
  plot_df_632plus,
  file.path(
    dir_results_intermediate_sccore_gep,
    "corr_632plus_CI_final_all_variables.csv"
  )
)


