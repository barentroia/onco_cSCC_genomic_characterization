#===============================================================================
# Aim:      Combine the outer bootstrap batches of SCCore-GEP, remove
#           duplicated bootstrap repetitions, compute the median and 95% CI of
#           the .632+ performance estimates and merge them with the
#           performance of the final SCCore-GEP model.
# Author:   R. Ruiter
# Input:    - Outer bootstrap batches (dir_results_intermediate_sccore_gep/CI/*)
#           - Final SCCore-GEP model (coxnet_Late_feature_integration_robj.rds)
# Output:   - corr_632plus_CI_final_all_variables.csv in
#             dir_results_intermediate_sccore_gep: .632+ performance of
#             SCCore-GEP with median and 95% CI per (sub)group
# Requires: dir_results_intermediate_sccore_gep, experiment_sccore_gep (main.R)
# Usage:    Sourced from main.R after the outer bootstrap (step 5.3)
#===============================================================================

# 1. Settings
#-------------------------------------------------------------------------------
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflicts_prefer(dplyr::filter)
#-------------------------------------------------------------------------------

# 2. Load all CI batches
#-------------------------------------------------------------------------------
batches_list <- Sys.glob(
  file.path(dir_results_intermediate_sccore_gep, "CI/*")
)

#' Name of a batch, taken from its folder name
#'
#' @param path character. Folder path of a batch.
#'
#' @return character. Folder name.
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


#-------------------------------------------------------------------------------

# 3. Remove potential duplicates
#-------------------------------------------------------------------------------
#
# Each bootstrap replicate (n_rep) should only occur once
# across all input models.
#
# The order of models in models_list determines which model
# keeps a duplicated n_rep.
#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------

# 4. Combine cleaned 0.632+ outer bootstrap results
#-------------------------------------------------------------------------------
combined_632plus <- bind_rows(
  models_list_632plus_clean,
  .id = "source_model"
)

#-------------------------------------------------------------------------------

# 5. Calculate median and 95% CI
#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------

# 6. Load whole SCCore-GEP
#-------------------------------------------------------------------------------
GEP23 <- readRDS(file.path(dir_results_intermediate_sccore_gep, experiment_sccore_gep, "coxnet_Late_feature_integration_robj.rds"))

#-------------------------------------------------------------------------------

# 7. Extract .632+ performance of SCCore-GEP
#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------

# 8. Merge CI with SCCore-GEP
#-------------------------------------------------------------------------------
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

# 9. Write final CV
#---------------------------------------------------
write_csv(
  plot_df_632plus,
  file.path(
    dir_results_intermediate_sccore_gep,
    "corr_632plus_CI_final_all_variables.csv"
  )
)
#-------------------------------------------------------------------------------


