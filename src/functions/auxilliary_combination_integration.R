# Auxilliary functions for model evaluation:
# ------------------------------------------------------------------------------
# Bootstrap_validation_inner function: This function performs bootstrap validation of a model fitted on data ds_x. Only inner loop, no confidence intervals are provided.
# @param model_method: character indicating which type of model should be applied to the data
# @param ds_x: dataset with all of the features, columns should be features with colnames identifying the features
# @param ds_additional: dataset with outcomes of interest and matching informations
# @param ds_x_counts: counts dataset with all of the features, columns should be features with colnames identifying the features
# @param outcome_name: character vector corresponding to the column name of the outcome of interest in the ds_additional dataframe
# @param fup_name: character vector corresponding to the column name of the follow-up time of the outcome of interest in the ds_additional dataframe
# @param weights_name: character vector corresponding to the column name of the weights in the ds_additional dataframe
# @param tp: number indicating the timepoint of interest
# @param B: number of bootstrap validations
# @param seed: number indicating starting seed
# @param verbose: logical indicating if current bootstrap run should be indicated
# @param model_params: vector of model parameters
# @param bad_samples: logical indicating if bootstrap cross-validation is performed
# @param remove_bad_samples: logical indicating if samples with bad QCs should be
# removed during the training, if TRUE, separate performances will be computed for these samples
bootstrap_validation_inner <- function(feat_sel_method,topfeats,model_method,ds_x_list, ds_additional, ds_x_counts_list,ds_mutation,outcome_name,fup_name,weights_name,tp,B,seed,verbose=TRUE,model.params=NULL,bad_samples=NA,remove_bad_samples=F,genes_oi=NULL){
  
  # Set seed to ensure reproducibility
  set.seed(seed)
  # Make sure the correct alpha is used for the mutation and combined mutation model
  model.params_mut <- model.params
  mut_digits <- model.params$mnumbers
  model.params_mut$mut_alphas <- mut_digits
  
  # Model fitting and performance in the entire dataset
  # Get the Bailey and remaining model separately for the training split
  # save training data for each dataset
  fit_train_list <- list()
  predictions_list <- list()
  perf_list <- list()
  ds_x_train_list <- list()
  ds_y_train_list <- list()
  ds_x_counts_train_list <- list()
  
  #-------------------------------------------------------------------------------
  # 1. For SCCore-GEP go through DvP and non-DvP to get selected genes
  #-------------------------------------------------------------------------------
  for (k in seq_along(ds_x_list)) {
    
    ds_x <- ds_x_list[[k]]
    ds_x_counts <- ds_x_counts_list[[k]]
   
  # ----------------------------------------------------------------------------
  # 1.1 If required, remove samples with bad QCs.
  if (!all(is.na(bad_samples)) & remove_bad_samples){
    # Subset clinical dataset
    ds_y_train <- ds_additional[!(rownames(ds_additional) %in% bad_samples),]
    # Subset gene expression matrices
    ds_x_train <- ds_x[!(rownames(ds_x) %in% bad_samples),]
    ds_x_counts_train <- ds_x_counts[, !(colnames(ds_x_counts) %in% bad_samples)]
    # Create vector where to save results
    predicted_probs_badQCsamples <- performance_badQCsamples <- vector("list", B)

    ds_m <- ds_mutation[!(rownames(ds_mutation) %in% bad_samples),]
  } else {
    ds_y_train <- ds_additional
    ds_x_train <- ds_x
    ds_x_counts_train <- ds_x_counts
  }
    ds_x_train_list[[k]] <- ds_x_train
    ds_y_train_list[[k]] <- ds_y_train
    ds_x_counts_train_list[[k]] <- ds_x_counts_train
  # ----------------------------------------------------------------------------
  # 1.2 Build model on entire dataset of DvP or non-DvP
    fit_train_list[[k]] <- build_model(feat_selection_method = feat_sel_method,
                           topfeats = topfeats,
                           model_method = model_method,
                           ds_x = ds_x_train,
                           ds_y = ds_y_train,
                           ds_x_counts = ds_x_counts_train,
                           outcome_name = outcome_name,
                           fup_name = fup_name,
                           tp = tp,
                           model.params = model.params,
                           genes_oi = genes_oi,
                           seed = seed)
    
    set.seed(seed)
    
    predictions_list[[k]] <- fit_train_list[[k]]$predictions
    perf_list[[k]] <- get_performance_measures(predictions_list[[k]], ds_y_train, outcome_name, fup_name, weights_name, tp, fit_train_list[[k]]$recalibration_cox_lp)
    
  }
  #-------------------------------------------------------------------------------
  # 2. Combine selected genes or DvP+non DvP predictions
  #-------------------------------------------------------------------------------
  # 2.1.1 Combine the scores from Bailey and remaining
  meta_input <- data.frame(
    bailey = predictions_list[[1]],
    datadriven = predictions_list[[2]]
  )
  # 2.1.2 Combine the selected genes from Bailey and remaining
  meta_input_features_comb <- unique(c(names(fit_train_list[[1]]$coefficients),names(fit_train_list[[2]]$coefficients)))
  # Remove elements with value "b"
  meta_input_features_comb <- meta_input_features_comb[meta_input_features_comb != "(Intercept)"]
  # ----------------------------------------------------------------------------
  # 2.2 Make sure that the training data only has the combined selected genes when running the combined model
  ds_x_train_list_filtered <- lapply(ds_x_train_list, function(ds) {
    genes_to_keep <- intersect(colnames(ds), meta_input_features_comb)
    ds[, genes_to_keep, drop = FALSE]
  })
  # Merge datasets by columns (genes) — samples are identical
  ds_x_train_stackedGEP <- do.call(cbind, ds_x_train_list_filtered)
  # For each counts dataset, keep only genes that exist and are in combined_selected_genes_bs
  ds_x_counts_train_list_filtered <- lapply(ds_x_counts_train_list, function(ds) {
    genes_to_keep <- intersect(rownames(ds), meta_input_features_comb)
    ds[genes_to_keep, , drop = FALSE]
  })
  # Merge datasets by rows (genes) — samples are identical
  ds_x_train_counts_stackedGEP <- do.call(rbind, ds_x_counts_train_list_filtered)
  
  
  #-------------------------------------------------------------------------------
  # 3. Run the mutation model
  #-------------------------------------------------------------------------------
  fit_mutation <- build_model(feat_selection_method = "none",
                                     topfeats = topfeats,
                                     model_method = "coxnet",
                                     ds_x = ds_m,
                                     ds_y = ds_y_train,
                                     ds_x_counts = NULL,
                                     outcome_name = outcome_name,
                                     fup_name = fup_name,
                                     tp = tp,
                                     model.params = model.params_mut,
                                     genes_oi = genes_oi,
                                     seed = seed)
  mut_feats <- names(fit_mutation$coefficients)
  rm(model.params_mut)
  #-------------------------------------------------------------------------------
  # 3.1 Remove cellularity so it is not included in the WES+GEP model
  if (!is.null(model.params$cellularity_check)) {
    mut_feats <- mut_feats[!mut_feats %in% c("(Intercept)", "Tumor_cellularity")] #we dont want tumor cellularity combined with genes
  } else {
    mut_feats <- mut_feats[mut_feats != "(Intercept)"]
  }
  
  #-------------------------------------------------------------------------------
  # 4. Combine WES AND GEP features or keep GEP separated
  #-------------------------------------------------------------------------------
  # If we want to combine GEP and WES features in a model, the GEP and WES
  # features will be combined and put in the 4.1 stackedGEP model as input
  ds_m_aligned <- ds_m[rownames(ds_x_train_stackedGEP), , drop = FALSE]
  ds_m_filtered <- ds_m_aligned[, intersect(colnames(ds_m_aligned), mut_feats), drop = FALSE]
  if (model.params$combi_mut_type == "late_feature_mut"){
    ds_x_final <- cbind(
      ds_x_train_stackedGEP,
      ds_m_filtered
    )
    stopifnot(nrow(ds_x_final) == nrow(ds_x_train_stackedGEP))
  }
  
  
  #-------------------------------------------------------------------------------
  # 4.1 Make  DvP+non-DvP+WES feature model
 if (model.params$combi_type == "stackedGEP"){
    # Put the combined selected genes in a un-regularized coxnet model
  fit_stacked_train <- build_model(
    feat_selection_method = "none",
    topfeats = NA,
    model_method = model.params$model_final,          # force coxnet for stacking
    ds_x = ds_x_final,
    ds_y = ds_y_train,
    ds_x_counts = ds_x_train_counts_stackedGEP,
    outcome_name = outcome_name,
    fup_name = fup_name,
    tp = tp,
    model.params = model.params,
    seed = seed)
  }
  
  
  #-------------------------------------------------------------------------------
  # 5. Get model predictions and performance in the training cohort
  #-------------------------------------------------------------------------------
  predictions_model_orig <- fit_stacked_train$predictions
  # ds_y_train should be the same for bailey and remaining
  perf_orig_on_orig <- get_performance_measures(predictions_model_orig, ds_y_train, outcome_name, fup_name, weights_name, tp, fit_stacked_train$recalibration_cox_lp)
  feats_to_keep <- fit_stacked_train$preselected_feats
  # To keep average lambda and alpha for final model
  if (grepl("coxnet|glmnet|lasso|ridge", model_method)){
    lambda1se_list <- vector("list", B)
    if (grepl("coxnet|glmnet", model_method)) {
      alpha_list <- vector("list", B)
    }
  }
  #To keep average lambda and alpha for final model 2
  if (grepl("model", model.params$model_final)) {
    lambda1se_list2 <- vector("list", B)
    alpha_list2 <- vector("list", B)
    
  }
  lambda1se_listmut <- vector("list", B)
  base_model_method <- model_method
  
  # 5.2 Save intermediate steps:
  optimism_bootstraps <- apparent_bootstraps <- oob_bootstraps <- predicted_probs <-  bootstrap_coeffs <- lambda_list <- alpha_list <-
    dca_optimism_bootstraps <- dca_apparent_bootstraps <- dca_oob_bootstraps <- dca_bsmodel_on_bootstrap <- dca_bsmodel_on_oob <- dca_bsmodel_on_orig <-
    threshold_metrics_optimism_bootstraps <- threshold_metrics_apparent_bootstraps <- threshold_metrics_oob_bootstraps  <- stacked_bootstrap_models <- vector("list", B)
  
  # The following seed will be used as input if the bootstrap sample does not
  # have enough cases or controls
  s <- 1
  
  #-------------------------------------------------------------------------------
  # 6 Bootstrap repetitions:
  #-------------------------------------------------------------------------------
  # To save bootstrap models per dataset
  fit_bootstrap_all <- vector("list", length(ds_x_list))
  selected_feats_mutations <- list()
  selected_feats_stackedGEP <- list ()
  for (k in seq_along(ds_x_list)) {
    fit_bootstrap_all[[k]] <- vector("list", B)
  }
  for (brep in c(1:B)){
    if (grepl("coxnet|glmnet|lasso|ridge", model_method)){
      lambda1se_list[[brep]] <- list() 
      
      if (grepl("coxnet|glmnet", model_method)) {
        alpha_list[[brep]] <- list() 
        
      }
    }
    #save results for each dataset
    meta_input_bootstrap <- list()
    meta_input_bootstrap_GEP <- list()
    ds_x_train_oob_list <- list()
    ds_y_train_oob_list <- list()
    ds_x_counts_train_bootstrap_list <- list()
    ds_x_train_removed_badQCsamples_list <- list()
    ds_y_train_removed_badQCsamples_list <- list()
    ds_x_train_bootstrap_list <- list()
    ds_y_train_bootstrap_list <- list()
    fit_bootstrap_list <- list()
    bootstrap_coeffs[[brep]] <- list()
    #-------------------------------------------------------------------------------
    # 6.1 For SCCore-GEP go through DvP and non-DvP to get selected genes per bs
    for (k in seq_along(ds_x_list)) {
      ds_x_train <- ds_x_train_list[[k]]
      ds_y_train <- ds_y_train_list[[k]]
      ds_x_counts_train <- ds_x_counts_train_list[[k]]
      
    if(verbose){
      print(paste0("Calculating bootstrap repetition: ", brep))
    }
    set.seed(brep) # fixed issue set seed glmnet
    #-------------------------------------------------------------------------------
    # 6.2 Draw an UNMATCHED bootstrap sample (ROW-LEVEL, NO Set_id logic)
    #-------------------------------------------------------------------------------
    n_rows <- nrow(ds_y_train)
    # sample rows with replacement (true bootstrap, no pairing constraint)
    j <- sample(seq_len(n_rows), size = n_rows, replace = TRUE)
    # OOB is everything not selected at row level
    oob_idx <- setdiff(seq_len(n_rows), unique(j))
    s <- s + brep
    used <- "No"
    # Extract bootstrapping datasets
    ds_y_train_bootstrap <- ds_y_train[j, , drop = FALSE]
    ds_x_train_bootstrap <- as.data.frame(ds_x_train[j, ])
    ds_x_counts_train_bootstrap <- ds_x_counts_train[, j, drop = FALSE]
    
    # Extract true OOB datasets
    ds_x_train_oob <- ds_x_train[oob_idx, , drop = FALSE]
    ds_y_train_oob <- ds_y_train[oob_idx, , drop = FALSE]
    # NOTE: All Set_id-based balancing / pairing logic REMOVED intentionally
    # (no filtering, no reindexing, no pairing checks)
    ds_y_train_bootstrap[rownames(ds_y_train_bootstrap) %in% rownames(ds_x_train_bootstrap), ]

    # save the RNG seed state after sampling
    rng_before_gep <- .Random.seed
    #-------------------------------------------------------------------------------
    # 6.3 Build a model in the bootstrap samples of DvP or Non-DvP
    fit_bootstrap_list[[k]] <- build_model(feat_selection_method = feat_sel_method,
                                 topfeats = topfeats,
                                 model_method = model_method,
                                 ds_x = ds_x_train_bootstrap,
                                 ds_y = ds_y_train_bootstrap,
                                 ds_x_counts = ds_x_counts_train_bootstrap,
                                 outcome_name = outcome_name,
                                 fup_name = fup_name,
                                 tp = tp,
                                 model.params = model.params,
                                 genes_oi = genes_oi,
                                 seed = brep)
    # rerun bootstrap with new unused seed (not a perfect solution)
    new_seed = 100
    while (!is.null(fit_bootstrap_list[[k]][["warning_method"]])){
      rerun_sampling <- rerun_sampling_bootstrap(ds_y_train=ds_y_train,
                                                 outcome_name=outcome_name,
                                                 ds_x_train=ds_x_train,
                                                 ds_x_counts_train=ds_x_counts_train,
                                                 new_seed=new_seed,
                                                 brep=bseed,
                                                 verbose=verbose)
      
      ds_x_train_bootstrap <- rerun_sampling$ds_x_train_bootstrap
      ds_y_train_bootstrap <- rerun_sampling$ds_y_train_bootstrap
      ds_x_counts_train_bootstrap <- rerun_sampling$ds_x_counts_train_bootstrap
      ds_x_train_oob <- rerun_sampling$ds_x_train_oob
      ds_y_train_oob <- rerun_sampling$ds_y_train_oob
      set.seed(brep)
      # 6.3 Build a model in the bootstrap samples:
      fit_bootstrap_list[[k]] <- build_model(feat_selection_method = feat_sel_method,
                                   topfeats = topfeats,
                                   model_method = model_method,
                                   ds_x = ds_x_train_bootstrap,
                                   ds_y = ds_y_train_bootstrap,
                                   ds_x_counts = ds_x_counts_train_bootstrap,
                                   outcome_name = outcome_name,
                                   fup_name = fup_name,
                                   tp = tp,
                                   model.params = model.params,
                                   genes_oi = genes_oi,
                                   seed = brep)
      new_seed = new_seed+1
    }
    #save bootstrap models
    fit_bootstrap_all[[k]][[brep]] <- fit_bootstrap_list[[k]]
    #save bootstrap data for each dataset
    ds_x_train_bootstrap_list[[k]] <- ds_x_train_bootstrap
    ds_y_train_bootstrap_list[[k]] <- ds_y_train_bootstrap
    ds_x_train_oob_list[[k]] <- ds_x_train_oob
    ds_y_train_oob_list[[k]] <- ds_y_train_oob
    ds_x_counts_train_bootstrap_list[[k]] <- ds_x_counts_train_bootstrap
    
    #-------------------------------------------------------------------------------
    # 7. Combine selected genes and DvP+non DvP predictions per bootstrap
    #-------------------------------------------------------------------------------
    if (is.null(names(ds_x_list))) {
      names(ds_x_list) <- paste0("dataset_", seq_along(ds_x_list))
    }
    #save the predictions and coefficients from each dataset at a specific bootstrap
    meta_input_bootstrap[[names(ds_x_list)[k]]] <- fit_bootstrap_list[[k]]$predictions
    meta_input_bootstrap_GEP[[names(ds_x_list)[k]]] <- names(fit_bootstrap_list[[k]]$coefficients)
    # maintain consistency of seed
    set.seed(brep)
    # Bootstrap results
    bootstrap_coeffs[[brep]][[names(ds_x_list)[k]]] <- fit_bootstrap_list[[k]][["coefficients"]]
    # 3. Apply model in the bootstrap sample and estimate performance:
    predictions_model <- fit_bootstrap_list[[k]]$predictions
    #-------------------------------------------------------------------------------
    # 7.2 Save lambda and alpha  
    # Take lambda and alpha
    if (grepl("cox|glmnet|lasso|ridge", model_method)) {
      lambda1se_list[[brep]][[names(ds_x_list)[k]]] <- fit_bootstrap_list[[k]]$fitted_model$lambda.1se
      alpha_list[[brep]][[names(ds_x_list)[k]]] <- fit_bootstrap_list[[k]]$alpha
    }
   
    
    }
    # ----------------------------------------------------------------------------
    # 8 Outside DvP+non-DvP loop: GEP data will have only combined selected genes 
    #-----------------------------------------------------------------------------
    meta_input_bootstrap <- as.data.frame(meta_input_bootstrap)
    combined_selected_genes_bs <- unique(unlist(meta_input_bootstrap_GEP))
    # 4. Prepare the bootstrap x to only have combined selected genes
    # For each dataset, keep only genes that exist and are in combined_selected_genes_bs
    ds_x_boot_list_filtered <- lapply(ds_x_train_bootstrap_list, function(ds) {
      genes_to_keep <- intersect(colnames(ds), combined_selected_genes_bs)
      ds[, genes_to_keep, drop = FALSE]
    })
    # Merge datasets by columns (genes) — samples are identical
    ds_x_boot_stackedGEP <- do.call(cbind, ds_x_boot_list_filtered)
    # For each counts dataset, keep only genes that exist and are in combined_selected_genes_bs
    ds_x_counts_boot_list_filtered <- lapply(ds_x_counts_train_bootstrap_list, function(ds) {
      genes_to_keep <- intersect(rownames(ds), combined_selected_genes_bs)
      ds[genes_to_keep, , drop = FALSE]
    })
    # Merge datasets by rows (genes) — samples are identical
    ds_x_boot_counts_stackedGEP <- do.call(rbind, ds_x_counts_boot_list_filtered)
    set.seed(brep)
    
    # ----------------------------------------------------------------------------
    # 9 Run mutation data per bootstrap
    #-----------------------------------------------------------------------------
    # Make sure mutation data match the bootstrap data samples
    boot_ids <- rownames(ds_y_train_bootstrap)
    # strip suffixes to recover original IDs
    base_ids <- sub("\\.\\d+$", "", boot_ids)
    # match to ds_m
    idx <- match(base_ids, rownames(ds_m))
    # Safety check
    stopifnot(all(!is.na(idx)))
    # subset ds_m
    ds_m_bootstrap <- ds_m[idx, , drop = FALSE]
    # restore original (duplicated) rownames
    rownames(ds_m_bootstrap) <- boot_ids
    # Final check
    stopifnot(identical(rownames(ds_m_bootstrap), rownames(ds_y_train_bootstrap)))
    ### RUn for mutations
    model.params_mut <- model.params
    model.params_mut$mut_alphas <- mut_digits
    
    fit_mutation_bootstrap <- build_model(feat_selection_method = "none",#"mutation"
                                topfeats = topfeats,
                                model_method = "coxnet",
                                ds_x = ds_m_bootstrap,
                                ds_y = ds_y_train_bootstrap,
                                ds_x_counts = ds_x_counts_train_bootstrap,
                                outcome_name = outcome_name,
                                fup_name = fup_name,
                                tp = tp,
                                model.params = model.params_mut,
                                genes_oi = genes_oi,
                                seed = seed)
    selected_feats_mutations[[brep]] <- fit_mutation_bootstrap$coefficients
    mut_feats_bs <- names(fit_mutation_bootstrap$coefficients)
    rm(model.params_mut)
    
    # Make sure to also save the average mutation lambda
    lambda1se_listmut[[brep]]<- fit_mutation_bootstrap$fitted_model$lambda.1se
    #-------------------------------------------------------------------------------
    # 9.1 Remove cellularity so it is not included in the WES+GEP model
    if (!is.null(model.params$cellularity_check)) {
      mut_feats_bs <- mut_feats_bs[!mut_feats_bs %in% c("(Intercept)", "Tumor_cellularity")] #we dont want tumor cellularity combined with genes
    } else {
      mut_feats_bs <- mut_feats_bs[mut_feats_bs != "(Intercept)"]
    }
    
    #-------------------------------------------------------------------------------
    # 10. Combine WES AND GEP features or keep GEP separated
    #-------------------------------------------------------------------------------
    # If we want to combine GEP and WES features in a model, the GEP and WES
    # features will be combined and put in the 10.2 stackedGEP model as input
    # if we want to combine GEP and WES scores in a model, only GEP features
    # will put in the 10.2 stackedGEP model as input
    ds_m_aligned_bs <- ds_m_bootstrap[rownames(ds_x_boot_stackedGEP), , drop = FALSE]
    ds_m_filtered_bs <- ds_m_aligned_bs[, intersect(colnames(ds_m_aligned_bs), mut_feats_bs), drop = FALSE]
    if (model.params$combi_mut_type == "late_feature_mut"){
      ds_x_final_bs <- cbind(
        ds_x_boot_stackedGEP,
        ds_m_filtered_bs
      )
      stopifnot(nrow(ds_x_final_bs) == nrow(ds_x_boot_stackedGEP))
    }
    
    ## Random seed state of the "preb" seed is the same as used in the GEP model
    set.seed(brep)
    .Random.seed <- rng_before_gep
    #-------------------------------------------------------------------------------
    # 10.2 Make DvP+non-DvP+WES feature model
    if (model.params$combi_type == "stackedGEP"){
    fit_stacked_boot <- build_model(
      feat_selection_method = "none",
      topfeats = NA,
      model_method = model.params$model_final,          # force coxnet for stacking
      ds_x = ds_x_final_bs,
      ds_y = ds_y_train_bootstrap,
      ds_x_counts = NULL,
      outcome_name = outcome_name,
      fup_name = fup_name,
      tp = tp,
      model.params = model.params,
      seed = brep)
    selected_feats_stackedGEP[[brep]] <- fit_stacked_boot$coefficients
    }

    if (grepl("model", model.params$model_final)) {
      lambda1se_list2[[brep]]<- fit_stacked_boot$fitted_model$lambda.1se
      alpha_list2[[brep]] <- fit_stacked_boot$alpha
    }
    
    #Save stacked bootstrap models
    stacked_bootstrap_models[[brep]] <- fit_stacked_boot
    
    #-------------------------------------------------------------------------------
    # 11. Get model predictions and performance
    #-------------------------------------------------------------------------------
    # 11.1 When applying predictions, make sure the function knows which model was used
    if (model.params$model_final == "model2") {
      final_apply_model <- "coxnet"
    } 
    
    #-------------------------------------------------------------------------------
    # 11.2 Apply bootstrap model in the entire, OOB datasets, and estimate performances:
 if(model.params$combi_type == "stackedGEP"){
      #-------------------------------------------------------------------------------
      # 12. For GEP+WES features or scores model, compure performance
      #-------------------------------------------------------------------------------
      # 12.1 only keep genes rows from combined genes for OOB and apparent data
      ds_x_oob_list_filtered <- lapply(ds_x_train_oob_list, function(ds) {
        genes_to_keep <- intersect(colnames(ds), combined_selected_genes_bs)
        ds[, genes_to_keep, drop = FALSE]
      })
      # Combine them by columns (genes), since rows = samples
      ds_x_oob_combined <- do.call(cbind, ds_x_oob_list_filtered)
      ds_x_orig_list_filtered <- lapply(ds_x_train_list, function(ds) {
        genes_to_keep <- intersect(colnames(ds), combined_selected_genes_bs)
        ds[, genes_to_keep, drop = FALSE]
      })
      # Combine them by columns (genes)
      ds_x_orig_combined <- do.call(cbind, ds_x_orig_list_filtered)
      
      #-------------------------------------------------------------------------------
      # 12.2 Make mutation data x OOB or apparent samples
      all(rownames(ds_x_oob_combined) %in% rownames(ds_m))
      ds_m_aligned_OOB <- ds_m[rownames(ds_x_oob_combined), , drop = FALSE]
      all(rownames(ds_x_orig_combined) %in% rownames(ds_m))
      ds_m_aligned_APP <- ds_m[rownames(ds_x_orig_combined), , drop = FALSE]
      
      #-------------------------------------------------------------------------------
      # 13. If mutations and gep features are combined
      if (model.params$combi_mut_type == "late_feature_mut"){

        #-------------------------------------------------------------------------------
        # 13.1 Align mutation data to OOB samples and get predictions
    ds_m_filtered_OOB <- ds_m_aligned_OOB[, intersect(colnames(ds_m_aligned_OOB), mut_feats_bs), drop = FALSE]
    ds_x_m_OOB <- cbind(
      ds_x_oob_combined,
      ds_m_filtered_OOB
    )
    
    # Apply model on OOB to get predictions
    res_bs_model_on_oob <- apply_model(
      final_apply_model,
      fit_stacked_boot,
      ds_x_m_OOB,
      ds_y_train_oob,  # combined y
      tp,
      fit_stacked_boot,
      model.params
    )
    
    # Get predictions and metrics
    predictions_bs_model_on_oob <- res_bs_model_on_oob$predictions
    
    perf_bsmodel_on_oob <- get_performance_measures(
      predictions_bs_model_on_oob, ds_y_train_oob,
      outcome_name, fup_name, weights_name, tp,
      fit_stacked_boot$recalibration_cox_lp
    )
    #-------------------------------------------------------------------------------
    # 13.2 Align mutation data to App samples and get predictions
    ds_m_filtered_APP <- ds_m_aligned_APP[, intersect(colnames(ds_m_aligned_APP), mut_feats_bs), drop = FALSE]
    ds_x_m_APP <- cbind(
      ds_x_orig_combined,
      ds_m_filtered_APP
    )
    
    # Get predictions and metrics
    # Apply model
    res_bs_model_on_orig <- apply_model(
      final_apply_model,
      fit_stacked_boot,
      ds_x_m_APP,
      ds_y_train,   # combined y
      tp,
      fit_stacked_boot,
      model.params
    )
    
    predictions_bs_model_on_orig <- res_bs_model_on_orig$predictions
    perf_bsmodel_on_orig <- get_performance_measures(
      predictions_bs_model_on_orig, ds_y_train,
      outcome_name, fup_name, weights_name, tp,
      fit_stacked_boot$recalibration_cox_lp
    )
      }
   
    }
    # Apparent: test performance of stacked bootstrap model on bootstrap sample predictions
    perf_bsmodel_on_bs <- get_performance_measures(
      fit_stacked_boot$predictions,
      ds_y_train_bootstrap,
      outcome_name, fup_name, weights_name, tp,
      fit_stacked_boot$recalibration_cox_lp
    )
    #-------------------------------------------------------------------------------
    # 12. Compute .632+
    #-------------------------------------------------------------------------------
    
    # 12.1 Save probabilities
    predicted_probs[[brep]] <- predictions_bs_model_on_oob %>%
      as.data.frame() %>%
      `colnames<-`("predicted_probs") %>%
      mutate(boots_sample = brep)
    if ("SkylineDx.ID" %in% colnames(ds_y_train)){
      predicted_probs[[brep]] <- predicted_probs[[brep]] %>%
        mutate(SkylineDx.ID = ds_y_train[-j,]$SkylineDx.ID) %>%
        relocate(SkylineDx.ID, predicted_probs, boots_sample)
    }
    
    # 12.2 Compute optimism-corrected metrics
    optimism <- full_join(perf_bsmodel_on_bs, perf_bsmodel_on_orig,
                          by = c("variable", "group", "metric")) %>%
      mutate(value = value.x - value.y) %>%
      dplyr::select(-value.x, -value.y)
    optimism_bootstraps[[brep]] <- optimism %>% mutate(boots_sample = brep)
    apparent_bootstraps[[brep]] <- perf_bsmodel_on_bs %>% mutate(boots_sample = brep)
    oob_bootstraps[[brep]] <- perf_bsmodel_on_oob %>% mutate(boots_sample = brep)
    
  }
  #-------------------------------------------------------------------------------
  # 13. Save results of first apparent and bootstraps
  #-------------------------------------------------------------------------------
optimism_bootstraps <- optimism_bootstraps %>%
  bind_rows()
apparent_bootstraps <- apparent_bootstraps %>%
  bind_rows()
oob_bootstraps <- oob_bootstraps %>%
  bind_rows()
optimism_estimates <- optimism_bootstraps %>%
  group_by(variable, group, metric) %>%
  summarise(value = mean(value, na.rm = T), .groups = "keep") %>%
  ungroup()
oob_estimates <- oob_bootstraps %>%
  group_by(variable, group, metric) %>%
  summarise(value = mean(value, na.rm = T), .groups = "keep") %>%
  ungroup()
apparent_performance <- perf_orig_on_orig
int_val_estimates <- full_join(apparent_performance, optimism_estimates,
                               by = c("variable", "group", "metric")) %>%
  mutate(value = value.x - value.y) %>%
  dplyr::select(-value.x, -value.y)
## .632 bias correction method
int_val_estimates0632 <- full_join(apparent_performance, oob_estimates %>%
                                     dplyr::filter(!grepl('pval', metric)), #pval with single qoates
                                   by = c("variable", "group", "metric")) %>%
  mutate(value = 0.368*value.x + 0.632*value.y) %>%
  dplyr::select(-value.x, -value.y)
## .632+ bias correction method
random_vals <- ds_y_train %>%
  mutate(All = "All") %>%
  dplyr::select(all_of(c(outcome_name, fup_name, as.character(unique(apparent_performance$variable)), weights_name))) %>%
  rownames_to_column(var = "Id") %>%
  gather(variable, group, -all_of(c("Id", outcome_name, fup_name, weights_name))) %>%
  group_by(variable, group) %>%
  summarise(# Random 0/E should be 2 * incidence
    wfit_intercept = ifelse(length(unique(.data[[outcome_name]])) == 1, NA,
                            2 * (1 - summary(survfit(Surv(.data[[fup_name]], .data[[outcome_name]]) ~ 1, weights = .data[[weights_name]]), times = tp)$surv)),
    # Random 0/E should be 2 * incidence. This value is not completly correct for survival model as it has a slight bias.
    fit_intercept = ifelse(length(unique(.data[[outcome_name]])) == 1, NA,
                           2 * (1 - summary(survfit(Surv(.data[[fup_name]], .data[[outcome_name]]) ~ 1),  times = tp)$surv)),
    .groups = "keep") %>%
  ungroup() %>%
  mutate(fit_intercept_bin = 0,
         wfit_intercept_bin = 0,
         fit_slope = 0,
         fit_slope_bin = 0,
         wfit_slope = 0,
         wfit_slope_bin = 0,
         auc_model = 0.5,
         wauc_model = 0.5,
         cind_model = 0.5,
         wcind_model = 0.5,
         median_pair_diff = 0) %>%
  gather(metric, random_value, -variable, -group)
int_val_estimates0632plus <- oob_estimates %>%
  dplyr::filter(!grepl('pval', metric)) %>%
  full_join(apparent_performance,
            by = c("variable", "group", "metric")) %>%
  full_join(random_vals,
            by = c("variable", "group", "metric")) %>%
  mutate(R = (value.x - value.y)/(random_value - value.y),
         R = ifelse(R > 1, 1, ifelse(R < 0, 0, R)), # R should range between 0 and 1
         weight = 0.632 / (1 - 0.368 * R),
         value = (1 - weight) * value.y + weight * value.x) %>%
  dplyr::select(-value.x, -value.y, -R, -weight, -random_value)


#-------------------------------------------------------------------------------
# 13. Make final model with optimal lambda and alpha the same way as training model
#-------------------------------------------------------------------------------
fit_final_list <- list()
predictions_list_final <- list()
if (grepl("coxnet|glmnet|lasso|ridge", base_model_method)) {
    for (k in seq_along(ds_x_list)) {
      ds_x <- ds_x_list[[k]]
    set.seed(seed)
    vals <- sapply(lambda1se_list, function(b) b[[k]])
    model.params$avg_lambda1se <- mean(unlist(vals), na.rm = TRUE)
    vals_alpha <- sapply(alpha_list, function(b) b[[k]])
    model.params$avg_alpha <- mean(unlist(vals_alpha), na.rm = TRUE)
    
    # 1. Build a model in the entire dataset for Bailey and Remaining with optimal lambda and alpha
    set.seed(seed)
    fit_final_list[[k]]  <- build_model(feat_selection_method = feat_sel_method,
                                        topfeats = topfeats,
                                        model_method = base_model_method,
                                        ds_x = ds_x_train_list[[k]] ,
                                        ds_y = ds_y_train_list[[k]] ,
                                        ds_x_counts = ds_x_counts_train_list[[k]],
                                        outcome_name = outcome_name,
                                        fup_name = fup_name,
                                        tp = tp,
                                        model.params = model.params,
                                        genes_oi = genes_oi,
                                        seed = seed)
    
    predictions_list_final[[k]] <- fit_final_list[[k]]$predictions
    
     }
    
    meta_input_features_final <- list(
      bailey     = fit_final_list[[1]]$coefficients,
      datadriven = fit_final_list[[2]]$coefficients
    )
    meta_input_final <- data.frame(
      bailey = predictions_list_final[[1]],
      datadriven = predictions_list_final[[2]]
    )

    # 1.1 Combine the selected genes from final model and subset x
    meta_input_features_comb_final <- unique(c(names(fit_final_list[[1]]$coefficients),names(fit_final_list[[2]]$coefficients)))
    samps_to_keep <- rownames(ds_x_train)
    meta_input_features_comb_final <- meta_input_features_comb_final[meta_input_features_comb_final != "(Intercept)"]
    
    ds_x_train_list_filtered_final <- lapply(ds_x_train_list, function(ds) {
      genes_to_keep <- intersect(colnames(ds), meta_input_features_comb_final)
      ds[, genes_to_keep, drop = FALSE]
    })
    
    # Merge datasets by columns (genes) — samples are identical
    ds_x_train_stackedGEP_final <- do.call(cbind, ds_x_train_list_filtered_final)
    # For each counts dataset, keep only genes that exist and are in combined_selected_genes_bs
    ds_x_counts_train_list_filtered_final <- lapply(ds_x_counts_train_list, function(ds) {
      genes_to_keep <- intersect(rownames(ds), meta_input_features_comb_final)
      ds[genes_to_keep, , drop = FALSE]
    })
    # Merge datasets by rows (genes) — samples are identical
    ds_x_train_counts_stackedGEP_final <- do.call(rbind, ds_x_counts_train_list_filtered_final)
    set.seed(seed)
    model.params_mut <- model.params
    model.params_mut$mut_alphas <- mut_digits
    model.params_mut$avg_lambda1se <- mean(unlist(lambda1se_listmut), na.rm = TRUE)
    model.params_mut$avg_alpha <- mean(mut_digits)
    fit_mutation_opt <- build_model(feat_selection_method = "none",
                                topfeats = topfeats,
                                model_method = "coxnet",
                                ds_x = ds_m,
                                ds_y = ds_y_train,
                                ds_x_counts = ds_x_counts_train,
                                outcome_name = outcome_name,
                                fup_name = fup_name,
                                tp = tp,
                                model.params = model.params_mut,
                                genes_oi = genes_oi,
                                seed = seed)
    rm(model.params_mut)
    mut_feats_opt <- names(fit_mutation_opt$coefficients)
    mut_feats_opt <- mut_feats_opt[mut_feats_opt != "(Intercept)"]
    if (!is.null(model.params$cellularity_check)) {
      mut_feats_opt <- mut_feats_opt[!mut_feats_opt %in% c("(Intercept)", "Tumor_cellularity")] #we dont want tumor cellularity combined with genes
    } else {
      mut_feats_opt <- mut_feats_opt[mut_feats_opt != "(Intercept)"]
    }
    
    ds_m_aligned_opt <- ds_m[rownames(ds_x_train_stackedGEP_final), , drop = FALSE]
    ds_m_filtered_opt <- ds_m_aligned_opt[, intersect(colnames(ds_m_aligned_opt), mut_feats_opt), drop = FALSE]
    
    if (model.params$combi_mut_type == "late_feature_mut"){
      ds_x_final_opt <- cbind(
        ds_x_train_stackedGEP_final,
        ds_m_filtered_opt
      )
      stopifnot(nrow(ds_x_final_opt) == nrow(ds_x_train_stackedGEP_final))
      head(ds_x_final_opt)
    }
    # 1.2 Run final coxnet model on combined selected genes
    if ( model.params$model_final == "model2" |  model.params$model_final == "model3"){
      model.params$avg_lambda1se2 <- mean(unlist(lambda1se_list2), na.rm = TRUE)
      model.params$avg_alpha2     <- mean(unlist(alpha_list2), na.rm = TRUE)
    }
    #2 Run the second coxnet model
    if (model.params$combi_type  == "stackedGEP"){
      combined <- meta_input_features_final
    fit_final2 <- build_model(
      feat_selection_method = "none",
      topfeats = NA,
      model_method = model.params$model_final,          # force coxnet for stacking
      ds_x = ds_x_final_opt,
      ds_y = ds_y_train,
      ds_x_counts = ds_x_train_counts_stackedGEP_final,
      outcome_name = outcome_name,
      fup_name = fup_name,
      tp = tp,
      model.params = model.params,
      seed = seed)
    }
    
    # 2. Get model predictions and performance in the training cohort
    predictions_model_orig_opt <- fit_final2$predictions
    
    }

# 10. Results
res <- list("fit_bootstrap_all" = fit_bootstrap_all,
            "Apparent_performance"= perf_orig_on_orig,
            
            "Optimism-Corrected" = int_val_estimates,
            "Optimism-Corrected0632" = int_val_estimates0632,
            "Optimism-Corrected0632plus" = int_val_estimates0632plus,
            "Optimism" = optimism_estimates,
            "All bootstrap optimism" = optimism_bootstraps,
            "All apparent bootstrap"= apparent_bootstraps,
            "All out of bag"= oob_bootstraps,
            "Model_fit" = fit_stacked_train,
            "Model_fit_pre" = fit_train_list,
            "OOB_predicted_probs" = predicted_probs %>% bind_rows(),
            "Bootstrap_coefficients" = bootstrap_coeffs,
            "stacked_bootstrap_models" = stacked_bootstrap_models,
            "selected_feats_mutations" = selected_feats_mutations,
            "selected_feats_stackedGEP" = selected_feats_stackedGEP)
# If removed bad QCs samples, save performances in this subset of samples
'if (!all(is.na(bad_samples)) & remove_bad_samples){
    res <- c(res,
             list("Bad_QCs_samples" = performance_badQCsamples %>% dplyr::select(-value, -boots_sample) %>% unique(),
                  "All_bad_QCs_samples" = performance_badQCsamples %>% dplyr::select(-avg_value),
                  "Bad_QCs_samples_predicted_probs" = predicted_probs_badQCsamples %>% bind_rows(),
                  "Training_model_on_bad_QCs_samples" = performance_train_badQCsamples))
  }'
# If fitted, save final model
if (grepl("coxnet|glmnet|lasso|ridge", base_model_method)) {
  res$final_selected_genes_WES <-  fit_mutation_opt$coefficients
  res$lambda_list2 <- lambda1se_list2
  res$lambda_listmut <- lambda1se_listmut
  res$lambda1se_list <- lambda1se_list
  res$alpha_list2 <- alpha_list2
  res$alpha_list <- alpha_list
  res$Final_model2 <- fit_final2
  res$Final_model_list <- fit_final_list
  res$final_selected_genes_GEP <- fit_final2$coefficients

}

  # Return results
  return(res)
  
}


# Bootstrap_validation function: This function performs bootstrap validation of a model fitted on data ds_x. It can produce confidence intervals for optimism-corrected estimates, if cis = TRUE and the number of outer loops are provided.
# @param model_method: character indicating which type of model should be applied to the data
# @param ds_x: dataset with all of the features, columns should be features with colnames identifying the features
# @param ds_additional: dataset with outcomes of interest and matching informations
# @param ds_x_counts: counts dataset with all of the features, columns should be features with colnames identifying the features
# @param outcome_name: character vector corresponding to the column name of the outcome of interest in the ds_additional dataframe
# @param fup_name: character vector corresponding to the column name of the follow-up time of the outcome of interest in the ds_additional dataframe
# @param weights_name: character vector corresponding to the column name of the weights in the ds_additional dataframe
# @param tp: number indicating the timepoint of interest
# @param B_inner: number of bootstrap repetitions in the inner loop
# @param seed: number indicating starting seed
# @param cis: logical indicating whether confidence intervals should be computed (TRUE) or not (FALSE)
# @param B_outer: number of bootstrap repetitions in the outer loop, for computation of confidence intervals
# @param model_params: vector of model parameters
# @param bad_samples: logical indicating if bootstrap cross-validation is performed
# @param remove_bad_samples: logical indicating if samples with bad QCs should be removed during the training,
# if TRUE, separate performances will be computed for these samples
# @param genes_oi: list of genes of interest
bootstrap_validation <- function(feat_sel_method,topfeats,model_method,ds_x_list, ds_additional, ds_x_counts_list,ds_mutation,outcome_name,fup_name,weights_name,tp,B_inner,seed,cis=FALSE,B_outer=NULL,model.params=NULL,bad_samples=NA,remove_bad_samples=F,genes_oi=NULL,dir_results = NULL){
  argg <- as.list(environment())
  argg <- argg[setdiff(names(argg),c("ds_x_list", "ds_additional", "ds_x_counts_list"))]  # Set seed to ensure reproducibility
  set.seed(seed)
  if (is.null(B_outer)){
    input_arg_check(ds_x_list=ds_x_list,ds_y=ds_additional,ds_x_counts_list=ds_x_counts_list,
                    outcome_name,fup_name,weights_name)
    
    bootstrap_results = bootstrap_validation_inner(feat_sel_method = feat_sel_method,
                                                   topfeats = topfeats,
                                                   model_method = model_method,
                                                   ds_x_list = ds_x_list,
                                                   ds_additional = ds_additional,
                                                   ds_x_counts_list = ds_x_counts_list,
                                                   ds_mutation = ds_mutation,
                                                   outcome_name = outcome_name,
                                                   fup_name = fup_name,
                                                   weights_name = weights_name,
                                                   tp = tp,
                                                   B = B_inner,
                                                   seed = seed,
                                                   model.params = model.params,
                                                   bad_samples = bad_samples,
                                                   remove_bad_samples = remove_bad_samples)
  }
  return(c(bootstrap_results, "Bootstrap_parameters" = argg))
}


# input_arg_check function: This function checks the input argument of the user
# @param ds_x_list: list of bailey and remaining datasets with all of the features, columns should be features with colnames identifying the features
# @param ds_additional: dataset with outcomes of interest and matching informations
# @param ds_x_counts_list: list of bailey and remaining counts datasets with all of the features, columns should be features with colnames identifying the features
# @param outcome_name: character vector corresponding to the column name of the outcome of interest in the ds_additional dataframe
# @param fup_name: character vector corresponding to the column name of the follow-up time of the outcome of interest in the ds_additional dataframe
# @param weights_name: column name of the weights in ds_y
input_arg_check <- function(ds_x_list, ds_y, ds_x_counts_list,
                            outcome_name, fup_name, weights_name) {
  
  # ─────────────────────────────────────────────
  # 1️⃣ Clinical / outcome data checks
  # ─────────────────────────────────────────────
  if (!(outcome_name %in% colnames(ds_y))) {
    stop(paste0("Outcome_name: ", outcome_name, " not found in clinical data (ds_y)."))
  }
  if (!(fup_name %in% colnames(ds_y))) {
    stop(paste0("Fup_name: ", fup_name, " not found in clinical data (ds_y)."))
  }
  if (!(weights_name %in% colnames(ds_y))) {
    stop(paste0("Weights_name: ", weights_name, " not found in clinical data (ds_y)."))
  }
  
  # ─────────────────────────────────────────────
  # 2️⃣ Check that lists are consistent
  # ─────────────────────────────────────────────
  if (length(ds_x_list) != length(ds_x_counts_list)) {
    stop("ds_x_list and ds_x_counts_list must have the same number of elements.")
  }
  
  # ─────────────────────────────────────────────
  # 3️⃣ Iterate over each dataset pair
  # ─────────────────────────────────────────────
  for (i in seq_along(ds_x_list)) {
    message(paste0("Checking dataset ", i, " of ", length(ds_x_list), "..."))
    ds_x <- ds_x_list[[i]]
    ds_x_counts <- ds_x_counts_list[[i]]
    
    if (!is.null(ds_x_counts)){
      if(!(all(colnames(ds_x_counts) == ds_y$SkylineDx.ID)) & identical(colnames(ds_x_counts), ds_y$SkylineDx.ID)){
        stop(paste0("Dataset ", i, ": sample mismatch between ds_x (rows) and ds_y (SkylineDx.ID)."))
      }
      if(!(all(colnames(ds_x) == rownames(ds_x_counts)))){
        stop(paste0("Dataset ", i, ": feature mismatch between ds_x (cols) and ds_x_counts (rows)."))
      }
    }
    if(!(all(rownames(ds_x) == ds_y$SkylineDx.ID)) & identical(rownames(ds_x), ds_y$SkylineDx.ID)){
      stop(paste0("Dataset ", i, ": sample mismatch in ds_x_counts (columns) vs ds_y (SkylineDx.ID)."))
    }
   
    
    message(paste0("Dataset ", i, " passed input checks."))
  }
  
  message("All datasets passed input checks.")
}

# get_performance_measures function: This function computes performance metrics for predictions given in model_predictions
# @param model_predictions: numeric vector with risk predictions from the model (higher prediction = higher risk)
# @param ds_additional: dataset with outcomes of interest and matching informations
# @param outcome: character vector corresponding to the column name of the outcome of interest in the ds dataframe
# @param outcome_FU: character vector corresponding to the column name of the follow-up time of the outcome of interest in the ds dataframe
# @param weights_name: character vector corresponding to the column name of the weights in the ds dataframe
# @param tp: number indicating the timepoint of interest
# @param lps: numeric vector with linear predictors from the model
get_performance_measures <- function(model_predictions,
                                     ds,
                                     outcome,
                                     outcome_FU,
                                     weights_name,
                                     tp,
                                     lps) {
  
  stopifnot(length(model_predictions) == nrow(ds))
  
  # Stratification variables (UNCHANGED)
  strat_vars <- c(
    "Biopsy_excision_2f", "AJCC_8", "BWH", "AJCC_bin", "BWH_bin",
    "AJCC_bin2", "BWH_bin2", "Sex", "Differentiation", "PNI_bin",
    "Lymphovascular_invasion_bin", "PNI_or_LVI", "Tissue_involvement",
    "Tumor_location", "HM_at_cSCC", "OTR_at_cSCC",
    "Number_of_cSCC_before_culprit_bin",
    "AJCC_and_Number_of_prior_cSCC", "BWH_and_Number_of_prior_cSCC",
    "Resection_margin_cat", "Immunosuppressed",
    "Percent_trimmed_quartiles", "Seq_depth_quartiles",
    "Procedure_number", "ESTIMATE_tumorpurity_quartiles"
  )
  
  # Predictions (UNCHANGED)
  if (!is.null(lps)) {
    stopifnot(length(lps) == nrow(ds))
    model_predictions <- data.frame(pred_score = model_predictions, lp = lps)
  } else {
    model_predictions <- data.frame(pred_score = model_predictions)
  }
  
  # Build dataset (UNCHANGED except naming safety)
  base_cols <- c(outcome, outcome_FU, weights_name)
  base_cols <- base_cols[base_cols %in% colnames(ds)]
  
  pred_ds <- cbind(
    model_predictions,
    ds[, base_cols, drop = FALSE]
  )
  
  colnames(pred_ds)[(ncol(model_predictions) + 1):ncol(pred_ds)] <-
    c("outcome", "outcome_FU", "weight")[seq_along(base_cols)]
  
  if ("SkylineDx.ID" %in% colnames(ds)) {
    pred_ds$SkylineDx.ID <- ds$SkylineDx.ID
  }
  
  # Stratification (UNCHANGED logic)
  if (any(strat_vars %in% colnames(ds))) {
    
    strat_vars <- intersect(strat_vars, colnames(ds))
    
    pred_ds <- pred_ds %>%
      cbind(ds[, strat_vars, drop = FALSE]) %>%
      tidyr::gather(variable, group,
                    -pred_score, -outcome, -outcome_FU,
                    -SkylineDx.ID, -weight) %>%
      dplyr::bind_rows(
        pred_ds %>%
          dplyr::mutate(variable = "All", group = "All")
      )
    
    strat_vars_levels <- c("All", strat_vars)
    
  } else {
    
    pred_ds <- pred_ds %>%
      dplyr::mutate(variable = "All", group = "All")
    
    strat_vars_levels <- c("All")
  }
  
  # Paired metric removed (Set_id removed)
  median_dif_ds <- NA
  
  # Metrics
  perf_ds <- pred_ds %>%
    dplyr::group_by(variable, group) %>%
    dplyr::summarise(
      
      auc_model = ifelse(length(unique(outcome)) == 1, NA,
                         as.numeric(pROC::auc(
                           pROC::roc(outcome, pred_score,
                                     quiet = TRUE,
                                     direction = "<")
                         ))),
      
      # FIXED: safe weights inside group (NO tryCatch, NO silent masking)
      
      wauc_model = {
        
        w <- weight
        valid <- !is.na(w) & w > 0
        
        if (length(unique(outcome)) == 1 || sum(valid) < 2 || length(unique(outcome[valid])) < 2) {
          NA
        } else {
          WeightedAUC(
            WeightedROC(
              pred_score[valid],
              outcome[valid],
              weight = w[valid]
            )
          )
        }
      },
      
      cind_model = ifelse(length(unique(outcome)) == 1,
                          NA,
                          1 - rcorr.cens(pred_score,
                                         survival::Surv(outcome_FU, outcome))[1]),
      
      wcind_model = ifelse(length(unique(outcome)) == 1 | all(is.na(weight)),
                           NA,
                           intsurv::cIndex(outcome_FU,
                                           event = outcome,
                                           pred_score,
                                           weight = weight)[[1]]),
      
      .groups = "keep"
    ) %>%
    dplyr::ungroup()
  
  # Final reshape
  perf_ds <- perf_ds %>%
    tidyr::gather(metric, value, -variable, -group) %>%
    dplyr::mutate(variable = factor(variable, levels = strat_vars_levels))
  
  perf_ds
}