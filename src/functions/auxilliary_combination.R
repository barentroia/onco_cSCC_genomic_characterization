# Auxilliary functions for model evaluation:
# ------------------------------------------------------------------------------
# Bootstrap_validation_inner function: This function performs bootstrap validation of a model fitted on data ds_x. Only inner loop, no confidence intervals are provided.
# @param feat_sel_method: character indicating which type of model should be applied to the data
# @param topfeats: number to indicate the amount of top features in feature selection, default NA
# @param ds_x_list: list of datasets of DvP and non_DvPwith all of the features, columns should be features with colnames identifying the features
# @param ds_additional: dataset with outcomes of interest and matching informations
# @param ds_x_counts_list: list of count datasets of DvP and non_DvPwith all of the features, columns should be features with colnames identifying the features
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
# @param genes_oi: genes of interest to keep in feature selection, on default not used
# removed during the training, if TRUE, separate performances will be computed for these samples
bootstrap_validation_inner <- function(feat_sel_method,topfeats,model_method,ds_x_list, ds_additional, ds_x_counts_list,outcome_name,fup_name,weights_name,tp,B,seed,verbose=TRUE,model.params=NULL,bad_samples=NA,remove_bad_samples=F,genes_oi=NULL){
  
  # Set seed to ensure reproducibility
  set.seed(seed)
  #-------------------------------------------------------------------------------
  
  # 1. Model fitting and performance in the entire dataset to calculate apparent
  #-------------------------------------------------------------------------------
  # Model fitting and performance in the entire dataset
  # save training data for each dataset
  fit_train_list <- list()
  predictions_list <- list()
  perf_list <- list()
  ds_x_train_list <- list()
  ds_y_train_list <- list()
  ds_x_counts_train_list <- list()
  # 1.1. Get the DvP and non_DvPmodel separately for the training split
  for (k in seq_along(ds_x_list)) {
    
    ds_x <- ds_x_list[[k]]
    ds_x_counts <- ds_x_counts_list[[k]]
   
  # 1.2 Remove samples with bad QCs.
  # Subset clinical dataset
    ds_y_train <- ds_additional[!(rownames(ds_additional) %in% bad_samples),]
    # Subset gene expression matrices
    ds_x_train <- ds_x[!(rownames(ds_x) %in% bad_samples),]
    ds_x_counts_train <- ds_x_counts[, !(colnames(ds_x_counts) %in% bad_samples)]
    # Create vector where to save results
    predicted_probs_badQCsamples <- performance_badQCsamples <- vector("list", B)

    ds_x_train_list[[k]] <- ds_x_train
    ds_y_train_list[[k]] <- ds_y_train
    ds_x_counts_train_list[[k]] <- ds_x_counts_train
  # 1.3. Build a model in the entire dataset
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
    perf_list[[k]] <- get_performance_measures(predictions_list[[k]], ds_y_train, outcome_name, fup_name, weights_name, tp, fit_train_list[[k]]$recalibration_cox_lp,model.params)
    
  }
  # 1.4.1 Combine the predictions from DvP and non_DvP
  meta_input <- data.frame(
    DvP = predictions_list[[1]],
    non_DvP= predictions_list[[2]]
  )

  # 1.4.2 Combine the features from DvP and non_DvP
  meta_input_features_comb <- unique(c(names(fit_train_list[[1]]$coefficients),names(fit_train_list[[2]]$coefficients)))
  meta_input_features_comb <- meta_input_features_comb[meta_input_features_comb != "(Intercept)"]
  # Make sure that the training data only has the combined selected genes 
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
  
  
  # 1.5 run the late prediction integration or late feature integration model
  if (model.params$combi_type == "Stacked") {
    # Put the predictions in an un-regularized model
    fit_stacked_train <- build_model(
      feat_selection_method = "none",
      topfeats = NA,
      model_method = model.params$model_final,
      ds_x = meta_input,
      ds_y = ds_y_train,
      ds_x_counts = NULL,
      outcome_name = outcome_name,
      fup_name = fup_name,
      tp = tp,
      model.params = model.params,
      seed = seed) 
  }else if (model.params$combi_type == "stackedGEP"){
    # Put the combined selected genes in an un-regularized model
  fit_stacked_train <- build_model(
    feat_selection_method = "none",
    topfeats = NA,
    model_method = model.params$model_final, 
    ds_x = ds_x_train_stackedGEP,
    ds_y = ds_y_train,
    ds_x_counts = ds_x_train_counts_stackedGEP,
    outcome_name = outcome_name,
    fup_name = fup_name,
    tp = tp,
    model.params = model.params,
    seed = seed)
  }
  # 1.6 Get model predictions and performance in the training cohort
  predictions_model_orig <- fit_stacked_train$predictions
  perf_orig_on_orig <- get_performance_measures(predictions_model_orig, ds_y_train, outcome_name, fup_name, weights_name, tp, fit_stacked_train$recalibration_cox_lp,model.params)
  feats_to_keep <- fit_stacked_train$preselected_feats
  
  # 1.7 Save model name and initiate alpha and lambda list
  #  To keep average lambda and alpha for final model
  if (grepl("coxnet|glmnet|lasso|ridge", model_method)){
    lambda1se_list <- vector("list", B)
    if (grepl("coxnet|glmnet", model_method)) {
      alpha_list <- vector("list", B)
    }}
  #To keep average lambda and alpha for the stacked final model
  if (grepl("model", model.params$model_final)) {
    lambda1se_list2 <- vector("list", B)
    alpha_list2 <- vector("list", B)
  }
  # Save the model name for the final model for  the correct model name
  base_model_method <- model_method
  
  #-------------------------------------------------------------------------------
  
  # 2. Model fitting and performance on bootstrapping
  #-------------------------------------------------------------------------------
  # 2.1. Save intermediate steps:
  optimism_bootstraps <- apparent_bootstraps <- oob_bootstraps <- predicted_probs <-  bootstrap_coeffs <- lambda_list <- alpha_list <- stacked_bootstrap_models <- vector("list", B)
  # The following seed will be used as input if the bootstrap sample does not
  # have enough cases or controls
  s <- 1
  
  # To save bootstrap models per dataset
  fit_bootstrap_all <- vector("list", length(ds_x_list))
  for (k in seq_along(ds_x_list)) {
    fit_bootstrap_all[[k]] <- vector("list", B)
  }
  
  # 2.2 Bootstrap repetitions:
  # ----------------------------------------------------------------------------
  for (brep in c(1:B)){
    if (grepl("coxnet|glmnet|lasso|ridge", model_method)){
      lambda1se_list[[brep]] <- list() 
      if (grepl("coxnet|glmnet", model_method)) {
        alpha_list[[brep]] <- list() 
      }}
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
    
    # 2.3 Run the bootstrap model training for each data set: DvP and non-DvP
    # --------------------------------------------------------------------------
    for (k in seq_along(ds_x_list)) {
      ds_x_train <- ds_x_train_list[[k]]
      ds_y_train <- ds_y_train_list[[k]]
      ds_x_counts_train <- ds_x_counts_train_list[[k]]
    if(verbose){
      print(paste0("Calculating bootstrap repetition: ", brep))
    }
    set.seed(brep) 

    # 2.4. Draw a bootstrap sample with replacement from the sample:
    # 2.4.1 with matched samples
    if (is.null( model.params$unmatched)){
    # Sample set ids:
    set_ids_sampled <- sample(unique(as.numeric(as.character(ds_y_train[,"Set_id"]))),replace=T)
    # Find indexes of sampled pairs:
    j <- unlist(purrr::map(set_ids_sampled, function(x) which(ds_y_train[,"Set_id"]==x)))
    # OOB is all rows with Set_id not in sampled Set_ids
    oob_idx <- which(!(ds_y_train$Set_id %in% set_ids_sampled))
    # Check that there are enough cases and controls in the bootstrap sample
    # and at least 1 case and 1 control in the OOB sample
    s <- s + brep
    used <- "No"
    while(any(table(ds_y_train[j,][[outcome_name]]) < 2) |
          length(unique(ds_y_train[j,][[outcome_name]])) == 1 |
          length(unique(ds_y_train[-j,][[outcome_name]])) == 1){

      # Sample set ids:
      set_ids_sampled <- sample(unique(as.numeric(as.character(ds_y_train[,"Set_id"]))),replace=T)
      # Find indexes of sampled pairs:
      j <- unlist(purrr::map(set_ids_sampled, function(x) which(ds_y_train[,"Set_id"]==x)))
      s <- s + 1
      oob_idx <- which(!(ds_y_train$Set_id %in% set_ids_sampled))
      
    }
    # Extract bootstrapping datasets
    ds_y_train_bootstrap <- ds_y_train[j,]
    ds_y_train_bootstrap$Set_id <- as.numeric(as.character(ds_y_train_bootstrap$Set_id))
    ds_x_train_bootstrap <- as.data.frame(ds_x_train[j,])
    ds_x_counts_train_bootstrap <- ds_x_counts_train[,j]
    # Extract true OOB datasets
    ds_x_train_oob <- ds_x_train[oob_idx, , drop = FALSE]
    ds_y_train_oob <- ds_y_train[oob_idx, , drop = FALSE]

    ## Exclude unmatched samples
    if (exists("ds_y_train_unmatched")){
      # Create index variable which will be used later to re-align the datasets
      ds_y_train_bootstrap$index <- seq_len(nrow(ds_y_train_bootstrap))
      ds_y_train_bootstrap_unmatched <- ds_y_train_bootstrap %>%
        dplyr::filter(Set_id %in% as.numeric(as.character(ds_y_train_unmatched$Set_id)))
      ds_y_train_bootstrap <- ds_y_train_bootstrap %>%
        dplyr::filter(!(Set_id %in% as.numeric(as.character(ds_y_train_unmatched$Set_id))))
    }
    ## Check that individuals in a pair are next to each other
    for (i in seq(1, nrow(ds_y_train_bootstrap)/2, 2)){
      stopifnot(ds_y_train_bootstrap$Set_id[[i]] == ds_y_train_bootstrap$Set_id[[i+1]])
    }
    ds_y_train_bootstrap <- ds_y_train_bootstrap %>%
      mutate(Set_id = as.factor(rep(1:(nrow(.)/2), each = 2)))
    ## Re-add unmatched samples and adjust their set ID
    if (exists("ds_y_train_unmatched")){
      ds_y_train_bootstrap <- ds_y_train_bootstrap %>%
        mutate(Set_id = as.numeric(as.character(Set_id))) %>%
        rbind(ds_y_train_bootstrap_unmatched %>%
                mutate(Set_id = max(as.numeric(as.character(ds_y_train_bootstrap$Set_id))) +
                         seq_len(nrow(ds_y_train_bootstrap_unmatched)))) %>%
        arrange(desc(index))
      # Align gene expression and clinical datasets
      ds_x_train_bootstrap <- ds_x_train_bootstrap[ds_y_train_bootstrap$index,]
      ds_x_counts_train_bootstrap <- ds_x_counts_train_bootstrap[, ds_y_train_bootstrap$index]
      stopifnot(all(rownames(ds_x_train_bootstrap) == gsub("\\..*", "", rownames(ds_y_train_bootstrap))))
      stopifnot(all(colnames(ds_x_counts_train_bootstrap) == rownames(ds_y_train_bootstrap)))
    }
    ds_y_train_bootstrap[rownames(ds_y_train_bootstrap) %in% rownames(ds_x_train_bootstrap), ]
    # 2.4.2 with unmatched samples
    }else{
      # All Set_id logic, pairing checks, filtering, and re-indexing removed intentionally
      n_rows <- nrow(ds_y_train)
      # bootstrap sample (row-level, replacement)
      j <- sample(seq_len(n_rows), size = n_rows, replace = TRUE)
      # OOB = rows not sampled at least once
      oob_idx <- setdiff(seq_len(n_rows), unique(j))
      s <- s + brep
      used <- "No"
      # bootstrap datasets
      ds_y_train_bootstrap <- ds_y_train[j, , drop = FALSE]
      ds_x_train_bootstrap <- as.data.frame(ds_x_train[j, ])
      ds_x_counts_train_bootstrap <- ds_x_counts_train[, j, drop = FALSE]
      # OOB datasets
      ds_x_train_oob <- ds_x_train[oob_idx, , drop = FALSE]
      ds_y_train_oob <- ds_y_train[oob_idx, , drop = FALSE]
      
      
    }
    
    # save the RNG seed state after sampling
    rng_before_gep <- .Random.seed
    
    # 2.5. Build a model in the bootstrap samples:
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
    
    # 2.6 save bootstrap results
    # Save bootstrap models
    fit_bootstrap_all[[k]][[brep]] <- fit_bootstrap_list[[k]]
    #save bootstrap data for each dataset
    ds_x_train_bootstrap_list[[k]] <- ds_x_train_bootstrap
    ds_y_train_bootstrap_list[[k]] <- ds_y_train_bootstrap
    ds_x_train_oob_list[[k]] <- ds_x_train_oob
    ds_y_train_oob_list[[k]] <- ds_y_train_oob
    ds_x_counts_train_bootstrap_list[[k]] <- ds_x_counts_train_bootstrap
    if (is.null(names(ds_x_list))) {
      names(ds_x_list) <- paste0("dataset_", seq_along(ds_x_list))
    }
    #save the predictions and coefficients from each dataset at a specific bootstrap
    meta_input_bootstrap[[names(ds_x_list)[k]]] <- fit_bootstrap_list[[k]]$predictions
    meta_input_bootstrap_GEP[[names(ds_x_list)[k]]] <- names(fit_bootstrap_list[[k]]$coefficients)
    # Bootstrap results
    bootstrap_coeffs[[brep]][[names(ds_x_list)[k]]] <- fit_bootstrap_list[[k]][["coefficients"]]
    feats_to_keep_bootstrap <- fit_bootstrap_list[[k]]$preselected_feats
    # Take lambda and alpha
    if (grepl("cox|glmnet|lasso|ridge", model_method)) {
      lambda1se_list[[brep]][[names(ds_x_list)[k]]] <- fit_bootstrap_list[[k]]$fitted_model$lambda.1se
      alpha_list[[brep]][[names(ds_x_list)[k]]] <- fit_bootstrap_list[[k]]$alpha
    }
    # maintain consistency of seed
    set.seed(brep)
    }
    
    # 2.7 Completed bootstrapping for each dataset. 
    #     Prepare the input for the stacked model.
    # The combined predicitons
    meta_input_bootstrap <- as.data.frame(meta_input_bootstrap)
    # The combined selected genes
    combined_selected_genes_bs <- unique(unlist(meta_input_bootstrap_GEP))
    # Prepare the bootstrap training data to only have the combined selected genes
    ds_x_boot_list_filtered <- lapply(ds_x_train_bootstrap_list, function(ds) {
      genes_to_keep <- intersect(colnames(ds), combined_selected_genes_bs)
      ds[, genes_to_keep, drop = FALSE]
    })
    ds_x_boot_stackedGEP <- do.call(cbind, ds_x_boot_list_filtered)
    ds_x_counts_boot_list_filtered <- lapply(ds_x_counts_train_bootstrap_list, function(ds) {
      genes_to_keep <- intersect(rownames(ds), combined_selected_genes_bs)
      ds[genes_to_keep, , drop = FALSE]
    })
    ds_x_boot_counts_stackedGEP <- do.call(rbind, ds_x_counts_boot_list_filtered)
    
    ## Random seed state of the "preb" seed is the same as used in the GEP model
    set.seed(brep)
    .Random.seed <- rng_before_gep

    # 2.8. Run combined model
    if (model.params$combi_type == "Stacked") {
      # Put the predictions in a un-regularized model
      fit_stacked_boot <- build_model(
        feat_selection_method = "none",
        topfeats = NA,
        model_method = model.params$model_final,
        ds_x = meta_input_bootstrap,
        ds_y = ds_y_train_bootstrap,  
        ds_x_counts = NULL,
        outcome_name = outcome_name,
        fup_name = fup_name,
        tp = tp,
        model.params = model.params,
        seed = brep)
    }else if (model.params$combi_type == "stackedGEP"){
      # Put the selected genes in a un-regularized model
    fit_stacked_boot <- build_model(
      feat_selection_method = "none",
      topfeats = NA,
      model_method = model.params$model_final,
      ds_x = ds_x_boot_stackedGEP,
      ds_y = ds_y_train_bootstrap,
      ds_x_counts = ds_x_boot_counts_stackedGEP,
      outcome_name = outcome_name,
      fup_name = fup_name,
      tp = tp,
      model.params = model.params,
      seed = brep)
}
    # 2.9. Save combined model results
    # lambda and alpha for the combined model
    if (grepl("model", model.params$model_final)) {
      lambda1se_list2[[brep]]<-fit_stacked_boot$fitted_model$lambda.1se
      alpha_list2[[brep]] <- fit_stacked_boot$alpha
    }
    # Save combined bootstrap models
    stacked_bootstrap_models[[brep]] <- fit_stacked_boot
    # Save model name used for the final model for signature
    base_model <-  model.params$model_final
    if (base_model == "model2") {
      final_apply_model <- "coxnet"
    } else if (base_model == "RSF") {
      final_apply_model <- "RSF"
      } 
    
    # 2.10. Apply combined bootstrap model in the entire, OOB datasets and, and estimate performances:
    if (model.params$combi_type == "Stacked") {
      # 2.10.1 You cannot directly use stacked model on OOB y, because model expects predictions, not original X, as input
      # Get predictions of OOB and apparent 
      pred_bailey <- apply_model(model_method, fit_bootstrap_list[[1]], 
                                            ds_x_train_oob_list[[1]], ds_y_train_oob_list[[1]],  
                                            tp, fit_bootstrap_list[[1]], model.params)$predictions
      pred_datadriven <- apply_model(model_method, fit_bootstrap_list[[2]], 
                                                ds_x_train_oob_list[[2]], ds_y_train_oob_list[[2]], 
                                                tp, fit_bootstrap_list[[2]], model.params)$predictions
      # stacked OOB predictions
      # With OOB predictions together, we can use it on the stacked model
      meta_input_oob <- data.frame(
        DvP = pred_bailey,
        non_DvP= pred_datadriven)
      # Run base GEP model training samples and get stacked training predictions
      meta_input_orig <- data.frame(
        DvP = apply_model(
          model_method, fit_bootstrap_list[[1]],
          ds_x_train_list[[1]],
          ds_y_train_list[[1]], tp, fit_bootstrap_list[[1]], model.params
        )$predictions,
        non_DvP= apply_model(
          model_method, fit_bootstrap_list[[2]],
          ds_x_train_list[[2]],
          ds_y_train_list[[2]], tp, fit_bootstrap_list[[2]], model.params
        )$predictions)
      # 2.10.2: Put predictions in stacked inbag model
      # OOB performance
      # Test performance of stacked model on combined oob predictions
      res_bs_model_on_oob <- apply_model(
        model_method = final_apply_model, fit_stacked_boot,
        meta_input_oob,
        ds_y_train_oob, tp, fit_stacked_boot, model.params)
      predictions_bs_model_on_oob <- res_bs_model_on_oob$predictions
      perf_bsmodel_on_oob <- get_performance_measures(
        predictions_bs_model_on_oob, ds_y_train_oob,
        outcome_name, fup_name, weights_name, tp,
        fit_stacked_boot$recalibration_cox_lp,model.params)
      
      # Original dataset (optimism estimation)
      # Test performance of stacked model on combined training predictions
      perf_bs_model_on_orig <- apply_model(
        model_method = final_apply_model, fit_stacked_boot,
        meta_input_orig,
        ds_y_train, tp, fit_stacked_boot, model.params)
      predictions_bs_model_on_orig <- perf_bs_model_on_orig$predictions
      perf_bsmodel_on_orig <- get_performance_measures(
        predictions_bs_model_on_orig, ds_y_train,
        outcome_name, fup_name, weights_name, tp,
        fit_stacked_boot$recalibration_cox_lp,model.params)
    }else if(model.params$combi_type == "stackedGEP"){
    # 2.11.1 Process OOB data for late feature integration prediction
    #  only keep OOB gene rows from combined genes
    ds_x_oob_list_filtered <- lapply(ds_x_train_oob_list, function(ds) {
      genes_to_keep <- intersect(colnames(ds), combined_selected_genes_bs)
      ds[, genes_to_keep, drop = FALSE]})
    # Combine them by columns (genes), since rows = samples
    ds_x_oob_combined <- do.call(cbind, ds_x_oob_list_filtered)
    
    # Apply model on OOB to get predictions
    res_bs_model_on_oob <- apply_model(
      final_apply_model,
      fit_stacked_boot,
      ds_x_oob_combined,
      ds_y_train_oob,  # combined y
      tp,
      fit_stacked_boot,
      model.params)
    predictions_bs_model_on_oob <- res_bs_model_on_oob$predictions
    # OOB performance
    perf_bsmodel_on_oob <- get_performance_measures(
      predictions_bs_model_on_oob, ds_y_train_oob,
      outcome_name, fup_name, weights_name, tp,
      fit_stacked_boot$recalibration_cox_lp,model.params)
    
    # 2.11.2 Bootstrap apparent performance
    # only keep Orig gene rows from combined genes
    ds_x_orig_list_filtered <- lapply(ds_x_train_list, function(ds) {
      genes_to_keep <- intersect(colnames(ds), combined_selected_genes_bs)
      ds[, genes_to_keep, drop = FALSE]})
    # Combine them by columns (genes)
    ds_x_orig_combined <- do.call(cbind, ds_x_orig_list_filtered)
    # Get predictions and metrics
    res_bs_model_on_orig <- apply_model(
      final_apply_model,
      fit_stacked_boot,
      ds_x_orig_combined,
      ds_y_train,   # combined y
      tp,
      fit_stacked_boot,
      model.params)
    predictions_bs_model_on_orig <- res_bs_model_on_orig$predictions
    perf_bsmodel_on_orig <- get_performance_measures(
      predictions_bs_model_on_orig, ds_y_train,
      outcome_name, fup_name, weights_name, tp,
      fit_stacked_boot$recalibration_cox_lp,model.params
    )}
    
    # 2.12. Apparent performance
    #: test performance of stacked bootstrap model on bootstrap sample predictions
    perf_bsmodel_on_bs <- get_performance_measures(
      fit_stacked_boot$predictions,
      ds_y_train_bootstrap,
      outcome_name, fup_name, weights_name, tp,
      fit_stacked_boot$recalibration_cox_lp,model.params
    )

    # 2.13. Save probabilities and corrected performance
    predicted_probs[[brep]] <- predictions_bs_model_on_oob %>%
      as.data.frame() %>%
      `colnames<-`("predicted_probs") %>%
      mutate(boots_sample = brep)
    if ("SkylineDx.ID" %in% colnames(ds_y_train)){
      predicted_probs[[brep]] <- predicted_probs[[brep]] %>%
        mutate(SkylineDx.ID = ds_y_train[-j,]$SkylineDx.ID) %>%
        relocate(SkylineDx.ID, predicted_probs, boots_sample)
    }
    # Compute optimism-corrected metrics
    optimism <- full_join(perf_bsmodel_on_bs, perf_bsmodel_on_orig,
                          by = c("variable", "group", "metric")) %>%
      mutate(value = value.x - value.y) %>%
      dplyr::select(-value.x, -value.y)
    # Save performance for each bootstrap
    optimism_bootstraps[[brep]] <- optimism %>% mutate(boots_sample = brep)
    apparent_bootstraps[[brep]] <- perf_bsmodel_on_bs %>% mutate(boots_sample = brep)
    oob_bootstraps[[brep]] <- perf_bsmodel_on_oob %>% mutate(boots_sample = brep)
    
}
  #-------------------------------------------------------------------------------
  
  # 3. Combine results and get corrected performance
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
## 3.1 .632 bias correction method
int_val_estimates0632 <- full_join(apparent_performance, oob_estimates %>%
                                     dplyr::filter(!grepl('pval', metric)), #pval with single qoates
                                   by = c("variable", "group", "metric")) %>%
  mutate(value = 0.368*value.x + 0.632*value.y) %>%
  dplyr::select(-value.x, -value.y)
## 3.2 .632+ bias correction method
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

# 4. Make regularized cox with average lambda and alpha for signature
#-------------------------------------------------------------------------------
fit_final_list <- list()
predictions_list_final <- list()
if (grepl("coxnet|glmnet|lasso|ridge", base_model_method)) {
  # 4.1. Make a final model with average lambda and alpha for each dataset
    for (k in seq_along(ds_x_list)) {
      ds_x <- ds_x_list[[k]]
    set.seed(seed)
    vals <- sapply(lambda1se_list, function(b) b[[k]])
    model.params$avg_lambda1se <- mean(unlist(vals), na.rm = TRUE)
    vals_alpha <- sapply(alpha_list, function(b) b[[k]])
    model.params$avg_alpha <- mean(unlist(vals_alpha), na.rm = TRUE)
    # 4.2. Build a model in the entire dataset for DvP and non_DvPwith optimal lambda and alpha
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
    # 4.3. Save predictions and selected features
    predictions_list_final[[k]] <- fit_final_list[[k]]$predictions}
    meta_input_features_final <- list(
      DvP     = fit_final_list[[1]]$coefficients,
      non_DvP= fit_final_list[[2]]$coefficients
    )
    meta_input_final <- data.frame(
      DvP = predictions_list_final[[1]],
      non_DvP= predictions_list_final[[2]]
    )
    # 4.4. Combine the selected genes from final model and subset x
    meta_input_features_comb_final <- unique(c(names(fit_final_list[[1]]$coefficients),names(fit_final_list[[2]]$coefficients)))
    samps_to_keep <- rownames(ds_x_train)
    meta_input_features_comb_final <- meta_input_features_comb_final[meta_input_features_comb_final != "(Intercept)"]
    # For each TPM dataset, keep only genes that exist and are in combined_selected_genes_bs
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
    
    # 4.5. Run final coxnet model on combined selected genes or predictions
    # Get average alpha and lambda from stacked models
    if ( model.params$model_final == "model2"){
      model.params$avg_lambda1se2 <- mean(unlist(lambda1se_list2), na.rm = TRUE)
      model.params$avg_alpha2     <- mean(unlist(alpha_list2), na.rm = TRUE)
      
    }
    set.seed(seed)
    
    # 4.6. Run the combined coxnet model
    if (model.params$combi_type == "Stacked") {
      combined <- meta_input_final
      # Put the predictions in an un-regularized  model
      fit_final2 <- build_model(
        feat_selection_method = "none",
        topfeats = NA,
        model_method = model.params$model_final,  
        ds_x = meta_input_final,
        ds_y = ds_y_train,
        ds_x_counts = NULL,
        outcome_name = outcome_name,
        fup_name = fup_name,
        tp = tp,
        model.params = model.params,
        seed = seed) 
    }else if (model.params$combi_type  == "stackedGEP"){
      combined <- meta_input_features_final
      # Put the selected features in an un-regularized  model
    fit_final2 <- build_model(
      feat_selection_method = "none",
      topfeats = NA,
      model_method = model.params$model_final,        
      ds_x = ds_x_train_stackedGEP_final,
      ds_y = ds_y_train,
      ds_x_counts = ds_x_train_counts_stackedGEP_final,
      outcome_name = outcome_name,
      fup_name = fup_name,
      tp = tp,
      model.params = model.params,
      seed = seed)
    } }

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
            "stacked_bootstrap_models" = stacked_bootstrap_models)
# If fitted, save final model
if (grepl("coxnet|glmnet|lasso|ridge", base_model_method)) {
  res$lambda_list2 <- lambda1se_list2
  res$lambda1se_list <- lambda1se_list
  res$alpha_list2 <- alpha_list2
  res$alpha_list <- alpha_list
  res$Final_model2 <- fit_final2
  res$Final_model_list <- fit_final_list
  res$final_selected_genes = combined
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
bootstrap_validation <- function(feat_sel_method,topfeats,model_method,ds_x_list, ds_additional, ds_x_counts_list,outcome_name,fup_name,weights_name,tp,B_inner,seed,cis=FALSE,B_outer=NULL,model.params=NULL,bad_samples=NA,remove_bad_samples=F,genes_oi=NULL,dir_results = NULL){
  argg <- as.list(environment())
  argg <- argg[setdiff(names(argg),c("ds_x_list", "ds_additional", "ds_x_counts_list"))]  # Set seed to ensure reproducibility
  set.seed(seed)

  if (cis & !is.null(B_outer)){
    bootstrap_results_all_outer = bootstrap_oc_estimates = bootstrap_oc0632_estimates = bootstrap_oc0632plus_estimates = vector(mode = "list", length = B_outer)
    # assuming batches is done in 10
    B_batch <- B_outer - 9
    count <- 1
    for (out_rep in c(B_batch:B_outer)){
      print(paste0("Outer bootstrap number ",out_rep, "/",B_outer))
      # Draw a bootstrap sample with replacement from the original sample:
      # ------------------------------
      # Step 0: remove bad samples
      # ------------------------------
      ds_additional <- ds_additional[!(rownames(ds_additional) %in% bad_samples), ]
      ds_x_list <- lapply(ds_x_list, function(df) df[!(rownames(df) %in% bad_samples), , drop = FALSE])
      ds_x_counts_list <- lapply(ds_x_counts_list, function(df) df[, !(colnames(df) %in% bad_samples), drop = FALSE])
      # ------------------------------
      # Step 1: bootstrap sampling
      # ------------------------------
      clean_names <- function(x) gsub("\\.\\d+$", "", x)
      # Sample
      # set seed again. If you dont set seed than each sampling will take a random state of seed 42
      set.seed(seed+out_rep) # make sure each sampling is different per outer rep, but independent and reproducible
      set_ids_sampled <- sample(unique(as.numeric(as.character(ds_additional[,"Set_id"]))), replace = TRUE)
      j <- unlist(purrr::map(set_ids_sampled, function(x) which(ds_additional[,"Set_id"] == x)))
      # 1. Basic bootstrap sampling
      ds_additional_bootstrap <- ds_additional[j, , drop = FALSE]
      ds_additional_bootstrap$Set_id <- as.numeric(as.character(ds_additional_bootstrap$Set_id))
      ds_additional_bootstrap$index <- seq_len(nrow(ds_additional_bootstrap))  # <-- FIX
      ds_x_list_bootstrap <- lapply(ds_x_list, function(df) as.data.frame(df[j, , drop = FALSE]))
      ds_x_counts_list_bootstrap <- lapply(ds_x_counts_list, function(df) df[, j, drop = FALSE])
      # for CI calculation
      
      invisible(capture.output(
        bootstrap_results_inner <- suppressMessages(
          bootstrap_validation_inner(
            feat_sel_method = feat_sel_method,
            topfeats = topfeats,
            model_method = model_method,
            ds_x_list = ds_x_list_bootstrap,
            ds_additional = ds_additional_bootstrap,
            ds_x_counts_list = ds_x_counts_list_bootstrap,
            outcome_name = outcome_name,
            fup_name = fup_name,
            weights_name = weights_name,
            tp = tp,
            B = B_inner,
            seed = seed,
            verbose = F,
            model.params = model.params,
            bad_samples = bad_samples,
            remove_bad_samples = remove_bad_samples
          )
        )
      ))
      
      
      
      bootstrap_oc0632plus_estimates[[out_rep]] = bootstrap_results_inner[["Optimism-Corrected0632plus"]] %>% mutate(n_rep = out_rep)
  
      # save results in a folder so results are already there
      bootstrap_oc0632plus_estimates_int = bootstrap_oc0632plus_estimates %>% bind_rows()
      # Clean up large objects to save memory if needed
      bootstrap_results_inner[["Model_fit"]] <- NULL
    }
    
    # Combine outer bootstrap loop optimism corrected performances
    bootstrap_oc0632plus_estimates = bootstrap_oc0632plus_estimates %>% bind_rows()

    # Compute confidence intervals
    bootstrap_oc0632plus_summary = bootstrap_oc0632plus_estimates %>%
      group_by(variable, group, metric) %>%
      summarise(mean = mean(value, na.rm = T),
                ci_025 = quantile(value, probs = 0.025, na.rm = T),
                ci_975 = quantile(value, probs = 0.975, na.rm = T),
                .groups = "keep") %>%
      ungroup()
    
    bootstrap_results <- list(
      "CIs Optimism-Corrected0632plus" = bootstrap_oc0632plus_summary,# if not run in batches
      "Estimates Optimism-Corrected0632plus" = bootstrap_oc0632plus_estimates # to combine with other batches
    )
    
  } else{
    input_arg_check(ds_x_list=ds_x_list,ds_y=ds_additional,ds_x_counts_list=ds_x_counts_list,
                    outcome_name,fup_name,weights_name)
    
    bootstrap_results = bootstrap_validation_inner(feat_sel_method = feat_sel_method,
                                                   topfeats = topfeats,
                                                   model_method = model_method,
                                                   ds_x_list = ds_x_list,
                                                   ds_additional = ds_additional,
                                                   ds_x_counts_list = ds_x_counts_list,
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
  #  Clinical / outcome data checks
  if (!(outcome_name %in% colnames(ds_y))) {
    stop(paste0("Outcome_name: ", outcome_name, " not found in clinical data (ds_y)."))
  }
  if (!(fup_name %in% colnames(ds_y))) {
    stop(paste0("Fup_name: ", fup_name, " not found in clinical data (ds_y)."))
  }
  if (!(weights_name %in% colnames(ds_y))) {
    stop(paste0("Weights_name: ", weights_name, " not found in clinical data (ds_y)."))
  }
  #  Check that lists are consistent
  if (length(ds_x_list) != length(ds_x_counts_list)) {
    stop("ds_x_list and ds_x_counts_list must have the same number of elements.")
  }
  #  Iterate over each dataset pair
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
get_performance_measures <- function(model_predictions, ds, outcome, outcome_FU, weights_name, tp, lps,model.params){
  
  stopifnot(length(model_predictions) == nrow(ds))
  stopifnot(names(model_predictions) == rownames(ds))
  
  
  # Identify stratification variables
  strat_vars <- c("Biopsy_excision_2f", "AJCC_8", "BWH", "AJCC_bin", "BWH_bin", "AJCC_bin2", "BWH_bin2",
                  "Sex", "Differentiation", "PNI_bin", "Lymphovascular_invasion_bin",
                  "PNI_or_LVI", "Tissue_involvement", "Tumor_location", "HM_at_cSCC",
                  "OTR_at_cSCC", "Number_of_cSCC_before_culprit_bin",
                  "AJCC_and_Number_of_prior_cSCC", "BWH_and_Number_of_prior_cSCC",
                  "Resection_margin_cat", "Immunosuppressed", "Percent_trimmed_quartiles",
                  "Seq_depth_quartiles", #"Year_of_obtained_material_bin",
                  "Procedure_number", "ESTIMATE_tumorpurity_quartiles")
  

  
  
  # If present, include linear predictors
  if (!is.null(lps)){
    stopifnot(length(lps) == nrow(ds))
    stopifnot(names(lps) == rownames(ds))
    model_predictions <- data.frame(pred_score = model_predictions,
                                    lp = lps)
  } else {
    model_predictions <- as.data.frame(model_predictions)
    colnames(model_predictions) <- "pred_score"
  }
  
  
  # Dataframe with predictions
  if (is.null( model.params$unmatched)){
  if ("SkylineDx.ID" %in% colnames(ds)){ # random bootstrap model does not have a skylinedx.id
    pred_ds <- model_predictions %>%
      cbind(ds %>%
              dplyr::select(all_of(c("Set_id", "SkylineDx.ID", outcome, outcome_FU, weights_name))) %>%
              `colnames<-`(c("Set_id", "SkylineDx.ID", "outcome", "outcome_FU", "weight")))
  } else {
    pred_ds <- model_predictions %>%
      cbind(ds %>%
              dplyr::select(all_of(c("Set_id", outcome, outcome_FU))) %>%
              `colnames<-`(c("Set_id", "outcome", "outcome_FU")))
  }}else{
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
    
  }
  
  if (any(strat_vars %in% colnames(ds))){
    # In the case stratification variables are available, performances are
    # computed for each stratification variable available and the entire dataset
    strat_vars <- intersect(strat_vars, colnames(ds))
    if (is.null( model.params$unmatched)){
    pred_ds <- pred_ds %>%
      cbind(ds %>%
              dplyr::select(all_of(strat_vars))) %>%
      gather(variable, group, -Set_id, -pred_score, -outcome, -outcome_FU, -SkylineDx.ID, -weight) %>%
      rbind(pred_ds  %>%
              mutate(variable = "All",
                     group = "All"))
    }else{
      pred_ds <- pred_ds %>%
      cbind(ds[, strat_vars, drop = FALSE]) %>%
        tidyr::gather(variable, group,
                      -pred_score, -outcome, -outcome_FU,
                      -SkylineDx.ID, -weight) %>%
        dplyr::bind_rows(
          pred_ds %>%
            dplyr::mutate(variable = "All", group = "All")
        )
      
    }
    strat_vars_levels <- c("All", strat_vars)
  } else {
    # In the case none of stratification variables are available, performances are
    # computed only for the entire dataset
    pred_ds <- pred_ds %>%
      dplyr::mutate(variable = "All",
             group = "All")
    strat_vars_levels <- c("All")
  }
  if (is.null( model.params$unmatched)){
  # Compute difference in predictions cases-controls per pair
  median_dif_ds <- tryCatch({pred_ds %>%
      dplyr::select(Set_id, variable, group, outcome, pred_score) %>%
      spread(outcome, pred_score) %>%
      group_by(variable, group, Set_id) %>%
      summarise(diff = `1` - `0`, .groups = "keep") %>%
      ungroup() %>%
      group_by(variable, group) %>%
      summarise(median_pair_diff = median(diff, na.rm = T), .groups = "keep") %>%
      ungroup()}, error = function(e){NA})
  }else{
    # Paired metric removed (Set_id removed)
    median_dif_ds <- NA
  }
  
  # Performance matched samples
  if (is.null( model.params$unmatched)){
  bad_idx <- which(!(pred_ds$weight > 0))
  if (length(bad_idx) > 0) {
    print("Some weights are not > 0 (or are NA):")
    print(pred_ds[bad_idx, ])
  }

  # Compute metrics:
  perf_ds <- pred_ds %>%
    group_by(variable, group) %>%
    summarise(auc_model = ifelse(length(unique(outcome)) == 1, NA,
                                 as.numeric(gsub(".*: ", "", pROC::roc(outcome, pred_score, quiet = TRUE, direction = "<")$auc))),
              wauc_model = ifelse(length(unique(outcome)) == 1 | all(is.na(weight)), NA, WeightedAUC(WeightedROC(pred_score, outcome, weight = weight))),
              cind_model = ifelse(length(unique(outcome)) == 1, NA,
                                  ifelse(is.nan(1 - rcorr.cens(pred_score, Surv(outcome_FU, outcome))[1] %>% `names<-`(NULL)), NA,
                                         1 - rcorr.cens(pred_score, Surv(outcome_FU, outcome))[1] %>% `names<-`(NULL))),
              wcind_model = ifelse(length(unique(outcome)) == 1 | all(is.na(weight)), NA, intsurv::cIndex(outcome_FU, event = outcome, pred_score, weight = weight)[[1]]),
              .groups = "keep") %>%
    ungroup()
  
  # Add difference in predictions cases-controls per pair
  if (all(is.na(median_dif_ds))){
    perf_ds$median_pair_diff <- as.numeric(NA)
  } else {
    perf_ds <- perf_ds %>%
      full_join(median_dif_ds, by = c("variable", "group"))
  }
  # Performance for unmatched samples
  }else{
    perf_ds <- pred_ds %>%
      dplyr::group_by(variable, group) %>%
      dplyr::summarise(
        auc_model = ifelse(length(unique(outcome)) == 1, NA,
                           as.numeric(pROC::auc(
                             pROC::roc(outcome, pred_score,
                                       quiet = TRUE,
                                       direction = "<")
                           ))),
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
                weight = w[valid]) ) } },
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
    
  }
  # Gather metrics
  perf_ds <- perf_ds %>%
    tidyr::gather(metric, value, -variable, -group) %>%
    dplyr::mutate(variable = factor(variable, levels = strat_vars_levels))
  # Return results
  perf_ds
}
