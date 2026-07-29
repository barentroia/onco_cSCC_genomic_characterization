# Auxilliary functions for model evaluation:
# ------------------------------------------------------------------------------
# Bootstrap_validation_inner function: This function performs bootstrap validation of a model fitted on data ds_x. Only inner loop, no confidence intervals are provided.
# @param feat_sel_method: character indicating which type of model should be applied to the data
# @param topfeats: number to indicate the amount of top features in feature selection, default NA
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
# @param genes_oi: genes of interest to keep in feature selection, on default not used
# removed during the training, if TRUE, separate performances will be computed for these samples
bootstrap_validation_inner <- function(feat_sel_method,topfeats=NA,model_method,ds_x,ds_additional,ds_x_counts,outcome_name,fup_name,weights_name,tp,B,seed,verbose=TRUE,model.params=NULL,bad_samples=NA,remove_bad_samples=F,genes_oi=NULL){

  # Set seed to ensure reproducibility
  set.seed(seed)
  #-------------------------------------------------------------------------------
  
  # 1. Model fitting and performance in the entire dataset to calculate apparent
  #-------------------------------------------------------------------------------
  # 1.0. Remove samples with bad QCs.
  # Subset clinical dataset
  ds_y_train <- ds_additional[!(rownames(ds_additional) %in% bad_samples),]
  # Subset gene expression matrices
  ds_x_train <- ds_x[!(rownames(ds_x) %in% bad_samples),]
  ds_x_counts_train <- ds_x_counts[, !(colnames(ds_x_counts) %in% bad_samples)]
  # Create vector where to save results
  predicted_probs_badQCsamples <- performance_badQCsamples <- vector("list", B)
  # 1.1. Build a model in the entire dataset
  fit_train <- build_model(feat_selection_method = feat_sel_method,
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

  # 1.2. Get model predictions and performance in the training cohort
  predictions_model_orig <- fit_train$predictions
  perf_orig_on_orig <- get_performance_measures(predictions_model_orig, ds_y_train, outcome_name, fup_name, weights_name, tp, fit_train$recalibration_cox_lp,model.params)
  ds_x_mod <- ds_x
  ds_x_train_mod <- ds_x_train
  feats_to_keep <- fit_train$preselected_feats

  #-------------------------------------------------------------------------------
  
  # 2. Model fitting and performance on bootstrapping
  #-------------------------------------------------------------------------------
  # 2.0. Save intermediate steps:
  optimism_bootstraps <- apparent_bootstraps <- oob_bootstraps <- predicted_probs <- selected_feats <- bootstrap_coeffs <- lambda_list <- alpha_list <- vector("list", B)

  # The following seed will be used as input if the bootstrap sample does not
  # have enough cases or controls
  s <- 1
  
  # Bootstrap repetitions:
  # ----------------------------------------------------------------------------
  for (brep in c(1:B)){
    if(verbose){
      print(paste0("Calculating bootstrap repetition: ", brep))
    }

    set.seed(brep) # set seed

    # 2.1. Draw a bootstrap sample with replacement from the sample:
    # 2.1.1 sampling for matches samples
    if (is.null( model.params$unmatched)){
    # Sample set ids:
    set_ids_sampled <- sample(unique(as.numeric(as.character(ds_y_train[,"Set_id"]))),replace=T)
    # Find indexes of sampled pairs:
    j <- unlist(purrr::map(set_ids_sampled, function(x) which(ds_y_train[,"Set_id"]==x)))

    # Check that there are enough cases and controls in the bootstrap sample
    # and at least 1 case and 1 control in the OOB sample
    s <- s + brep
    while(any(table(ds_y_train[j,][[outcome_name]]) < 2) |
          length(unique(ds_y_train[j,][[outcome_name]])) == 1 |
          length(unique(ds_y_train[-j,][[outcome_name]])) == 1){
      set.seed(s)
      # Sample set ids:
      set_ids_sampled <- sample(unique(as.numeric(as.character(ds_y_train[,"Set_id"]))),replace=T)
      # Find indexes of sampled pairs:
      j <- unlist(purrr::map(set_ids_sampled, function(x) which(ds_y_train[,"Set_id"]==x)))
      s <- s + 1
    }

    # Extract bootstrapping datasets
    ds_y_train_bootstrap <- ds_y_train[j,]
    ds_y_train_bootstrap$Set_id <- as.numeric(as.character(ds_y_train_bootstrap$Set_id))
    ds_x_train_bootstrap <- as.data.frame(ds_x_train[j,])
    ds_x_counts_train_bootstrap <- ds_x_counts_train[,j]

    # Because some set IDs are selected multiple times, they need to be renamed
    # otherwise it is not possible to compute difference in predictions
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
    # 2.1.2 unmatched sampling
    }else{
      n <- nrow(ds_y_train)
      j <- sample.int(n, size = n, replace = TRUE)
      
      # Ensure at least 2 classes in bootstrap sample + OOB sanity
      s <- s + brep
      while (
        length(unique(ds_y_train[j, ][[outcome_name]])) == 1 |
        any(table(ds_y_train[j, ][[outcome_name]]) < 2) |
        length(unique(ds_y_train[-unique(j), ][[outcome_name]])) == 1
      ) {
        set.seed(s)
        j <- sample.int(n, size = n, replace = TRUE)
        s <- s + 1
      }
      # Build bootstraps samples
      ds_y_train_bootstrap <- ds_y_train[j, , drop = FALSE]
      ds_x_train_bootstrap <- as.data.frame(ds_x_train[j, , drop = FALSE])
      ds_x_counts_train_bootstrap <- ds_x_counts_train[, j, drop = FALSE]

      # NO pairing assumption anymore → assign independent IDs
      ds_y_train_bootstrap$Set_id <- seq_len(nrow(ds_y_train_bootstrap))
      # index for downstream alignment (keep your logic safe)
      ds_y_train_bootstrap$index <- seq_len(nrow(ds_y_train_bootstrap))
    }
    # 2.2. Build a model in the bootstrap samples:
    fit_bootstrap <- build_model(feat_selection_method = feat_sel_method,
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
                                 seed = seed)

    # maintain consistency of seed
    set.seed(brep)

    # Bootstrap results
    bootstrap_coeffs[[brep]] <- fit_bootstrap[["coefficients"]]
    selected_feats[[brep]] <- data.frame(feats = unique(names(fit_bootstrap[["coefficients"]])),
                                         brep = 1) %>%
      `colnames<-`(c("feats", brep))

    # 2.3. Apply model in the bootstrap sample and estimate performance:
    predictions_model <- fit_bootstrap$predictions
    perf_bsmodel_on_bs <- get_performance_measures(predictions_model, ds_y_train_bootstrap, outcome_name, fup_name, weights_name, tp, fit_bootstrap$recalibration_cox_lp,model.params)
    feats_to_keep_bootstrap <- fit_bootstrap$preselected_feats

    # 2.4. Save lambda and alpha
    if (!grepl("cox|glmnet|lasso|ridge", model_method)) {
      lambda_list[[brep]] <- NULL
      alpha_list[[brep]] <- NULL
    } else {
      lambda_list[[brep]] <- fit_bootstrap[["fitted_model"]]$lambda.1se
      alpha_list[[brep]] <- fit_bootstrap[["alpha"]]
    }
   
    # 2.5. Apply bootstrap model in the entire, OOB datasets and bad QCs samples, and estimate performances:
    if (feat_sel_method != "none") {
    res_bs_model_on_orig <- apply_model(model_method,fit_bootstrap,ds_x_train_mod[,feats_to_keep_bootstrap],ds_y_train,tp,fit_bootstrap, model.params)
    res_bs_model_on_oob <- apply_model(model_method,fit_bootstrap,ds_x_train_mod[-j,feats_to_keep_bootstrap],ds_y_train[-j,],tp,fit_bootstrap,model.params)
    }else{
      res_bs_model_on_orig <- apply_model(model_method,fit_bootstrap,ds_x_train_mod,ds_y_train,tp,fit_bootstrap, model.params)
      res_bs_model_on_oob <- apply_model(model_method,fit_bootstrap,ds_x_train_mod[-j,],ds_y_train[-j,],tp,fit_bootstrap,model.params)
    }
    predictions_bs_model_on_orig <- res_bs_model_on_orig$predictions
    predictions_bs_model_on_oob <- res_bs_model_on_oob$predictions
    ### Bootstrap model on entire dataset
    perf_bsmodel_on_orig <- get_performance_measures(predictions_bs_model_on_orig, ds_y_train, outcome_name, fup_name, weights_name, tp, res_bs_model_on_orig$recalibration_cox_lp,model.params)
    ### Bootstrap model on out of bag
    perf_bsmodel_on_oob <- get_performance_measures(predictions_bs_model_on_oob, ds_y_train[-j,], outcome_name, fup_name, weights_name, tp, res_bs_model_on_oob$recalibration_cox_lp,model.params)

    # 2.6. Save probabilities
    predicted_probs[[brep]] <- predictions_bs_model_on_oob %>%
          as.data.frame() %>%
          `colnames<-`("predicted_probs") %>%
          mutate(boots_sample = brep)
    if ("SkylineDx.ID" %in% colnames(ds_y_train)){
        predicted_probs[[brep]] <- predicted_probs[[brep]] %>%
          mutate(SkylineDx.ID = ds_y_train[-j,]$SkylineDx.ID) %>%
          relocate(SkylineDx.ID, predicted_probs, boots_sample)
    }

    # 2.7. Compute optimism-corrected metrics
    optimism <- full_join(perf_bsmodel_on_bs, perf_bsmodel_on_orig,
                         by = c("variable", "group", "metric")) %>%
      mutate(value = value.x - value.y) %>%
      dplyr::select(-value.x, -value.y)
    optimism_bootstraps[[brep]] <- optimism %>% mutate(boots_sample = brep)
    apparent_bootstraps[[brep]] <- perf_bsmodel_on_bs %>% mutate(boots_sample = brep)
    oob_bootstraps[[brep]] <- perf_bsmodel_on_oob %>% mutate(boots_sample = brep)
  }
  #-------------------------------------------------------------------------------
  
  # 3. Make regularized cox with average lambda and alpha for signature
  #-------------------------------------------------------------------------------
  if (!grepl("cox", model_method)) {
    lambda_list[[brep]] <- NULL
    alpha_list[[brep]] <- NULL
    fit_train_opt <- NULL
  } else {

    set.seed(seed)
    mean_lambda <- mean(unlist(lambda_list))
    mean_alpha <- mean(unlist(alpha_list))
    model.params$avg_lambda1se <- mean_lambda
    model.params$avg_alpha <- mean_alpha
    # 3.1. Build a model in the entire dataset with optimal lambda and alpha
    fit_train_opt <- build_model(feat_selection_method = feat_sel_method,
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
  }

  #-------------------------------------------------------------------------------
  
  # 4. Combine results and get corrected performance
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

  # 4.1. Corrected metrics
  ## Harrel's bias correction method
  int_val_estimates <- full_join(apparent_performance, optimism_estimates,
                                by = c("variable", "group", "metric")) %>%
    mutate(value = value.x - value.y) %>%
    dplyr::select(-value.x, -value.y)
  ## .632 bias correction method
  int_val_estimates0632 <- full_join(apparent_performance, oob_estimates %>%
                                      filter(!grepl('pval', metric)),
                                    by = c("variable", "group", "metric")) %>%
    mutate(value = 0.368*value.x + 0.632*value.y) %>%
    dplyr::select(-value.x, -value.y)
  ## 4.2 .632+ bias correction method
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
    filter(!grepl('pval', metric)) %>%
    full_join(apparent_performance,
              by = c("variable", "group", "metric")) %>%
    full_join(random_vals,
              by = c("variable", "group", "metric")) %>%
    mutate(R = (value.x - value.y)/(random_value - value.y),
           R = ifelse(R > 1, 1, ifelse(R < 0, 0, R)), # R should range between 0 and 1
           weight = 0.632 / (1 - 0.368 * R),
           value = (1 - weight) * value.y + weight * value.x) %>%
    dplyr::select(-value.x, -value.y, -R, -weight, -random_value)

  # 4.3. Selected features
  selected_feats <- selected_feats %>%
    purrr::reduce(full_join, by = "feats") %>%
    replace(is.na(.), 0) %>%
    column_to_rownames(var = "feats")

  # 4.5. Results
  # Save alpha and lambda
  parametersperboot <- data.frame(lambda = lambda_list, alpha = alpha_list)
  res <- list("Optimism-Corrected" = int_val_estimates,
              "Optimism-Corrected0632" = int_val_estimates0632,
              "Optimism-Corrected0632plus" = int_val_estimates0632plus,
              "Optimism" = optimism_estimates,
              "All bootstrap optimism" = optimism_bootstraps,
              "All apparent bootstrap"= apparent_bootstraps,
              "All out of bag"= oob_bootstraps,
              "Model_fit" = fit_train,
              "Final_model2" = fit_train_opt,
              "Selected_Features"= selected_feats,
              "Bootstrap_coefficients" = bootstrap_coeffs,
              "Bootstrap_alphalambda" = parametersperboot,
              "Apparent_performance"= perf_orig_on_orig
  )
  
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
bootstrap_validation <- function(feat_sel_method,topfeats,model_method,ds_x,ds_additional,ds_x_counts,outcome_name,fup_name,weights_name,tp,B_inner,seed,cis=FALSE,B_outer=NULL,model.params=NULL,bad_samples=NA,remove_bad_samples=F,genes_oi=NULL){
  argg <- as.list(environment())
  argg <- argg[setdiff(names(argg),c("ds_x", "ds_additional", "ds_x_counts"))]
  # Set seed to ensure reproducibility
  set.seed(seed)

    input_arg_check(ds_x=ds_x,ds_y=ds_additional,ds_x_counts=ds_x_counts,
                      outcome_name,fup_name,weights_name)

    bootstrap_results = bootstrap_validation_inner(feat_sel_method = feat_sel_method,
                                                   topfeats = topfeats,
                                                   model_method = model_method,
                                                   ds_x = ds_x,
                                                   ds_additional = ds_additional,
                                                   ds_x_counts = ds_x_counts,
                                                   outcome_name = outcome_name,
                                                   fup_name = fup_name,
                                                   weights_name = weights_name,
                                                   tp = tp,
                                                   B = B_inner,
                                                   seed = seed,
                                                   model.params = model.params,
                                                   bad_samples = bad_samples,
                                                   remove_bad_samples = remove_bad_samples)
  return(c(bootstrap_results, "Bootstrap_parameters" = argg))
}


# input_arg_check function: This function checks the input argument of the user
# @param ds_x: dataset with all of the features, columns should be features with colnames identifying the features
# @param ds_additional: dataset with outcomes of interest and matching informations
# @param ds_x_counts: counts dataset with all of the features, columns should be features with colnames identifying the features
# @param outcome_name: character vector corresponding to the column name of the outcome of interest in the ds_additional dataframe
# @param fup_name: character vector corresponding to the column name of the follow-up time of the outcome of interest in the ds_additional dataframe
# @param weights_name: column name of the weights in ds_y
input_arg_check <- function(ds_x,ds_y,ds_x_counts,outcome_name,fup_name,weights_name){
  if(!(outcome_name %in% colnames(ds_y))){
    stop(paste0("Outcome_name: ",outcome_name," not in clinical data (ds_additional/ds_y)"))
  }
  if(!(fup_name %in% colnames(ds_y))){
    stop(paste0("Fup_name: ",fup_name," not in clinical data (ds_additional/ds_y)"))
  }
  if(!(weights_name %in% colnames(ds_y))){
    stop(paste0("Weights_name: ",weights_name," not in clinical data (ds_additional/ds_y)"))
  }
  if (!is.null(ds_x_counts)){
    if(!(all(colnames(ds_x_counts) == ds_y$SkylineDx.ID)) & identical(colnames(ds_x_counts), ds_y$SkylineDx.ID)){
      stop("Missing samples | not the same order of samples in counts data (samples=columns)")
    }
    if(!(all(colnames(ds_x) == rownames(ds_x_counts)))){
      stop("Features not the same between ds_x and count data")
    }
  }
  if(!(all(rownames(ds_x) == ds_y$SkylineDx.ID)) & identical(rownames(ds_x), ds_y$SkylineDx.ID)){
    stop("Missing samples | not the same order of samples in ds_x (samples=rows)")
  }

  print("Input check DONE.")
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
  }
  }else{
  # Build dataset 
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
                     group = "All"))}
    else{
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
      mutate(variable = "All",
             group = "All")
    strat_vars_levels <- c("All")
  }

  # Compute difference in predictions cases-controls per pair
  if (is.null( model.params$unmatched)){
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

  if (is.null( model.params$unmatched)){
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
  }else{
    # Metrics
    perf_ds <- pred_ds %>%
      dplyr::group_by(variable, group) %>%
      dplyr::summarise(
        auc_model = ifelse(length(unique(outcome)) == 1, NA,as.numeric(pROC::auc(pROC::roc(outcome, pred_score, quiet = TRUE, direction = "<")))),
        wauc_model = {
          w <- weight
          valid <- !is.na(w) & w > 0
          if (length(unique(outcome)) == 1 || sum(valid) < 2 || length(unique(outcome[valid])) < 2) {
            NA
          } else {
            WeightedAUC(WeightedROC(pred_score[valid],outcome[valid], weight = w[valid] )  )} },
        cind_model = ifelse(length(unique(outcome)) == 1, NA, 1 - rcorr.cens(pred_score, survival::Surv(outcome_FU, outcome))[1]),
        wcind_model = ifelse(length(unique(outcome)) == 1 | all(is.na(weight)),
                             NA, intsurv::cIndex(outcome_FU,event = outcome,pred_score,weight = weight)[[1]]),
        .groups = "keep"
      ) %>%
      dplyr::ungroup()
    
  }
  if (is.null( model.params$unmatched)){
  # Add difference in predictions cases-controls per pair
  if (all(is.na(median_dif_ds))){
    perf_ds$median_pair_diff <- as.numeric(NA)
  } else {
    perf_ds <- perf_ds %>%
      full_join(median_dif_ds, by = c("variable", "group"))
  }
  }
  
  # Gather metrics
  perf_ds <- perf_ds %>%
    gather(metric, value, -variable, -group) %>%
    mutate(variable = factor(variable, levels = strat_vars_levels))
  # Return results
  perf_ds
}
