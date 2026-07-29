source(file.path(dir_scripts_functions, "features_selection.r"))

# build_model function: This function builds the model specified by "model_method"
# @param feat_selection_method: character indicating the feature selection method to use
# @param topfeats: integer indicating the number of features to select
# @param model_method: character indicating which type of model should be applied to the data
# @param ds_x: dataset with all of the features, columns should be features with colnames identifying the features
# @param ds_y: dataset with outcomes of interest and matching information
# @param ds_x_counts: counts dataset with all of the features, columns should be features with colnames identifying the features
# @param outcome_name: character corresponding to the column name of the outcome of interest in the ds_y dataframe
# @param fup_name: character corresponding to the column name of the follow-up time of the outcome of interest in the ds_y dataframe
# @param pairs_name: character corresponding to the column name of the pairs in the ds_y dataframe
# @param sampleid_name: character corresponding to the column name of the sample IDs in the ds_y dataframe
# @param weights_name: character corresponding to the column name of the weights in the ds_y dataframe
# @param tp: number indicating the timepoint of interest
# @param model_params: vector of model parameters
# @param genes_oi: list of genes of interest
build_model <- function(feat_selection_method, topfeats, model_method, ds_x, ds_y, ds_x_counts = NULL,
                       outcome_name, fup_name, weights_name = "Weight_rescaled", pairs_name = "Set_id",
                       sampleid_name = "SkylineDx.ID", tp, model.params = NULL, genes_oi = NULL,seed){
  
  set.seed(seed)
    # Check that datasets are aligned
    stopifnot(all(rownames(ds_y) == rownames(ds_x)))
    # Ensure ds_x_counts is genes x samples
    if (!all(colnames(ds_x_counts) == rownames(ds_x))) {
      print("ds_counts has incorrect dim")
      # Try transpose
      if (all(rownames(ds_x_counts) == rownames(ds_x))) {
        ds_x_counts <- t(ds_x_counts)
      } else {
        stop("ds_x_counts cannot be aligned to ds_x")
      }
    }

    # Gene expression dataset before features selection
    ds_orig <- ds_x

    # 1. Features selection
    ordered_feats <- select_features(ds_x, ds_y, feat_selection_method, topfeats, 
                                     outcome_name, fup_name, weights_name, 
                                     pairs_name, sampleid_name, ds_x_counts, 
                                     genes_oi, model.params)

    # Subset ds_x only if ordered_feats is not NULL
    if (!is.null(ordered_feats)){
      ds_x <- ds_x[, ordered_feats]
    }

    # 2. Model building
    # 2.1 Make sure unweighted models are used
    w <- NULL
    # 2.2 Cox regression
    if (grepl("(?!.*coxnet)(?=.*cox)", model_method, perl = T)){
        output_model <- build_cox(ds_x, ds_y, outcome_name, fup_name, tp, w)}
    # 2.3 Cox and GLM elastic-net/lasso/ridge
    else if (grepl("coxnet", model_method)){
     # Define link function and alphas
      # Cox models
      model_family <- "cox"
      alphas <- c(0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)
      
      # For WES Models, run on alpha = 0.5
      if (!is.null( model.params$mut_alphas)) {
        alphas <- model.params$mut_alphas
      }
      
      # Check if mean_alpha exist in model.params
      if (!is.null( model.params$avg_alpha)) {
        alphas <-  model.params$avg_alpha  # Set alpha to mean_alpha
      }
      output_model <- build_cv.glmnet(ds_x, ds_y, outcome_name, fup_name, tp, pairs_name, model_family, alphas, w, model.params$avg_lambda1se,seed,model.params)
    }
    # 2.3 Combined cox model
    else if (grepl("model2", model_method)){
      model_family <- "cox"
       if (!is.null(model.params$avg_alpha2)){
         ##  average alpha from combined bootstrap models
        alphas <- model.params$avg_alpha2
      } else {
        ## alpha for combined bootstrap models, which is 0
        alphas <- model.params$alphas2
      }
      output_model <- build_cv.glmnet(ds_x, ds_y, outcome_name, fup_name, tp, pairs_name, model_family, alphas,w, model.params$avg_lambda1se2,seed,model.params)
    }
    # 2.4 Random survival forest
    else if (grepl("RSF", model_method)){
        output_model <- build_rsf(ds_x, ds_y, outcome_name, fup_name, tp, w)}
    
    #3 Save model
    # Add preselected features to model's output
    output_model$preselected_feats <- ordered_feats
    # Save parameters used in building model
    params <- list("weights_name" = weights_name,
                   "pairs_name" = pairs_name,
                   "sampleid_name" = sampleid_name,
                   "w" = w)
    output_model$parameters <- params
    output_model
}


# build_cv.glmnet function: This function builds models with the cv.glmnet function
# @param x: dataset with preselected features, columns should be features with colnames identifying the features
# @param y: dataset with outcomes of interest and matching information
# @param outcome_n: character  corresponding to the column name of the outcome of interest in the y dataframe
# @param fup_n: character  corresponding to the column name of the follow-up time of the outcome of interest in the y dataframe
# @param tp: number indicating the timepoint of interest
# @param pairs_n: character corresponding to the column name of the pairs in the y dataframe
# @param family: model family type corresponding to family parameter in cv.glmnet function
# @param alphas: numeric vector of alpha values
# @param weights: numeric vectors with weights or NULL
build_cv.glmnet <- function(x, y, outcome_n, fup_n, tp, pairs_n, family, alphas, weights, lambda = NULL,seed = NULL,model.params){
    # 0. Initialize warning variable and object variable
    warning_glmnet <- NULL
    out_obj <- Surv(y[, fup_n], y[, outcome_n])
    type_measure <- "C"

    # 1. Prepare data
    # cv.glmnet does not work with categorical data and these need to be
    # encoded into dummy variables. Check if there are categorical variables
    # and encode them into dummy ones
    if (any(sapply(colnames(x), function(a) is.character(x[,a])|is.factor(x[,a])))){
        # Drop categorical variables with a single value
        for (col in colnames(x)){
            if (length(unique(x[, col])) == 1){
                x[, col] <- NULL
            }
        }
        mat_x <- model.matrix(~ ., x)
        # Drop intercept
        mat_x <- mat_x[, setdiff(colnames(mat_x), "(Intercept)")]
    } else {
        mat_x <- as.matrix(x)
    }

    # Predefine folds so that pairing is kept
    if (is.null( model.params$unmatched)){
        pred_folds <- predefine_folds(y, outcome_n, pairs_n, 0,seed)
        stopifnot(rownames(pred_folds) == rownames(y))
        foldid <- pred_folds$fold_n
    } else {
        foldid <- NULL
    }
    penalty_factor <- rep(1, ncol(mat_x))
    # 2. Optimize alpha for elastic-net or use the average alpha if given
    if (length(alphas) > 1){
      # Elastic-net: optimize alpha
      
      # Initialize vector for c-index for alpha optimization
      cind <- vector(mode = "numeric", length(alphas))
      for (i in 1:length(alphas)){
        model_fit <- cv.glmnet(mat_x, out_obj, alpha = alphas[i], family = family, type.measure = type_measure,
                               weights = weights, foldid = foldid)
        model_predictions <- predict(model_fit, s = "lambda.1se", newx = mat_x, type = "response")
        cind[i] <- intsurv::cIndex(y[, fup_n], event = y[, outcome_n], model_predictions, weight = weights)[[1]]
      }
      # Pick alpha giving the highest c-index
      max_cind <- cind[order(cind, decreasing = TRUE)][1]
      opt_alpha <- alphas[which.max(cind)]
    } else {
        # Ridge/lasso: optimal alpha is the only one given as input
        opt_alpha <- alphas[[1]]
    }
    set.seed(seed)
    # 3. Optimize lambda with cross-validation if `lambda` (average lambda) is not provided
    # 3.1 Measure type C gives convergence issues if amount of genes is low and alpha = 0
    if (is.null(lambda)) {
      if (opt_alpha == 0){
        # use deviance if alpha =0 to not get convergence issues
        type_measure <- "deviance"
        model_fit <- cv.glmnet(mat_x, out_obj, alpha = 0, family = family, type.measure = type_measure,
                               weights = weights, foldid = y$fold_n, penalty.factor = penalty_factor)
        set.seed(seed)
      }else{
        model_fit <- cv.glmnet(mat_x, out_obj, alpha = opt_alpha, family = family, type.measure = type_measure,
                               weights = weights, foldid = y$fold_n, penalty.factor = penalty_factor)
        set.seed(seed)
      }
    }else{
      if (opt_alpha == 0){
        type_measure <- "deviance"
        model_fit <- glmnet(mat_x, out_obj, alpha = 0, family = family, type.measure = type_measure,lambda = lambda,
                               weights = weights, foldid = y$fold_n, penalty.factor = penalty_factor)
        set.seed(seed)
      }else{
      model_fit <- glmnet(mat_x, out_obj, alpha = opt_alpha, family = family,
                          lambda = lambda, type.measure = type_measure,
                          weights = weights, foldid = y$fold_n, penalty.factor = penalty_factor)
      set.seed(seed)
      } }

    # 4. Check convergence
    ns = 1
    # It can happen that for some reason, the glmnet falls into a local optimal where nothing is selected.
    # In that case, we repeat the procedure with another seed.
    if (opt_alpha != 0) {
    while(length(unique(predict(model_fit, s = "lambda.1se", newx = mat_x, type = "response"))) == 1){
        set.seed(1 + ns)
        print(paste0("Repeating the procedure with another seed:", 1 + ns))
        # Change predifined folds so that pairing is kept
        if (is.null(model.params$unmatched)){
            pred_folds <- predefine_folds(y, outcome_n, pairs_n, 0,seed)
            stopifnot(rownames(pred_folds) == rownames(y))
            foldid <- pred_folds$fold_n
        } else {
            foldid <- NULL
        }
        model_fit <- cv.glmnet(mat_x, out_obj, alpha = opt_alpha, family = family, type.measure = type_measure,
                               weights = weights, foldid = foldid)
        set.seed(seed)
        ns <- ns + 1
        if(ns == 20){
            print("WARNING: Did not converge!")
            warning_glmnet <- "no_converge"
            break()
        }
    }
}
  
    # 5. Extract predictions and coefficients 
    model_predictions <- as.vector(1 - summary(survfit(model_fit, s = "lambda.1se", x = mat_x, y = out_obj, weights = weights, newx = mat_x), times = tp)$surv)
    names(model_predictions) <- rownames(mat_x)
   
    df_fit <- as.matrix(coef(model_fit, s = "lambda.1se"))
    coeffs <- df_fit[df_fit != 0, ]
    names(coeffs) <- rownames(df_fit)[df_fit != 0]

    # 6. Output results
    res <- list(fitted_model = model_fit,
                coefficients = coeffs,
                predictions = model_predictions,
                warning_method = warning_glmnet,
                alpha = opt_alpha,
                x_fitted_model = mat_x,
                y_fitted_model = out_obj,
                w_fitted_model = weights)
    res
}


# build_rsf function: This function builds random survival forest models
# @param x: dataset with preselected features, columns should be features with colnames identifying the features
# @param y: dataset with outcomes of interest and matching information
# @param outcome_n: character  corresponding to the column name of the outcome of interest in the y dataframe
# @param fup_n: character corresponding to the column name of the follow-up time of the outcome of interest in the y dataframe
# @param tp: number indicating the timepoint of interest
# @param weights: numeric vectors with weights or NULL
build_rsf <- function(x, y, outcome_n, fup_n, tp, weights){
    # 0. Prepare settings and data to build model
    # Combine outcome and features data
    stopifnot(rownames(x) == rownames(y))
    comb_ds <- cbind(y[, c(outcome_n, fup_n)], x) %>% as.data.frame()

    # 1. Fit model
    # Model tuning
    tuned_model <- suppressMessages(tune.rfsrc(formula = as.formula(paste0("Surv(", fup_n, ",", outcome_n, ") ~ .")),
                              data = as.data.frame(comb_ds),
                              mtryStart = max(1, ceiling(sqrt(ncol(x)))),       # Starting point for mtry
                              nodesizeTry = seq(20, 30, by = 5),                 # Nodesize values to try
                              ntreeTry = 500,                                   # Initial number of trees
                              nsplit = 1,                                       # Default number of splits
                              stepFactor = 1.25,
                              improve = 1e-3,
                              strikeout = 3,
                              maxIter = 25,
                              doBest = TRUE,
                              nodedepth = 6,
                              trace = FALSE # to not print tuning
                              ))
    best_nodesize <- tuned_model$optimal[names(tuned_model$optimal) == "nodesize"]
    best_mtry <- tuned_model$optimal[names(tuned_model$optimal) == "mtry"]
    best_ntree <- tuned_model$rf$ntree  # Get number of trees from the fitted model
    # Fit best model
    model_fit <- rfsrc(formula = as.formula(paste0("Surv(", fup_n, ",", outcome_n, ") ~ .")),
                       data = comb_ds,
                       ntree = best_ntree,
                       mtry = best_mtry,
                       nodesize = best_nodesize,
                       importance = "permute",
                       case.wt = weights,
                       forest = TRUE,
                       nodedepth = 6)

    # 2. Extract model info and output
    model_predictions <- 1 - pec::predictSurvProb(model_fit, newdata = comb_ds, times = tp) %>% as.vector()
    names(model_predictions) <- rownames(comb_ds)
    # Variable importance
    imp_named_vector <- vimp(model_fit)$importance[which(vimp(model_fit)$importance != 0)]

    # 3. Output results
    list(fitted_model = model_fit,
         coefficients = imp_named_vector,
         predictions = model_predictions,
         best_params = list(ntree = best_ntree,
                            mtry = best_mtry,
                            nodesize = best_nodesize))

}

# predifine_folds function: This function predifines the folds to be used with the cv.glmnet function
# @param ds: dataset with outcomes of interest and matching information
# @param outcome_n: character  corresponding to the column name of the outcome of interest in the y dataframe
# @param pairs_n: character corresponding to the column name of the pairs in the y dataframe
# @param ctrl_value: value to indicate controls
predefine_folds <- function(ds, outcome_n, pairs_n, ctrl_value, seed){
  #set.seed(seed)
    ix_controls <- which(ds[[outcome_n]] == ctrl_value)
    folds <- createFolds(ix_controls, k = 10)

    ds$fold_n <- NA
    for (fold_ix in 1:length(folds)){
        pairs_fold <- ds[which(ds[[outcome_n]] == 0),][[pairs_n]][folds[[fold_ix]]]
        ds[ds[[pairs_n]] %in% pairs_fold,]$fold_n <- fold_ix
    }
    ds %>%
        dplyr::select(fold_n)
}
