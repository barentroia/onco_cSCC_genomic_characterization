#===============================================================================
# Aim:      Feature selection and fitting of survival prediction models: Cox,
#           regularized Cox (elastic-net) and Random Survival Forest.
# Contents: - build_model(): feature selection + model fitting (entry point)
#           - build_cv.glmnet(): regularized Cox with alpha/lambda tuning
#           - build_rsf(): Random Survival Forest with hyperparameter tuning
#           - predefine_folds(): CV folds that keep case-control sets together
# Requires: dir_scripts_functions; select_features() (features_selection.r)
# Usage:    Sourced by Model_input_processing.R
# Notes:    build_cox(), used when model_method contains "cox" but not
#           "coxnet", is not defined in this repository.
#===============================================================================

source(file.path(dir_scripts_functions, "features_selection.r"))


#' Select features and build a survival model
#'
#' Runs feature selection (select_features()) and fits the requested model:
#' Cox regression, regularized Cox (Coxnet, elastic-net with alpha tuning),
#' the combined "model2" Coxnet used in late integration, or a Random
#' Survival Forest. Models are fitted unweighted (w = NULL).
#'
#' @param feat_selection_method character. Feature selection method
#'   ("deseq2" or "none").
#' @param topfeats integer. Number of features to select.
#' @param model_method character. Model type; matched with grepl() against
#'   "cox", "coxnet", "model2" and "RSF".
#' @param ds_x data.frame. Features (samples x features).
#' @param ds_y data.frame. Outcomes and matching information, rows aligned
#'   with ds_x.
#' @param ds_x_counts data.frame or NULL. Raw counts (genes x samples) for
#'   DESeq2 feature selection; transposed automatically if needed.
#' @param outcome_name character. Outcome column in ds_y.
#' @param fup_name character. Follow-up time column in ds_y.
#' @param weights_name character. Weights column in ds_y (stored only).
#' @param pairs_name character. Case-control set column in ds_y.
#' @param sampleid_name character. Sample ID column in ds_y.
#' @param tp numeric. Time point (years) for the predicted probabilities.
#' @param model.params list or NULL. Model parameters, e.g. mut_alphas,
#'   avg_alpha, avg_lambda1se (final model with averaged hyperparameters),
#'   alphas2/avg_alpha2/avg_lambda1se2 (model2) and unmatched.
#' @param genes_oi character vector or NULL. Genes always kept.
#' @param seed integer. Seed for reproducibility.
#'
#' @return list from the model builder (fitted_model, coefficients,
#'   predictions, ...) plus preselected_feats and parameters.
build_model <- function(feat_selection_method, topfeats, model_method, ds_x, ds_y, ds_x_counts = NULL,
                       outcome_name, fup_name, weights_name = "Weight_rescaled", pairs_name = "Set_id",
                       sampleid_name = "SkylineDx.ID", tp, model.params = NULL, genes_oi = NULL,seed){
  
  #set.seed(seed)
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


#' Build a regularized Cox model with (cv.)glmnet
#'
#' Tunes alpha (if several are given) by the apparent C-index at lambda.1se,
#' then tunes lambda by 10-fold cross-validation (case-control sets kept in
#' the same fold in the matched design), unless a fixed lambda is given. If
#' the model selects no features, the fit is repeated with another seed (max
#' 20 times).
#'
#' @param x data.frame. Preselected features (samples x features);
#'   categorical variables are dummy-encoded.
#' @param y data.frame. Outcomes and matching information.
#' @param outcome_n character. Outcome column in y.
#' @param fup_n character. Follow-up time column in y.
#' @param tp numeric. Time point (years) for the predicted probabilities.
#' @param pairs_n character. Case-control set column in y.
#' @param family character. glmnet family (here "cox").
#' @param alphas numeric vector. Candidate alpha values.
#' @param weights numeric vector or NULL. Observation weights.
#' @param lambda numeric or NULL. Fixed lambda (e.g. average lambda.1se over
#'   bootstraps); if NULL, lambda is tuned by cross-validation.
#' @param seed integer or NULL. Seed for reproducibility.
#' @param model.params list. Model parameters; unmatched disables set-based folds.
#'
#' @return list with fitted_model, coefficients (non-zero), predictions
#'   (event probability at tp), warning_method, alpha, and the x/y/weights
#'   used for fitting (x_fitted_model, y_fitted_model, w_fitted_model).
build_cv.glmnet <- function(x, y, outcome_n, fup_n, tp, pairs_n, family, alphas, weights,
                            lambda = NULL, seed = NULL, model.params){
  # 0. Initialize warning variable, outcome object and seed switch
  warning_glmnet <- NULL
  out_obj <- Surv(y[, fup_n], y[, outcome_n])
  type_measure <- "C"
  use_seed <- !is.null(model.params$unmatched)
  
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
  
  # Predefine folds so that pairing is kept (matched design only)
  if (is.null(model.params$unmatched)){
    pred_folds <- predefine_folds(y, outcome_n, pairs_n, 0, seed)
    stopifnot(rownames(pred_folds) == rownames(y))
    foldid <- pred_folds$fold_n
  } else {
    foldid <- NULL
  }
  penalty_factor <- rep(1, ncol(mat_x))
  
  ## 2. Optimize alpha for elastic-net or use the average alpha if given
  if (length(alphas) > 1){
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
  
  # 3. Optimize lambda with cross-validation if `lambda` (average lambda) is not provided
  # NOTE: foldid = y$fold_n is NULL (y has no fold_n column), so cv.glmnet
  # draws random folds here; with use_seed = TRUE they are fixed by the seed
  if (use_seed) set.seed(seed)
  if (is.null(lambda)) {
    if (opt_alpha == 0){
      # Deviance avoids convergence issues of the C-index when alpha = 0
      type_measure <- "deviance"
      model_fit <- cv.glmnet(mat_x, out_obj, alpha = 0, family = family, type.measure = type_measure,
                             weights = weights, foldid = y$fold_n, penalty.factor = penalty_factor)
    } else {
      model_fit <- cv.glmnet(mat_x, out_obj, alpha = opt_alpha, family = family, type.measure = type_measure,
                             weights = weights, foldid = y$fold_n, penalty.factor = penalty_factor)
    }
  } else {
    print(paste0("The optimal lambda:", lambda))
    if (opt_alpha == 0){
      type_measure <- "deviance"
      model_fit <- glmnet(mat_x, out_obj, alpha = 0, family = family, type.measure = type_measure, lambda = lambda,
                          weights = weights, foldid = y$fold_n, penalty.factor = penalty_factor)
    } else {
      model_fit <- glmnet(mat_x, out_obj, alpha = opt_alpha, family = family,
                          lambda = lambda, type.measure = type_measure,
                          weights = weights, foldid = y$fold_n, penalty.factor = penalty_factor)
    }
  }
  if (use_seed) set.seed(seed)
  
  # 4. Check convergence: if nothing is selected, repeat with another seed
  ns = 1
  if (opt_alpha != 0) {
    while(length(unique(predict(model_fit, s = "lambda.1se", newx = mat_x, type = "response"))) == 1){
      set.seed(1 + ns)
      print(paste0("Repeating the procedure with another seed:", 1 + ns))
      # Change predefined folds so that pairing is kept
      if (is.null(model.params$unmatched)){
        pred_folds <- predefine_folds(y, outcome_n, pairs_n, 0, seed)
        stopifnot(rownames(pred_folds) == rownames(y))
        foldid <- pred_folds$fold_n
      } else {
        foldid <- NULL
      }
      model_fit <- cv.glmnet(mat_x, out_obj, alpha = opt_alpha, family = family, type.measure = type_measure,
                             weights = weights, foldid = foldid)
      if (use_seed) set.seed(seed)
      ns <- ns + 1
      if (ns == 20){
        print("WARNING: Did not converge!")
        warning_glmnet <- "no_converge"
        break()
      }
    }
  }
  
  # 5. Extract predictions and coefficients
  model_predictions <- as.vector(1 - summary(survfit(model_fit, s = "lambda.1se", x = mat_x, y = out_obj,
                                                     weights = weights, newx = mat_x), times = tp)$surv)
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


#' Build a Random Survival Forest
#'
#' Tunes nodesize and mtry with randomForestSRC::tune.rfsrc() and fits the
#' final forest (maximum node depth 6) with permutation importance.
#'
#' @param x data.frame. Preselected features (samples x features).
#' @param y data.frame. Outcomes, rows aligned with x.
#' @param outcome_n character. Outcome column in y.
#' @param fup_n character. Follow-up time column in y.
#' @param tp numeric. Time point (years) for the predicted probabilities.
#' @param weights numeric vector or NULL. Case weights.
#'
#' @return list with fitted_model, coefficients (non-zero variable
#'   importances), predictions (event probability at tp) and best_params.
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

#' Predefine cross-validation folds that keep case-control sets together
#'
#' Controls are split into 10 folds and each case is assigned to the fold of
#' its matched control(s).
#'
#' @param ds data.frame. Outcomes and matching information.
#' @param outcome_n character. Outcome column in ds.
#' @param pairs_n character. Case-control set column in ds.
#' @param ctrl_value value of the outcome that indicates controls (0).
#' @param seed integer. Not used (set.seed is commented out).
#'
#' @return data.frame with one column, fold_n, rows aligned with ds.
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
