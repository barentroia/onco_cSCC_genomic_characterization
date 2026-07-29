# apply_model function: This function applies the model provided in fitted model to the dataset ds_new
# @param model_method: character indicating which type of model should be applied to the data
# @param model_res: list with model results
# @param ds_new: dataset where the model should be applied and which contains all of the features, columns should be features with colnames identifying the features
# @param ds_new_y: dataset with outcomes of interest and matching informations
# @param tp: number indicating the timepoint of interest
# @param extra_model_info: additional information for modelling
# @param model_params: vector of model parameters
apply_model <- function (model_method, model_res, ds_new, ds_new_y, tp, extra_model_info, model.params=NULL){
  # Model
  fitted_model <- model_res$fitted_model

 # Apply model and get preditions
  if (grepl("coxnet", model_method)){
    # If late feature integration or late prediction integration is chosen
    if (!is.null(model.params$combi_type)) {
      # If combination type is to combine predictions but this model is not the combined prediction model
      if (model.params$combi_type == "Stacked" & length(colnames(model_res$x_fitted_model)) > 2) {
        mat_ds_new_all <- as.matrix(ds_new)
        mat_ds_new_all <- mat_ds_new_all[, colnames(model_res$x_fitted_model)]
        # If combination type is to combine predictions and this model is the combined prediction model
      } else if (model.params$combi_type == "Stacked" & length(colnames(model_res$x_fitted_model)) == 2) {
        mat_ds_new_all <- as.matrix(ds_new)
      } else {
        # If combination type is to combine features
        mat_ds_new_all <- as.matrix(ds_new)
        mat_ds_new_all <- mat_ds_new_all[, colnames(model_res$x_fitted_model)]
      }
    } else {
      # If DvP, non-DvP or early integration is chosen
      mat_ds_new_all <- as.matrix(ds_new)
      mat_ds_new_all <- mat_ds_new_all[, colnames(model_res$x_fitted_model)]
    }
    model_predictions <- as.vector(1 - summary(survfit(fitted_model, s = "lambda.1se", x = model_res$x_fitted_model, y = model_res$y_fitted_model, weights = model_res$w_fitted_model, newx = mat_ds_new_all), times = tp)$surv)
    names(model_predictions) <- rownames(ds_new)
    
    
  } else if (grepl("RSF", model_method)){
    # Apply Random survival forest
    # Make sure it works on prediction input data
    if (!is.null(model.params$combi_type) && model.params$combi_type == "Stacked") {
      train_vars <- fitted_model$forest$xvar.names
      new_vars   <- colnames(ds_new)
      if (!identical(new_vars, train_vars)) {
        # same variables, wrong order
        if (setequal(new_vars, train_vars)) {
          ds_new <- ds_new[, train_vars, drop = FALSE]
          # unnamed / same length case
        } else if (ncol(ds_new) == length(train_vars)) {
          colnames(ds_new) <- train_vars
          ds_new <- ds_new[, train_vars, drop = FALSE]
        }
      }
      model_predictions <- 1 - as.vector(
        pec::predictSurvProb(
          fitted_model,
          newdata = as.data.frame(ds_new),
          times = tp))
    } else {
      model_predictions <- 1 - as.vector(
        pec::predictSurvProb(
          fitted_model,
          newdata = as.data.frame(ds_new),
          times = tp
        )) }

    names(model_predictions) <- rownames(ds_new)
  } 

  # Save results
  res <- list(predictions = model_predictions)
  return(res)
}
