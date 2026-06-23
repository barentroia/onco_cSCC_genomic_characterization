#---------------------------------------------------
# Aim: Util functions for model building
# Author: B. Rentroia Pacheco
# Input: Statistical tests for genes and pathways
# Output: Plots summarizing hazard ratios
#---------------------------------------------------

# Genetic algorithm:
# Fitness function: Picks samples that lead to highest c-index, with or without validation:
fitness_fun_cox <- function(chromosome, x, y,validation_split = FALSE) {
  
  # chromosome is a binary vector indicating selected features
  selected <- which(chromosome == 1)
  
  # limit max features to 10
  if(length(selected) == 0 || length(selected) > 10) return(0)  # penalize invalid solutions
  
  x_sub <- x[, selected, drop = FALSE]
  
  # fit Cox model
  if(validation_split){
    n = nrow(x)
    i.train <- sample(seq_len(n), size = floor(0.8 * n))
      
    smp_ids_train = unique(rownames(x)[i.train])
    smp_ids_test = setdiff(rownames(x),smp_ids_train)
  }else{
    smp_ids_train = rownames(x)
    smp_ids_test = rownames(x)
  }
  i.train.ids =  which(rownames(x)%in%smp_ids_train)
  i.test.ids =  which(rownames(x)%in%smp_ids_test)
  
  fit <- try(coxph(Surv(y[i.train.ids, "time"], y[i.train.ids, "status"]) ~ ., data = as.data.frame(x_sub)[smp_ids_train,]), silent = TRUE)
  
  if(inherits(fit, "try-error")) return(0)
  
  # compute C-index
  lp <- predict(fit, type = "lp",newdata=as.data.frame(x_sub)[i.test.ids,])
  c_idx <- 1-concordance(Surv(y[i.test.ids, "time"], y[i.test.ids, "status"]) ~ lp)$concordance
  
  return(c_idx)
}

# Fitness function: Picks features that lead to highest separation between cases and controls:
fitness_fun_oncoscore <- function(chromosome, x, y,validation_split = FALSE,penalty_feats = TRUE) {
  
  # chromosome is a binary vector indicating selected features
  selected <- which(chromosome == 1)
  
  # limit max features to 10
  if(length(selected) == 0 || length(selected) > 10) return(0)  # penalize invalid solutions
  
  x_sub <- x[, selected, drop = FALSE]
  
  # fit Cox model
  if(validation_split){
    n = nrow(x)
    i.train <- sample(seq_len(n), size = floor(0.66 * n))
    
    smp_ids_train = unique(rownames(x)[i.train])
    smp_ids_test = setdiff(rownames(x),smp_ids_train)
  }else{
    smp_ids_train = rownames(x)
    smp_ids_test = rownames(x)
  }
  i.train.ids =  which(rownames(x)%in%smp_ids_train)
  i.test.ids =  which(rownames(x)%in%smp_ids_test)
  
  #fit <- try(coxph(Surv(y[i.train.ids, "time"], y[i.train.ids, "status"]) ~ ., data = as.data.frame(x_sub)[smp_ids_train,]), silent = TRUE)
  
  #if(inherits(fit, "try-error")) return(0)
  
  sign_gene <-c()
  for (g in colnames(x_sub)) {
    cox_fit <- coxph(Surv(y[i.train.ids, "time"], y[i.train.ids, "status"]) ~  as.data.frame(x_sub)[i.train.ids, g,drop=TRUE])
    gene_coef = summary(cox_fit)$coefficients[,"coef"]
    if(!is.na(gene_coef)){
      sign_gene[g] = ifelse(gene_coef>0,1,-1)
    }else{
      sign_gene[g] = 0
    }
    
  }
  # print(sign_gene)
  # Compute C-index
  # sign_gene = as.numeric(ifelse(fit$coefficients>0,1,-1))
  # names(sign_gene) = names(fit$coefficients)
  if(length(sign_gene)>1){
    surv_prob = 1-as.numeric((as.matrix(x_sub)[i.test.ids,, drop = FALSE]) %*%sign_gene)
  }else{
    surv_prob =  1-as.numeric((as.matrix(x_sub)[i.test.ids,, drop = FALSE]) %*%sign_gene)
  }
  #lp <- predict(fit, type = "lp",newdata=as.data.frame(x_sub)[i.test.ids,])
  #print(sign_gene)
  #print(length(surv_prob))
  #print(length(i.test.ids))
  c_idx <- concordance(Surv(y[i.test.ids, "time"], y[i.test.ids, "status"]) ~ surv_prob)$concordance
  
 
  if(penalty_feats){
    # favor fewer features
    penalty_factor <- 0.01
    num_features <- length(selected)
    c_idx <- c_idx * (1 - penalty_factor * (num_features - 1))
  }
 
  return(c_idx)
}

build_model = function(x,y,model_method,survival_method = "survival",seed=123,params = list()){
  set.seed(seed)
  feat_sel = params[["fsel"]]
  if (!is.null(feat_sel)){
    if(feat_sel =="cox"){
      fsel_genes <- c()
      sign_gene <-c()
      for (g in colnames(x)) {
        cox_fit <- coxph(Surv(y[,"time"], y[,"status"]) ~ x[, g, drop=TRUE])
        
        p_val <- summary(cox_fit)$coefficients[,"Pr(>|z|)"]
        p_thresh = params[["p_thresh"]]
        if (!is.na(p_val) && p_val < p_thresh) {
          fsel_genes <- c(fsel_genes, g)
          sign_gene[g] = ifelse(summary(cox_fit)$coefficients[,"coef"]>0,1,-1)
        }else{
          sign_gene[g] =0
        }
      }
      x <- x[, fsel_genes, drop = FALSE]  # keep only significant features
    }else if(feat_sel=="GA"){
      #Pre-selection, if indicated:
      if(!is.null(params[["p_thresh"]])){
        fsel_genes <- c()
        for (g in colnames(x)) {
          cox_fit <- coxph(Surv(y[,"time"], y[,"status"]) ~ x[, g, drop=TRUE])
          
          p_val <- summary(cox_fit)$coefficients[,"Pr(>|z|)"]
          p_thresh = params[["p_thresh"]]
          if (!is.na(p_val) && p_val < p_thresh) {
            fsel_genes <- c(fsel_genes, g)
            #sign_gene[g] = ifelse(summary(cox_fit)$coefficients[,"coef"]>0,1,-1)
          }else{
            #sign_gene[g] =0
          }
        }
        x <- x[, fsel_genes, drop = FALSE]  # keep only significant features
      }
      
      n_features <- ncol(x)
      GA_res <- ga(
        type = "binary",
        fitness = function(chrom) fitness_fun_oncoscore(chrom, x, y,params[["GA.val"]],params[["GA.feat_penalty"]]),
        nBits = n_features,
        popSize = 100,         # population size
        maxiter = 100,        # max generations
        run = 50,             # stop if no improvement after 50 generations
        pmutation = 0.1      # mutation probability 
      )
      fsel_genes = colnames(x)[GA_res@solution[1,]==1]
      sign_gene <-c()
      for (g in fsel_genes) {
        cox_fit <- coxph(Surv(y[,"time"], y[,"status"]) ~ x[, g, drop=TRUE])
        sign_gene[g] = ifelse(summary(cox_fit)$coefficients[,"coef"]>0,1,-1)
      }
      x <- x[, fsel_genes, drop = FALSE]  # keep only significant features
    }
    
    if(feat_sel =="NMF"){
      #fsel_genes <- c()
      
      nmf_results = nmf_model(x, y[,"status"], n_rank = params[["n_rank_nmf"]], cnv_cols =params[["cnv_cols"]],pseudo_count = params[["pseudocount"]])
      
      x = nmf_results$new_df
      
    }
  }
  
  x_transf_info=NULL
  base_surv = NULL
 
  if(model_method =="coxnet"){
    if(ncol(x)>0){
      model_fit <- cv.glmnet(x = as.matrix(x), y = y, family = "cox", alpha = params[["alpha"]],type.measure="C")
      coef_full <- coef(model_fit, s = params[["s"]])
      selected_coefs <- rownames(coef_full)[as.numeric(coef_full) != 0]
      coef_active_full <- coef_full[as.numeric(coef_full) != 0]
      
      surv_prob <- 1-as.numeric(as.matrix(x[, selected_coefs, drop = FALSE]) %*% as.matrix(coef_active_full)) # this is the same as 1-predict(full_cv,s="lambda.min",newx=x)[,1]. So it is monotically relate with the survival probability, but it cannot be interepreted as a probability!
      # For computing survival:
      lp_train = predict(model_fit,s=params[["s"]],newx=as.matrix(x))[,1]
      # Fit a Cox model using the linear predictor as offset
      cox_fit <- coxph(y ~ offset(lp_train))
      
      # Get baseline survival
      base_surv <- list("Survfit"= survfit(cox_fit),"LP_avg" = mean(lp_train))
      
    }else{
      model_fit<-cv.glmnet(x = as.matrix(rep(0.5,nrow(x)),ncol=1), y = y, family = "cox", alpha = params[["alpha"]],type.measure="C")
      selected_coefs <-colnames(x)
      surv_prob = rep(0.5,nrow(x))
      # For computing survival:
      lp_train = rep(1,nrow(x))
      # Fit a Cox model using the linear predictor as offset
      cox_fit <- coxph(y ~ offset(lp_train))
      
      # Get baseline survival
      base_surv <- list("Survfit"= survfit(cox_fit),"LP_avg" = mean(lp_train))
      
    }
    
  }else if(model_method =="rsf"){
    df_rsf <- as.data.frame(x)
    df_rsf$FU <- y[, "time"]
    df_rsf$Event <- y[, "status"]
    
    model_fit <- ranger(
      formula = Surv(FU, Event) ~ .,
      data = df_rsf,
      num.trees = params[["num.trees"]],
      importance = "permutation",
      max.depth = params[["max.depth"]],
      min.node.size=params[["min.node.size"]],
      seed = seed
    )
    selected_coefs <-colnames(x)
    surv_prob <- predict(model_fit, data = df_rsf)$survival[,5]
  }else if(model_method == "Onco_score"){
    model_fit<-sign_gene
    surv_prob = 1-as.numeric(x %*%sign_gene[sign_gene!=0])
    selected_coefs <-colnames(x)
  }else if(model_method =="late_integration_GEP_other"){
    GEP_var="GEP_Pred"
    x_GEP =  data.frame(Intercept=1,GEP_model=x[,"GEP_Pred"])
    #print(x_GEP)
    # Fit a model with GEP only:
    model_fit_GEP <- cv.glmnet(x = as.matrix(x_GEP), y = y, family = "cox", alpha = params[["alpha"]],type.measure="C")
    coef_full <- coef(model_fit_GEP, s = params[["s"]])
    selected_coefs <- rownames(coef_full)[as.numeric(coef_full) != 0]
    coef_active_full <- coef_full[as.numeric(coef_full) != 0]
    
    surv_prob_GEP <- 1-as.numeric(as.matrix(x_GEP[, selected_coefs, drop = FALSE]) %*% as.matrix(coef_active_full)) # this is the same as predict(full_cv,s="lambda.min",newx=x)[,1]
    #print( concordance(y ~ surv_prob_GEP)$concordance)
    
    #print(model_fit)
    x_remaining = x[,setdiff(colnames(x),"GEP_Pred")]
    
    #print(x_remaining)
    # Fit a model with the remaining genes:
    set.seed(seed)
    #model_fit_remaining <- cv.glmnet(x = as.matrix(x_remaining), y = y, family = "cox", alpha = params[["alpha"]],type.measure="C")
    #coef_full <- coef(model_fit_remaining, s = params[["s"]])
    #selected_coefs <- rownames(coef_full)[as.numeric(coef_full) != 0]
    #coef_active_full <- coef_full[as.numeric(coef_full) != 0]
    
    #surv_prob_other <- 1-as.numeric(as.matrix(x_remaining[, selected_coefs, drop = FALSE]) %*% as.matrix(coef_active_full)) # this is the same as predict(full_cv,s="lambda.min",newx=x)[,1]
    df_rsf <- as.data.frame( x_remaining )
    df_rsf$FU <- y[, "time"]
    df_rsf$Event <- y[, "status"]
    
    model_fit_remaining <- ranger(
      formula = Surv(FU, Event) ~ .,
      data = df_rsf,
      num.trees = params[["num.trees"]],
      importance = "permutation",
      max.depth = params[["max.depth"]],
      min.node.size=params[["min.node.size"]],
      seed = seed
    )
    selected_coefs <-colnames( x_remaining )
    surv_prob_other <- predict(model_fit_remaining, data = df_rsf)$survival[,5]
    #print( concordance(y ~ surv_prob_other)$concordance)
    
    late_int_matrix = data.frame(GEP_score = surv_prob_GEP , Other_alt_score=surv_prob_other)
    set.seed(seed)
    model_fit = cv.glmnet(x = as.matrix(late_int_matrix), y = y, family = "cox", alpha = 0,type.measure="C")
    
    coef_full <- coef(model_fit, s = params[["s"]])
    
    selected_coefs <- rownames(coef_full)[as.numeric(coef_full) != 0]
    coef_active_full <- coef_full[as.numeric(coef_full) != 0]
    
    surv_prob <- 1-as.numeric(as.matrix(late_int_matrix[, selected_coefs, drop = FALSE]) %*% as.matrix(coef_active_full)) # this is the same as 1-predict(full_cv,s="lambda.min",newx=x)[,1]. So it is monotically relate with the survival probability, but it cannot be interepreted as a probability!
    #print( concordance(y ~ surv_prob)$concordance)
    x_transf_info = list(GEP_model = model_fit_GEP,Other_model = model_fit_remaining)
    
  }
  
  if(feat_sel =="NMF"){
    #model_fit$NMF_basis = nmf_results$nmf.w
    print("This function is not completed.")
  }
  #print(surv_prob)
  return(list(model_fit = model_fit,selected_coefs = selected_coefs,surv_prob = surv_prob,x_transf_info = x_transf_info,base_surv = base_surv ))
}

apply_model = function(x, built_model,model_method,params=list(),y){
  model_fit = built_model$model_fit
  if(model_method =="coxnet"){
    coef_full <- coef(model_fit, s = params[["s"]])
    coef_active_full <- coef_full[as.numeric(coef_full) != 0]
    selected_coefs <- rownames(coef_full)[as.numeric(coef_full) != 0]
    #surv_prob <- 1-as.numeric(as.matrix(x[, selected_coefs, drop = FALSE]) %*% as.matrix(coef_active_full)) # this is the same as predict(full_cv,s="lambda.min",newx=x)[,1]
    
    lp_new = predict(model_fit,s=params[["s"]],newx=as.matrix(x[,rownames(model_fit$glmnet.fit$beta)]))[,1]
    
    surv_prob = get_survival_prob(built_model$base_surv$Survfit,built_model$base_surv$LP_avg,lp_new,t=5)
    #print(surv_prob)
  }else if (model_method == "rsf") {
    df_rsf <- as.data.frame(x)
    surv_prob <- predict(model_fit, data = df_rsf)$survival[,5]
  }else if (model_method =="Onco_score"){
    surv_prob = 1 - as.numeric(x[,names(model_fit[model_fit!=0])]%*%model_fit[model_fit!=0])
  }else if(model_method =="late_integration_GEP_other"){
    # separate the variables:
    x_GEP =  data.frame(Intercept=1,GEP_model=x[,"GEP_Pred"])
    x_remaining = x[,setdiff(colnames(x),"GEP_Pred")]
    #print(dim(x_GEP))
    #print(dim(x_remaining))
    # apply the two models:
    surv_prob_GEP <- 1-predict(built_model$x_transf_info$GEP_model,s=params[["s"]],newx=as.matrix(x_GEP))[,1]
    surv_prob_other  <- predict(built_model$x_transf_info$Other_model, data = as.data.frame(x_remaining))$survival[,5] # if rsf
    #surv_prob_other <- 1-predict(built_model$x_transf_info$Other_model,s=params[["s"]],newx=as.matrix(x_remaining))[,1] #if coxnet
    late_int_matrix = data.frame(GEP_score = surv_prob_GEP , Other_alt_score=surv_prob_other)
    surv_prob <- 1-predict(model_fit,s=params[["s"]],newx=as.matrix(late_int_matrix))[,1]
    
  }
  return(surv_prob)
}

get_survival_prob = function(base_surv,lp_train_avg,lp_new,t=5){
  # Example: survival probability at time t 
  S_t <- summary(base_surv,times = t)$surv ^ exp(lp_new-lp_train_avg)
  
  return(S_t)
}
bootstrap_model <- function(x, y, model_method, n_boot = 100, params = list(),save_path=NULL) {
  
  # Initialize containers:
  n <- nrow(x)
  oob_cindex <- numeric(n_boot)
  optimism_vec <- numeric(n_boot)
  gene_selection_freq <- setNames(rep(0, ncol(x)), colnames(x))
  nr_genes_model <-numeric(n_boot)
  models_list <- vector("list", n_boot)
  C_632b <- numeric(n_boot)
  C_632plus_b <- numeric(n_boot)
  oob_predictions = matrix(NA,nrow = nrow(x),ncol=n_boot)
  rownames(oob_predictions)=rownames(x)
  bs_predictions= matrix(NA,nrow = nrow(x),ncol=n_boot)
  rownames(bs_predictions)=rownames(x)
  
  fsel = params[["fsel"]]
  # Apparent C-index:
  model_training <- build_model(x, y, model_method,params = params)
  app_seed = 1
  while(length(model_training$selected_coefs)==0 & app_seed <20){
    model_training <- build_model(x, y, model_method,seed = app_seed+1000, params = params)
    app_seed=app_seed+1
  }
  
  C_apparent <- survival::concordance(y ~ model_training$surv_prob)$concordance
  cat("Apparent C-index:", round(C_apparent, 3), "\n")
  
  # Calcualte C_no_info:
  n_perm <- 100
  C_noinfo_vec <- numeric(n_perm)
  for (i in 1:n_perm) {
    random_pred <- sample(model_training$surv_prob)
    C_noinfo_vec[i] <- concordance(y ~ random_pred)$concordance
  }
  C_noinfo <- mean(C_noinfo_vec)
  cat("No-information C-index:", round(C_noinfo, 3), "\n")
  
  #C_noinfo <- 0.5
  print(C_noinfo)
  for (b in 1:n_boot) {
    cat(sprintf("[%s] Bootstrap iteration: %d\n", model_method, b))
    set.seed(b)
    
    boot_idx <- sample(1:n, replace = TRUE)
    oob_idx <- setdiff(1:n, boot_idx)
    stopifnot(length(oob_idx)!=0)
    
    x_boot <- x[boot_idx, ]
    y_boot <- y[boot_idx]
    x_oob  <- x[oob_idx, ]
    y_oob  <- y[oob_idx]
    
    fit_bootstrap <- build_model(x_boot, y_boot, model_method,seed = b, params = params)
    bs_seed = 1
    while(length(fit_bootstrap$selected_coefs)==0 & bs_seed <20){
      fit_bootstrap <- build_model(x_boot, y_boot, model_method,seed = bs_seed+1000, params = params)
      bs_seed=bs_seed+1
    }
    active_genes <- fit_bootstrap$selected_coefs
    nr_genes_model[b]=length(active_genes)
    print(active_genes)
    if (!is.null(active_genes) && length(active_genes) > 0) {
      gene_selection_freq[active_genes] <- gene_selection_freq[active_genes] + 1
    }
    surv_boot <- apply_model(x_boot, fit_bootstrap, model_method, params = params)
    surv_oob  <- apply_model(x_oob, fit_bootstrap, model_method, params = params)
    
    C_boot <- concordance(y_boot ~ surv_boot)$concordance
    #print(C_boot)
    #print(concordance(y_boot ~ fit_bootstrap$surv_prob)$concordance)
    C_oob  <- concordance(y_oob ~ surv_oob)$concordance
    # Save OOB:
    oob_predictions[rownames(x_oob ),b] = surv_oob
    bs_predictions[rownames(x_boot),b]= surv_boot
    
    oob_cindex[b] <- C_oob
    optimism_vec[b] <- C_boot - C_oob
    models_list[[b]] <- fit_bootstrap$model_fit
    C_632b[b] <- 0.632 * C_oob + 0.368 * C_apparent
    
    R_b <- (C_apparent - C_oob) / (C_apparent - C_noinfo)
    R_b <-ifelse(is.na(R_b),0,R_b)
    R_b <- pmin(pmax(R_b, 0), 1)
    w_b <- 0.632 / (1 - 0.368 * R_b)
    C_632plus_b[b] <- (1 - w_b) * C_apparent + w_b * C_oob
  }
  
  # Harrel:
  mean_oob_cindex <- mean(oob_cindex, na.rm = TRUE)
  mean_optimism <- mean(optimism_vec, na.rm = TRUE)
  C_corrected <- C_apparent - mean_optimism
  
  # C_632:
  C_632 = 0.632*mean_oob_cindex +0.368* C_apparent
  
  #.632+:
  
  #R <- (C_apparent - mean_oob_cindex ) / (C_apparent-C_noinfo)
  #print(R)
  #R <- min(max(R, 0), 1)  # constrain between 0 and 1
  
  #w <- 0.632 / (1 - 0.368 * R)
  #C_632plus <- (1 - w) * C_apparent + w * mean_oob_cindex 
  C_632plus = mean(C_632plus_b)
  # confidence intervals:
  ci_632 <- quantile(C_632b, c(0.025, 0.975), na.rm=TRUE)
  ci_632plus <- quantile(C_632plus_b, c(0.025, 0.975), na.rm=TRUE)
  CIs <- list(
    C632 = ci_632,
    C632plus = ci_632plus
  )
  
  results = list(
    model_method = model_method,
    model_training=model_training,
    nr_genes_bs_model = nr_genes_model,
    C_apparent = C_apparent,
    mean_oob_cindex = mean_oob_cindex,
    mean_optimism = mean_optimism,
    C_corrected = C_corrected,
    C_632 = C_632,
    C_632plus = C_632plus,
    oob_cindex = oob_cindex,
    optimism_vec = optimism_vec,
    gene_selection_freq = gene_selection_freq,
    #models = models_list,
    C_632plus_all = C_632plus_b,
    C_632_all = C_632b,
    C_harrel_all = C_apparent-optimism_vec,
    C_oob_all = oob_cindex,
    oob_predictions=oob_predictions,
    bs_predictions=bs_predictions
  )
  
  # Save detailed results as .rds (optional)
  if (!is.null(save_path)) {
    saveRDS(results, file = file.path(save_path, paste0("b03_02_06_bootstrap_results_", model_method, ".rds")))
  }
  
  
  return(results)
}

run_bootstrap_all_models <- function(x_list, y, model_list, n_boot = 100, outdir = "bootstrap_results",out_models = NULL) {
  dir.create(outdir, showWarnings = FALSE)
  
  results_all <- list()
  summary_df <- data.frame()
  
  
  for(feats_set in names (x_list)){
    x = x_list[[feats_set]]
    for (m in names(model_list)) {
      
      if (m %in% c("random_coxnet", "random_rsf")) {
        results_random = list()
        for (i in 1:20) {
          set.seed(i)
          # fresh permutation every time
          y_random <- Surv(time=y[,"time"][sample(1:length(y),replace=FALSE)],
                           event = y[,"status"][sample(1:length(y),replace=FALSE)])
          
          
          if(m =="random_coxnet"){
            #y_random <- y[sample(seq_along(y))]  # full random permutation
            res <- bootstrap_model(x, y_random, model_method = "coxnet", n_boot = n_boot, params = model_list[[m]], save_path = out_models)
            
          }else if(m =="random_rsf"){
            #y_random <- y[sample(seq_along(y))]  # full random permutation
            res <- bootstrap_model(x, y_random, model_method = "rsf", n_boot = n_boot, params = model_list[[m]], save_path = out_models)
          }
          # unique id for this run
          run_id <- paste0(feats_set, "_", m, "_rep", i)
          results_random[[run_id]] <- res
        } 
        
        res$C_apparent = mean(sapply(results_random, function(x) x$C_apparent), na.rm = TRUE)
        res$C_632plus = mean(sapply(results_random, function(x) x$C_632plus), na.rm = TRUE)
        res$C_632 = mean(sapply(results_random, function(x) x$C_632), na.rm = TRUE)
        res$C_corrected = mean(sapply(results_random, function(x) x$C_corrected ), na.rm = TRUE)
        
      }else{
        res <- bootstrap_model(x, y, model_method = model_list[[m]][["method"]], n_boot = n_boot, params = model_list[[m]], save_path = out_models)
        
      }
      results_all[[paste0(feats_set,"_",m)]] <- res
      
      summary_df <- rbind(summary_df, data.frame(
        feats = feats_set,
        model = m,
        pipeline_id = paste0(feats_set,"_",m),
        C_apparent = round(res$C_apparent, 3),
        C_632plus = round(res$C_632plus, 3),
        C_632 = round(res$C_632, 3),
        C_corrected = round(res$C_corrected, 3),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  
  summary_path <- file.path(outdir, paste0("b03_02_06_bootstrap_summary_",n_boot,".csv"))
  write.csv(summary_df, summary_path, row.names = FALSE)
  
  cat("\n Summary saved to:", summary_path, "\n")
  
  for (var_c in c("C_632plus_all","C_632_all","C_harrel_all","C_oob_all","nr_genes_bs_model")){
    df_all <- get_all_c_estimates(results_all,var_c)
    
    plot_c_estimates(
      df_all,
      c_var = var_c,
      outpath = file.path(outdir, paste0("boxplot_",gsub("_all","",var_c), n_boot, ".pdf")),
      title_prefix = "Comparison of "
    )
    # Define colors: gray for random, blue for real models
    #colors <- c("Random" = "gray70", "Real" = "#1f77b4")
    
    # Plot
    #p_boxplot = ggplot(df_all, aes(x = model, y = C_632plus, fill = model_type)) +
    #  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    #  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
    #  scale_fill_manual(values = colors) +
    #  theme_minimal(base_size = 14) +
    #  labs(
    #    title = "Comparison of .632+ Bootstrap C-index Across Models",
    #    x = "Model",
    #    y = ".632+ C-index",
    #    fill = "Model Type"
    #  ) +
    #  theme(
    #    axis.text.x = element_text(angle = 45, hjust = 1),
    #    legend.position = "top"
    #  )
    
    #ggsave(file.path(outdir,paste0("boxplot_comparison_C632p_",n_boot,"bs.pdf")),p_boxplot,width=7,height=7)
    
  }
  
  # Show most important variables:
  plot_top_features_facet(results_all, n_boot, outdir, cutoff = 0)
  
  return(list(results = results_all, summary = summary_df))
}

get_all_c_estimates <- function(results,c_estimate_var) {
  df_list <- lapply(names(results), function(model_name){
    res <- results[[model_name]]
    data.frame(
      model = model_name,
      value= res[[c_estimate_var]]
    )
  })
  df_all <- bind_rows(df_list)
  
  # Add a column to distinguish random models
  print(colnames(df_all))
  df_all <- df_all %>%
    dplyr::mutate(model_type = ifelse(grepl("random", model), "Random", "Real"))
  
  # Keep the factor order the same as in final_results
  df_all$model <- factor(df_all$model, levels = names(results))
  
  return(df_all)
}

plot_c_estimates <- function(df_all, c_var, outpath, title_prefix = "") {
  
  # Define colors: gray for random, blue for real models
  colors <- c("Random" = "gray70", "Real" = "#1f77b4")
  
  p <- ggplot(df_all, aes(x = model, y = value, fill = model_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
    scale_fill_manual(values = colors) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste0(title_prefix, c_var, " Across Models"),
      x = "Model",
      y = gsub("_all","",c_var),
      fill = "Model Type"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )
  
  ggsave(outpath, p, width = 7, height = 7)
  
  return(p)
}

plot_top_features_facet <- function(results_all, n_boot, outdir, cutoff = 0) {
  
  df_list <- lapply(names(results_all), function(m) {
    freq <- results_all[[m]]$gene_selection_freq
    df <- data.frame(
      model = m,
      feature = names(freq),
      freq = as.numeric(freq),
      prop = as.numeric(freq) / n_boot
    )
    df[df$prop >= cutoff, ]
  })
  
  df_all <- do.call(rbind, df_list)
  
  if (nrow(df_all) == 0) {
    message("No features exceed the cutoff across any model.")
    return(NULL)
  }
  
  p <- ggplot(df_all, aes(x = reorder_within(feature, freq, model), y = prop*100)) +
    geom_col(fill = "#1f77b4") +
    coord_flip() +
    facet_wrap(~ model, scales = "free_y") +
    scale_x_reordered() +
    labs(
      title = paste0("Features Selected in ≥", cutoff*100, "% of Bootstraps"),
      x = "Feature",
      y = paste0("Percentage of bootstraps in which feature was selected (%)")
    ) +
    theme_minimal(base_size = 14)+
    theme(strip.text.x = element_text(size = 8),  # column facets
                  strip.text.y = element_text(size = 8) )
  
  ggsave(
    file.path(outdir, paste0("top_features_faceted_", n_boot, "bs.pdf")),
    p,
    width = 10, height = 8
  )
  
  return(p)
}

# This function is an internal function from the c060 package, and allows to compute the baseline hazard/survival
basesurv_c060 <- function (response, lp, times.eval = NULL, centered = FALSE) 
{
  if (is.null(times.eval)) 
    times.eval <- sort(unique(response[, 1]))
  t.unique <- sort(unique(response[, 1][response[, 2] == 1]))
  alpha <- length(t.unique)
  for (i in 1:length(t.unique)) {
    alpha[i] <- sum(response[, 1][response[, 2] == 1] == 
                      t.unique[i])/sum(exp(lp[response[, 1] >= t.unique[i]]))
  }
  obj <- stats::approx(t.unique, cumsum(alpha), yleft = 0, xout = times.eval, 
                       rule = 2)
  if (centered) 
    obj$y <- obj$y * exp(mean(lp))
  obj$z <- exp(-obj$y)
  names(obj) <- c("times", "cumBaseHaz", "BaseSurv")
  return(obj)
}

predict_survival_coxnet = function(fit_glmnet, s,sur_tr,tr_data,newdata,times){
  lp1_tr <- as.numeric(predict(fit_glmnet, newx = as.matrix(tr_data), s = s, type = "link")) # get the linear predictor on the old data
  lp2_newdata <- as.numeric(predict(fit_glmnet, newx = as.matrix(newdata), s = s, type = "link")) # get the linear predictor on the new data
  coef_full <- coef(fit_glmnet, s =  s)
  selected_coefs <- rownames(coef_full)[as.numeric(coef_full) != 0]
  lp2_newdata_check =  coef_full[selected_coefs,] %*% t(as.matrix(newdata)[,selected_coefs])
  
  print(lp2_newdata[1:10])
  print(lp2_newdata_check[1:10])
  basesur <- basesurv_c060(response=sur_tr, lp=lp1_tr, times.eval=times)
  probs <- exp(exp(lp2_newdata) %*% -t(basesur$cumBaseHaz))    # same as ePCR package
  
  # checked that it is the same result as with baseline survival:
  S0 <- exp(-basesur$cumBaseHaz)
  probs = S0 ^ exp(lp2_newdata)
  return(list("baseline_survival"=S0,"survival_prob"=probs))
  
}
