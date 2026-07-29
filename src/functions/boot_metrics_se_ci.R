library(boot)

boot_hr_se_ci <- function(df, outcome, time, weight, preds, pair, n_boot, metric) {
  #' Compute standard error of log HR using bootstrap
  #'
  #' @param df dataframe. Dataframe
  #' @param outcome string. String indicating column with outcomes
  #' @param time string. String indicating column with follow-up time
  #' @param weight string. String indicating column with weights
  #' @param preds string vector. String vector indicating column with predictors to use when fitting Cox model
  #' @param pair string or NULL. String indicating column with pairing info (set IDs)
  #' @param n_boot integer. Integer indicating number of bootstrap runs
  #' @param metric string. String indicating metric to compute (whr, hr)
  #' 
  #' @return standard error and 95% CI of log HR

  # Define bootstrap function
  boot_whr <- function(d, i) {
      di <- d[i, ]
      fit <- coxph(as.formula(paste0("Surv(", time, ", ", outcome, ") ~ ", paste(preds, collapse = " + ") )),
                   data = di,
                   weights = di[[weight]])
      log_whr <- coef(fit)[group]
      log_whr
  }
  boot_hr <- function(d, i) {
      di <- d[i, ]
      fit <- coxph(as.formula(paste0("Surv(", time, ", ", outcome, ") ~ ", paste(preds, collapse = " + ") )),
                   data = di)
      log_hr <- coef(fit)[group]
      log_hr
  }

  # Run (paired) bootstrap
  df <- as.data.frame(df)
  groups <- setdiff(colnames(model.matrix(~ . , df %>% select(all_of(preds)))), "(Intercept)")
  boot_res_groups <- vector("list", length(groups))
  names(boot_res_groups) <- groups
  if (is.null(pair)){
    set.seed(1234)
    if (metric == "whr"){
      boot_fun <- boot_whr
    } else if (metric == "hr"){
      boot_fun <- boot_hr
    }
    for (group in groups){
      boot_res_groups[[group]] <- boot(data = df, statistic = boot_fun, R = n_boot)
    }
    boot_res_df <- sapply(boot_res_groups, "[[", "t")
    i <- 1
    seed <- 1
    while (any(is.na(boot_res_df))){
      seed <- seed + i
      set.seed(seed)
      for (group in groups){
        boot_res_groups[[group]] <- boot(data = df, statistic = boot_fun, R = n_boot)
      }
      boot_res_df <- sapply(boot_res_groups, "[[", "t")
    }
    avg <- apply(boot_res_df, 2, function(x) mean(as.numeric(x)))
    se <- apply(boot_res_df, 2, function(x) sd(as.numeric(x)))
    ci95low <- apply(boot_res_df, 2, function(x) quantile(as.numeric(x), probs = 0.025))
    ci95up <- apply(boot_res_df, 2, function(x) quantile(as.numeric(x), probs = 0.975))
    ci95 <- paste(round(ci95low, 4), round(ci95up, 4), sep = "-")
    ci95 <- sapply(ci95, function(x) paste0("(", x, ")"))
    names(ci95) <- names(avg)
    est <- sapply(boot_res_groups, "[[", "t0")
    names(est) <- names(avg)
  } else {
    boot_res <- vector("list", n_boot)
    set.seed(1234)
    for (i in 1:n_boot){
      set_ids_sampled <- sample(unique(as.numeric(as.character(df[[pair]]))), replace = T)
      j <- unlist(purrr::map(set_ids_sampled, function(x) which(df[[pair]] == x)))
      boot_res[[i]] <- vector("list", length(groups))
      names(boot_res[[i]]) <- groups
      for(group in groups){
        if (metric == "whr"){
          boot_res[[i]][[group]] <- boot_whr(df, j)
        } else if (metric == "hr"){
          boot_res[[i]][[group]] <- boot_hr(df, j)
        }
      }
      s <- 1
      seeds <- round(runif(100)*10000)
      while (any(is.na(boot_res[[i]]))){
        seed <- seeds[[s]]
        set.seed(seed)
        set_ids_sampled <- sample(unique(as.numeric(as.character(df[[pair]]))), replace = T)
        j <- unlist(purrr::map(set_ids_sampled, function(x) which(df[[pair]] == x)))
        for(group in groups){
          if (metric == "whr"){
            boot_res[[i]][[group]] <- boot_whr(df, j)
          } else if (metric == "hr"){
            boot_res[[i]][[group]] <- boot_hr(df, j)
          }
        }
        s <- s + 1
      }
    }
    boot_res_df <- as.matrix(do.call(rbind, boot_res))
    avg <- apply(boot_res_df, 2, function(x) mean(as.numeric(x)))
    se <- apply(boot_res_df, 2, function(x) sd(as.numeric(x)))
    ci95low <- apply(boot_res_df, 2, function(x) quantile(as.numeric(x), probs = 0.025))
    ci95up <- apply(boot_res_df, 2, function(x) quantile(as.numeric(x), probs = 0.975))
    ci95 <- paste(round(ci95low, 4), round(ci95up, 4), sep = "-")
    ci95 <- sapply(ci95, function(x) paste0("(", x, ")"))
    names(ci95) <- names(avg)
    if (metric == "whr"){
      group <- groups
      est <- boot_whr(df, 1:nrow(df))
    } else if (metric == "hr"){
      group <- groups
      est <- boot_hr(df, 1:nrow(df))
    }
    names(est) <- names(avg)
  }  
  # Compute p-value: null hypothesis HR = 1 --> log(HR) = 0, double side test
  tstar <- 0
  ## Left side
  pl <- apply(boot_res_df, 2, function(x) sum(tstar < as.numeric(x)) / n_boot)
  ## Right side
  pr <- apply(boot_res_df, 2, function(x) sum(tstar > as.numeric(x)) / n_boot)
  p_comb <- apply(rbind(pl, pr), 2, function(x) max(2 * 1 / n_boot, 2 * min(x)))
  # print(p_comb)
  # Extract SE and 95% CI
  list(se = se, ci95 = ci95, ci95low = ci95low, ci95up = ci95up, avg = avg, est = est, pvalue = p_comb, boot_ests = boot_res_df)
}

boot_metrics_se_ci <- function(df, score, outcome, time, weight, pair, n_boot, metric) {
  #' Compute standard error of (w)Cindex/(w)AUC using bootstrap
  #'
  #' @param df dataframe. Dataframe
  #' @param score string. String indicating column with scores
  #' @param outcome string. String indicating column with outcomes
  #' @param time string. String indicating column with follow-up time
  #' @param weight string. String indicating column with weights
  #' @param pair string or NULL. String indicating column with pairing info (set IDs)
  #' @param n_boot integer. Integer indicating number of bootstrap runs
  #' @param metric string. String indicating metric to compute (wauc, wcindex, auc, cindex)
  #' 
  #' @return standard error and 95% CI of (w)Cindex/(w)AUC

  # Define bootstrap functions
  boot_wauc <- function(d, i) {
      di <- d[i, ]
      if (length(unique(di[[outcome]])) == 1){
        wauc <- NA
      } else {
        wauc <- WeightedAUC(WeightedROC(di[[score]], di[[outcome]], weight = di[[weight]]))
      }
      wauc
  }
  boot_wcindex <- function(d, i) {
      di <- d[i, ]
      if (length(unique(di[[outcome]])) == 1){
        wcindex <- NA
      } else {
        wcindex <- cIndex(di[[time]], event = di[[outcome]], di[[score]], weight = di[[weight]])[[1]]
      }
      wcindex
  }
  boot_auc <- function(d, i) {
      di <- d[i, ]
      if (length(unique(di[[outcome]])) == 1){
        auc <- NA
      } else {
        auc <- as.numeric(gsub(".*: ", "", pROC::roc(di[[outcome]], di[[score]], quiet = TRUE, direction = "<")$auc))
      }
      auc
  }
  boot_cindex <- function(d, i) {
      di <- d[i, ]
      if (length(unique(di[[outcome]])) == 1){
        cindex <- NA
      } else {
        cindex <- 1 - rcorr.cens(di[[score]], Surv(di[[time]], di[[outcome]]))["C Index"] %>% `names<-`(NULL)
      }
      cindex
  }

  # Run (paired) bootstrap
  df <- as.data.frame(df)
  if (is.null(pair)){
    if (metric == "wauc"){
      boot_fun <- boot_wauc
    } else if (metric == "wcindex"){
      boot_fun <- boot_wcindex
    } else if (metric == "auc"){
      boot_fun <- boot_auc
    } else if (metric == "cindex"){
      boot_fun <- boot_cindex
    }
    set.seed(123)
    boot_res <- boot(data = df, statistic = boot_fun, R = n_boot)
    i <- 1
    seed <- 123
    while (sum(is.na(boot_res$t)) > 0){
      seed <- seed + i
      set.seed(seed)
      boot_res <- boot(data = df, statistic = boot_fun, R = n_boot)
      i <- i + 1
      if (i > 200){
        boot_res$t0 <- na.omit(boot_res$t0)
        boot_res$t <- na.omit(boot_res$t)
        boot_res$R <- length(boot_res$t[,1])
        boot_res$call[4] <- boot_res$R
        break
        }
    }
    se <- sd(boot_res$t)
    avg <- mean(boot_res$t)
    est <- boot_res$t0
    ci95_res <- boot.ci(boot_res, conf = 0.95, type = "perc")
    #ci95_res <- boot.ci(boot_res, conf = 0.95, type = "perc") 
    if (is.null(ci95_res)){
        ci95 <- "(NA-NA)"
    } else {
        ci95 <- paste0("(", round(ci95_res$percent[1, 4], 4), "-", round(ci95_res$percent[1, 5], 4), ")")
    }
  } else {
    boot_res <- vector("list", n_boot)
    set.seed(123)
    for (i in 1:n_boot){
        set_ids_sampled <- sample(unique(as.numeric(as.character(df[[pair]]))), replace = T)
        # Find indexes of sampled pairs
        j <- unlist(purrr::map(set_ids_sampled, function(x) which(df[[pair]] == x)))
        if (metric == "wauc"){
            boot_res[[i]] <- boot_wauc(df, j)
        } else if (metric == "wcindex"){
            boot_res[[i]] <- boot_wcindex(df, j)
        } else if (metric == "auc"){
            boot_res[[i]] <- boot_auc(df, j)
        } else if (metric == "cindex"){
            boot_res[[i]] <- boot_cindex(df, j)
        }
    }
    se <- sd(unlist(boot_res))
    avg <- mean(unlist(boot_res))
    ci95low <- quantile(unlist(boot_res), probs = 0.025, na.rm = T)
    ci95up <- quantile(unlist(boot_res), probs = 0.975, na.rm = T)
    ci95 <- paste0("(", round(ci95low, 4), "-", round(ci95up, 4), ")")
    if (metric == "wauc"){
      est <- boot_wauc(df, 1:nrow(df))
    } else if (metric == "wcindex"){
      est <- boot_wcindex(df, 1:nrow(df))
    } else if (metric == "auc"){
      est <- boot_auc(df, 1:nrow(df))
    } else if (metric == "cindex"){
      est <- boot_cindex(df, 1:nrow(df))
    }
  }
  # Extract SE and 95% CI
  list(se = se, ci95 = ci95, avg = avg, est = est)
}