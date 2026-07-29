run_pca <- function(object, ntop = 500){
  #' Plot PCA giving a DESeq2 object
  #'
  #' @param object matrix or dataframe. Gene expression matrix or dataframe
  #' @param ntop numeric. Number of genes (with highest variance) to compute PCs
  #'
  #' @return list with PCA prcomp results
  #'
  
  # Gene variance
  rv <- matrixStats::rowVars(object)
  
  # Select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # Perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(object[select,]), scale = T)
  
  pca
  
}


prepare_pca_dfs <- function(object, clnc){
  #' Prepare datasets with PCA results
  #'
  #' @param object prcomp results
  #' @param clnc dataframe. Dataframe with clinical info
  #'
  #' @return list dataframes used for plotting
  #'
  
  # Combine PCA info with clinical info
  # (datafrme with PC values for each sample)
  df <- object$x %>%
    as.data.frame() %>%
    mutate(samples = rownames(.)) %>%
    left_join(clnc %>%
                mutate(samples = rownames(.)), by = "samples")
  
  # Contribution to the total variance for each component
  # (dataframe with % variance explained by each PC)
  percent_var_df <- (object$sdev^2/sum(object$sdev^2)) %>% `names<-`(colnames(object$x)) %>%
    sapply(., function(x) round(100*x, 2)) %>%
    t() %>%
    as.data.frame()
  
  list("df_res" = df,
       "df_var_exp" = percent_var_df)
}

plot_pca <- function(df, percent_var_df, color, shape = "NA"){
  #' Plot PCA giving a prcomp result object
  #'
  #' @param object dataframe. Dataframe with PCA results
  #' @param percent_var_df dataframe. Dataframe with percentage of explained variances
  #' @param color character. Variable for color aesthetics
  #' @param shape character. Variable for shape aesthetics
  #'
  #' @return ggplot object
  #'
  
  if (shape != "NA"){
    p1 <- ggplot(data = df,
                 mapping = aes(x = PC1,
                               y = PC2,
                               color = .data[[color]],
                               shape = .data[[shape]]))+
      geom_text(hjust = 1, vjust = 1)
    p2 <- ggplot(data = df,
                 mapping = aes(x = PC3,
                               y = PC4,
                               color = .data[[color]],
                               shape = .data[[shape]]))+
      geom_text(hjust = 1, vjust = 1)
    p3 <- ggplot(data = df,
                 mapping = aes(x = PC1,
                               y = PC3,
                               color = .data[[color]],
                               shape = .data[[shape]]))+
      geom_text(hjust = 1, vjust = 1)
    p4 <- ggplot(data = df,
                 mapping = aes(x = PC1,
                               y = PC4,
                               color = .data[[color]],
                               shape = .data[[shape]]))+
      geom_text(hjust = 1, vjust = 1)
    
  } else {
    p1 <- ggplot(data = df,
                 mapping = aes(x = PC1,
                               y = PC2,
                               color = .data[[color]]))
    p2 <- ggplot(data = df,
                 mapping = aes(x = PC3,
                               y = PC4,
                               color = .data[[color]]))
    p3 <- ggplot(data = df,
                 mapping = aes(x = PC1,
                               y = PC3,
                               color = .data[[color]]))
    p4 <- ggplot(data = df,
                 mapping = aes(x = PC1,
                               y = PC4,
                               color = .data[[color]]))
    
  }
  
  p1 <- p1 + 
    geom_point(size=3) + 
    xlab(paste0("PC1: ", percent_var_df$PC1, "% variance")) +
    ylab(paste0("PC2: ", percent_var_df$PC2, "% variance")) +
    theme_bw()
  p2 <- p2+
    geom_point(size=3) + 
    xlab(paste0("PC3: ", percent_var_df$PC3, "% variance")) +
    ylab(paste0("PC4: ", percent_var_df$PC4, "% variance")) +
    theme_bw()
  p3 <- p3 + 
    geom_point(size=3) + 
    xlab(paste0("PC1: ", percent_var_df$PC1, "% variance")) +
    ylab(paste0("PC3: ", percent_var_df$PC3, "% variance")) +
    theme_bw()
  p4 <- p4+
    geom_point(size=3) + 
    xlab(paste0("PC1: ", percent_var_df$PC1, "% variance")) +
    ylab(paste0("PC4: ", percent_var_df$PC4, "% variance")) +
    theme_bw()
  
  if (is.character(df[[color]]) & length(unique(df[[color]])) > 10){
    p1 <- p1 +
      theme(legend.position = "none")
    p2 <- p2 +
      theme(legend.position = "none")
    p3 <- p3 +
      theme(legend.position = "none")
    p4 <- p4 +
      theme(legend.position = "none")
  }
  
  grid.arrange(p1, p2, p3, p4, nrow = 2)
  
}