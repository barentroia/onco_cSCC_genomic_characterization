
default_seurat_processing_PCA <- function(seurat_object) {
  #' Default processing of seurat object, normalization + scaling + PCA
  #' 
  #' @param seurat_object seurat object 
  #'
  #' @return processed seurat object
  #'  
  
  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(seurat_object)
  seurat_object <- ScaleData(object = seurat_object, features =all.genes)
  
  top.genes<-VariableFeatures(seurat_object)
  seurat_object <- RunPCA(seurat_object, npcs = 25, verbose = FALSE,features=top.genes)
  
  return(seurat_object)
}

default_seurat_UMAP <- function(seurat_object, reduction="pca", dims=15) {
  #' Default processing of seurat object, UMAP
  #' 
  #' @param seurat_object seurat object 
  #' @param reduction character. Reduction, default is pca 
  #' @param dims numeric. Number of dimensisions, default is 15. 
  #'
  #' @return processed seurat object
  #'  
  
  seurat_object <- FindNeighbors(seurat_object, reduction = reduction, dims = 1:dims)
  seurat_object <- FindClusters(seurat_object)
  seurat_object <- RunUMAP(seurat_object, reduction = reduction, dims = 1:dims)
  
  return(seurat_object)
}

custom_dotPlot <- function(seurat_obj, features, group_by){
  #' Custom dot plot based on the dotPlot function in Seurat
  #' 
  #' @param seurat_object seurat object 
  #' @param features list. List of features (genes) to be plotted 
  #' @param group_by character. Variable by which the cells are grouped
  #'
  #' @return data frame for plotting the dotplot with ggplot2
  #'  
  
  # extract expression and generate dot plot of model genes
  data.features <- FetchData(object = seurat_obj, vars = features, cells = rownames(seurat_obj@meta.data))
  data.features$ident_custom <- seurat_obj@meta.data[match(rownames(data.features),rownames(seurat_obj@meta.data)),][[group_by]]
  
  data.plot <- data.features %>%
    group_split(ident_custom) %>%
    map(function(df) {
      df_use <- df %>% select(-ident_custom)
      # average expression per cell type
      avg.exp <- df_use %>%
        summarise(across(everything(), ~ mean(expm1(.x)))) %>%
        unlist()
      # percentage of cells expression gene per cell type
      pct.exp <- df_use %>%
        summarise(across(everything(), ~ PercentAbove(.x, threshold = 0))) %>%
        unlist()
      
      list(ident_custom = unique(df$ident_custom),avg.exp = avg.exp, pct.exp = pct.exp)
    })
  
  # turn into data frame for plotting
  data.plot.df <- map_dfr(data.plot, function(x) {
    tibble(ident_custom = x$ident_custom, feature = names(x$avg.exp), avg.exp = unname(x$avg.exp), 
           pct.exp = unname(x$pct.exp))
  }) %>% as.data.frame()
  
  # scale data per gene
  data.plot.df$avg.exp.scaled <- log1p(data.plot.df$avg.exp+1) # similar to Seurat
  for (feature_x in unique(data.plot.df$feature)){
    data.plot.df[data.plot.df$feature == feature_x,]$avg.exp.scaled <- scale(data.plot.df[data.plot.df$feature == feature_x,]$avg.exp.scaled)
  }
  
  return(data.plot.df)
}


custom_dotPlot_sum <- function(seurat_obj, features, group_by){
  #' Custom dot plot based on the dotPlot function in Seurat
  #' 
  #' @param seurat_object seurat object 
  #' @param features list. List of features (genes) to be plotted 
  #' @param group_by character. Variable by which the cells are grouped
  #'
  #' @return data frame for plotting the dotplot with ggplot2
  #'  
  
  # extract expression and generate dot plot of model genes
  data.features <- FetchData(object = seurat_obj, vars = features, cells = rownames(seurat_obj@meta.data))
  data.features$ident_custom <- seurat_obj@meta.data[match(rownames(data.features),rownames(seurat_obj@meta.data)),][[group_by]]
  
  data.plot <- data.features %>%
    group_split(ident_custom) %>%
    map(function(df) {
      df_use <- df %>% select(-ident_custom)
      # average expression per cell type
      avg.exp <- df_use %>%
        summarise(across(everything(), ~ sum(.x))) %>%
        unlist()
      # percentage of cells expression gene per cell type
      pct.exp <- df_use %>%
        summarise(across(everything(), ~ PercentAbove(.x, threshold = 0))) %>%
        unlist()
      
      list(ident_custom = unique(df$ident_custom),avg.exp = avg.exp, pct.exp = pct.exp)
    })
  
  # turn into data frame for plotting
  data.plot.df <- map_dfr(data.plot, function(x) {
    tibble(ident_custom = x$ident_custom, feature = names(x$avg.exp), avg.exp = unname(x$avg.exp), 
           pct.exp = unname(x$pct.exp))
  }) %>% as.data.frame()
  
  # scale data per gene
  data.plot.df$avg.exp.scaled <- log1p(data.plot.df$avg.exp+1) # similar to Seurat
  for (feature_x in unique(data.plot.df$feature)){
    data.plot.df[data.plot.df$feature == feature_x,]$avg.exp.scaled <- scale(data.plot.df[data.plot.df$feature == feature_x,]$avg.exp.scaled)
  }
  
  return(data.plot.df)
}
