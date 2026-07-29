# select_features function: Wrapper function for different features selection method
# @param ds_x: dataset with all of the features, columns should be features with colnames identifying the features
# @param ds_y: dataset with outcomes of interest and matching information
# @param sel_method: character indicating the feature selection method to use
# @param N: integer indicating the number of features to select
# @param outcome_name: character vector corresponding to the column name of the outcome of interest in the ds_y dataframe
# @param fup_name: character vector corresponding to the column name of the follow-up time of the outcome of interest in the ds_y dataframe
# @param weights_name: character vector corresponding to the column name of the weights in the ds_y dataframe
# @param pairs_name: character corresponding to the column name of the pairs in the ds_y dataframe
# @param sampleid_name: character corresponding to the column name of the sample IDs in the ds_y dataframe
# @param ds_x_counts: counts dataset with all of the features, columns should be features with colnames identifying the features
# @param model_params: vector of model parameters
# @param genes_oi: list of genes of interest
select_features <- function(ds_x, ds_y, sel_method, N, outcome_name, fup_name, weights_name, pairs_name, sampleid_name, ds_x_counts, genes_oi, model_params){

    # Select features based on DGEA with DESeq2
     if (sel_method == "deseq2"){ranked_feats <- select_features_deseq2(ds_x, ds_y, N, outcome_name, ds_x_counts, genes_oi, model_params)
     }
  # No features selection
  else if (sel_method == "none"){
    # But genes of interest are given as input
    if (!is.null(genes_oi)){
      ranked_feats <- genes_oi
    } else {
      ranked_feats <- NULL
    }}
    # Wrong features selection method
    else {
        stop("Wrong features selection method provided.")
    }

    ranked_feats
}


# select_features_deseq2 function: Function to select features based on DGEA with DESeq2
# @param x: dataset with all of the features, columns should be features with colnames identifying the features
# @param y: dataset with outcomes of interest and matching information
# @param n: integer indicating the number of features to select
# @param outcome_n: character vector corresponding to the column name of the outcome of interest in the y dataframe
# @param x_counts: dataset with all of the features, columns should be features with colnames identifying the features (counts data)
# @param goi: list of genes of interest
# @param model_params: vector of model parameters
select_features_deseq2 <- function(x, y, n, outcome_n, x_counts, goi, model_params){
    if (is.null(x_counts)){stop("No counts data provided to run DGEA with DESeq2.")}
    if (is.null(model_params$deseq2_design)){stop("No DESeq2 design provided or provided with wrong variable name.")}
    # Extract design
    design <- model_params$deseq2_design
    # During exploratory analyses, some samples were problematic and
    # we excluded them, we are going to do the same thing here
    samples_to_exclude <- c("S7032", "S7886", "S7099", "S8254")
    y_deseq <- y %>% dplyr::filter(!(rownames(.) %in% samples_to_exclude))
    y_deseq[, outcome_n] <- factor(y_deseq[, outcome_n], levels = c(0, 1))
    x_counts_deseq <- x_counts[, rownames(y_deseq)]
    # Check rows/columns order
    stopifnot(nrow(y_deseq) == ncol(x_counts_deseq))
    if (!all(rownames(y_deseq) == colnames(x_counts_deseq))){
        x_counts_deseq <- x_counts_deseq[, rownames(y_deseq)]
    }
    stopifnot(all(rownames(y_deseq) == colnames(x_counts_deseq)))
    # Run DGEA
    dds <- suppressMessages(DESeqDataSetFromMatrix(countData = round(x_counts_deseq),
                                  colData = y_deseq,
                                  design = as.formula(design)))
    dds <- suppressMessages(DESeq(dds, quiet = TRUE))
    res <- results(dds, alpha = 0.05) %>% as.data.frame() %>% arrange(pvalue)

    deseq_nr <- rownames(res)
    # Filter results based on significance, logFC value and average expression
    padj_thresh <- ifelse(is.null(model_params$deseq2_padj_thresh), 1, model_params$deseq2_padj_thresh)
    logfc_thresh <- ifelse(is.null(model_params$deseq2_logfc_thresh), 0, model_params$deseq2_logfc_thresh)
    baseMean_thresh <- ifelse(is.null(model_params$deseq2_baseMean_thresh), 0, model_params$deseq2_baseMean_thresh)
    res <- res %>%
        dplyr::filter(padj < padj_thresh &
                      abs(log2FoldChange) > logfc_thresh)
    deseq_afterfilt_nr <- rownames(res)
    # If required, iteratively filter results until the number of features is below the
    # max number of features requested by iteratively increasing the logFC threshold of 0.1
    if (!is.null(model_params$deseq2_feats_thresh_by_logfc)){
      initial_logfc_thresh <- logfc_thresh
        delta <- 0.10
        while(nrow(res) > model_params$deseq2_feats_thresh_by_logfc){
            res <- res %>% dplyr::filter(abs(log2FoldChange) > initial_logfc_thresh)
            initial_logfc_thresh <- initial_logfc_thresh + delta
        }
    }
    deseq_afterlarafilt_nr <- rownames(res)
    # Return filtered features
    ordered_feats <- rownames(res)
    # If required,, select only genes of interest
    if (!is.null(goi)){
        ordered_feats <- rownames(res)[rownames(res) %in% goi]
    }
    # If required, select top N features
    if (!is.na(n)){
      n <- min(n, length(ordered_feats))
        ordered_feats <- ordered_feats[1:n]
    }
    
    
    return(
      features = ordered_feats
)
}
