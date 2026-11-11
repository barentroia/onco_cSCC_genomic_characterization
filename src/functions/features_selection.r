#===============================================================================
# Aim:      Feature selection for the survival models, based on differential
#           gene expression analysis with DESeq2.
# Contents: - select_features(): wrapper for the feature selection methods
#           - select_features_deseq2(): DESeq2-based gene selection
# Usage:    Sourced by build_models.r
#===============================================================================

#' Select features with the requested method
#'
#' Wrapper around the available feature selection methods: "deseq2"
#' (select_features_deseq2()) or "none" (no selection; genes_oi is returned
#' if given).
#'
#' @param ds_x data.frame. Features (samples x features).
#' @param ds_y data.frame. Outcomes and matching information.
#' @param sel_method character. Feature selection method ("deseq2" or "none").
#' @param N integer. Maximum number of features to select (NA = no limit).
#' @param outcome_name character. Outcome column in ds_y.
#' @param fup_name character. Follow-up time column in ds_y (not used).
#' @param weights_name character. Weights column in ds_y (not used).
#' @param pairs_name character. Case-control set column in ds_y (not used).
#' @param sampleid_name character. Sample ID column in ds_y (not used).
#' @param ds_x_counts data.frame. Raw counts (genes x samples).
#' @param genes_oi character vector or NULL. Genes of interest to keep.
#' @param model_params list. Model parameters with the DESeq2 settings.
#'
#' @return Character vector of selected features ordered by p-value, or NULL
#'   when no selection is performed.
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


#' Select differentially expressed genes with DESeq2
#'
#' Runs DESeq2 with model_params$deseq2_design and keeps genes passing the
#' adjusted p-value (deseq2_padj_thresh, default 1) and |log2FC|
#' (deseq2_logfc_thresh, default 0) thresholds. If deseq2_feats_thresh_by_logfc
#' is set, the |log2FC| threshold is raised in steps of 0.1 until at most that
#' many genes remain. Four samples flagged during exploratory analyses are
#' excluded.
#'
#' @param x data.frame. Features (samples x genes); not used directly.
#' @param y data.frame. Outcomes, rownames = sample IDs.
#' @param n integer. Maximum number of genes to return (NA = no limit).
#' @param outcome_n character. Outcome column in y (0/1).
#' @param x_counts data.frame. Raw counts (genes x samples).
#' @param goi character vector or NULL. Genes of interest; if given, only
#'   these genes are kept.
#' @param model_params list. DESeq2 settings: deseq2_design (required),
#'   deseq2_padj_thresh, deseq2_logfc_thresh, deseq2_feats_thresh_by_logfc.
#'
#' @return Character vector of selected genes ordered by p-value.
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
