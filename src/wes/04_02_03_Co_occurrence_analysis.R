#---------------------------------------------------
# Aim: Co-occurence analysis
# Author: B. Rentroia Pacheco
# Input: Driver mutations
# Output: Overview of mutation co-occurrence and exclusivity
#---------------------------------------------------

gene_subset = genes_to_show
generate_co_occur_mat = function(df,gene_subset,cutoff_genes,mega_MML = FALSE,all_samples=NULL,sample_info_metrics = NA, Fisher_LogReg = "LogReg"){
  
  if(mega_MML){
    
    if(is.null(all_samples)){
      # Save ALL samples before mutation filtering
      all_samples <- unique(df$Tumor_Sample_Barcode)
    }
    
    # Apply some prefiltering to the MML table:
    df <-df %>%
      mutate(Metastasis = ifelse(grepl("_1", Tumor_Sample_Barcode), "Case", "Control")) %>%
      filter(Hugo_Symbol%in%gene_subset)%>%
      filter(Pathogenic %in% c("Probably"))%>%as.data.frame()
    
    # Apply cutoff:
    tab_genes = df%>%distinct(Hugo_Symbol, Tumor_Sample_Barcode, .keep_all = TRUE)%>%pull(Hugo_Symbol)%>%table()
    genes_with_significant_muts_mcutoff =names(tab_genes)[which(tab_genes>cutoff_genes)]
    
    
    # Prepare data for the heatmap:
    # 1. Find all mutations
    heatmap_data <-  df%>%
      filter(Hugo_Symbol %in% genes_with_significant_muts_mcutoff) %>%
      distinct(Hugo_Symbol, Tumor_Sample_Barcode, .keep_all = TRUE) %>%
      mutate(has_mutation = 1)
    
    # Create binary matrix from existing mutations, filling all other combinations with 0:
    mat <- heatmap_data %>%
      select(Hugo_Symbol, Tumor_Sample_Barcode, has_mutation) %>%
      complete(
        Hugo_Symbol,
        Tumor_Sample_Barcode = all_samples,
        fill = list(has_mutation = 0)
      ) %>%
      pivot_wider(names_from = Tumor_Sample_Barcode, values_from = has_mutation, values_fill = 0) %>%
      textshape::column_to_rownames("Hugo_Symbol") %>%
      as.matrix()
    
  }else{
    # Use the input matrix directly, after filtering:
    genes_with_significant_muts_mcutoff =setdiff(colnames(df)[colSums(df>0)>cutoff_genes],"Sample")
    
    mat <- df %>%
      textshape::column_to_rownames("Sample") %>%
      select(genes_with_significant_muts_mcutoff )%>%
      t()%>%
      as.matrix()
    mat[abs(mat) > 0] <- 1 #binarize
  }
  
  print(dim(mat))
  # Initialize matrices
  genes_local <- rownames(mat)
  n <- length(genes_local)
  logOR_mat <- matrix(NA_real_, n, n, dimnames = list(genes_local, genes_local))
  FDR_mat   <- matrix(NA_real_, n, n, dimnames = list(genes_local, genes_local))
  
  # Compute pairwise Fisher's tests
  gene_pairs <- combn(rownames(mat), 2, simplify = FALSE)
  self_pairs <- lapply(genes_local, function(g) c(g, g))
  
  # Compute all pairwise tests
  results <- lapply(gene_pairs, function(pair) {
    a <- mat[pair[1], ]
    b <- mat[pair[2], ]
    
    if(Fisher_LogReg =="Fisher"){
      # Fisher exact test:
      tbl <- table(a, b)
      if (all(dim(tbl) == c(2,2))) {
        test <- suppressWarnings(fisher.test(tbl))
        data.frame(
          Gene1 = pair[1],
          Gene2 = pair[2],
          OddsRatio = ifelse(is.finite(test$estimate), test$estimate, NA),
          p.value = test$p.value
        )
      }else {
        data.frame(Gene1 = pair[1], Gene2 = pair[2], OddsRatio = NA, p.value = NA)
      }
    }else if (Fisher_LogReg=="LogReg"){
        # Logistic regression:
        df_model <- data.frame(
          a = as.numeric(a),
          b = as.numeric(b),
          callable_cov =  sample_info_metrics[match(names(a),sample_info_metrics$Sample.ID),c("Percent_Bait_with_callable_coverage")],
          tmb = sample_info_metrics[match(names(a),sample_info_metrics$Sample.ID),c("Burden.(Mutations/megabase).(MB)")]

        )
        fit <-glm(a ~ b + callable_cov,
                  data = df_model,
                  family = binomial)
        coef_summary <- summary(fit)$coefficients
        
        beta <- coef_summary["b", "Estimate"]
        pval <- coef_summary["b", "Pr(>|z|)"]
        
        data.frame(
          Gene1 = pair[1],
          Gene2 = pair[2],
          OddsRatio = exp(beta),
          p.value = pval
        )
    }
  })
  results_df <- bind_rows(results)
  
  # Apply multiple-testing correction (Benjamini–Hochberg)
  results_df <- results_df %>%
    mutate(
      FDR = p.adjust(p.value, method = "BH")
    )
  
  # For the heatmap, we need the diagonal:
  # Compute all pairwise tests
  results_df_self <- lapply(self_pairs, function(pair) {
    a <- mat[pair[1], ]
    b <- mat[pair[2], ]
    tbl <- table(a, b)
    if (all(dim(tbl) == c(2,2))) {
      test <- suppressWarnings(fisher.test(tbl))
      data.frame(
        Gene1 = pair[1],
        Gene2 = pair[2],
        OddsRatio = ifelse(is.finite(test$estimate), test$estimate, NA),
        p.value = test$p.value
      )
    } else {
      data.frame(Gene1 = pair[1], Gene2 = pair[2], OddsRatio = NA, p.value = NA)
    }
  })
  results_df_self <- bind_rows(results_df_self )
  
  # Apply multiple-testing correction (Benjamini–Hochberg)
  results_df_self <- results_df_self %>%
    mutate(
      p.value=NA,
      FDR = NA
    )
  results_df_wself = rbind(results_df, results_df_self)%>%
    mutate(
      logOR = log2(OddsRatio),
      logOR = ifelse(is.infinite(logOR),5,logOR),
      logOR = ifelse(is.nan(logOR), NA, logOR))
  
  # Create symmetric logOR and FDR matrices
  logOR_mat <- reshape2::acast(results_df_wself, Gene1 ~ Gene2, value.var = "logOR", fill = NA)
  logOR_mat[lower.tri(logOR_mat)] <- t(logOR_mat)[lower.tri(logOR_mat)]
  diag(logOR_mat) <- 0
  logOR_mat[is.na(logOR_mat)] <- 0
  
  FDR_mat <- reshape2::acast(results_df_wself, Gene1 ~ Gene2, value.var = "FDR", fill = NA)
  FDR_mat[lower.tri(FDR_mat)] <- t(FDR_mat)[lower.tri(FDR_mat)]
  diag(FDR_mat) <- 1
  FDR_mat[is.na(FDR_mat)] <- 1
  
  # Significance symbols (FDR < 0.05)
  sig_mat <- ifelse(FDR_mat < 0.05, "*", "")
  
  # Define color scale
  col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  
  # Plot heatmap
  ht <- ComplexHeatmap::Heatmap(
    logOR_mat,
    name = "log2(OR)",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    border = TRUE,
    column_title = "Gene Co-occurrence (FDR < 0.05 = *)",
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (sig_mat[i, j] == "*") {
        grid.text("*", x, y, gp = gpar(fontsize = 10, col = "black", fontface = "bold"))
      }
    },
    row_names_side = "left",
    show_column_names = TRUE,
    show_row_names = TRUE
  )
  
  rownames(results_df)=NULL
  print(ht)
  return(results_df)
}
# Save ALL samples before mutation filtering, to make sure all samples are taken into account
all_samples <- unique(mega_MML$Tumor_Sample_Barcode)

pdf(file.path(dir_results_script,"Co_occurence_fisher_tests.pdf"))
res_all_Fisher = mega_MML_clonal_visualization%>%generate_co_occur_mat(gene_subset,9,mega_MML = TRUE,all_samples=all_samples,sample_info_metrics = sample_summary, Fisher_LogReg = "Fisher")
dev.off()
write.xlsx(res_all_Fisher,file.path(dir_results_script,"Co_occurence_fisher_tests_Fisher.xlsx"))

pdf(file.path(dir_results_script,"Co_occurence_LogReg_tests.pdf"))
res_all_LogReg = mega_MML_clonal_visualization%>%generate_co_occur_mat(gene_subset,9,mega_MML = TRUE,all_samples=all_samples,sample_info_metrics = sample_summary, Fisher_LogReg = "LogReg")
dev.off()
write.xlsx(res_all_LogReg,file.path(dir_results_script,"Co_occurence_fisher_tests_LogReg.xlsx"))

# Now the same but for pathways:
pathways_co_occur_mat = genomic_df_summary[,c("Sample",pathways_with_enough_mutations)]
pdf(file.path(dir_results_script,"Co_occurence_fisher_tests_pathway.pdf"))
res_pathway_fisher = pathways_co_occur_mat%>%generate_co_occur_mat(gene_subset,9,mega_MML = FALSE,sample_info_metrics = sample_summary, Fisher_LogReg = "Fisher")
dev.off()
write.xlsx(res_pathway_fisher ,file.path(dir_results_script,"Co_occurence_fisher_tests_pathway.xlsx"))

pdf(file.path(dir_results_script,"Co_occurence_fisher_tests_logreg_pathway.pdf"))
res_pathway_LogReg = pathways_co_occur_mat%>%generate_co_occur_mat(gene_subset,9,mega_MML = FALSE,sample_info_metrics = sample_summary, Fisher_LogReg = "LogReg")
dev.off()
write.xlsx(res_pathway_LogReg,file.path(dir_results_script,"Co_occurence_LogReg_tests_pathway.xlsx"))

#res_pathway_LogReg_cases = pathways_co_occur_mat%>%filter(grepl("_1",Sample))%>%generate_co_occur_mat(gene_subset,19,mega_MML = FALSE,sample_info_metrics = sample_summary, Fisher_LogReg = "LogReg")
#res_pathway_LogReg_controls = pathways_co_occur_mat%>%filter(grepl("_0",Sample))%>%generate_co_occur_mat(gene_subset,19,mega_MML = FALSE,sample_info_metrics = sample_summary, Fisher_LogReg = "LogReg")