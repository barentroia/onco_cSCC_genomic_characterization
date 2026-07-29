# Functions used to interpret the sample clusters


# Function to compute overrepresentation analysis (ORA) for genes within each cluster, using output from cola
# 
# Args:
#   annot             - string specifying one or more MSigDB categories, separated by "_"
#                        (e.g. "H_C2" for Hallmark + C2 gene sets)
#   sign_res          - data frame that is the output of cola's function get_signatures
#   output_folder_cola - path to the folder where the output Excel file will be saved
#   universe_genes    - character vector of background/universe gene symbols used for
#                        the enrichment test (denominator of the hypergeometric test)
#
# Output:
#   Writes an .xlsx file with one worksheet per cluster, each containing the ORA  results
ORA_gene_clusters = function(annot,sign_res,output_folder_cola,universe_genes){
  # --- Step 1: Split genes by cluster ---
  # Group the significant genes by their cluster label ("km" column)
  genes_by_cluster_named <- split(sign_res$Gene, sign_res$km)
  
  # Rename list elements to "Cluster_<number>" for clarity in downstream outputs
  names(genes_by_cluster_named) <- paste0("Cluster_", names(genes_by_cluster_named))
  
  # --- Step 2: Retrieve MSigDB gene sets ---
  # annot may contain multiple categories separated by "_", e.g. "H_C2" -> c("H","C2")
  msig_categories <- unlist(strsplit(annot, "_"))
  # Download/retrieve gene sets for each requested MSigDB category
  msig_list <- lapply(msig_categories, function(cat) {
    msigdbr(species = "Homo sapiens", category = cat)
  })
  
  # Combine all category gene sets into a single data frame
  msig_all <- do.call(rbind, msig_list)
  
  # Keep only the columns needed by enricher(): gene set name and gene symbol
  # (this forms the TERM2GENE mapping required by clusterProfiler::enricher)
  term2gene_all <- msig_all[, c("gs_name", "gene_symbol")]
  
  # --- Step 3: Convert cluster gene IDs (Ensembl) to gene symbols ---
  genes_by_cluster_symbols <- lapply(genes_by_cluster_named, function(genes) {
    # Strip version suffix from Ensembl IDs (e.g. "ENSG000001.5" -> "ENSG000001")
    ensembl_ids <- do.call(rbind, strsplit(genes, "[.]"))[,1]
    # Map Ensembl IDs to gene symbols using the external annotation table "out"
    # (assumed to be defined earlier in the script/environment,
    #  with columns "ensembl_gene_id" and "external_gene_name")
    symbols <- out$external_gene_name[
      match(ensembl_ids, out$ensembl_gene_id)
    ]
    # Remove unmapped (NA) or empty symbols, then deduplicate
    symbols <- symbols[!is.na(symbols) & symbols != ""]
    unique(symbols)
  })
  
  # --- Step 4: Run overrepresentation analysis (ORA) per cluster ---
  # enricher() performs a hypergeometric test for each cluster's gene list
  # against the combined MSigDB gene sets, using universe_genes as background
  ora_results <- lapply(genes_by_cluster_symbols, function(genes) {
    enricher(
      gene          = genes,
      universe      = universe_genes,
      TERM2GENE     = term2gene_all,
      pAdjustMethod = "BH",
      qvalueCutoff  = 0.05
    )
  })
  
  # --- Step 5: Extract and tidy up ORA results into data frames ---
  ora_df_list <- lapply(seq_along(ora_results), function(i) {
    res <- ora_results[[i]]
    # Handle clusters with no enrichment result (NULL object or empty result table)
    if (is.null(res) || nrow(res@result) == 0) {
      return(data.frame(Message = "No enriched pathways"))  # handle empty clusters
    }
    # Keep only the relevant columns from the enrichment result
    df <- res@result[, c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count")]
    df$Cluster <- i
    return(df)
  })
  
  # Restore cluster names (lost when using seq_along above)
  names(ora_df_list) <- names( genes_by_cluster_named )
  
  # --- Step 6: Reorder clusters numerically ---
  # Since cluster names are strings ("Cluster_1", "Cluster_10", "Cluster_2", ...),
  # sort them by their numeric suffix to avoid lexicographic ("1,10,2,...") ordering
  cluster_nums <- as.integer(sub("Cluster_", "", names(ora_df_list)))
  ora_df_list <- ora_df_list[order(cluster_nums)]
  
  # --- Step 7: Write results to an Excel workbook (one sheet per cluster) ---
  wb <- createWorkbook()
  
  for (cluster_name in names(ora_df_list)) {
    addWorksheet(wb, cluster_name)
    writeData(wb, sheet = cluster_name, ora_df_list[[cluster_name]])
  }
  
  # --- Step 8: Save the workbook to disk ---
  xlsx_file_name = file.path(output_folder_cola,paste0("ORA_results_clusters_",annot,".xlsx"))
  
  saveWorkbook(wb, file = xlsx_file_name, overwrite = TRUE)
  
}


# Function to generate Kaplan-Meier survival curves comparing metastasis-free survival across
# clusters (e.g. from a heatmap/hierarchical clustering). Returns a ggsurvplot object.
#
# Arguments:
#   ph                  - hierarchical clustering object (e.g. hclust/pheatmap) passed to
#                         get_samples_in_each_cluster() to derive cluster membership
#   nr_clusters         - number of clusters to cut the tree into
#   clinical_tumors_df  - clinical data frame; rownames are sample IDs (Skyline_ID),
#                         must contain "Metastasis" and "FU_metastasis_years" columns
#   i.subset            - optional row indices/names to restrict clinical_tumors_df to a subset
#   cluster_labels      - optional custom names to assign to the clusters
#   def_sample_clusters - optional pre-defined list of samples per cluster (bypasses ph clustering)
#   color_scheme        - optional named palette mapping cluster names to colors
#   max_time            - optional follow-up time cap (censors events beyond this horizon)
#   title_plot          - optional plot title

# Returns:
#   A ggsurvplot object (list containing the KM plot and risk table)
generate_surv_curv_hm_clusters<-function(ph,nr_clusters,clinical_tumors_df, i.subset=NULL, cluster_labels=NULL,def_sample_clusters=NULL,color_scheme=NULL,max_time=NULL,title_plot=""){
  
  # Restrict the clinical data to the requested subset of samples, if provided
  if(!is.null(i.subset)){
    clinical_tumors_df = clinical_tumors_df[i.subset,]
  }
  
  # Determine cluster membership, if not provided, and assign color palette
  if(is.null(def_sample_clusters)){
    # No predefined clusters: derive them from the clustering object.
    # Drop the color scheme since it may not match the auto-generated cluster names.
    samples_clustered=get_samples_in_each_cluster(ph,nr_clusters)
    color_scheme=NULL
  }else{
    # Use the supplied predefined clusters
    samples_clustered=def_sample_clusters
    # Keep the color scheme only if every cluster name has a matching color; otherwise drop it
    if(all(names(samples_clustered)%in%names(color_scheme))){
      color_scheme=color_scheme
    }else{
      color_scheme=NULL
    }
  }
  
  # Optionally overwrite the cluster names with user-supplied labels
  if(!is.null(cluster_labels)){
    names(samples_clustered)=cluster_labels
  }
  
  # Convert the per-cluster list of samples into a long data frame (sample, cluster).
  # cluster is a factor with levels ordered as in samples_clustered.
  cluster_df <- do.call(rbind, lapply(names(samples_clustered), function(cl) {
    data.frame(
      sample  = samples_clustered[[cl]],
      cluster = factor(cl,levels = names(samples_clustered))
    )
  }))
  # Expose the rownames as an explicit ID column so they can be used as a merge key
  clinical_tumors_df$Skyline_ID = rownames(clinical_tumors_df)
  
  # Join cluster assignments with the survival-relevant clinical columns
  surv_df <- merge(
    cluster_df,
    clinical_tumors_df[,c("Skyline_ID","Metastasis","FU_metastasis_years")],
    by.x = "sample",
    by.y = "Skyline_ID"
  )
  # Encode the event: 1 = metastasis ("Case"), 0 = censored
  surv_df$Metastasis_numeric = ifelse(surv_df$Metastasis=="Case",1,0)
  
  # Optionally apply an administrative censoring horizon at max_time:
  # cap follow-up time and treat any event occurring beyond max_time as censored
  if(!is.null(max_time)){
    surv_df$Metastasis_numeric <- ifelse(surv_df$FU_metastasis_years > max_time, 0, surv_df$Metastasis_numeric)
    surv_df$FU_metastasis_years <- pmin(surv_df$FU_metastasis_years, max_time)
  }
  # Build the survival object (time-to-event with censoring indicator)
  surv_obj <- Surv(
    time  = surv_df$FU_metastasis_years,
    event = surv_df$Metastasis_numeric
  )
  # Fit Kaplan-Meier survival curves stratified by cluster
  fit <- survminer::surv_fit(surv_obj ~ cluster, data = surv_df)
  
  # Build the survival plot. Two branches only differ in how legend labels are set.
  if(is.null(cluster_labels)){
    g_plot <-ggsurvplot(
      fit,
      data        = surv_df,
      risk.table  = TRUE,
      pval        = TRUE,
      conf.int    = FALSE,
      xlab        = "Years",
      ylab        = "Metastasis-free survival",
      legend.title = "Cluster",
      legend.labs = levels(surv_df$cluster),
      palette = color_scheme,
      title=title_plot
    )
  }else{
    g_plot <-ggsurvplot(
      fit,
      data        = surv_df,
      risk.table  = TRUE,
      pval        = TRUE,
      conf.int    = FALSE,
      xlab        = "Years",
      ylab        = "Metastasis-free survival",
      legend.labs = cluster_labels, 
      palette = color_scheme,
      title=title_plot
    )
  }
  # Return the ggsurvplot object (plot + risk table)
  return(g_plot)
}
