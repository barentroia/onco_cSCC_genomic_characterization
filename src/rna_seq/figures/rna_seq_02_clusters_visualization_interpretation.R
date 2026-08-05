#-------------------------------------------------------------------------------
# Aim: Visualize clusters found on RNA-seq data
# Author: B. Rentroia Pacheco
# Input: Cluster solution found in rna_seq_01
# Output: Several clustering solutions performed by the cola R package
#-------------------------------------------------------------------------------

# Setup - Load libraries & define output directory
#-------------------------------------------------------------------------------
#library(cola)
#library(openxlsx)
#library(plyr)
#library(dplyr)
#library(DESeq2)
#library(clusterProfiler)
#library(msigdbr)
#library(tibble)
#library(rtracklayer)
#library(fgsea)
#library(biomaRt)
#library(pheatmap)
#library(ggplot2)
#library(ggalluvial)
#library(officer)
#library(patchwork)
#library(survival)
#library(survminer)
#library(gtsummary)

# If running this script by itself, run config.R to define file paths
# Load useful functions:

# Output directory:
output_dir <- file.path(results_dir,"rna_seq","clustering_cola")

# Cluster discovery step:
source(file.path(code_dir,"functions","rnaseq_cluster_interpretation.R"))
source(file.path(code_dir,"rna_seq","analyses","rna_seq_01_clustering_analyses.R"))

# Use data from WES:
genomic_df = read.xlsx("/home/cscc_project_230717/brentroia/data/S04_02_03_genomic_df_summary.xlsx")

# Fetch final cola solution and also important intermediate data for clustering:
#-------------------------------------------------------------------------------
output_folder_cola = file.path(output_dir,"scaled","VST","adj","Excision")
rl = readRDS(file = file.path(output_folder_cola,"/rds_file.RDS"))
res = rl["MAD", "skmeans"]
chosen_k_cola = 3

# Fetch corresponding input matrix:
vt_mat_model_input = readRDS(file =  file.path(output_dir,"scaled","VST","VST_input.RDS"))

# Preprocess output from cola and add sample cluster names:
#-------------------------------------------------------------------------------
# Convert cola output into dataframe:
df_cola_cl = data.frame(Sample_id = rownames(res@object_list[[as.character(chosen_k_cola)]]$class_df),
                        Cluster=res@object_list[[as.character(chosen_k_cola)]]$class_df$class)
rownames(df_cola_cl) = df_cola_cl$Sample_id

# Add sample cluster names:
df_clusters<- df_cola_cl %>%
  mutate(Cluster_name = ifelse(Cluster==3,"Basal-like",ifelse(Cluster==1,"Differentiated","Mesenchymal-like")))%>%
  mutate(Cluster_name = factor(Cluster_name, levels = c(
    "Differentiated",
    "Basal-like",
    "Mesenchymal-like"
  )))

# Preprocess clinical and rnaseq datasets: 
#-------------------------------------------------------------------------------
# This can be omitted if running this script after rna_seq_01_clustering_analyses.R
source(file.path(code_dir,"rna_seq","analyses","rna_seq_a1_preprocessing.R"))
#-------------------------------------------------------------------------------

# Cluster annotations:
#-------------------------------------------------------------------------------
# Process the annotations defined in rna_seq_a1_preprocessing.R script:
# Select rows and columns of interest:
annotation_cols <- annotation_df[colnames(vt_mat_model_input), c("Metastasis","Differentiation","Biopsy_excision_2f","Tumor_location","Tumor_diameter","Morphology_subtype","PNI_or_LVI","Tissue_involvement","Tumor_purity","Sex","Age","Depth_of_Invasion","BWH","IS.at.CSCC")]
# Check rownames match the correlation matrix:
stopifnot(rownames(annotation_cols)==colnames(vt_mat_model_input))

# Add the sample clusters to the annotation dataframe:
annotation_cols$Cluster_cola = df_clusters[rownames(annotation_cols),"Cluster_name"]
annotation_cols$Cluster_cola <- factor(
  annotation_cols$Cluster_cola,
  levels = c(
    "Differentiated",
    "Basal-like",
    "Mesenchymal-like"
  ),
  ordered = TRUE
)
annotation_cols$Differentiation = ifelse(is.na(annotation_cols$Differentiation),"NA",as.character(annotation_cols$Differentiation))

# Stratifications by sample type:
#-------------------------------------------------------------------------------
# We want to check the clusters in biopsies and excisions separately, so we need to stratify the input matrix by sample type:
inp_mat = as.matrix(vt_mat_model_input)
inp_mat_adjusted = t(apply(inp_mat ,1,adjust_outlier)) # Note: we dont use adjust_matrix because we dont need to remove columns, it's just for visualization

i.mat_type_exc = which(annotation_cols$Biopsy_excision_2f=="Excision")
inp_mat_exc = inp_mat[,i.mat_type_exc]
inp_mat_exc_adjusted = adjust_matrix(inp_mat_exc) # Note: we dont use adjust_matrix because we dont need to remove columns, it's just for visualization

i.mat_type_biop = which(annotation_cols$Biopsy_excision_2f=="Biopsy")
inp_mat_biop = inp_mat[,i.mat_type_biop]
inp_mat_biop_adjusted = t(apply(inp_mat_biop,1,adjust_outlier)) # Note: we dont use adjust_matrix because we dont need to remove columns, it's just for visualization

# Visualize clustering solution:
#-------------------------------------------------------------------------------
# Cluster visualized with cola:
top_n = 1000
sign_res = get_signatures(res, k = chosen_k_cola,,anno = annotation_cols[i.mat_type_exc,],anno_col = annotation_colors,row_km = 7,top_signatures=top_n)
stopifnot(dim(inp_mat_exc_adjusted)[1]==length(res@row_index))
sign_res$Gene = rownames(inp_mat_exc_adjusted[sign_res$which_row,])

# Recreate heatmap for publication:

## Gene visualization (rows):
### Select genes used to visualize clustering:
heatmap_df = inp_mat_exc_adjusted[sign_res$Gene[order(sign_res$fdr,decreasing=FALSE)[1:top_n]],]
sign_res = sign_res[sign_res$Gene%in%rownames(heatmap_df),]
### Change the order of the genes in the heatmap, to make the row clusters more evident:
df_genes_heatmap = data.frame(Gene =sign_res$Gene,Cluster=sign_res$km)
desired_order <- c(7,2,3,1,4,5,6)
df_genes_heatmap$Clusters_reordered <- factor(df_genes_heatmap$Cluster,
                                              levels = desired_order)
rownames(df_genes_heatmap) = df_genes_heatmap$Gene

###  Compute where each cluster ends, so that the row clusters are visualized with breaks:
cluster_ends <- cumsum(table(df_genes_heatmap$Clusters_reordered[order(df_genes_heatmap$Clusters_reordered)]))
###  Remove separation of the last two clusters (no need for a gap after the last):
gaps_row <- cluster_ends[-length(cluster_ends)]
gaps_row <- gaps_row[-length(gaps_row)] # We merge the last two clusters into one, since the annotation is the same: cell-cycle proliferation
### Reorder heatmap rows:
heatmap_df =heatmap_df[df_genes_heatmap$Gene[order(df_genes_heatmap$Clusters_reordered)],intersect(rownames(annotation_cols)[order(annotation_cols$Cluster_cola)],colnames(heatmap_df))]

### Manually define colors for each row cluster, to make it easier to identify them:
cluster_colors <- c(
  "1" = "#1f77b4",
  "2" = "#ff7f0e",
  "3" = "#2ca02c",
  "4" = "#d62728",
  "5" = "#9467bd",
  "6" = "#8c564b",
  "7" = "#e377c2",
  "8" = "#7f7f7f"
)

annotation_colors_row <- list(Cluster = cluster_colors)
df_genes_heatmap$Cluster= as.character(df_genes_heatmap[,c("Cluster")])
row_annotations_hm = data.frame("Cluster"=df_genes_heatmap[,"Cluster"],row.names =rownames(df_genes_heatmap),Gene_names = rownames(df_genes_heatmap))

## Sample visualization (columns):
### Cluster samples within each cola cluster:
col_order <- unlist(
  lapply(levels(annotation_cols$Cluster_cola), function(clust) {
    # samples in this cluster
    samples_in_clust <- colnames(heatmap_df)[annotation_cols[colnames(heatmap_df), "Cluster_cola"] == clust]
    
    if(length(samples_in_clust) > 1) {
      # compute correlation between samples
      cor_mat <- cor(heatmap_df[, samples_in_clust], use = "pairwise.complete.obs",method="spearman")
      # hierarchical clustering using correlation distance
      hc <- hclust(as.dist(1 - cor_mat), method = "average")
      samples_in_clust[hc$order]
    } else {
      samples_in_clust
    }
  })
)
heatmap_df <- heatmap_df[, col_order]

## Visualization of genomic variables:
### Add 2 important genomic variables, that are significantly associated with the clusters: 
annotation_cols=annotation_cols %>% rownames_to_column(var = "SkylineDx.ID")%>%left_join(
  genomic_df[, c("SkylineDx.ID", "Notch", "pct_altered")],
  by = "SkylineDx.ID"
)

### Fraction of genome altered legend
# Categorizing Fraction of genome altered because pheatmap does not accept NAs in continuous variables:
annotation_cols$pct_altered_cat <- cut(
  annotation_cols$pct_altered,
  breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.000001),
  labels = c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1"),
  include.lowest = TRUE,
  right = FALSE
)
# Color legend for the fraction of the genome altered:
annotation_colors$pct_altered_cat <- c(
  "NA"       = "gray99",  
  "0–0.2"   = "#c6dbef",
  "0.2–0.4" = "#9ecae1",
  "0.4–0.6" = "#4292c6",
  "0.6–0.8" = "#2171b5",
  "0.8–1"   = "#084594"
)
annotation_cols$pct_altered_cat = as.character(annotation_cols$pct_altered_cat)
annotation_cols$pct_altered_cat[is.na(annotation_cols$pct_altered)]="NA"

### Notch alterations:
annotation_cols$Notch_cat = as.character(annotation_cols$Notch)
annotation_cols$Notch_cat[is.na(annotation_cols$Notch)]="NA"
annotation_colors$Notch_cat <- c(
  "NA"       = "gray99",  
  "0"   = "#fff7bc",
  "1" = "#fe9929",
  "2" = "#ec7014",
  "3" = "#cc4c02",
  "4"   = "#993404"
)

### Rownames:
rownames(annotation_cols) = annotation_cols$SkylineDx.ID

### Scale for gene expression values:
bk <- seq(-2, 2, length.out = 101)

### Plot heatmap for excisions:
p_excisions = pheatmap(heatmap_df,scale = "row",
                       cluster_rows = FALSE,show_rownames = FALSE,show_colnames = FALSE,cluster_cols  = FALSE,
                       annotation_col = annotation_cols[col_order,c("Metastasis","Differentiation","Notch_cat","pct_altered_cat","IS.at.CSCC","Cluster_cola")],
                       annotation_colors = c(annotation_colors,annotation_colors_row),
                       gaps_row = gaps_row,
                       annotation_row = row_annotations_hm[rownames(heatmap_df),"Cluster",drop=FALSE],
                       breaks=bk,
                       na_col="white")

### Save heatmap for publication (Fig 4A):
### Image: 
pdf(file.path(output_dir,"/rna_seq_02_heatmap_excisions.pdf"),width=8,height=10)
print(p_excisions)
dev.off()

### Source data:
write.xlsx(heatmap_df,file.path(output_dir,"/rna_seq_02_heatmap_excisions.xlsx"))

# Pathway enrichment analysis:
#-------------------------------------------------------------------------------
## Before performing enrichment analysis, we need to define the universe genes:
ensembl_genes <-rownames(vt_mat_model_input) 

# Convert Ensembl to gene symbols:
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
out <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","entrezgene_id"), 
             filters = "ensembl_gene_id",
             values = do.call(rbind,strsplit(ensembl_genes,"[.]"))[,1],
             mart = mart)

gene_symbols <- out[match(do.call(rbind,strsplit(rownames(inp_mat),"[.]"))[,1],out$ensembl_gene_id),]$external_gene_name
gene_symbols <-gene_symbols[which(gene_symbols!="")]
universe_genes = gene_symbols

# Perform pathway enrichment analysis:
for(annot in c("H","H_C2_C7_C8","H_C7_C8","C7","C6")){
  ORA_gene_clusters(annot,sign_res,output_folder_cola,universe_genes)
}

# Survival curves:
#-------------------------------------------------------------------------------
# On the excisions first:
stopifnot(all(rownames(annotation_cols)==rownames(clinical)))
new_samples_by_cluster_excisions <- split(
  rownames(df_clusters),
  df_clusters$Cluster_name
)

p_surv_exc = generate_surv_curv_hm_clusters(NULL,NULL,clinical, i.subset=i.mat_type_exc, cluster_labels=NULL,def_sample_clusters=new_samples_by_cluster_excisions,color_scheme=annotation_colors[["Cluster_cola"]],max_time=5,title_plot="Excisions")


# Now on biopsies. 
# Goal: build a PCA space from the excision (training) samples using the top
# discriminating genes, then project the biopsy (new) samples into that same space.

# First have cluster information in a vector:
cluster_vec <- rep(NA, ncol(inp_mat_exc))
names(cluster_vec) <- colnames(inp_mat_exc)
# Assign a cluster label to each sample
for (i in 1:nrow(df_clusters)) {
  cluster_vec[df_clusters$Sample_id[i]] <- as.character(df_clusters$Cluster_name)[i]
}

## PCA on training data (excisions), using only the chosen top genes from the excisions:
expr_t <- t(as.matrix(inp_mat_exc[sign_res$Gene[order(sign_res$fdr,decreasing=FALSE)[1:top_n]],])) #sign_res$fdr<0.00000001
pca <- prcomp(expr_t,scale. = TRUE,center = TRUE)

# Separate the two sample types:
df_train_sub <-  as.matrix(vt_mat_model_input)[colnames(expr_t),which(annotation_cols$Biopsy_excision_2f=="Excision")]
df_new_sub   <-  as.matrix(vt_mat_model_input)[colnames(expr_t),which(annotation_cols$Biopsy_excision_2f=="Biopsy")]

# Extract the PCA loadings (gene -> principal component weights) learned on excisions.
rotation <- pca$rotation

# transpose + scale biopsy samples:
expr_new_t <- t(as.matrix(df_new_sub))
expr_new_t <-scale(expr_new_t)

# Project biopsies into the excision space:
pca_new_scores <- expr_new_t %*% rotation

# Assemble a tidy data frame of the excision (training) PCA scores: one row per
# sample with its first six principal component coordinates, its cluster label
# (looked up by sample ID), and a type tag marking these as excisions.
pca_old_df <- data.frame(
  sample  = rownames(pca$x),
  PC1     = pca$x[,1],
  PC2     = pca$x[,2],
  PC3     = pca$x[,3],
  PC4    = pca$x[,4],
  PC5    = pca$x[,5],
  PC6    = pca$x[,6],
  cluster = cluster_vec[rownames(pca$x)],
  type    = "Excision"
)
# Same structure for the projected biopsies. These samples have no cluster
# assignment yet, so cluster is set to "New" and type marks them as biopsies.
pca_new_df <- data.frame(
  sample  = rownames(pca_new_scores),
  PC1     = pca_new_scores[,1],
  PC2     = pca_new_scores[,2],
  PC3     = pca_new_scores[,3],
  PC4     = pca_new_scores[,4],
  PC5     = pca_new_scores[,5],
  PC6     = pca_new_scores[,6],
  cluster = "New",
  type    = "Biopsy"
)

# Stack excisions and biopsies into a single data frame for combined plotting.
pca_all <- rbind(pca_old_df, pca_new_df)

# Plot projection:
# Plot 1: excisions only, coloured by their cluster (PC1 vs PC2).
pca_excisions = ggplot(pca_old_df, aes(PC1, PC2, color = cluster)) +
  geom_point(size = 4, alpha = 1)  +
  scale_color_manual(values = annotation_colors$Cluster_cola,
                     breaks = c("Differentiated",
                                "Basal-like",
                                "Mesenchymal-like"))+
  theme_classic(base_size = 16)+ theme(legend.position = "bottom")+  ylim(c(min(pca_all$PC2)-1),max(pca_all$PC2)+1)+
  xlim(c(min(pca_all$PC1)-1),max(pca_all$PC1)+1)

# Plot 2: excisions + biopsies together (before biopsies are assigned to clusters).
# Colour still encodes cluster ("New" for biopsies); point shape encodes sample type.
pca_bipsies_excisions_before_ass =ggplot(pca_all, aes(PC1, PC2, color = cluster, shape = type)) +
  geom_point(size = 4, alpha = 1)  +
  scale_color_manual(values = annotation_colors$Cluster_cola,
                     breaks = c("Differentiated",
                                "Basal-like",
                                "Mesenchymal-like"))+
  theme_classic(base_size = 16)+ theme(legend.position = "bottom")+
  scale_shape_manual(
    values = c(
      "Excision" = 16,
      "Biopsy"      = 17
    )
  )+  ylim(c(min(pca_all$PC2)-1),max(pca_all$PC2)+1)+
  xlim(c(min(pca_all$PC1)-1),max(pca_all$PC1)+1)

# Compute the centroids:
# Work in the first three PCs only. Group the excision samples by their cluster
# and average each PC within a cluster, giving one centroid coordinate per cluster.
pcs <- c("PC1", "PC2","PC3")
centroids_df <- pca_old_df %>%
  group_by(cluster) %>%
  summarise(across(all_of(pcs), mean), .groups = "drop")

# Assign clusters based on euclidean distance in the excision cohort:
# For a single sample's PC coordinates (x), compute its distance to every centroid
# and return the name of the closest cluster.
assign_cluster <- function(x, centroids,pcs,dist) {
  if(dist=="euclidean"){
    d <- apply(centroids[,pcs], 1, function(c) {
      sum((x - c)^2)
    })
    
  }
  centroids$cluster[which.min(d)]
  
}

# Predict clusters based on the nearest centroid:
# Apply the assignment to every excision sample (row-wise over its PC coordinates).
# This is essentially a self-consistency check: re-classifying the training samples
# by nearest centroid and comparing against their original cluster.
pred_euclid = apply(pca_old_df[,pcs],1,assign_cluster,centroids = centroids_df,pcs=pcs,dist="euclidean")
# Store predictions back on the data frame, indexed by sample name to keep alignment.
pca_old_df$pred_euclid <- pred_euclid[pca_old_df$sample]

#Check performance: 95% good!
mean(pca_old_df$pred_euclid == pca_old_df$cluster)

# Apply nearest centroid classifier to the biopsy samples:
# Same row-wise nearest-centroid assignment as for excisions, but now on the
# projected biopsy scores
pca_new_df$predicted_cluster <- apply(
  pca_new_df[, pcs],
  1,
  assign_cluster,
  centroids = centroids_df,
  pcs=pcs,
  dist="euclidean"
)

# Check how biopsy samples are classified:
# Convert the predicted labels to a factor with a fixed level order
pca_new_df$cluster <- factor(pca_new_df$predicted_cluster,levels=c("Differentiated",
                                                                   "Basal-like",
                                                                   "Mesenchymal-like"))
pca_new_df$predicted_cluster = NULL

# Rebuild the combined data frame, this time with biopsies carrying their PREDICTED
# cluster (previously they were all tagged "New"). Keep only cluster + PCs and re-tag
# sample type for each cohort.
pca_all2 <- rbind(
  cbind(pca_old_df[,c("cluster",pcs)],type="Excision"),
  cbind(pca_new_df[,c("cluster",pcs)],type="Biopsy")
)
# Plot excisions + biopsies coloured by (assigned) cluster, shape by sample type,
# so you can visually check that projected biopsies land near the expected clusters.
pca_bipsies_excisions_after_ass =ggplot(pca_all2, aes(PC1, PC2, color = cluster, shape = type)) +
  geom_point(size = 4, alpha = 1) +
  theme_classic() +
  scale_color_manual(values = annotation_colors$Cluster_cola,
                     breaks = c("Differentiated",
                                "Basal-like",
                                "Mesenchymal-like")) +
  labs(color = "Cluster", shape = "Sample type")+  
  theme_classic(base_size = 16)+ theme(legend.position = "bottom")+
  scale_shape_manual(
    values = c(
      "Excision" = 16,
      "Biopsy"      = 17
    )
  )+  ylim(c(min(pca_all$PC2)-1),max(pca_all$PC2)+1)+
  xlim(c(min(pca_all$PC1)-1),max(pca_all$PC1)+1)

# Find min/max across your data
x_limits <- range(pca_old_df$PC1)
y_limits <- range(pca_old_df$PC2)

# Biopsy-only version of the projection plot (excisions omitted), coloured by the
# cluster each biopsy was assigned to.
pca_bipsies_after_ass =ggplot(pca_new_df, aes(PC1, PC2, color = cluster, shape = type)) +
  geom_point(size = 4, alpha = 1) +
  theme_classic() +
  scale_color_manual(values = annotation_colors$Cluster_cola,
                     breaks = c("Differentiated",
                                "Basal-like",
                                "Mesenchymal-like")) +
  labs(color = "Cluster", shape = "Sample type")+  
  theme_classic(base_size = 16)+ theme(legend.position = "bottom")+
  scale_shape_manual(
    values = c(
      "Excision" = 16,
      "Biopsy"      = 17
    )
  )+  ylim(c(min(pca_all$PC2)-1),max(pca_all$PC2)+1)+
  xlim(c(min(pca_all$PC1)-1),max(pca_all$PC1)+1)

# Save transformations
# Combine in one row
combined <- pca_excisions+ pca_bipsies_excisions_before_ass + pca_bipsies_after_ass + plot_layout(nrow = 1)

# Save figure
ggsave(
  file.path(output_dir,"rna_seq_02_PCA_biopsies.pdf"),
  combined,
  width = 18,      # inches
  height = 6,      # inches (1 row)
  dpi = 600
)

# Group biopsy sample IDs by their assigned cluster -> named list (cluster -> samples),
# in the shape expected by generate_surv_curv_hm_clusters(def_sample_clusters=...).
new_samples_by_cluster <- split(
  pca_new_df$sample,
  pca_new_df$cluster
)


# Check survival on the biopsies:
i.biop = which(annotation_cols$Biopsy_excision_2f=="Biopsy")
p_surv_biop=generate_surv_curv_hm_clusters(NULL,NULL,clinical, i.subset=i.biop , cluster_labels=NULL,def_sample_clusters=new_samples_by_cluster,color_scheme=annotation_colors[["Cluster_cola"]],max_time=5,title_plot="Biopsies")

# Check cluster assignment:
# Add sample IDs as an explicit column and export the combined PCA + cluster table.
pca_all2$Sample_id = rownames(pca_all2)
write.xlsx(pca_all2,file.path(output_dir,"rna_seq_02_PCA_biopsy_excision.xlsx"))

# Save two curves together:
p_combined <- arrange_ggsurvplots(list(p_surv_biop, p_surv_exc),ncol=2,nrow=1)
ggsave(file.path(output_dir,"rna_seq_02_scurve_biopsy_vs_excision.pdf"),p_combined, width = 15,height =8)

# Save source data:
write.xlsx(p_surv_biop$data.survplot,file.path(output_dir,"rna_seq_02_scurve_biopsy_source_data.xlsx"))
write.xlsx(p_surv_biop$data.survtable,file.path(output_dir,"rna_seq_02_stable_biopsy_source_data.xlsx"))
write.xlsx(p_surv_exc$data.survplot,file.path(output_dir,"rna_seq_02_scurve_excision_source_data.xlsx"))
write.xlsx(p_surv_exc$data.survtable,file.path(output_dir,"rna_seq_02_stable_excision_source_data.xlsx"))

# Check heatmap of biopsies:
# Fill in annotations:
# For each predicted biopsy cluster, write that cluster name into the shared
# annotation table — but only for samples whose Cluster_cola is still NA, so
# existing (excision) assignments are never overwritten.
for (cl in names(new_samples_by_cluster)) {
  samples <- new_samples_by_cluster[[cl]]
  
  to_fill <- samples[is.na(annotation_cols[samples, "Cluster_cola"])]
  
  annotation_cols[to_fill, "Cluster_cola"] <- cl
}

# Heatmap rows for biopsies:
# Subset the biopsy expression matrix to the same genes (rows) and row order as the
# excision heatmap, so both heatmaps are directly comparable.
heatmap_df_biopsies = inp_mat_biop_adjusted[rownames(heatmap_df),]

# Build the column (sample) ordering: group samples by cluster, and within each
# cluster order them by similarity with hierarchical clustering, like what it was done for the excision cohort.
col_order <- unlist(
  lapply(levels(annotation_cols$Cluster_cola), function(clust) {
    # samples in this cluster
    samples_in_clust <- colnames(heatmap_df_biopsies)[annotation_cols[colnames(heatmap_df_biopsies), "Cluster_cola"] == clust]
    
    if(length(samples_in_clust) > 1) {
      # compute correlation between samples
      cor_mat <- cor(heatmap_df_biopsies[, samples_in_clust], use = "pairwise.complete.obs",method="spearman")
      # hierarchical clustering using correlation distance
      hc <- hclust(as.dist(1 - cor_mat), method = "average")  # you can change method
      samples_in_clust[hc$order]
    } else {
      samples_in_clust
    }
  })
)
stopifnot(length(col_order) == ncol(heatmap_df_biopsies))
# Draw the biopsy heatmap: columns in the cluster-grouped order computed above,
# rows fixed (no row clustering) to match the excision heatmap layout.
p_biopsies = pheatmap(heatmap_df_biopsies[,col_order] ,scale = "row",
                      cluster_rows = FALSE,show_rownames = FALSE,show_colnames = FALSE,cluster_cols  = FALSE,
                      annotation_col = annotation_cols[col_order,c("Metastasis","Differentiation","IS.at.CSCC","Cluster_cola")],
                      annotation_colors = c(annotation_colors,annotation_colors_row),
                      gaps_row = gaps_row,
                      annotation_row = as.data.frame(df_genes_heatmap[,"Cluster"],row.names =rownames(df_genes_heatmap)),
                      breaks=bk)

pdf(file.path(output_dir,"rna_seq_02_Heatmap_biopsies.pdf"),width=8,height=10)
print(p_biopsies)
dev.off()

png(
  filename =file.path(output_dir,"rna_seq_02_Heatmap_biopsies.png"),
  width  = 8000,   # pixels
  height = 10000,   # pixels
  res    = 600      # DPI (important for print)
)
print(p_biopsies)
dev.off()
# Export the plotted matrix (in display column order) as source data.
write.xlsx(heatmap_df_biopsies[,col_order],file.path(output_dir,"rna_seq_02_Heatmap_biopsies.xlsx"))

# Repeat the same (heatmap + survival curves) but for all samples:
heatmap_df_all= inp_mat_adjusted[rownames(heatmap_df),]
col_order <- unlist(
  lapply(levels(annotation_cols$Cluster_cola), function(clust) {
    # samples in this cluster
    samples_in_clust <- colnames(heatmap_df_all)[annotation_cols[colnames(heatmap_df_all), "Cluster_cola"] == clust]
    
    if(length(samples_in_clust) > 1) {
      # compute correlation between samples
      cor_mat <- cor(heatmap_df_all[, samples_in_clust], use = "pairwise.complete.obs",method="spearman")
      # hierarchical clustering using correlation distance
      hc <- hclust(as.dist(1 - cor_mat), method = "average")  # you can change method
      samples_in_clust[hc$order]
    } else {
      samples_in_clust
    }
  })
)

pheatmap(heatmap_df_all[,col_order] ,scale = "row",
         cluster_rows = FALSE,show_rownames = FALSE,show_colnames = FALSE,cluster_cols  = FALSE,
         annotation_col = annotation_cols[,c("Metastasis","Differentiation","Cluster_cola","Biopsy_excision_2f")],
         annotation_colors = c(annotation_colors,annotation_colors_row),
         gaps_row = gaps_row,
         annotation_row = as.data.frame(df_genes_heatmap[,"Cluster"],row.names =rownames(df_genes_heatmap)),
         breaks=bk)


new_samples_by_cluster_all =  Map(c, new_samples_by_cluster_excisions, new_samples_by_cluster)
generate_surv_curv_hm_clusters(NULL,NULL,clinical, i.subset=1:nrow(clinical), cluster_labels=NULL,def_sample_clusters=new_samples_by_cluster_all,color_scheme=annotation_colors[["Cluster_cola"]],max_time=5,title_plot="All")

# Association tests between clinico-pathological variables and the clusters:
#-------------------------------------------------------------------------------
# Tidy up clinical data:
association_df_not_imputed = annotation_df%>%
  dplyr::select(Metastasis, Sex,Age,IS.at.CSCC, Tumor_location, BWH,AJCC_8, Number_of_cSCC_before_culprit,Tumor_diameter,Tissue_involvement,Differentiation,PNI_or_LVI,Depth_of_Invasion,Breslow_thickness,Peritumoral_infiltration,Solar_elastosis,Tumor_budding, Mitotic_rate,Morphology_subtype)%>%
  dplyr::mutate(Metastasis = factor(Metastasis,levels=c("Control","Case")), 
                Sex = factor(Sex,levels=c("Female","Male")),
                Tumor_location=factor(Tumor_location,levels=c("Face","Scalp/neck","Trunk/Extremities")),
                BWH = factor(BWH,levels=c("T1","T2a","T2b","T3")), 
                AJCC_8 = factor(AJCC_8,levels=c("T1","T2","T3","T4")),
                Tissue_involvement=factor(Tissue_involvement,levels=c("Dermis","Subcutaneous fat","Beyond subcutaneous fat")),
                PNI_or_LVI = factor(PNI_or_LVI,levels=c("No","Yes")),
                Peritumoral_infiltration=factor(Peritumoral_infiltration,levels=c("Absent/moderate","Abundant")),
                Solar_elastosis=factor(Solar_elastosis,levels=c("Absent/moderate","Extensive")),
                Morphology_subtype = factor(Morphology_subtype,levels=c("Not otherwise specified",setdiff(names(table(annotation_df$Morphology_subtype)),"Not otherwise specified"))))
# Add genomic clusters to the clinical dataset:
association_df_not_imputed = merge(association_df_not_imputed, pca_all2[,c("cluster","Sample_id")],by = "row.names")
association_df_not_imputed = association_df_not_imputed%>%
  dplyr::select(-c(Sample_id,Row.names))%>%
  dplyr::mutate(cluster = factor(cluster,levels=c("Differentiated", "Basal-like", "Mesenchymal-like")))

#Perform tests:
tbl <- association_df_not_imputed %>%
  tbl_summary(
    by = cluster,
    missing = "ifany",
    statistic = list(
      all_continuous() ~ "{median} ({IQR})",
      all_categorical() ~ "{n} ({p}%)"
    )
  ) %>%
  add_p(test.args = list(
    Morphology_subtype ~ list(simulate.p.value = TRUE, B = 10000) # Morphology subtype has too many categories
  ))

# Save tests:
tbl %>%  as_tibble( col_labels = TRUE)%>%
  write.xlsx(file.path(output_dir,"rna_seq_02_Association_with_gclusters_notimputed.xlsx"))
