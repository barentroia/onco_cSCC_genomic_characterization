#===============================================================================
# Aim:      Preprocess the clinical data and the RNA-seq counts/TPMs of the
#           D-SQUAME discovery dataset for the clustering analyses: filter
#           samples and genes, prepare annotations and fit a PCA model.
# Author:   B. Rentroia Pacheco, based on J. Traets code
# Input:    - Clinical dataset (p3_clinical), bad QC list (p7_badQC)
#           - TPMs and counts (p1_tpms, p2_counts), gene biotypes (p12_genes_to_keep)
# Output:   Objects in the calling environment: clinical, counts,
#           counts_model_input, log2tpm_model_input, feats_to_exlude,
#           annotation_cols, annotation_colors. Files in results_dir/rna_seq:
#           - rna_seq_a1_PCA_discovery.csv (PC1/PC2 per sample)
#           - rna_seq_a1_pca_model_top1000.rds (PCA projection model)
# Requires: results_dir (main.R), config.R paths, openxlsx/dplyr/DESeq2 loaded
# Usage:    Sourced by rna_seq_01_clustering_analyses.R and
#           rna_seq_02_clusters_visualization_interpretation.R
#===============================================================================

# 1. Read data
#-------------------------------------------------------------------------------
clinical_samples <- read.xlsx(p3_clinical, sheet = "Request Lara - Samples (values)")
clinical_tumors <- read.xlsx(p3_clinical, sheet = "Request Lara - Tumor (values)")
bad_qc_samples <- read.csv(p7_badQC)
tpm <- read.csv(p1_tpms)
counts <- read.csv(p2_counts)
genes_to_keep <- readRDS(p12_genes_to_keep)
# Output directory:
output_dir <- file.path(results_dir,"rna_seq")
#-------------------------------------------------------------------------------

# 2. Data processing
#-------------------------------------------------------------------------------
# 2.1 Clinical dataset (bad QC samples removed):
# PNI_bin column is repeated, drop one (column 28)
clinical_tumors <- clinical_tumors[, -28]
clinical <- clinical_samples %>%
  dplyr::rename(Sample_ID = Sample_id) %>%
  dplyr::full_join(clinical_tumors,
                   suffix = c("_samplessheet", "_tumorssheet"),
                   by = c("Skyline_ID", "Patient_ID_SKY", "Sample_ID")) %>%
  dplyr::filter(!Skyline_ID %in% bad_qc_samples$Skyline_ID) %>%
  mutate(`Sample type` = ifelse(Biopsy_excision == "Biopsy", "Biopsy", "Excision/Curettage/Excochleation")) %>%
  tibble::column_to_rownames(var = "Skyline_ID")
#-------------------------------------------------------------------------------

# 2.2 Counts
counts <- counts %>%
  tibble::column_to_rownames(var = "ID") %>%
  dplyr::select(all_of(rownames(clinical)))
#-------------------------------------------------------------------------------

# 2.3 log2(TPM + 1)
log2tpm <- log2(tpm %>%
                  tibble::column_to_rownames(var = "ID") %>%
                  dplyr::select(all_of(rownames(clinical))) + 1)

# 2.4 Gene filtering:
## Remove features with low expression and keep protein coding + lncRNA
## Zero variance
feats_zero_var <- rownames(log2tpm)[apply(log2tpm, 1, stats::var) == 0]
# Lowly expressed in half of the cohort
feats_high_exp_half_cohort <- rownames(log2tpm)[apply(log2tpm, 1, function(x) sum(x > 1)) > (ncol(log2tpm) * 0.5)]
feats_low_exp_half_cohort <- base::setdiff(rownames(log2tpm), feats_high_exp_half_cohort)
## Low average expression
feats_low_mean_exp <- rownames(log2tpm)[apply(log2tpm, 1, function(x) mean(x) < 1)]
## Non protein coding/lncRNA
feats_not_protcoding_lncRNA <- rownames(log2tpm)[!gsub("\\..*", "", rownames(log2tpm)) %in% genes_to_keep$ensembl_gene_id]
## Exclude features:
feats_to_exlude <- unique(c(feats_zero_var, feats_low_exp_half_cohort, feats_low_mean_exp, feats_not_protcoding_lncRNA))
counts_model_input <- counts[which(!rownames(counts) %in% feats_to_exlude), ]
log2tpm_model_input <- log2tpm[which(!rownames(log2tpm) %in% feats_to_exlude), ]
#-------------------------------------------------------------------------------

# 2.5 Sanity check:
## Check that colnames of counts are in the same order of rownames in clinical
if (!all(rownames(clinical) == colnames(counts))){
  counts <- counts[, rownames(clinical)]
}
stopifnot(all(rownames(clinical) == colnames(counts)))
#-------------------------------------------------------------------------------

# 2.6 Clinical dataset in the form of annotation data:
## Prepare annotation data
### Match samples in vt_mat with clinical data:
annotation_df <- clinical[match(colnames(counts), rownames(clinical)), ]
###  Simplify clinical variables:
annotation_df = annotation_df%>%
  mutate(Biopsy_excision_2f = factor(ifelse(Biopsy_excision%in%c("Curettage","Excochleation"),"Excision",Biopsy_excision)),
         Tumor_purity = as.numeric(as.character(Tumor_purity)),
         IS.at.CSCC = factor(ifelse(OTR_at_cSCC=="Yes"|HM_at_cSCC=="Yes","Yes","No"),levels=c("No","Yes")),
         Differentiation = factor(Differentiation,levels=c("Good/moderate","Poor/undifferentiated")))
### Select specific columns to annotate:
annotation_cols <- annotation_df[, c("Metastasis","Differentiation","Biopsy_excision_2f","IS.at.CSCC")]
###  Set rownames to match the correlation matrix:
rownames(annotation_cols) <- colnames(counts)

# Annotation colors:
annotation_colors <- list(
  Metastasis = c("Case" = "red", "Control" = "blue"),
  Differentiation = c("Good/moderate"="yellow","Poor/undifferentiated"="brown","NA"="white"),
  Biopsy_excision_2f = c("Biopsy"="gray","Excision"="black"),
  PNI_or_LVI = c("Yes"="black","No"="gray"),
  Tumor_location = c("Face"="lightblue","Scalp/neck"="darkblue","Trunk/Extremities"="darkgreen"),
  Morphology_subtype = c("Acantholytic"="brown","Clear cell"="peru","dd Keratoacanthoma"="magenta2","dd Keratoacanthoom or folliculair"="magenta4","Desmoplastic"="black","Not classifiable"="gray","Not otherwise specified"="grey75","Spindle cell"="darkorange") ,
  Tissue_involvement = c("Dermis"= "wheat1","Subcutaneous fat"="tan1","Beyond subcutaneous fat"="tan4"),
  Sex=c("Male"="lightblue","Female"= "lightpink"),
  BWH=c("T1"="lightblue","T2a"="turquoise3","T2b"="dodgerblue3","T3"="darkblue"),
  Cluster_cola = c("Differentiated"="wheat2","Basal-like"="tan1","Mesenchymal-like"= "tan4"),
  IS.at.CSCC = c("Yes"="darkred" ,"No"="gray"),
  pct_altered = colorRampPalette(c("white", "darkorchid4"))(100),
  Notch = colorRampPalette(c("white", "goldenrod3"))(100))

# 2.7 PCA on the 1000 most variable genes (used to project other datasets):
top_pca = 1000
highest_variance <- sort(apply(log2tpm_model_input,1,function(x)var(x)),decreasing=TRUE)[1:top_pca]
expr_t <- t(as.matrix(log2tpm_model_input[names(highest_variance),])) 
pca <- prcomp(expr_t, scale = T)
# % variance explained per PC
pv <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 2)
df <-data.frame(PC1 = pca$x[,"PC1"],PC2 = pca$x[,"PC2"], Metastasis = clinical$Metastasis[match(colnames(log2tpm_model_input), rownames(clinical))])
# Data check: 
#df%>%ggplot(aes(x=PC1, y=PC2,colour =  Metastasis)) +
#  geom_point(size = 4, alpha = 1)
# Add variance to the column names:
write.csv(df[,c(1,2)],file.path(output_dir,"rna_seq_a1_PCA_discovery.csv"))

# Save minimal projection object:

#' Create a minimal PCA projection model
#'
#' Keeps only what is needed to project new samples with predict().
#'
#' @param pca prcomp object.
#' @param n_pc integer. Number of principal components to keep.
#'
#' @return Object of class "prcomp" with genes, center, scale, rot
make_pca_model <- function(pca, n_pc = 2) {
  n_pc <- min(n_pc, ncol(pca$rotation))
  structure(
    list(
      genes    = rownames(pca$rotation),           # genes in training order
      center   = pca$center,                       # training gene means
      scale    = pca$scale,                        # training gene SDs
      rotation = pca$rotation[, 1:n_pc, drop = FALSE],
      pct_var  = round(100 * pca$sdev^2 / sum(pca$sdev^2), 2)[1:n_pc]
    ),
    class = "prcomp"                               # lets predict() work
  )
}
pca_model <- make_pca_model(pca, n_pc = 2)
saveRDS(pca_model, file.path(output_dir,paste0("rna_seq_a1_pca_model_top",top_pca,".rds")))
#-------------------------------------------------------------------------------