#===============================================================================
# Aim:      Unsupervised clustering of the RNA-seq data of the D-SQUAME
#           discovery dataset with multiple consensus clustering methods
#           (cola R package), for all samples and for biopsies/excisions
#           separately.
# Author:   B. Rentroia Pacheco
# Input:    - Clinical and RNA-seq count data, preprocessed by
#             rna_seq_a1_preprocessing.R
# Output:   - cola reports and result objects (rds_file.RDS) in
#             results_dir/rna_seq/clustering_cola/scaled/VST/adj/<subset>
# Requires: code_dir, results_dir (main.R) and config.R paths
# Usage:    Sourced by rna_seq_02_clusters_visualization_interpretation.R
#===============================================================================

# 1. Setup - Load libraries & define output directory
#-------------------------------------------------------------------------------
# Libraries (other packages are loaded by rna_seq_a1_preprocessing.R)
library(cola)
# If running this script by itself, run config.R to define file paths

# Output folder:
output_dir <- file.path(results_dir,"rna_seq","clustering_cola")
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}
#-------------------------------------------------------------------------------

# 2. Preprocess clinical and rnaseq datasets
#-------------------------------------------------------------------------------
source(file.path(code_dir,"rna_seq","analyses","rna_seq_a1_preprocessing.R"))
#-------------------------------------------------------------------------------

# 3. Generate input matrices for clustering and annotations
#-------------------------------------------------------------------------------
# Variance stabilizing transformation (VST) of the counts:
## Create dds object:
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = clinical,
                              design = ~1)
## Estimate size factor
dds <- estimateSizeFactors(dds)

## Variance stabilization
vt <- vst(dds, blind = T)
vt_mat <-assay(vt)
# Remove genes flagged during preprocessing
vt_mat_model_input <-vt_mat[which(!rownames(vt_mat) %in% feats_to_exlude), ]

## Include other run of hierarchical clustering with a different setting:
register_partition_methods(
  hclust_wd = function(mat, k) cutree(hclust(dist(t(mat), method = "euclidean"),method="ward.D2"), k)
)
#-------------------------------------------------------------------------------

# 4. Clustering using the cola package
#-------------------------------------------------------------------------------
# Settings to run cola:
scale_df = "scaled" # Scale input matrix
df_in = "VST" # Use VST transformation
adj_mat = TRUE # Use adjustment function from the cola package
dir.create(file.path(output_dir,scale_df,df_in),recursive = TRUE)

# Loop allows to apply cola on all data or just biopsies/excisions separately:
for(subset_smp in c("All","Excision","Biopsy")){
  print(subset_smp)
  
  # Save input matrix:
  inp_mat = as.matrix(vt_mat_model_input)
  saveRDS(inp_mat,file.path(output_dir,scale_df,df_in,"VST_input.RDS"))
  
  # Use the desired sample subset:
  if(subset_smp =="Excision"){
    i.mat_type = which(annotation_cols$Biopsy_excision_2f=="Excision")
  }else if(subset_smp =="Biopsy"){
    i.mat_type = which(annotation_cols$Biopsy_excision_2f=="Biopsy")
  }else{
    i.mat_type = 1:ncol(inp_mat)
  }

  inp_mat = inp_mat[,i.mat_type]
  
  # Cola adjustment to the input matrix:
  if(adj_mat){
    inp_mat = adjust_matrix(inp_mat)  # optional
    adj_txt ="adj"
  }else{
    adj_txt ="not_adj"
  }
  
  # Perform scaling or not:
  if(scale_df=="scaled"){
    scale_df_comm = TRUE
  }else if(scale_df=="not_scaled"){
    scale_df_comm = FALSE
  }
  
  # Create folder according to options chosen:
  output_folder_cola = file.path(output_dir,scale_df,df_in,adj_txt,subset_smp)
  dir.create(output_folder_cola, recursive = TRUE)
  
  # Run cola:
  print(dim(inp_mat))
  rl = run_all_consensus_partition_methods(inp_mat , cores = 4,anno =annotation_cols[i.mat_type,],anno_col = annotation_colors[colnames(annotation_cols)],scale_rows = scale_df_comm,top_value_method=c("SD","MAD","CV"))
  
  # Save report:
  cola_report(rl, output_dir = output_folder_cola, cores = 2)
  
  # Save rds:
  saveRDS(rl,file = file.path(output_folder_cola,"rds_file.RDS"))
}
      
#-------------------------------------------------------------------------------
