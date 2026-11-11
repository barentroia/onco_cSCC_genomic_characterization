#===============================================================================
# Aim:      Central configuration: defines the file paths of all input data
#           used by the analysis and figure scripts.
# Requires: project_dir (defined in main.R)
# Usage:    Sourced from main.R. Paths are assigned with <<- so they are
#           available globally to every script sourced afterwards.
# Notes:    Prefixes indicate the dataset: p = D-SQUAME discovery,
#           e = D-SQUAME validation, n = Nassir et al.
#===============================================================================

# 1. D-SQUAME discovery dataset
#-------------------------------------------------------------------------------
# Gene-level TPMs and raw counts (Salmon, txi)
p1_tpms <<- file.path(project_dir,"fixedData","seqData","foxglove_v2.05","ABE_Salmon","txi.gene.tpm.csv")
p2_counts <<- file.path(project_dir,"fixedData","seqData","foxglove_v2.05","ABE_Salmon","txi.gene.counts.csv")
# Raw clinical dataset (Excel) and parsed clinical dataset (csv)
p3_clinical <<- file.path(project_dir, "lpozza", "data", "clinical" , "20231019_Clinical_Dataset_cSCC_final_v6.xlsx")
p4_clinical_parsed <<- file.path(project_dir, "lpozza", "step-ident-2023-eda-and-comparison-other-signatures", "output", "parsed_clinical_data", "parsed_data_clinical_discovery_v6.csv")
# Sample IDs (and pairs) with bad RNA-seq QC
p7_badQC <<- file.path(project_dir,"lpozza","data","qc","sample_ids_bad_QC.csv")
p8_badQC_pairs <<- file.path(project_dir, "lpozza", "data", "qc", "sample_ids_bad_QC_pairs.csv") 
# GENCODE v38 gene annotation
p9_gtf <<- file.path(project_dir,"fixedData","seqData","hawkweed_v2.05","GTF","gencode.v38.annotation.gtf.gz")
# Sampling weights of the nested case-control (NCC) design: 183 sets (main
# analyses) and 195 sets (unmatched integration analyses)
p10_weights <<- file.path(project_dir, "lpozza", "data", "weights" , "WC_01_rec_based_samp_prob_weights183_NCC_cohort_noADMIN.csv")
p11_weights_unmatched <- weights_ncc_fn <- file.path(project_dir, "lpozza", "data", "weights", paste0("WC_01_rec_based_samp_prob_weights", "195", "_NCC_cohort_noADMIN.csv"))
# Gene IDs with biotype "protein_coding" or "lncRNA" (retrieved from biomaRt on 2024-10-14 15:40:04 CEST)
p12_genes_to_keep <<- file.path(project_dir,"rruiter", "gep_publication_data/Genes_biotype.rds")
# WES mutation summary, curated genes/pathways for the tile plot, and the dataframe of
# good QC WES samples (integration analyses)
p13_mutation_data <<- file.path(project_dir, "jtraets/_OLD/gep_publication/data/S04_02_03_genomic_df_summary.xlsx")
p14_all_pathways <<- file.path(project_dir,"rruiter", "gep_publication_data", "2025_10_20_Tileplot_genes_and_pathways.xlsx")
p15_input138 <<- file.path(project_dir,"rruiter","gep_publication_data","genomic_df_input.csv") 
#-------------------------------------------------------------------------------

# 2. D-SQUAME validation dataset
#-------------------------------------------------------------------------------
# Output folder of the RNA-seq (Salmon) pipeline: gene counts and TPMs
dir_output_pipeline <<- file.path(project_dir, "jtraets", "cscc_validation_250214","seqData", "merged_output_set2","salmon")
e1_counts <<- file.path(dir_output_pipeline,"salmon.merged.gene_counts.all.tsv")
e2_tpms <<- file.path(dir_output_pipeline,"salmon.merged.gene_tpm.all.tsv")
# Parsed clinical data of the complete case-control sets
e3_clinical <<- file.path(project_dir,"lpozza","stepident-2025-validation-random-controls","output","clinical_data","parsed_clinical_data_v3lp_valset_complete_sets.csv")
# Sample IDs with bad RNA-seq QC and NCC sampling weights
e6_badQC <<- file.path(project_dir,"jtraets","cscc_validation_250214","BadQCs_set1.txt")
e7_weights <<- file.path(project_dir, "lpozza", "data", "weights_emc_validation", "WC_LR_NCC_cohort_valset_noADMIN.csv")
#-------------------------------------------------------------------------------

# 3. Nassir et al. dataset (external validation)
#-------------------------------------------------------------------------------
# Gene TPMs, gene counts and clinical metadata (Supplementary Table S3)
n1_tpms <<- file.path(project_dir,"jtraets","external_data","NassirBS_Mayo_2024","output_sky-rnaseq","merged_output","salmon","salmon.merged.gene_tpm.all.tsv")
n2_counts <<- file.path(project_dir,"jtraets","external_data","NassirBS_Mayo_2024","output_sky-rnaseq","merged_output","salmon","salmon.merged.gene_counts.all.tsv")
n3_clinical <<- file.path(project_dir, "jtraets","external_data","NassirBS_Mayo_2024","Metadata_supp_tableS3_SRA.csv")
#-------------------------------------------------------------------------------

# 4. Spatial sequencing dataset (Seurat object, 14 samples)
#-------------------------------------------------------------------------------
data_ST_dir <<- file.path(project_dir,"jtraets","external_data","Vignesh_spatial_2026")
seurat_object_fn <<- file.path(data_ST_dir,"/combined_seurat_obj_normalized_14Samples_Oct2026.RDS")
#-------------------------------------------------------------------------------

# 5. Single-cell dataset of Andrew L. Ji et al. 2020 (GEO: GSE144236)
#-------------------------------------------------------------------------------
data_SC_dir    <<- file.path(project_dir, "jtraets", "external_data", "AndrewLJi_spatial_2020", "data_GSE144236")
SC_counts_fn   <<- file.path(data_SC_dir, "GSE144236_cSCC_counts.txt")
SC_metadata_fn <<- file.path(data_SC_dir, "GSE144236_patient_metadata_new.txt")
#-------------------------------------------------------------------------------

# 6. Gene signatures
#-------------------------------------------------------------------------------
p6_bailey <<- file.path(project_dir,"jtraets","data","signatures","BaileyP_UK_DvP_SuppDataFile5.csv")
p7_bailey_xsls <<- file.path(project_dir,"jtraets","data","signatures","41467_2023_40822_MOESM8_ESM.xlsx")
#-------------------------------------------------------------------------------



