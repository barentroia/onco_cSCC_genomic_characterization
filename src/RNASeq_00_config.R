#-------------------------------------------------------------------------------
# Config file with all the paths
#-------------------------------------------------------------------------------


# D-SQUAME discovery dataset
#-------------------------------------------------------------------------------
p1_tpms <<- file.path(project_dir,"fixedData","seqData","foxglove_v2.05","ABE_Salmon","txi.gene.tpm.csv")
p2_counts <<- file.path(project_dir,"fixedData","seqData","foxglove_v2.05","ABE_Salmon","txi.gene.counts.csv")
p3_clinical <<- file.path(project_dir, "lpozza", "data", "clinical" , "20231019_Clinical_Dataset_cSCC_final_v6.xlsx")
p4_clinical_parsed <<- file.path(project_dir, "lpozza", "step-ident-2023-eda-and-comparison-other-signatures", "output", "parsed_clinical_data", "parsed_data_clinical_discovery_v6.csv")
# TODO fix paths to multiple imputation data sets?
p7_badQC <<- file.path(project_dir,"lpozza","data","qc","sample_ids_bad_QC.csv")
p8_badQC_pairs <<- file.path(project_dir, "lpozza", "data", "qc", "sample_ids_bad_QC_pairs.csv") 
p9_gtf <<- file.path(project_dir,"fixedData","seqData","hawkweed_v2.05","GTF","gencode.v38.annotation.gtf.gz")
p10_weights <<- file.path(project_dir, "lpozza", "data", "weights" , "WC_01_rec_based_samp_prob_weights183_NCC_cohort_noADMIN.csv")
#p11_tp_data <- file.path(project_dir, "ychen", "data", "results", "tumorPurity", "cSCC_discovery_tumor_purity_ESTIMATE.csv")
#p12_qcs_data <- file.path(project_dir, "fixedData", "QC", "cSCC_discovery_QC.csv")
p11_weights_unmatched <- weights_ncc_fn <- file.path(project_dir, "lpozza", "data", "weights", paste0("WC_01_rec_based_samp_prob_weights", "195", "_NCC_cohort_noADMIN.csv"))
# Genes ids with biotype "protein_coding" or "lncRNA": Retrieved from biomart on 2024-10-14 15:40:04 CEST
p12_genes_to_keep <<- file.path(project_dir,"rruiter", "gep_publication_data/Genes_biotype.rds")
p13_mutation_data <<- file.path(project_dir, "jtraets/gep_publication/data/S04_02_03_genomic_df_summary.xlsx")
p14_all_pathways <<- file.path(project_dir,"rruiter", "gep_publication_data", "2025_10_20_Tileplot_genes_and_pathways.xlsx")
p15_input138 <<- file.path(project_dir,"rruiter","gep_publication_data","genomic_df_input.csv") # to know the good QC WES samples, can be removed when wes is put in

#-------------------------------------------------------------------------------

# D-SQUAME validation dataset
#-------------------------------------------------------------------------------
dir_output_pipeline <<- file.path(project_dir, "jtraets", "cscc_validation_250214","seqData", "merged_output_set2","salmon")
e1_counts <<- file.path(dir_output_pipeline,"salmon.merged.gene_counts.all.tsv")
e2_tpms <<- file.path(dir_output_pipeline,"salmon.merged.gene_tpm.all.tsv")
e3_clinical <<- file.path(project_dir,"lpozza","stepident-2025-validation-random-controls","output","clinical_data","parsed_clinical_data_v3lp_valset_complete_sets.csv")
# TODO fix paths to multiple imputation data sets?
e6_badQC <<- file.path(project_dir,"jtraets","cscc_validation_250214","BadQCs_set1.txt")
e7_weights <<- file.path(project_dir, "lpozza", "data", "weights_emc_validation", "WC_LR_NCC_cohort_valset_noADMIN.csv")
#-------------------------------------------------------------------------------


# Nassir et al. dataset
#-------------------------------------------------------------------------------
n1_tpms <<- file.path(project_dir,"jtraets","external_data","NassirBS_Mayo_2024","output_sky-rnaseq","merged_output","salmon","salmon.merged.gene_tpm.all.tsv")
n2_counts <<- file.path(project_dir,"jtraets","external_data","NassirBS_Mayo_2024","output_sky-rnaseq","merged_output","salmon","salmon.merged.gene_counts.all.tsv")
n3_clinical <<- file.path(project_dir, "jtraets","external_data","NassirBS_Mayo_2024","Metadata_supp_tableS3_SRA.csv")
#-------------------------------------------------------------------------------


# Spatial sequencing data set
#-------------------------------------------------------------------------------
data_ST_dir <<- file.path(project_dir,"jtraets","external_data","Vignesh_spatial_2026")
seurat_object_fn <<- file.path(data_ST_dir,"/combined_seurat_obj_normalized_genes_of_interest.RDS")
#-------------------------------------------------------------------------------


# Signatures
#-------------------------------------------------------------------------------
p6_bailey <<- file.path(project_dir,"jtraets","data","signatures","BaileyP_UK_DvP_SuppDataFile5.csv")
#-------------------------------------------------------------------------------