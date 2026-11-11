#===============================================================================
# Aim:      Main script of the GEP publication. Runs, in order, every analysis
#           and figure script in src/ to reproduce all results, figures and
#           tables of the publication.
# Input:    - Raw/processed data files listed in src/config.R
# Output:   Written to results_dir by the sourced scripts:
#           - data_rnaseq/                : processed RNA-seq data objects (.rds)
#           - intermediate/clinical/      : multiple imputation results
#           - intermediate/sccore_gep/    : SCCore-GEP models, performances, source data
#           - intermediate/integration/   : GEP, WES and GEP + WES models, source data
#           - intermediate/validation/    : validation predictions and results
#           - dgea/                       : DESeq2 results per design
#           - rna_seq/clustering_cola/    : clustering results (cola)
#           - publication/                : final figures, tables, source data
# Requires: The project directory structure on the analysis server (see
#           "Settings" below) and the R environment restored with renv.
# Usage:    Rscript main.R   (or source("main.R") from the project root)
# Notes:    Scripts are sourced into the global environment, so objects created
#           here (code_dir, results_dir, experiment_setup, subset_oi, ...) are
#           used as implicit inputs by the sourced scripts.
#===============================================================================

# 1. Settings
#-------------------------------------------------------------------------------
# Root folders; adapt project_dir/working_dir when running on another machine
project_dir <-  file.path("/home","cscc_project_230717")
working_dir <- file.path(project_dir, "rruiter", "gep_publication")
code_dir <- file.path(working_dir, "src")
results_dir <- file.path(working_dir, "results")
#-------------------------------------------------------------------------------

# 2. Load config file (global file paths to all input data)
#-------------------------------------------------------------------------------
source(file.path(code_dir, "RNASeq_00_config.R"))
#-------------------------------------------------------------------------------

# 3. Prepare clinical data: multiple imputation of missing clinical variables
#-------------------------------------------------------------------------------
# D-SQUAME discovery dataset
print("Multiple imputation in D-SQUAME discovery dataset...")
source(file.path(code_dir, "clinical", "multiple_imputation_clinical_data_d_squame_discovery.R"))
# D-SQUAME validation dataset
print("Multiple imputation in D-SQUAME validation dataset...")
source(file.path(code_dir, "clinical", "multiple_imputation_clinical_data_d_squame_validation.R"))
#-------------------------------------------------------------------------------

# 4. Load RNA-seq data and save it as data objects (.rds) for later scripts
#-------------------------------------------------------------------------------

source(file.path(code_dir,"functions","loadData.R"))
output_data_dir <<- file.path(results_dir,"data_rnaseq")
if (!dir.exists(output_data_dir)){dir.create(output_data_dir, recursive = T)}
# Discovery dataset with complete bad QC sets removed (N = 366)
load_discovery(p1_tpms,p2_counts,p4_clinical_parsed,p7_badQC,p9_gtf,output_data_dir)
# Discovery dataset with only the bad QC samples removed (N = 378)
load_discovery_all(p1_tpms,p2_counts,p4_clinical_parsed,p7_badQC,p9_gtf,output_data_dir)
# Load and save validation data set
load_validation_d_squame(e1_counts,e2_tpms,e3_clinical,e6_badQC,p9_gtf,output_data_dir)
# Load and save Nassir et al. data set
load_validation_nassir(n1_tpms,n2_counts,n3_clinical,p9_gtf,output_data_dir)
output_data_dir <- file.path(results_dir, "data_rnaseq")
#-------------------------------------------------------------------------------

# 5. SCCore-GEP development in the D-SQUAME discovery dataset
#-------------------------------------------------------------------------------
# Directories used by the SCCore-GEP scripts
dir_scripts_sccore_gep<- file.path(code_dir ,"sccore_gep","analyses")
dir_scripts_functions <- file.path(code_dir ,"functions")
dir_scripts_figures_sccore_gep <- file.path(code_dir,"sccore_gep", "figures")
dir_results_intermediate_sccore_gep <- file.path(results_dir, "intermediate", "sccore_gep")
dir_results_publication <- file.path(results_dir, "publication")
# Experiment identifiers read by Model_input_processing.R / the model scripts
experiment_setup <- "sccore_gep"
experiment_sccore_gep <- "Discovery_models"

# 5.1 Model development (inner bootstrap only; outer = NULL means no outer CI loop)
message("Model development sccore-gep")
# Late feature integration and late prediction integration models
source(file.path(dir_scripts_sccore_gep, "Late_integration.R"))
message("Model development DvP, non-DvP and early integration")

# DvP, non-DvP and early integration models
source(file.path(dir_scripts_sccore_gep, "DvP_early_integration.R"))
message("Training discovery models complete.")

# 5.2 Collect model performances and draw Supplementary Figure 10A
source(file.path(dir_scripts_sccore_gep, "process_output_models.R"))
message("Performance data processing complete.")
source(file.path(dir_scripts_figures_sccore_gep, "supp_figure_10A_performance_boxplot.R"))
message("Supplementary figure 10a complete.")

# 5.3 Outer bootstrap for confidence intervals of SCCore-GEP performance.
# Late_integration.R is run as 10 parallel Rscript processes, each covering 10
# outer bootstrap iterations (i-9 .. i). All paths are passed as command-line
# arguments because the child processes do not share this R session.
message("running outer CI on SCCore-GEP")
outer_values <- seq(10, 100, 10)
mclapply(outer_values,function(i) { message(paste0("Running outer bootstrap on ", i - 9, "-", i))
    system2("Rscript",args = c( file.path(dir_scripts_sccore_gep, "Late_integration.R"),
        i, 
        dir_scripts_sccore_gep,
        dir_results_intermediate_sccore_gep,
        dir_scripts_functions,
        p1_tpms ,
        p2_counts,
        p7_badQC,
        project_dir,
        p11_weights_unmatched,
        p10_weights,
        p12_genes_to_keep,
        p6_bailey,
        p13_mutation_data,
        p14_all_pathways,
        p15_input138,
        experiment_setup,
        results_dir,
        p9_gtf))}, mc.cores = 10)

message("running outer CI complete")
# Combine the outer bootstrap batches into median and 95% CI per model
source(file.path(dir_scripts_sccore_gep, "process_CI_sccore_gep.R"))
message("processing CI complete")

message("Compute operating point")
# 5.4 Choose the SCCore-GEP operating point (risk threshold) on discovery data
source(file.path(dir_scripts_sccore_gep, "Discovery_operating_point.R"))


# 5.5 Performances of BWH, AJCC8 and EMC model in the D-SQUAME discovery dataset
# (SCCore-GEP performance comes from the internal validation above)
print("Compute performances of BWH, AJCC8, EMC model in D-SQUAME discovery dataset...")
source(file.path(code_dir, "sccore_gep", "analyses", "BWH_AJCC8_EMCmodel_performances.R"))

# 5.6 Figures: Supplementary Figure 10B
source(file.path(code_dir, "sccore_gep", "figures", "supp_figure_10B_performance_forestplot.R"))

#-------------------------------------------------------------------------------

# 6. SCCore-GEP integration with WES (126 unmatched samples)
#-------------------------------------------------------------------------------
message("SCCore-GEP integration with WES...")
# Directories for integration
dir_scripts_integration <- file.path(code_dir ,"integration","analyses")
# Directories used by the integration scripts
dir_scripts_figures_integration <- file.path(code_dir, "integration", "figures")
dir_results_intermediate_integration <- file.path(results_dir, "intermediate", "integration")
experiment_setup <- "integration"
# GEP-only model on the 126 unmatched samples
outer <- NULL
source(file.path(dir_scripts_integration, "unmatched_GEPmodel.R"))
message("GEP model complete.")
# WES-only model on the 126 unmatched samples
source(file.path(dir_scripts_integration, "WESmodel.R"))
message("WES model complete.")
# Combined GEP + WES model
source(file.path(dir_scripts_integration, "GEPWESmodel.R"))
message("Combined model complete.")
# Supplementary Figure 12: comparison of GEP, WES and GEP + WES models
source(file.path(dir_scripts_figures_integration, "supp_figure_12_WESGEP_integration.R"))
message("Supplementary figure 12 complete.")

#-------------------------------------------------------------------------------

# 7. SCCore-GEP validation (D-SQUAME validation and Nassir et al. datasets)
#-------------------------------------------------------------------------------
# 7.1 Apply SCCore-GEP to the D-SQUAME validation dataset
print("Apply SCCore-GEP in D-SQUAME validation dataset...")
source(file.path(code_dir, "validation", "d_squame", "analyses", "apply_SCCoreGEP.R"))

# 7.2 SCCore-GEP, BWH, AJCC8, EMC model performances in D-SQUAME validation.
# subset_oi selects the patient subset analysed by the sourced script.
## BWH T1-T2a subset
print("Compute performances of SCCore-GEP, BWH, AJCC8, EMC model in T1-T2a subset of D-SQUAME validation dataset...")
subset_oi <- "t1t2a"
source(file.path(code_dir, "validation", "d_squame", "analyses", "SCCoreGEP_BWH_AJCC8_EMCmodel_performances.R"))
## Entire dataset
print("Compute performances of SCCore-GEP, BWH, AJCC8, EMC model in entire D-SQUAME validation dataset...")
subset_oi <- "entiredataset"
source(file.path(code_dir, "validation", "d_squame", "analyses", "SCCoreGEP_BWH_AJCC8_EMCmodel_performances.R"))
## AJCC8 T1-T2 subset
print("Compute performances of SCCore-GEP, BWH, AJCC8, EMC model in T1-T2 subset of D-SQUAME validation dataset...")
subset_oi <- "t1t2"
source(file.path(code_dir, "validation", "d_squame", "analyses", "SCCoreGEP_BWH_AJCC8_EMCmodel_performances.R"))

# 7.3 SCCore-GEP and BWH performances in the Nassir et al. dataset
print("Compute performances of SCCore-GEP, BWH, in Nassir et al. validation dataset...") 
source(file.path(code_dir, "validation", "nassir", "analyses", "SCCoreGEP_BWH_performances.R"))

# 7.4 Multivariable analysis in D-SQUAME validation (requires recalibration first)
print("Before performing multivariable analysis, perform recalibration in D-SQUAME validation dataset...")
source(file.path(code_dir, "validation", "d_squame", "analyses", "recalibration.R"))
## BWH T1-T2a subset
print("Multivariable analysis in T1-T2a subset of D-SQUAME validation dataset...")
subset_oi <- "t1t2a"
source(file.path(code_dir, "validation", "d_squame", "analyses", "multivariable_analysis.R"))
## Entire dataset
print("Multivariable analysis in entire D-SQUAME validation dataset...")
subset_oi <- "entiredataset"
source(file.path(code_dir, "validation", "d_squame", "analyses", "multivariable_analysis.R"))

# 7.5 Threshold-based metrics (sensitivity, specificity, ...) in D-SQUAME validation
print("Compute threshold-based metrics in D-SQUAME validation dataset...")
source(file.path(code_dir, "validation", "d_squame", "analyses", "threshold_based_metrics.R"))

# 7.6 Figures
print("Publication figures of SCCore-GEP validation results...")
source(file.path(code_dir, "validation", "figures", "figure_5b_t1t2a_performance_forestplot.R"))
source(file.path(code_dir, "validation", "d_squame", "figures", "figure_5c_t1t2a_emcmodel_performance_comp_forestplot.R"))
# NOTE: Supplementary Figure 13D must run before Figure 5D, because 5D uses the
# sensitivity computed in 13D
source(file.path(code_dir, "validation", "d_squame", "figures", "supp_figure_13D_t1t2a_scatterplot.R"))
source(file.path(code_dir, "validation", "d_squame", "figures", "figure_5d_t1t2a_precision_recall_plot.R"))
source(file.path(code_dir, "validation", "d_squame", "figures", "supp_figure_13B_t1t2a_performance_ext_forestplot.R"))
source(file.path(code_dir, "validation", "d_squame", "figures", "supp_figure_13C_t1t2a_multivariable_analysis.R"))
source(file.path(code_dir, "validation", "d_squame", "figures", "supp_figure_14A_entiredataset_performance_forestplot.R"))
source(file.path(code_dir, "validation", "d_squame", "figures", "supp_figure_14B_entiredataset_multivariable_analysis.R"))
source(file.path(code_dir, "validation", "figures", "supp_figure_14C_entiredataset_performance_roc_curve.R"))
source(file.path(code_dir, "validation", "d_squame", "figures", "supp_figure_15A_t1t2_performance_forestplot.R"))
source(file.path(code_dir, "validation", "d_squame", "figures", "supp_figure_15B_t1t2_precision_recall_plot.R"))
#-------------------------------------------------------------------------------

# 8. RNA-seq, single-cell and spatial sequencing analyses
#-------------------------------------------------------------------------------
# 8.1 Differential gene expression analysis on the discovery dataset (N = 378)
source(file.path(code_dir,"rna_seq","analyses","DGEA_discovery.R"))
# 8.2 Unsupervised clustering: Figure 4A-B and Supplementary Figure 7A-C
source(file.path(code_dir,"rna_seq","figures","rna_seq_02_clusters_visualization_interpretation.R"))
# 8.3 Heatmaps of the discovery and validation datasets
source(file.path(code_dir,"rna_seq","figures","figure_5a_heatmap_discovery.R"))
source(file.path(code_dir,"rna_seq","figures","supp_figure_13A_heatmaps_validation.R"))
# 8.4 DGEA volcano plots and GSEA bar plots
source(file.path(code_dir,"rna_seq","figures","figure_4c_and_supp_figure_7d_DEGs.R"))
source(file.path(code_dir,"rna_seq","figures","supp_figure_8a_DEG.R"))
source(file.path(code_dir,"rna_seq","figures","figure_4d_GSEA.R"))
source(file.path(code_dir,"rna_seq","figures","supp_figure_8b_GSEA.R"))
# 8.5 Spatial sequencing: fraction of reads per cell type
source(file.path(code_dir,"single_cell_spatial","figure_5a_spatial_data_N14.R"))
# 8.6 Single-cell sequencing: fraction of reads per cell type
source(file.path(code_dir,"single_cell_spatial","supp_figure_11A_single_cell_Andrew.R"))

#-------------------------------------------------------------------------------
