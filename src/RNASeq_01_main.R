#-------------------------------------------------------------------------------
# Aim: Main script calling all scripts to run analyses and generate publication results
#-------------------------------------------------------------------------------

# Settings
#-------------------------------------------------------------------------------
project_dir <-  file.path("/home","cscc_project_230717")
working_dir <- file.path(project_dir, "jtraets", "tmp","gep_publication")
code_dir <- file.path(working_dir, "src")
results_dir <- file.path(working_dir, "results")
#-------------------------------------------------------------------------------

# Load config file (Global)
#-------------------------------------------------------------------------------
source(file.path(code_dir, "RNASeq_00_config.R"))
#-------------------------------------------------------------------------------

# Prepare clinical data running multiple imputation
#-------------------------------------------------------------------------------
# D-SQUAME discovery dataset
print("Multiple imputation in D-SQUAME discovery dataset...")
source(file.path(code_dir, "clinical", "multiple_imputation_clinical_data_d_squame_discovery.R"))
# D-SQUAME validation dataset
print("Multiple imputation in D-SQUAME validation dataset...")
source(file.path(code_dir, "clinical", "multiple_imputation_clinical_data_d_squame_validation.R"))
#-------------------------------------------------------------------------------

# Load RNAseq data
#-------------------------------------------------------------------------------
source(file.path(code_dir,"functions","loadData.R"))
output_data_dir <<- file.path(results_dir,"data_rnaseq")
if (!dir.exists(output_data_dir)){dir.create(output_data_dir, recursive = T)}
# Load and save discovery data set, N366
load_discovery(p1_tpms,p2_counts,p4_clinical_parsed,p7_badQC,p9_gtf,output_data_dir)
# Load and save discovery data set, N378
load_discovery_all(p1_tpms,p2_counts,p4_clinical_parsed,p7_badQC,p9_gtf,output_data_dir)
# Load and save validation data set
load_validation_d_squame(e1_counts,e2_tpms,e3_clinical,e6_badQC,p9_gtf,output_data_dir)
# Load and save Nassir et al. data set
load_validation_nassir(n1_tpms,n2_counts,n3_clinical,p9_gtf,output_data_dir)
#-------------------------------------------------------------------------------

# SCCore-GEP related analysis in D-SQUAME discovery dataset
#-------------------------------------------------------------------------------
# Directories for sccore-gep
dir_scripts_sccore_gep<- file.path(code_dir ,"sccore_gep","analyses")
dir_scripts_functions <- file.path(code_dir ,"functions")
dir_scripts_figures_sccore_gep <- file.path(code_dir,"sccore_gep", "figures")
dir_results_intermediate_sccore_gep <- file.path(results_dir, "intermediate", "sccore_gep")
dir_results_publication <- file.path(results_dir, "publication")
experiment_setup <- "sccore_gep"
experiment_sccore_gep <- "Discovery_models"
# Model development sccore-gep
message("Model development sccore-gep")
outer <- NULL
source(file.path(dir_scripts_sccore_gep, "Late_integration.R"))
message("Model development DvP, non-DvP and early integration")
source(file.path(dir_scripts_sccore_gep, "DvP_early_integration.R"))
message("Training discovery models complete.")

source(file.path(dir_scripts_sccore_gep, "process_output_models.R"))
message("Performance data processing complete.")

source(file.path(dir_scripts_figures_sccore_gep, "supp_figure_10a_performance_boxplot.R"))
message("Supplementary figure 10a complete.")

## running outer CI
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
source(file.path(dir_scripts_sccore_gep, "process_CI_sccore_gep.R"))
message("processing CI complete")

# BWH, AJCC8, EMC model performances (SCCore-GEP performance from internal validation) in D-SQUAME discovery dataset
print("Compute performances of BWH, AJCC8, EMC model in D-SQUAME discovery dataset...")
source(file.path(code_dir, "sccore_gep", "analyses", "BWH_AJCC8_EMCmodel_performances.R"))

# Figures
source(file.path(code_dir, "sccore_gep", "figures", "supp_figure_10b_performance_forestplot.R"))

#-------------------------------------------------------------------------------

# SCCore-GEP integration with WES
#-------------------------------------------------------------------------------
message("SCCore-GEP integration with WES...")
# Directories for integration
dir_scripts_integration <- file.path(code_dir ,"integration","analyses")
dir_scripts_figures_integration <- file.path(code_dir, "integration", "figures")
dir_results_intermediate_integration <- file.path(results_dir, "intermediate", "integration")
experiment_setup <- "integration"
#Run GEP model on 126 samples unmatched
outer <- NULL
source(file.path(dir_scripts_integration, "unmatched_GEPmodel.R"))
message("GEP model complete.")
#Run WES model on 126 samples unmatched
source(file.path(dir_scripts_integration, "WESmodel.R"))
message("WES model complete.")
# Run WES+GEP model on 136 samples unmatched
source(file.path(dir_scripts_integration, "GEPWESmodel.R"))
message("Combined model complete.")
source(file.path(dir_scripts_figures_integration, "supp_figure_11_WESGEP_integration.R"))
message("Supplementary figure 11 complete.")

#-------------------------------------------------------------------------------

# SCCore-GEP validation analyses
#-------------------------------------------------------------------------------
# SCCore-GEP, BWH, AJCC8, EMC model performances in D-SQUAME validation dataset
print("Apply SCCore-GEP in D-SQUAME validation dataset...")
source(file.path(code_dir, "validation", "d_squame", "analyses", "apply_SCCoreGEP.R"))
## T1-T2a
print("Compute performances of SCCore-GEP, BWH, AJCC8, EMC model in T1-T2a subset of D-SQUAME validation dataset...")
subset_oi <- "t1t2a"
source(file.path(code_dir, "validation", "d_squame", "analyses", "SCCoreGEP_BWH_AJCC8_EMCmodel_performances.R"))
## Entire dataset
print("Compute performances of SCCore-GEP, BWH, AJCC8, EMC model in entire D-SQUAME validation dataset...")
subset_oi <- "entiredataset"
source(file.path(code_dir, "validation", "d_squame", "analyses", "SCCoreGEP_BWH_AJCC8_EMCmodel_performances.R"))
## T1-T2
print("Compute performances of SCCore-GEP, BWH, AJCC8, EMC model in T1-T2 subset of D-SQUAME validation dataset...")
subset_oi <- "t1t2"
source(file.path(code_dir, "validation", "d_squame", "analyses", "SCCoreGEP_BWH_AJCC8_EMCmodel_performances.R"))

# SCCore-GEP, BWH performances in Nassir et al. dataset
print("Compute performances of SCCore-GEP, BWH, in Nassir et al. validation dataset...") 
source(file.path(code_dir, "validation", "nassir", "analyses", "SCCoreGEP_BWH_performances.R"))

# Multivariable analysis in D-SQUAME validation dataset
print("Before performing multivariable analysis, perform recalibration in D-SQUAME validation dataset...")
source(file.path(code_dir, "validation", "d_squame", "analyses", "recalibration.R"))
## T1-T2a
print("Multivariable analysis in T1-T2a subset of D-SQUAME validation dataset...")
subset_oi <- "t1t2a"
source(file.path(code_dir, "validation", "d_squame", "analyses", "multivariable_analysis.R"))
## Entire dataset
print("Multivariable analysis in entire D-SQUAME validation dataset...")
subset_oi <- "entiredataset"
source(file.path(code_dir, "validation", "d_squame", "analyses", "multivariable_analysis.R"))

# Threshold-based metrics in D-SQUAME validation dataset
print("Compute threshold-based metrics in D-SQUAME validation dataset...")
source(file.path(code_dir, "validation", "d_squame", "analyses", "threshold_based_metrics.R"))

# Figures
print("Publication figures of SCCore-GEP validation results...")
source(file.path(code_dir, "validation", "figures", "figure_5b_t1t2a_performance_forestplot.R"))
source(file.path(code_dir, "validation", "d_squame", "figures", "figure_5c_t1t2a_emcmodel_performance_comp_forestplot.R"))
source(file.path(code_dir, "validation", "d_squame", "figures", "figure_5d_t1_precision_recall_plot.R"))
source(file.path(code_dir, "validation", "d_squame", "figures", "supp_figure_12b_t1t2a_performance_ext_forestplot.R"))
source(file.path(code_dir, "validation", "d_squame", "figures", "supp_figure_12c_t1t2a_multivariable_analysis.R"))
source(file.path(code_dir, "validation", "d_squame", "figures", "supp_figure_12d_t2a_precision_recall_plot.R"))
source(file.path(code_dir, "validation", "d_squame", "figures", "supp_figure_13a_entiredataset_performance_forestplot.R"))
source(file.path(code_dir, "validation", "d_squame", "figures", "supp_figure_13b_entiredataset_multivariable_analysis.R"))
source(file.path(code_dir, "validation", "figures", "supp_figure_13c_entiredataset_performance_roc_curve.R"))
source(file.path(code_dir, "validation", "d_squame", "figures", "supp_figure_14a_t1t2_performance_forestplot.R"))
source(file.path(code_dir, "validation", "d_squame", "figures", "supp_figure_14b_t1t2_precision_recall_plot.R"))
#-------------------------------------------------------------------------------

# RNAseq and spatial sequencing related analysis
#-------------------------------------------------------------------------------
# Generate DEGs on the discovery data set (N=378)
source(file.path(code_dir,"rna_seq","analyses","DGEA_discovery.R"))
# Generate clustering results:
source(file.path(code_dir,"rna_seq","figures","rna_seq_02_clusters_visualization_interpretation.R"))
# Heatmap discovery and validation data sets
source(file.path(code_dir,"rna_seq","figures","figure_5a_heatmap_discovery.R"))
source(file.path(code_dir,"rna_seq","figures","supp_figure_12a_heatmaps_validation.R"))
# DGEA and GSEA
source(file.path(code_dir,"rna_seq","figures","figure_4c_and_supp_figure_7d_DEGs.R"))
source(file.path(code_dir,"rna_seq","figures","supp_figure_8a_DEG.R"))
source(file.path(code_dir,"rna_seq","figures","figure_4d_GSEA.R"))
source(file.path(code_dir,"rna_seq","figures","supp_figure_8b_GSEA.R"))
# Spatial sequencing, fraction of reads per cell type
source(file.path(code_dir,"single_cell_spatial","figure_5a_spatial_data.R"))

#-------------------------------------------------------------------------------

