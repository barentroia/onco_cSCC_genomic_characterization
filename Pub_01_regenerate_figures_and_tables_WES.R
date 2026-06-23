#---------------------------------------------------
# Aim: Master script to generate figures and tables
# Author: B. Rentroia Pacheco
# Input: Mutect2, CNV analyses from DNAnexus, compiled clinical set
# Output: Figures and Table for publication
#---------------------------------------------------

#---------------------------------------------------
# 0. Libraries
#---------------------------------------------------
dir_scripts = file.path("T:","Barbara","WES_analyses","Publication","UCSF_WES_CSCC")
source(file.path(dir_scripts,"Pub_00_required_libraries.R"))


#---------------------------------------------------
# 1. Directories
#---------------------------------------------------
dir_results =file.path("T:","Barbara","WES_analyses","Publication")
#---------------------------------------------------
# 2. Assembling clinical information and manual inspection information
#---------------------------------------------------
# Note - these need to be moved to the publication results folder:
sample_summary_total = read.xlsx("T:/Barbara/WES_analyses/Results/Meta_data/DA_04_Sample_summary_extended_version_incl_exclusion_v2.xlsx")
sample_file_names = read.xlsx("T:/Barbara/WES_analyses/Results/Meta_data/File_names_all_batches.xlsx")

# Get patient IDS:
df_map=read.xlsx("T:/Barbara/EMC_cSCC/clinical_data/Mapping_Sky_UCSF_EMC.xlsx")
# Get RNA-seq exclusion:
df_rnaseq_pairs = read.csv("T:/Barbara/EMC_cSCC/clinical_data/sample_ids_bad_QC_pairs.csv")
df_rnaseq_exc = read.csv("T:/Barbara/EMC_cSCC/clinical_data/sample_ids_bad_QC.csv")

# Focal copy number:
focal_copy_nr = read.xlsx("T:/Barbara/WES_analyses/Results/Manual_inspections/Focal_CNVS/Revised_focal_CN_v1.xlsx")
#copy_number_data = read.xlsx("T:/Barbara/WES_analyses/Results/CNV/Summary_CNV/Cutoff_investigation/Excluding_lower_0.268/b03_02_05_bin_matrix_cnv_armlevel_0.1.xlsx")
#frac_genome_altered = read.xlsx("T:/Barbara/WES_analyses/Results/CNV/Summary_CNV/Cutoff_investigation/Excluding_lower_0.268/b03_02_05_alt_genome_0.1.xlsx")

# HLA information
HLA_info = read.xlsx("T:/Harsh/HLA/Paper/Genetic_alteration_samples.xlsx")

#---------------------------------------------------
# 3. Settings
#---------------------------------------------------
exclusion_callable_coverage = TRUE
rescued_mutations = TRUE
IS_classification = "IS.at.cSCC"
exome_mut_analysis_paramter="TRUE"
if(exclusion_callable_coverage){
  sample_summary = sample_summary_total %>%filter(Exclude_sample=="No")%>%as.data.frame()
}

if(rescued_mutations){
  sample_summary = sample_summary %>%
    select(-contains("without.RESCUED"))%>%
    rename_with(~ gsub("\\.with\\.RESCUED", "", .x), contains("with.RESCUED"))%>%as.data.frame()
  resc_info="with_rescued_muts"
  mega_MML =  read.xlsx("T:/Barbara/WES_analyses/Results/MML/Tumor_only_analyses/03_Annotation_MMLists/b02_01_13_mega_MML_147_woncokb_05_revised.xlsx")
  gr_analysis = "Somatic_and_Rescued_Clonal" # For 04_01_02 script on mutational signatures
}else{
  resc_info="without_rescued_muts"
}

mega_MML_clonal = mega_MML %>%
  dplyr::filter(Clonal=="Yes")%>%
  as.data.frame()

# Threshold for cancer gene analysis:
threshold = 0.01

# Alterations:
# Mutation type: synonymous and nonsynonymous
nonsynonymous <- c("START_CODON_SNP", "Nonstop_Mutation", "DE_NOVO_START_OUT_FRAME", 
                   "DE_NOVO_START_IN_FRAME", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Splice_Site", "Nonsense_Mutation")
synonymous <- c("Silent", "Intron", "IGR", "RNA", "COULD_NOT_DETERMINE","3'UTR", "5'UTR","5'Flank")

# Load genes and pathways:
genes_and_pathways_correspondence = read.xlsx("T:/Barbara/WES_analyses/Methods/Pathogenic_mutations/2025_10_20_Tileplot_genes_and_pathways.xlsx")
gene_path_dictionary = list("genes"=genes_and_pathways_correspondence$Gene,"pathways"=genes_and_pathways_correspondence$Pathway)
oncogenes = c("HRAS","PIK3CA","EZH2","RRAS2","TERT","CCND1","KRAS","RAC1","GNA11","GNAS","RAP1B","ERBB3","MTOR","MYC","EGFR","YAP1","ERBB2","FGFR2","FGFR3","CARD11","JAK3","MAP2K1","SF3B1","PPM1D","BRAF","FOXA1","PTPN11")
tsps = c("TP53","CDKN2A","ZNF750","KMT2D","NOTCH1","FAT1","FAT2","ARID2","NOTCH2","CASP8","PTEN","CHUK","PBRM1","TGFBR2","USP9X","IGF2BP2","RB1", "CREBBP","FAS","PMS2","AJUBA","USP28","NFE2L2","ARID1A","EP300","RHOA","LAMC1")
pathogenic_mutations_shain_lab = read.xlsx("T:/Barbara/WES_analyses/Methods/Pathogenic_mutations/2025_09_08_Pathogenic_mutations_info_BR_working_version.xlsx")

# Color palette:
#group_colors_IS = c("Yes"="#8CCC74","No"="#1B5E20")
group_colors_IS = c("Yes"="darkgray","No"="darkgray")
#group_colors_Mets = c("Case"="#D95F5F","Control"="#5F9ED1")
group_colors_Mets = c("Case"="black","Control"="darkgray")
group_colors_Mets_IS = c("Case_IS_No"="#F8766D","Control_IS_No"="#00BFC4","Case_IS_Yes"="palevioletred3","Control_IS_Yes"="turquoise4")

# Results directory with these settings:
dir_results = file.path(dir_results,paste0("Tot_",nrow(sample_summary),"_",resc_info))
dir.create(dir_results)

# Connect to biomart:
mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  GRCh = 37
)

# CNV information:
chromosome_bands =read.table("T:/Barbara/WES_analyses/Methods/CNV/cytoBand.txt")
ct_cnv = 0.1 # This was the cutoff used for determining copy number alterations
include_CNV_samples="All_CC_but_tcel_0.135"
#---------------------------------------------------
# 5. Generate figures
#---------------------------------------------------

# Save dataframe with full sample IDS:
#stopifnot(sample_summary_total$Sample.ID==sample_file_names$Folder_names[match(sample_summary_total$Sample.ID,sample_file_names$Folder_names)])
#sample_summary_total$Full_Sample_ID =  gsub("_.*","",gsub("S","",gsub("^Hyb(?:[^_]*_){2}","",gsub("_DNA.*","",gsub("_fastp.*","",sample_file_names$Tumor_prefix[match(sample_summary_total$Sample.ID,sample_file_names$Folder_names)])))))
sample_summary_total$Full_Sample_ID = df_map$WES_sample_ID[match(sample_summary_total$Sample.ID,df_map$UCSF_short_ID)]
sample_summary_total$Pt_ID =  gsub("[A-Za-z]", "",sample_summary_total$Full_Sample_ID)

# Patient characteristics table:
dir_results_script = file.path(dir_results,"00_Pat_Tumor_characteristics")
#source(file.path(dir_scripts,"04_Publication_scripts","04_00_01_Pat_Tum_characteristics_table.R"))

# Preprocess meta-analysis cohort, so that it is compatible with our dataframes:
#source(file.path(dir_scripts,"04_Publication_scripts","04_00_03_Meta_analysis_cohort_preprocessing.R"))

# Figure 1:
# Mutation burden - plots:
dir_results_script = file.path(dir_results,"01_Mutations_summary",IS_classification)
#source(file.path(dir_scripts,"04_Publication_scripts","04_01_01_MML_summaries.R"))
source(file.path(dir_scripts,"04_Publication_scripts","04_01_A2_MML_tumor_in_normal_contamination_figures.R"))

# Mutational signature plot:
dir_results_script = file.path(dir_results,"02_Mutational_signatures")
dir.create(dir_results_script)
source(file.path(dir_scripts,"04_Publication_scripts","04_01_02_MML_Mutational_signatures.R"))

# New classification of immunosuppressed patients:
IS_classification = "IS.at.cSCC.updatedS32"
pts_with_SBS_32=read.csv(file.path(dir_results_script,"Data","SigProfilerAssignment",paste0("ExomeRN_",exome_mut_analysis_paramter),gr_analysis,"/S04_01_02_SBS32_patients.csv"))[,2]

# Repeat previous points with new classification of immunosuppressed patients:
dir_results_script = file.path(dir_results,"01_Mutations_summary",IS_classification)
sample_summary$IS.at.cSCC.updatedS32 = ifelse(sample_summary$Sample.ID%in%gsub("S","",pts_with_SBS_32),"Yes",sample_summary$IS.at.cSCC)
#source(file.path(dir_scripts,"04_Publication_scripts","04_01_01_MML_summaries.R"))
dir_results_script = file.path(dir_results,"02_Mutational_signatures")
#source(file.path(dir_scripts,"04_Publication_scripts","04_01_02_MML_Mutational_signatures.R"))

# Cancer Gene Analyses - Figure 2:
dir_results_script = file.path(dir_results,"03_Cancer_gene_analysis")
dir.create(dir_results_script)
dir_results_script = file.path(dir_results,"03_Cancer_gene_analysis",IS_classification)
dir.create(dir_results_script)
#source(file.path(dir_scripts,"04_Publication_scripts","04_02_01_Intogen.R")) # To run intogen in server

# Revision of intogen results for pathogenic calls:
if(exclusion_callable_coverage & rescued_mutations){
  driver_genes_oficial_list = read.xlsx(file.path(dir_results_script,"04_02_01_intogen_results_all_with_man_revision_v3.xlsx"),sheet="Intogen_results",startRow = 2)
  oncokb_muts =  read.xlsx(file.path(dir_results_script,"04_02_01_oncokb_v3.xlsx"),sheet = "OncoKB pathogenic")
}else{
  print("Revision missing")
}

# Copy number alterations:
dir_results_script = file.path(dir_results,"04_CNV")
dir.create(dir_results_script)
source(file.path(dir_scripts,"04_Publication_scripts","04_02_02_CNV.R")) 

# Tileplot:
# Add copy number alterations to the sample summary:
copy_number_data_tileplot = read.xlsx(file.path(dir_results_script,"cns","All_CC_but_tcel_0","S04_02_02_bin_matrix_cnv_armlevel_0.1.xlsx"))
dir_results_script = file.path(dir_results,"05_Tileplot")
dir.create(dir_results_script)
dir_results_script = file.path(dir_results,"05_Tileplot",IS_classification)
dir.create(dir_results_script)
source(file.path(dir_scripts,"04_Publication_scripts","04_02_03_Tileplot.R")) 

# Forest plots:
dir_results_script = file.path(dir_results,"06_Forest_plots")
dir.create(dir_results_script)
# Driver genes and copy number forest plots:
source(file.path(dir_scripts,"04_Publication_scripts","04_02_03_Driver_mutations_differences.R"))
# Supplementary Figure: IS patients:
source(file.path(dir_scripts,"04_Publication_scripts","04_02_04_Forest_plots_IS.R"))

# Co-occurence analysis: 
dir_results_script = file.path(dir_results,"05_Tileplot",IS_classification)
source(file.path(dir_scripts,"04_Publication_scripts","04_02_03_Co_occurrence_analysis.R")) 

# Analysis of all alterations together:
dir_results_script = file.path(dir_results,"07_Bootstrapping")
dir.create(dir_results_script)
recompute_BS = FALSE
source(file.path(dir_scripts,"04_Publication_scripts","04_03_02_Bootstrapping.R")) 

# Integrate driver mutations with genomic clusters:
dir_results_script = file.path(dir_results,"08_RNAseq_integration")
dir.create(dir_results_script)
source(file.path(dir_scripts,"04_Publication_scripts","04_04_Integration.R"))

# Save sample summary with extra information:
source(file.path(dir_scripts,"04_Publication_scripts","05_01_01_Sample_Summary.R"))
