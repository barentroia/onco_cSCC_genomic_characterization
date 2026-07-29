# Aim: Cancer gene analysis - Preprocessing
# Author: B. Rentroia Pacheco
# Input: MM lists
# Output: Genes under selection
#---------------------------------------------------

#-----------------------
# 0. Setup
#-----------------------
#sample_summary = read.xlsx("T:/Barbara/WES_analyses/Results/Meta_data/DA_04_Sample_summary_extended_version_incl_exclusion_v2.xlsx")
#dir_results = "T:/Barbara/WES_analyses/Results/Cancer_gene_discovery/"
#-----------------------
# 1. Load and preprocess mutation data
#-----------------------

########## EMC Cohort ###############

# We load the muation data for 2 settings: with and without rescued mutations:
#mega_MML = read.xlsx("T:/Barbara/WES_analyses/Results/MML/Tumor_only_analyses/02_Rescued_MMLists/b_02_01_10_mega_MML_147.xlsx")
#colnames(mega_MML)[which(colnames(mega_MML)=="Sample_id")] ="Tumor_Sample_Barcode"
# Code NA in ref and altered allele as follows:
stopifnot(length(mega_MML$Reference_Allele[is.na(mega_MML$Reference_Allele)])==0)
stopifnot(length(mega_MML$Tumor_Seq_Allele2[is.na(mega_MML$Tumor_Seq_Allele2)])==0)

# mega_MML for different subgroups:
mega_MML_IS <- mega_MML%>%filter(Tumor_Sample_Barcode%in%sample_summary$Sample.ID[which(sample_summary[,IS_classification]=="Yes")]) %>%as.data.frame()
mega_MML_IC <- mega_MML%>%filter(Tumor_Sample_Barcode%in%sample_summary$Sample.ID[which(sample_summary[,IS_classification]=="No")]) %>%as.data.frame()
mega_MML_Cases <-mega_MML%>%filter(grepl("_1",Tumor_Sample_Barcode)) %>%as.data.frame()
mega_MML_Controls <-mega_MML%>%filter(grepl("_0",Tumor_Sample_Barcode)) %>%as.data.frame()


setwd(dir_results_script)

########## Meta-analysis cohort ###############
# Preprocess Darwins MML:
mega_MML_chang = mega_MML_chang %>%
  rename(Tumor_Sample_Barcode = Sample,
         Ref_normal = Normal_Ref,
         Mut_normal = Normal_Mut,
         MAF_normal = Normal_MAF,
         Ref=Tumor_Ref,
         Mut=Tumor_Mut,
         MAF=Tumor_MAF,
         Adjusted_MAF = Normalized_MAF)%>%
  mutate(Rescued_SNP="No")%>%
  as.data.frame()

# Code NA in ref and altered allele as follows:
mega_MML_chang$Reference_Allele[is.na(mega_MML_chang$Reference_Allele)]= "-"
mega_MML_chang$Tumor_Seq_Allele2[is.na(mega_MML_chang$Tumor_Seq_Allele2)]= "-"

# Join mega MML lists:
cols_of_interest = c("Tumor_Sample_Barcode","Hugo_Symbol","Chromosome","Start_Position","End_Position","Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele2", "dbSNP_RS","dbSNP_Val_Status" ,           
                     "Genome_Change", "Protein_Change","COSMIC_overlapping_mutations", "UV","Ref","Mut","MAF", "Ref_normal","Mut_normal","MAF_normal", "Rescued_SNP","Clonal","Adjusted_MAF")
mega_MML_to_join = mega_MML[,setdiff(cols_of_interest,c("Indel","Pathogenic","TERT_promoter","cSCC_driver_mutations"))]
mega_MML_to_join =mega_MML_to_join %>%filter(Clonal=="Yes")%>%as.data.frame()
mega_MML_EMC_Chang = rbind(mega_MML_to_join [,setdiff(cols_of_interest,c("Indel","Pathogenic","TERT_promoter","cSCC_driver_mutations","Clonal"))],mega_MML_chang[,setdiff(cols_of_interest,c("Indel","Pathogenic","TERT_promoter","cSCC_driver_mutations","Clonal"))])

# Create the immunosuppressed and immunocompetent groups:
mega_MML_EMC_Chang_IC = mega_MML_EMC_Chang%>%filter(Tumor_Sample_Barcode%in%c(unique(mega_MML_IC$Tumor_Sample_Barcode),chang_IC)) %>%as.data.frame()
mega_MML_EMC_Chang_IS = mega_MML_EMC_Chang%>%filter(Tumor_Sample_Barcode%in%c(unique(mega_MML_IS$Tumor_Sample_Barcode),chang_IS)) %>%as.data.frame()

########## Bowens disease ###############
# MM lists preprocessed for each tool are available here: T:/Harsh/Bowens_disease_paper/cancer_gene_discovery/