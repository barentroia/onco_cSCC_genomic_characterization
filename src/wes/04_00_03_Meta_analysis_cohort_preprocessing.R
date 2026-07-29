#---------------------------------------------------
# Aim: Preprocessing of Chang and Shain 2021 cohort
# Author: B. Rentroia Pacheco
# Input: MM lists
# Output: Genes under selection
#---------------------------------------------------

#-----------------------
# 0. Load data
#-----------------------
# Original data from Chang and Shain 2021:
original_data_chang = read.xlsx("T:/Barbara/WES_analyses/Results/Chang_2021/Stab1_Chang_Shain.xlsx")
original_data_repeated_sample  = original_data_chang[which(original_data_chang $Tissue=="Primary2"),]
original_data_chang = original_data_chang[which(original_data_chang$Tissue=="Normal"),] # Because clinical information is in these rows
original_data_chang = rbind(original_data_chang,original_data_repeated_sample)
# Here there are 105 samples:
# We now exclude the ones marked as excluded from analysis:
original_data_chang = original_data_chang[which(original_data_chang $`Included.in.All,.Partial,.or.None.Analyses`=="All"),]
original_data_chang$`Tumor.Cellularity.(%)` = round(as.numeric(original_data_chang$`Tumor.Cellularity.(%)`),2)
original_data_chang$Tcel_char = as.character(original_data_chang$`Tumor.Cellularity.(%)`)

# MML from the meta-analysis of Chang and Shain:
#darwins_mml_72 = read.xlsx("T:/Darwin/Shain_Calls.xlsx")
#colnames(darwins_mml_72)= darwins_mml_72[1,]
#darwins_mml_72=darwins_mml_72%>% dplyr::slice(-1)%>%as.data.frame()
#stopifnot(length(unique(darwins_mml_72$Sample))==72)

#darwins_mml_other = read.xlsx("T:/Darwin/Shain_Calls_Yilmaz_Cho_Durinck.xlsx")
#colnames(darwins_mml_other)= darwins_mml_other[1,]
#darwins_mml_other=darwins_mml_other%>% slice(-1)%>%as.data.frame()

# Sample summary of the meta-analysis of Chang and Shain:
sample_summary_chang = read.xlsx("T:/Barbara/WES_analyses/Results/Chang_2021/map.xlsx")
sample_summary_chang= sample_summary_chang[1:83,]
sample_summary_chang = sample_summary_chang[,-1]
sample_summary_chang$Tumor_cellularity =round(as.numeric(sample_summary_chang$Tumor_cellularity),2)
sample_summary_chang$Tcel_char = as.character(sample_summary_chang$Tumor_cellularity)

# Combined the two via Tcel:
merged_data = merge(original_data_chang,sample_summary_chang[,c("Sample","Tcel_char","Subtype")],by="Tcel_char",all.x=TRUE) # Only a XP sample cannot be mapped this way.

# Combine all master mutation lists:
mega_MML_chang=data.frame()
for(i in 1:nrow(sample_summary_chang)){
  paper_folder = sample_summary_chang$Paper[i]
  sample_folder = sample_summary_chang$Sample[i]
  if(sample_folder =="AJ-01"){
    dir_mml_folder = "T:/Darwin/Ji_Cell/Removed_Samples/Low_Tumor_Cellularity/AJ-01/Mutations"
  }else{
    dir_mml_folder = file.path("T:/Darwin/",paper_folder,sample_folder,"Mutations")
  }
  dir_mml_files = list.files(dir_mml_folder,, pattern = "^[^~]")
  mml_file_path = dir_mml_files[grepl("Pass_Filter_MasterMutationList.xlsx",dir_mml_files)]
  

  if(length(mml_file_path)>0){
    mml_file = read.xlsx(file.path(dir_mml_folder,mml_file_path))
    #print(head(mml_file))
    colnames(mml_file)=mml_file[1,]
    mml_file = mml_file[-1,]
    mml_file$Sample = sample_folder
    
    if("Call"%in%colnames(mml_file)){
      stopifnot(sum(mml_file$Call=="PASS")==nrow(mml_file))
      mml_file$Call=NULL
    }
    
    if("Normalized"%in%colnames(mml_file)){
      mml_file$Normalized=NULL
    }
    if("Normalized_MAF"%in%colnames(mml_file)){
      mml_file$Normalized_MAF=NULL
    }
    
    if(c("Met_Ref")%in%colnames(mml_file)){
      mml_file[,c("Met_Ref","Met_Mut","Met_MAF")]=NULL
    }
    
    #Add Normalized MAF:
    tcel = sample_summary_chang$Tumor_cellularity[i]/100
    mml_file$Normalized_MAF = as.numeric(mml_file$Tumor_MAF)/tcel
    mega_MML_chang = rbind(mega_MML_chang,mml_file)
  }else{

  }
}
mega_MML_chang$Normalized_MAF = ifelse(mega_MML_chang$Normalized_MAF >1,1,mega_MML_chang$Normalized_MAF )
# Sanity check: Checked that data encompassess the Shain_calls file:
#for (smp in intersect(darwins_mml_72$Sample,mega_MML_chang$Sample)){
#  d_file = darwins_mml_72 %>%filter(Sample==smp)%>%nrow()
#  mml_mine = mega_MML_chang %>%filter(Sample==smp)%>%nrow()
#  stopifnot(d_file == mml_mine)
#}

write.xlsx(mega_MML_chang,"T:/Barbara/WES_analyses/Results/Chang_2021/04_00_03_mega_MML_Chang_and_Shain.xlsx")

# Obtain information of Immunosuppressed and immunocompetent patients:
sample_summary_extensive_chang = read.xlsx("T:/Barbara/WES_analyses/Results/Chang_2021/clinical_annotations_updated.xlsx")
sample_summary_extensive_chang = sample_summary_extensive_chang[1:139,]# the remaining ones are AKs
sample_summary_extensive_chang = sample_summary_extensive_chang[which(sample_summary_extensive_chang$Tissue=="Normal"),] # Because clinical information is in these rows
sample_summary_extensive_chang$Sample_ID = gsub("\\/N","",gsub("_N.*","",sample_summary_extensive_chang$X7))
#sample_summary_darwin_extended=merge(sample_summary_darwin,sample_summary_extensive_chang[,c("Sample_ID","Burden.(Mutations/.megabase)","Immune.Status","Fraction.of.Footprint.to.Target.Territory.for.Reference/Tumor.Pair.(%)")],by.x="Sample",by.y="Sample_ID",all.x=TRUE)

chang_IS = sample_summary_extensive_chang$Sample_ID[which(sample_summary_extensive_chang$Immune.Status=="Suppressed")] # 36, out of which 30 are kept in the analysis
chang_IC = sample_summary_extensive_chang$Sample_ID[which(sample_summary_extensive_chang$Immune.Status=="Competent")] # 31, out of which 28 are kept in the analysis
rdeb_samples = sample_summary_chang$Sample[sample_summary_chang$Subtype=="RDEB"]# 25 samples
# Note these are 36+31 samples, the remaining samples are all RDEB samples
stopifnot(unique(gsub("-.*","",setdiff(sample_summary_chang$Sample,c(chang_IS,chang_IC))))=="RC")

# Check that we have all immunosuppressed and immunocompetent:
length(unique(original_data_chang$Sample)) #83
length(unique(mega_MML_chang$Sample)) #83
length(intersect(chang_IS,mega_MML_chang$Sample)) # 30
length(intersect(chang_IC,mega_MML_chang$Sample)) # 28

# Preprocess original data so we can extract mutation burden and callable coverage correctly
original_data_chang = original_data_chang %>%
  mutate(`Fraction.of.Footprint.to.Target.Territory.for.Reference/Tumor.Pair.(%)`=as.numeric(`Fraction.of.Footprint.to.Target.Territory.for.Reference/Tumor.Pair.(%)`),
         `Burden.(Mutations/.megabase)`=as.numeric(`Burden.(Mutations/.megabase)`))%>%
  as.data.frame()
