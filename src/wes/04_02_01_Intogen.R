#---------------------------------------------------
# Aim: Cancer gene analysis - intogen data preparation
# Author: B. Rentroia Pacheco
# Input: MM lists
# Output: Genes under selection
#---------------------------------------------------

#-----------------------
# 0. Setup
#-----------------------

#sample_summary = read.xlsx("T:/Barbara/WES_analyses/Results/Meta_data/DA_04_Sample_summary_extended_version_incl_exclusion_v2.xlsx")
#dir_results = "T:/Barbara/WES_analyses/Results/Cancer_gene_discovery/"
#threshold=0.05
#-----------------------
# 1. Load and preprocess mutation data
#-----------------------

source("T:/Barbara/WES_analyses/Publication/UCSF_WES_CSCC/04_Publication_scripts/04_02_01_A1_Preprocessing_cancer_gene_analysis.R")

#-----------------------
# 2. Preprocess dataframes to a format that can be interpreted 
#-----------------------
# EMC_cohort:
options(scipen=10)

# Note: the strand should be positive, let's check if it is:
#sample_file_names = read.xlsx(paste0("T:/Barbara/WES_analyses/Results/Meta_data/File_names_all_batches.xlsx"))
#all_files = list.files("T:/Barbara/WES_analyses/Results/MML/Tumor_only_analyses/")
#funcotator_files = all_files[grepl("fastp_tonly\\.funcotator",all_files)]
#for(t_prefix in funcotator_files ){
#  if(grepl("N_",t_prefix)){
    
#  }else{
#    if(!t_prefix%in%c("Hyb-21_5_085-1-A1_DNA_S5_L002_fastp.funcotator.txt","Hyb-21_5_085-1-A1_DNA_S5_L002_fastp.funcotator_file-J19k2q80GGpVvKQyZBB02Bgf.txt")){
#        funcotator_output = read.delim(paste0("T:/Barbara/WES_analyses/Results/MML/Tumor_only_analyses/",t_prefix),          sep="\t", fill=TRUE,skip=156)
#        stopifnot(names(table(funcotator_output$Strand))=="+")
#    }
#
#  }

#}
sbatch_script_template = read_file(file.path(dir_scripts,"04_Publication_scripts/04_02_01_A1_sbatch_intogen.sh"))
annotation_template<- yaml.load_file("T:/Barbara/WES_analyses/Results/Cancer_gene_discovery/intogen_plus/annotation_mega_MML_147.yaml")
dir.create(file.path(dir_results_script,"EMC_cohort"))
for( group in c("Entire_cohort","Cases", "Controls","ISuppressed","ICompetent")){
  if(group == "Entire_cohort"){
    mega_MML_filtered = mega_MML%>%filter(Clonal=="Yes")%>%as.data.frame()
  }else if(group =="Cases"){
    mega_MML_filtered = mega_MML_Cases%>%filter(Clonal=="Yes")%>%as.data.frame()
  }else if(group =="Controls"){
    mega_MML_filtered = mega_MML_Controls%>%filter(Clonal=="Yes")%>%as.data.frame()
  }else if (group=="ISuppressed"){
    mega_MML_filtered = mega_MML_IS%>%filter(Clonal=="Yes")%>%as.data.frame()
  }else if (group=="ICompetent"){
    mega_MML_filtered = mega_MML_IC%>%filter(Clonal=="Yes")%>%as.data.frame()
  }
  print(length(unique(mega_MML_filtered$Tumor_Sample_Barcode)))
  options(scipen=10)
  dir.create(file.path(dir_results_script,"EMC_cohort",group))
  # Save the MM list
  mega_MML_filtered  %>%
    rename(Ref_tumor=Ref,
           Alt_tumor=Mut)%>%
    mutate(Strand = "+")%>%
    write.table(file.path(dir_results_script,"EMC_cohort",group,paste0("mega_MML_EMC_",group,".txt")),sep="\t" ,row.names=FALSE,quote=FALSE)
  
  # Save the associated metadata file:
  annotation_edited = annotation_template
  annotation_edited$pattern = paste0("mega_MML_EMC_",group,".txt")
  annotation_edited$annotation[[11]]$value = paste0("EMC_",group)
  write_yaml(annotation_edited, file.path(dir_results_script,"EMC_cohort",group,paste0("annotation_mega_MML_",group,".yaml")), fileEncoding = "UTF-8")
  
  # Adapt the sbatch file:
  sbatch_script = gsub("GROUP",group,sbatch_script_template)
  write_file(sbatch_script,file.path(dir_results_script,"EMC_cohort",group,paste0("sbatch_intogen_",group,".sh")))
}

# Now for the meta-analysis cohort:
dir.create(file.path(dir_results_script,"EMC_Chang"))
for( group in c("Chang_only","EMC_Chang","ISuppressed","ICompetent","EMC_Chang_noRDEB")){
  group_folder = group
  if(group == "EMC_Chang"){
    mega_MML_filtered = mega_MML_EMC_Chang
    group_folder = "Entire_cohort"
  }else if (group=="ISuppressed"){
    mega_MML_filtered = mega_MML_EMC_Chang_IS
    group="EMC_Chang_ISuppressed"
  }else if (group=="ICompetent"){
    mega_MML_filtered = mega_MML_EMC_Chang_IC
    group="EMC_Chang_ICompetent"
  }else if (group=="Chang_only"){
    mega_MML_filtered =mega_MML_chang
  }else if(group=="EMC_Chang_noRDEB"){
    mega_MML_filtered =mega_MML_EMC_Chang%>%filter(!Tumor_Sample_Barcode%in%rdeb_samples)
  }
  dir.create(file.path(dir_results_script,"EMC_Chang",group_folder))
  
  dir.create(file.path(dir_results_script,"EMC_Chang",group_folder,"input"))
  
  print(length(unique(mega_MML_filtered$Tumor_Sample_Barcode)))
  options(scipen=999)
  # Save the MM list
  mega_MML_filtered  %>%
    rename(Ref_tumor=Ref,
           Alt_tumor=Mut)%>%
    mutate(Strand = "+",
           Start_Position = as.numeric(Start_Position),
           End_Position = as.numeric(End_Position))%>%
    write.table(file.path(dir_results_script,"EMC_Chang",group_folder,"input",paste0("mega_MML_",group,".txt")),sep="\t" ,row.names=FALSE,quote=FALSE)
  
  # Save the associated metadata file:
  annotation_edited = annotation_template
  annotation_edited$pattern = paste0("mega_MML_",group,".txt")
  annotation_edited$annotation[[11]]$value = paste0(group)
  write_yaml(annotation_edited, file.path(dir_results_script,"EMC_Chang",group_folder,"input",paste0("annotation_mega_MML_",group,".yaml")), fileEncoding = "UTF-8")
  
  # Adapt the sbatch file:
  sbatch_script = gsub("GROUP",group,sbatch_script_template)
  sbatch_script = gsub("EMC_cohort","EMC_Chang",sbatch_script )
  write_file(sbatch_script,file.path(dir_results_script,"EMC_Chang",group_folder,"input",paste0("sbatch_intogen_",group,".sh")))
  
}

#-----------------------
# 3. Process immediate outputs from intogen
#-----------------------

# Plot of methods used to find genes
upset_from_intogen_drivers = function(df,file_path,df_name){
  print("To do: define methods already!")
  methods_list <- str_split(df$METHODS, ",")
  names(methods_list) <- df$SYMBOL
  all_methods <- unique(unlist(methods_list))
  binary_df <- data.frame(
    gene = df$SYMBOL,
    stringsAsFactors = FALSE
  )
  for(method in all_methods) {
    binary_df[[method]] <- sapply(methods_list, function(x) as.numeric(method %in% x))
  }
  
  p_upset = upset(binary_df, 
                  sets = all_methods,
                  order.by = "freq",
                  decreasing = TRUE,
                  mb.ratio = c(0.6, 0.4),
                  number.angles = 0,
                  point.size = 3.5,
                  line.size = 2,
                  mainbar.y.label = "Number of Genes",
                  sets.x.label = "Genes per Method",
                  text.scale = c(1.3, 1.3, 1, 1, 2, 1.2))
  pdf(paste0(file_path,"b_02_06_A4_sing_methods_intogen_",df_name,".pdf"),width=12, height=6)
  print(p_upset)
  dev.off()
  #ggsave(paste0(file_path,"b_02_06_A4_sing_methods_intogen_",df_name,".pdf"),p_upset,width=8, height=6)
}

upset_from_intogen_potential_drivers=function(df,file_path,df_name){
  warning_sets <- c("WARNING_EXPRESSION", "WARNING_GERMLINE", "OR_WARNING","SAMPLES_3MUTS","WARNING_ARTIFACT", "KNOWN_ARTIFACT", "DRIVER", "CGC_GENE","WARNING_ENSEMBL_TRANSCRIPTS")
  binary_df <- df %>%
    select(SYMBOL, all_of(warning_sets)) %>%
    mutate(SAMPLES_3MUTS = ifelse(SAMPLES_3MUTS ==0,"False","True"))%>%
    mutate(across(all_of(warning_sets), ~ ifelse(.x == "True", 1, 0)))
  
  summary_stats <- binary_df %>%
    select(-SYMBOL) %>%
    summarise_all(sum) %>%
    pivot_longer(everything(), names_to = "Category", values_to = "Count") %>%
    arrange(desc(Count))
  
  p_upset = upset(binary_df, 
                  sets = warning_sets,
                  order.by = "freq",
                  decreasing = TRUE,
                  mb.ratio = c(0.6, 0.4),
                  number.angles = 0,
                  point.size = 4,
                  line.size = 2,
                  mainbar.y.label = "Number of Genes",
                  sets.x.label = "Genes per Category",
                  text.scale = c(1.3, 1.3, 1, 1, 2, 1.2),
                  main.bar.color = "steelblue",
                  sets.bar.color = "darkblue",
                  matrix.color = "black")
  
  pdf(paste0(file_path,"b_02_06_A4_postprocessing_intogen_",df_name,".pdf"),width=12, height=6)
  print(p_upset)
  dev.off()
}

intogen_results =data.frame()
intogen_results_unfiltered = data.frame()
for(cht in c("EMC_cohort","EMC_Chang")){
  for( group in c("Entire_cohort","Cases", "Controls","ISuppressed","ICompetent","EMC_Chang_noRDEB")){
    dir_output_intogen = file.path(dir_results_script,cht,group,"output")
    if(file.exists(file.path(dir_output_intogen,"drivers.tsv"))){
      # Read all tables and save them in a  more convenient format:
      intogen_drivers =read.table(file.path(dir_output_intogen,"drivers.tsv"),sep="\t",header=TRUE)
      write.xlsx(intogen_drivers,file.path(dir_output_intogen,"drivers.xlsx"))
      intogen_mutations = read.table(file.path(dir_output_intogen,"mutations.tsv"),sep="\t",header=TRUE)
      #intogen_unique_drivers =read.table(paste0(dir_output_intogen,"unique_drivers.tsv"),sep="\t",header=TRUE) #Note: this does not add a lot of information to the previous drivers tsv.
      intogen_unfiltered_drivers = read.table(file.path(dir_output_intogen,"unfiltered_drivers.tsv"),sep="\t",header=TRUE) # Genes are sorted out here, depending on the level of evidence. 
      write.xlsx(intogen_unfiltered_drivers,file.path(dir_output_intogen,"unfiltered_drivers.xlsx"))
      
      # Generate diagnostic plots: 
      #upset_from_intogen_drivers(intogen_drivers,dir_output_intogen ,paste0(cht,"_",group)) 
      #upset_from_intogen_potential_drivers(intogen_unfiltered_drivers,dir_output_intogen ,paste0(cht,"_",group))  
      
      # Add results to summary table:
      intogen_drivers$Cohort = cht
      intogen_drivers$Group = group
      intogen_results=rbind(intogen_results,intogen_drivers)
      
      intogen_unfiltered_drivers$Cohort = cht
      intogen_unfiltered_drivers$Group = group
      intogen_results_unfiltered = rbind(intogen_results_unfiltered,intogen_unfiltered_drivers)
    }
  }
}

write.xlsx(intogen_results,file.path(dir_results_script,"04_02_01_intogen_results_all.xlsx"))
write.xlsx(intogen_results_unfiltered,file.path(dir_results_script,"04_02_01_intogen_results_all_unfiltered.xlsx"))

#-----------------------
# 4. Check concordance with other tools
#-----------------------
intogen_results_unfiltered_expression = intogen_results_unfiltered %>%
  filter(SIGNATURE10==0,SIGNATURE9==0,WARNING_GERMLINE=="False",SAMPLES_3MUTS==0,OR_WARNING=="False",WARNING_ARTIFACT=="False",KNOWN_ARTIFACT=="False",WARNING_ENSEMBL_TRANSCRIPTS=="False",WARNING_EXPRESSION=="True",DRIVER=="False")%>%
  mutate(Exclusion_reason ="Expression_intogen")%>%as.data.frame()
intogen_results_unfiltered_3muts = intogen_results_unfiltered %>%
  filter(SIGNATURE10==0,SIGNATURE9==0,WARNING_GERMLINE=="False",SAMPLES_3MUTS!=0,OR_WARNING=="False",WARNING_ARTIFACT=="False",KNOWN_ARTIFACT=="False",WARNING_ENSEMBL_TRANSCRIPTS=="False",WARNING_EXPRESSION=="False",DRIVER=="False")%>%
  mutate(Exclusion_reason ="3muts")%>%as.data.frame()
intogen_results_unfiltered_3muts_exp = intogen_results_unfiltered %>%
  filter(SIGNATURE10==0,SIGNATURE9==0,WARNING_GERMLINE=="False",SAMPLES_3MUTS!=0,OR_WARNING=="False",WARNING_ARTIFACT=="False",KNOWN_ARTIFACT=="False",WARNING_ENSEMBL_TRANSCRIPTS=="False",WARNING_EXPRESSION=="True",DRIVER=="False")%>%
  mutate(Exclusion_reason ="3muts_and_expr")%>%as.data.frame()
intogen_results_unfiltered_drivers = intogen_results_unfiltered %>%
  filter(DRIVER=="True")%>%
  mutate(Exclusion_reason ="DRIVER")%>%as.data.frame() # Note: the drivers table only has genes that are noted as Filter=Pass
intogen_results_unfiltered_not_enough_evidence_A_M = intogen_results_unfiltered %>%
  filter(SIGNATURE10==0,SIGNATURE9==0,WARNING_GERMLINE=="False",SAMPLES_3MUTS==0,OR_WARNING=="False",WARNING_ARTIFACT=="False",KNOWN_ARTIFACT=="False",WARNING_ENSEMBL_TRANSCRIPTS=="False",WARNING_EXPRESSION=="False",DRIVER=="False")%>%
  filter(str_detect(SYMBOL, "^[A-M]")) %>% # For plotiting purporses
  mutate(Exclusion_reason ="No_big_red_flag_but_not_enough_evidence_AtoM")%>%as.data.frame() 
intogen_results_unfiltered_not_enough_evidence_N_Z = intogen_results_unfiltered %>%
  filter(SIGNATURE10==0,SIGNATURE9==0,WARNING_GERMLINE=="False",SAMPLES_3MUTS==0,OR_WARNING=="False",WARNING_ARTIFACT=="False",KNOWN_ARTIFACT=="False",WARNING_ENSEMBL_TRANSCRIPTS=="False",WARNING_EXPRESSION=="False",DRIVER=="False")%>%
  filter(str_detect(SYMBOL, "^[N-Z]")) %>% # For plotiting purporses
  mutate(Exclusion_reason ="No_big_red_flag_but_not_enough_evidence_N_to_Z")%>%as.data.frame() 

intogen_results_unfiltered_drivers$Exclusion_reason[which(intogen_results_unfiltered_drivers$FILTER!="PASS")] ="Driver_but_FAIL"
all_intogen_results_merged=rbind(intogen_results_unfiltered_expression,intogen_results_unfiltered_3muts,intogen_results_unfiltered_3muts_exp,intogen_results_unfiltered_drivers,intogen_results_unfiltered_not_enough_evidence_A_M,intogen_results_unfiltered_not_enough_evidence_N_Z)

# Load results of other tools for reference:
df_all_results_all_methods = read.xlsx(paste0("T:/Barbara/WES_analyses/Results/Cancer_gene_discovery/b02_06_summary_all_tools.xlsx"))
all_intogen_results_merged$Found_in_other_tools = NA
for(cht in c("EMC_cohort","EMC_Chang")){
  for( gr in c("Entire_cohort","Cases", "Controls","ISuppressed","ICompetent","EMC_Chang_noRDEB")){
    df_all_results_3_prev_methods = df_all_results_all_methods %>%filter(Cohort==ifelse(cht=="EMC_Chang","EMC_Darwin",cht),Group==gr)%>%as.data.frame()
    
    genes_and_tools = table(df_all_results_3_prev_methods$Gene)
    i.gr_of_interest = intersect(which(all_intogen_results_merged$Cohort==cht),which(all_intogen_results_merged$Group==gr))
    
    for(nr in c(1,2,3)){
      if(length(names(genes_and_tools[genes_and_tools ==nr]))>0){
        i.gr_of_interest_Xtool = intersect(i.gr_of_interest,
                                           which((all_intogen_results_merged$SYMBOL%in% names(genes_and_tools[genes_and_tools ==nr])) ==TRUE) )
        all_intogen_results_merged$Found_in_other_tools[i.gr_of_interest_Xtool] = nr
      }
    }
    
  }
}

# Additional exclusion criterion:
# Check if oncogenic mutations are present for genes that are deemed as oncogenic:
genes_deemed_oncogenic = all_intogen_results_merged %>%
  filter(ROLE=="Act",DRIVER=="True",QVALUE_COMBINATION<threshold)%>%pull(SYMBOL)%>%unique()

# Check which ones should be excluded:
oncogenes_no_oncokb_evidence = mega_MML %>%filter(
  Hugo_Symbol%in%genes_deemed_oncogenic)%>%
  group_by(Hugo_Symbol)%>%
  summarise(Count = sum(ONCOGENIC != "Unknown"& 	
                          MUTATION_EFFECT!="Likely Loss-of-function", na.rm = TRUE))%>%
  filter(Count==0)%>%pull(Hugo_Symbol)%>%unique()

all_intogen_results_merged=all_intogen_results_merged %>%
  mutate(Exclusion_reason=ifelse(DRIVER=="True" & SYMBOL%in%oncogenes_no_oncokb_evidence,"Excluded_oncokb_revision",Exclusion_reason))%>%
  as.data.frame()

# Plot:
exclusion_reasons <- unique(all_intogen_results_merged$Exclusion_reason)
plot_list <- list()
colors_5 <- c("Expression_intogen" = "#DC143C",
              "Driver_but_FAIL" = "#4682B4", 
              "None" = "#228B22",
              "3muts_and_expr" = "#FF8C00",
              "3muts" = "#9370DB")
colors_nr_genes <- c("0" = "#DC143C",   # Red
                     "1" = "#ADD8E6",   # Light blue
                     "2" = "#1E90FF",   # Dodger blue
                     "3" = "#0047AB")   # Cobalt blue

for(reason in exclusion_reasons) {
  plot_list[[reason]] <- all_intogen_results_merged %>%
    filter(QVALUE_COMBINATION < threshold & Exclusion_reason == reason) %>%
    mutate(Found_in_other_tools = as.character(ifelse(is.na(Found_in_other_tools),0,Found_in_other_tools)))%>%
    ggplot(aes(x = SYMBOL, y = Cohort)) +
    geom_point(aes(fill = Found_in_other_tools, color = Exclusion_reason, 
                   size = -log10(QVALUE_COMBINATION)), 
               shape = 21, stroke = 0.5) +
    facet_wrap(~Group, ncol = 1, strip.position = "left") +
    scale_fill_manual(values = colors_nr_genes, name = "Found_in_other_tools") +
    scale_color_manual(values =  colors_nr_genes, name = "Found_in_other_tools") +
    scale_size_continuous(range = c(1, 6), name = "-log10(q-value)") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text.y.left = element_text(angle = 0),
      strip.placement = "outside",
      panel.grid.major = element_line(color = "lightgray", size = 0.3),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    ) +
    labs(x = "Genes", y = "Patient Groups", 
         title = paste("Exclusion Reason:", reason))
}

# Save all plots to PDF
pdf(file.path(dir_results_script,"04_02_01_intogen_results_overview_exclusion_other3tools.pdf"), width = 18, height = 8)
for(reason in exclusion_reasons) {
  print(plot_list[[reason]])
}
dev.off()

#-----------------------
# 5. Convert results into a table that is easy to annotate:
#-----------------------
final_table_intogen= read.xlsx(file.path(dir_results_script,"04_02_01_intogen_results_all.xlsx"))

# Let's reshape it in a way that is easier to parse:
final_table_wide <- final_table_intogen %>%
  mutate(Excluded_oncokb_criterion = ifelse(SYMBOL%in%oncogenes_no_oncokb_evidence,"Excluded","PASS"))%>%
  mutate(Excluded_combination = ifelse(METHODS=="combination","Excluded","PASS"))%>%
  filter(QVALUE_COMBINATION<threshold)%>%
  # Create a combined identifier for Cohort-Group
  mutate(Cohort_Group = paste0(Cohort, "_", Group)) %>%
  arrange(SYMBOL, Cohort_Group != "EMC_Chang_Entire_cohort") %>% # This ensures that the annotations that prevail are of the combined cohort (entire cohort)
  # For each SYMBOL, collapse the non-varying columns
  group_by(SYMBOL) %>%
  mutate(
    # Keep first value for annotation columns (assuming they're consistent per gene)
    METHODS_all = first(METHODS),
    ROLE_all = first(ROLE),
    CGC_GENE_all = first(CGC_GENE),
    CGC_CANCER_GENE_all = first(CGC_CANCER_GENE)
  ) %>%
  ungroup() %>%
  
  # Select columns for pivoting
  select(SYMBOL, Cohort_Group, QVALUE_COMBINATION, 
         METHODS_all, ROLE_all, CGC_GENE_all, CGC_CANCER_GENE_all,Excluded_oncokb_criterion,Excluded_combination) %>%
  
  # Pivot Q-values to wide format
  pivot_wider(
    names_from = Cohort_Group,
    values_from = QVALUE_COMBINATION,
    names_prefix = "QVALUE_"
  ) %>%
  
  # Rename annotation columns
  rename(
    METHODS = METHODS_all,
    ROLE = ROLE_all,
    CGC_GENE = CGC_GENE_all,
    CGC_CANCER_GENE = CGC_CANCER_GENE_all
  )
#write.xlsx(final_table_wide,file.path(dir_results_script,"04_02_01_intogen_results_all_with_man_revision.xlsx"))

#-----------------------
# 6. Loliplots for all genes:
#-----------------------
# Genes of interest:
driver_genes_oficial_list = final_table_intogen %>%
  mutate(Excluded_oncokb_criterion = ifelse(SYMBOL%in%oncogenes_no_oncokb_evidence,"Excluded","PASS"))%>%
  mutate(Excluded_combination = ifelse(METHODS=="combination","Excluded","PASS"))%>%
  dplyr::filter(QVALUE_COMBINATION<threshold,COHORT=="EMC_Chang",Excluded_oncokb_criterion=="PASS",Excluded_combination=="PASS")%>%
  pull(SYMBOL)%>%unique()
  
driver_genes=intogen_results_unfiltered%>%
  mutate(Excluded_oncokb_criterion = ifelse(SYMBOL%in%oncogenes_no_oncokb_evidence,"Excluded","PASS"))%>%
  mutate(Excluded_combination = ifelse(SIG_METHODS=="combination","Excluded","PASS"))%>%
  dplyr::filter(QVALUE_COMBINATION<threshold,Excluded_oncokb_criterion=="PASS",Excluded_combination=="PASS",DRIVER=="True",!FILTER%in%c("Warning expression","Samples with more than 2 mutations","Lack of literature evidence"))%>%pull(SYMBOL)%>%unique()
rescued_genes = c("FAT2") #We considered rescuing these genes as they were filtered out by intogen without any good reason
final_genes_of_interest =intogen_results_unfiltered%>%
  dplyr::filter(SYMBOL%in%c(driver_genes,rescued_genes),QVALUE_COMBINATION<threshold)%>%as.data.frame()

# Assign role to genes:
final_genes_of_interest = merge(final_genes_of_interest,pathogenic_mutations_shain_lab[,c("Gene","Type")],by.x="SYMBOL",by.y="Gene",all.x=TRUE)
final_genes_of_interest$Type[is.na(final_genes_of_interest$Type)] = final_genes_of_interest$ROLE[is.na(final_genes_of_interest$Type)]
final_genes_of_interest$Type[which(final_genes_of_interest$Type=="Act")]="Oncogene"
final_genes_of_interest$Type[which(final_genes_of_interest$Type=="LoF")]="TSP"

#Patient groups:
isup_patients = unique(mega_MML_IS$Tumor_Sample_Barcode)
icomp_patients = unique(mega_MML_IC$Tumor_Sample_Barcode)

cases_pts = unique(mega_MML_Cases$Tumor_Sample_Barcode)
control_pts = unique(mega_MML_Controls$Tumor_Sample_Barcode)

# Information needed for lolliplots:
edbx <- EnsDb.Hsapiens.v75
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Loliplots for each candidate gene:

make_lolliplot = function(gene_symbol,mega_MML){
  gene_id <- get(gene_symbol, org.Hs.egSYMBOL2EG)
  gene_ensb = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                    filters = "hgnc_symbol",
                    values = gene_symbol,
                    mart = mart)[,1]

  
  # sample track:
  mml_gene = mega_MML %>%dplyr::filter(Hugo_Symbol==gene_symbol,!Variant_Classification%in%c("Intron","3'UTR","5'UTR","5'Flank"))%>%arrange(Start_Position)%>%as.data.frame()
  mml_gene$aa_pos = as.numeric(str_extract(mml_gene$Protein_Change, "\\d+"))
  
  #TODO: revise with Hunter/Harsh:
  for(i in which(is.na(mml_gene$aa_pos))){
    if(i!=1 & i!=nrow(mml_gene)){
      prev <- max(mml_gene$aa_pos[1:(i-1)], na.rm = TRUE)      # last known position before NA
      next_pos <- min(mml_gene$aa_pos[(i+1):nrow(mml_gene)], na.rm = TRUE)  # next known position after NA
      mml_gene$aa_pos[i] <- round(median(c(prev, next_pos)))       # fill with median
      mml_gene$Protein_Change[i] <- paste0("SS_",round(median(c(prev, next_pos))))       # fill with median
      
    }else if(i==nrow(mml_gene)){
      prev <- max(mml_gene$aa_pos[1:(i-1)], na.rm = TRUE)      # last known position before NA
      mml_gene$aa_pos[i] <- prev+1
      mml_gene$Protein_Change[i] <- paste0("SS_last")       # fill with median
    }
    else if(i==1){
      next_pos <- min(mml_gene$aa_pos[(i+1):nrow(mml_gene)], na.rm = TRUE)  # next known position after NA
      mml_gene$aa_pos[i] <-next_pos-1
      mml_gene$Protein_Change[i] <- paste0("SS_first")       # fill with median
    }
  }
  
  # Add score based on the number of samples with such mutation:
  mutation_counts <- mml_gene %>%
    group_by(Protein_Change) %>%
    summarise(score = n()) %>%
    as.data.frame()
  
  mml_gene = mml_gene  %>%          # optional: remove rows with NA positions first
    distinct(Protein_Change, .keep_all = TRUE) %>%as.data.frame()
  sample.gr <- GRanges(
    seqnames = gene_symbol,
    IRanges(names = mml_gene$Protein_Change,start = mml_gene$aa_pos, width = 1),
    mutation = mml_gene$Protein_Change,
    type = mml_gene$Variant_Classification
  )
  
  # change color of mutations
  sample.gr$color <- unlist(lapply(1:length(sample.gr$type), function(z){
    if(sample.gr$type[z]== "Missense_Mutation") return("blue")
    if(sample.gr$type[z]== "Silent") return("yellow")
    if(sample.gr$type[z] == "DE_NOVO_START_OUT_FRAME") return("black")
    if(sample.gr$type[z] == "Nonsense_Mutation") return("#9400d3")
    if(sample.gr$type[z] == "Splice_Site") return("#9400d3")
    if(sample.gr$type[z] == "START_CODON_SNP") return("black")
    if(sample.gr$type[z] == "Frame_Shift_Ins") return("black")
    if(sample.gr$type[z] == "Frame_Shift_Del") return("black")
  }))
  
  
  
  sample.gr$score = mutation_counts$score[match(sample.gr$mutation, mutation_counts$Protein_Change)]
  
  #gnm <- GRanges(unique(mml_gene$Chromosome), IRanges(start = mml_gene$Start_Position))
  #gnm_prt <- genomeToProtein(gnm, edbx)
  #prt <- proteins(edbx, filter = ProteinIdFilter(names(gnm_prt[[3]])))
  #gnm$score <- sample.int(5, length(gnm), replace = TRUE)
  #gnm$feature.height = .5
  # gene track:
  #gene_positions <- geneTrack(get(gene_symbol, org.Hs.egSYMBOL2EG), TxDb.Hsapiens.UCSC.hg19.knownGene)[[1]]
  
  
  
  
  # Get protein data
  APIurl <- "https://www.ebi.ac.uk/proteins/api/" # base URL of the API
  taxid <- "9606" # human tax ID
  gene <- gene_symbol # target gene
  orgDB <- "org.Hs.eg.db" # org database to get the uniprot accession id
  eid <- mget(gene_symbol, get(sub(".db", "SYMBOL2EG", orgDB)))[[1]]
  chr <- mget(eid, get(sub(".db", "CHR", orgDB)))[[1]]
  accession <- unlist(lapply(eid, function(.ele){
    mget(.ele, get(sub(".db", "UNIPROT", orgDB)))
  }))
  #stopifnot(length(accession)<=20) # max number of accession is 20
  
  featureURL <- paste0(APIurl, 
                       "features?offset=0&size=-1&reviewed=true",
                       "&types=DNA_BIND%2CMOTIF%2CDOMAIN",
                       "&taxid=", taxid,
                       "&accession=", paste(accession, collapse = "%2C")
  )
  response <- GET(featureURL)
  if(!http_error(response)){
    content <- httr::content(response)
    if(length(content)!=0){
      content <- content[[1]]
      acc <- content$accession
      sequence <- content$sequence
      gr <- GRanges(gene_symbol, IRanges(1, nchar(sequence)))
      domains <- do.call(rbind, content$features)
      domains <- GRanges(gene_symbol, IRanges(as.numeric(domains[, "begin"]),
                                              as.numeric(domains[, "end"]))) # removed, names=domains[,"description"] to make plots less crowded
      #names(domains)[1] <- "DNA_BIND" ## this is hard coding.
      domains$fill <- 1+seq_along(domains)
      domains$height <- 0.04
      
      lolliplot(sample.gr, domains, ranges = gr,score = sample.gr$score,legend = NULL)
      grid.text(gene_symbol, x=.5, y=.9, just="top", 
                gp=gpar(cex=1.5, fontface="bold"))
    }else{
      message("Can not get variations. http error")
      # We get the data from ensembl
      prot_data <- getBM(attributes = c("ensembl_transcript_id", "ensembl_peptide_id", "peptide"),
                         filters = "ensembl_gene_id",
                         values = gene_ensb,
                         mart = mart)
      
      prot_data$protein_length <- nchar(prot_data$peptide)
      prot_data=prot_data[which(prot_data$ensembl_peptide_id!=""),]
      prot_data=prot_data[which.max(prot_data$protein_length), ]
      protein_length = prot_data$protein_length
      protein_id <- prot_data$ensembl_peptide_id
      # If domains cannot be fetched, we can still plot the lolliplot without the domains:
      features <- GRanges(gene_symbol,
                          IRanges(start = 1, end = protein_length),
                          feature = "protein")
      
      lolliplot(sample.gr, features = features,legend = NULL)
      grid.text(gene_symbol, x=.5, y=.9, just="top", 
                gp=gpar(cex=1.5, fontface="bold"))
    }
  }
}

stripchart_zygosity = function(gene_symbol,mml,pathogenic=NULL){
  
  mml_gene = mml %>%
    dplyr::filter(Hugo_Symbol==gene_symbol,!Variant_Classification%in%c("Intron","3'UTR","5'UTR","5'Flank","Silent"))%>%
    as.data.frame()
  if(!(is.null(pathogenic))){
    mml_gene = mml_gene %>%dplyr::filter(Pathogenic=="Probably")%>%as.data.frame()
    ylab_path = "(Pathogenic mutations)"
  }else{
    ylab_path = "(All mutations)"
  }
  mml_gene$Adjusted_MAF_hits = ifelse(mml_gene$Adjusted_MAF>1,1,mml_gene$Adjusted_MAF)*2
  
  sample_counts <- mml_gene %>%
    count(Tumor_Sample_Barcode, name = "mut_count")
  
  mml_gene <- mml_gene %>%
    left_join(sample_counts, by = "Tumor_Sample_Barcode") %>%
    mutate(shape_group = ifelse(mut_count > 1, "Multiple", "Single"))
  
  mml_gene <- mml_gene %>%
    mutate(
      group = ifelse(grepl("_1$", Tumor_Sample_Barcode), "Case", "Control"),
      shape_group = factor(shape_group, levels = c("Single", "Multiple"))
    )
  
  mml_gene <- mml_gene %>%
    arrange(group, shape_group, Tumor_Sample_Barcode) %>%
    mutate(Tumor_Sample_Barcode = factor(Tumor_Sample_Barcode, levels = unique(Tumor_Sample_Barcode)))
  
  # Find where to place separator (between last Control and first Case)
  group_positions <- mml_gene %>%
    distinct(Tumor_Sample_Barcode, group) %>%
    arrange(Tumor_Sample_Barcode) %>%
    pull(group)
  
  sep_pos <- sum(group_positions == group_positions[1]) + 0.5
  
  # Make stripchart
  p <- ggplot(mml_gene, aes(x = Tumor_Sample_Barcode, y = Adjusted_MAF_hits, shape = shape_group,col=group)) +
    geom_point(size = 3, width = 0.15, alpha = 0.8) +
    scale_shape_manual(values = c(Single = 16, Multiple = 17)) +
    geom_vline(xintercept = sep_pos, linetype = "dashed", color = "grey50")+
    theme_bw() +
    geom_hline(yintercept = c(0.75, 1.25), linetype = "dashed", color = "grey50", alpha = 0.6) +
    labs(
      title = paste0(gene_symbol, " – Adjusted MAF per sample ", ylab_path),
      x = "Sample",
      y = "Adjusted MAF",
      shape = "Mutations per sample"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", size = 14),
      plot.margin = margin(10, 10, 10, 10)
    )+ylim(c(0,2.1))
  
  
  return(p)
}

boxplot_freq <-function(gene_symbol,mml,s_summary,variable_of_interest,pathogenic=NULL){
  mml = merge(mml,s_summary[,c("Sample.ID",variable_of_interest)],by.x="Tumor_Sample_Barcode",by.y="Sample.ID",all.x=TRUE)
  
  if(!(is.null(gene_symbol))){
    mml_gene = mml %>%
      dplyr::mutate(Tumor_Sample_Barcode=factor(Tumor_Sample_Barcode))%>%
      dplyr::filter(Hugo_Symbol==gene_symbol,!Variant_Classification%in%c("Intron","3'UTR","5'UTR","5'Flank","Silent"))%>%
      as.data.frame()
  }else{
    mml_gene = mml %>%
      dplyr::mutate(Tumor_Sample_Barcode=factor(Tumor_Sample_Barcode))%>%
      dplyr::filter(!Variant_Classification%in%c("Intron","3'UTR","5'UTR","5'Flank","Silent"))%>%
      as.data.frame()
  }
  if(!(is.null(pathogenic))){
    mml_gene = mml_gene %>%dplyr::filter(Pathogenic=="Probably")%>%as.data.frame()
    ylab_path = "% of samples with pathogenic mutation"
  }else{
    ylab_path = "% of samples with mutation"
  }
  # Check whether each sample has any mutation in the gene
  samples_with_mutation <- unique(mml_gene$Tumor_Sample_Barcode)
  
  # Add mutation status to all samples
  s_summary_complete <- s_summary %>%
    dplyr::mutate(
      has_mutation = ifelse(Sample.ID %in% samples_with_mutation, "Mutated", "Wild-type")
    )
  
  freq_summary <- s_summary_complete %>%
    dplyr::group_by(!!sym(variable_of_interest)) %>%
    dplyr::summarise(
      n_total = n(),
      n_mutated = sum(has_mutation == "Mutated"),
      freq_mutated = (n_mutated / n_total) * 100,
      .groups = 'drop'
    )
  
  # define color palette conditionally
  if (variable_of_interest == "Metastasis") {
    my_colors <- c( "#F8766D","#00BFC4")  # ggplot default two colors
  } else if (variable_of_interest == "IS.at.cSCC") {
    my_colors <- c("gray30", "darkred")
  } else {
    my_colors <- c("gray70", "gray30")  # fallback
  }
  
  p <- ggplot(freq_summary, aes(x = !!sym(variable_of_interest), y = freq_mutated,fill=!!sym(variable_of_interest))) +
    geom_bar(stat = "identity", alpha = 0.7) +
    scale_fill_manual(values = my_colors) +
    geom_text(aes(label = paste0(round(freq_mutated, 1), "%\n(", n_mutated, "/", n_total, ")")),
              vjust = -0.5, size = 3.5) +
    theme_bw() +
    labs(
      title = paste0(gene_symbol, " mutation frequency"),
      x = variable_of_interest,
      y = ylab_path
    ) +
    ylim(0, max(freq_summary$freq_mutated) * 1.15)
  
  return(p)
}

dir.create(file.path(dir_results_script,"Gene_investigation"))
for(gn in c(driver_genes,rescued_genes)){
  dir.create(file.path(dir_results_script,"Gene_investigation",gn))
  print(gn)
  pdf(file.path(dir_results_script,"Gene_investigation",gn,paste0(gn,".pdf")),width=12,height=3)
  make_lolliplot(gn,mega_MML_clonal)
  dev.off()
  pdf(file.path(dir_results_script,"Gene_investigation",gn,paste0(gn,"_stripchart.pdf")),width=12,height=3)
  if(gn %in%c("LAMC1","MLH1")){
    stripchart_zygosity(gn, mega_MML_clonal)
  }else{
    stripchart_zygosity(gn, mega_MML_clonal,pathogenic = "Yes")
  }
  
  dev.off()
}

# Save all information on a slide deck:
# Make it into a powerpoint slide deck:
ppt <- read_pptx()
tsps_without_path = setdiff(unique(final_genes_of_interest$SYMBOL[which(final_genes_of_interest$Type=="TSP")]), mega_MML$Hugo_Symbol[which(mega_MML$Pathogenic=="Probably")])
for (gn in c(driver_genes,rescued_genes)) {
  print(gn)
  
  # Add slide
  ppt <- ppt %>%
    add_slide(layout = "Title and Content", master = "Office Theme") %>%
    ph_with(value = ifelse(is.null(gn),"All genes",gn), location = ph_location_type(type = "title"))
  
  if(gn %in% tsps_without_path ){
    mega_MML_clonal_altered = mega_MML_clonal %>%
      dplyr::filter(Hugo_Symbol==gn)%>%as.data.frame()
    mega_MML_clonal_altered$Pathogenic =
      ifelse( mega_MML_clonal_altered$Variant_Classification%in%c("Missense_Mutation","Nonsense_Mutation","Splice_Site"), "Probably","Unlikely")
    print("Pathogenic mutations were added to this gene")
    pathogenic_label= "Yes"
  }else{
    mega_MML_clonal_altered = mega_MML_clonal
    pathogenic_label = "Yes"
  }
  met_bar_obj<- rvg::dml(code = print(boxplot_freq(gn, mega_MML_clonal_altered ,sample_summary,"Metastasis",pathogenic = pathogenic_label)))
  is_bar_obj<- rvg::dml(code = print(boxplot_freq(gn, mega_MML_clonal_altered ,sample_summary,"IS.at.cSCC",pathogenic = pathogenic_label)))
  
  
  if(!is.null(gn)){
    # Safely try to create plot
    lolliplot_obj <- rvg::dml(code = make_lolliplot(gn, mega_MML_clonal_altered))
    stripchart_obj <- rvg::dml(code = print(stripchart_zygosity(gn, mega_MML_clonal_altered,pathogenic = pathogenic_label)))
    
    
    # Add lolliplot on top left
    ppt <- ppt %>%
      ph_with(value = lolliplot_obj, 
              location = ph_location(left = 0.5, top = 1.2, width = 6, height = 3))
    
    # Add stripchart on bottom left
    ppt <- ppt %>%
      ph_with(value = stripchart_obj, 
              location = ph_location(left = 0.5, top = 4.4, width = 6, height = 3))
    
    # Add barplot Metastasis on top right
    ppt <- ppt %>%
      ph_with(value = met_bar_obj, 
              location = ph_location(left = 7, top = 1.2, width = 3, height = 3))
    
    # Add barplot Immunosuppresion on bottom right
    ppt <- ppt %>%
      ph_with(value = is_bar_obj, 
              location = ph_location(left = 7, top = 4.4, width = 3, height = 3))
    
  }else{
    # Add barplot Metastasis on right
    ppt <- ppt %>%
      ph_with(value = met_bar_obj, 
              location = ph_location(left = 0.5, top = 1.2, width = 4, height = 3))
    
    # Add barplot Immunosuppresion on left
    ppt <- ppt %>%
      ph_with(value = is_bar_obj, 
              location = ph_location(left = 5, top = 1.2, width = 4, height = 3))
    
  }
  
}


# Save the PowerPoint
output_file <- file.path(dir_results_script,"Gene_investigation","04_02_01_all_genes_with_stripcharts_probable_pathogenic.pptx")
print(ppt, target = output_file)


# Record oncokb mutations that are not in genes identified by intogen:
mega_MML_filtered = mega_MML%>%filter(Clonal=="Yes")%>%as.data.frame()

final_selection = read.xlsx(file.path(dir_results_script,"04_02_01_intogen_results_all_with_man_revision_v2.xlsx"),sheet="Intogen_results")
final_selection=final_selection[-1,]
intogen_drivers = final_selection[,1][which(final_selection$X20=="Yes")]

# ONCOKB mutations:
oncokb_oncogenic = mega_MML_filtered %>%filter(ONCOGENIC%in%c("Oncogenic","Likely Oncogenic"))%>%
  mutate(Protein_Change = if_else(
    Hugo_Symbol == "TERT",
    as.character(Start_Position),
    Protein_Change
  ))%>%
  select(Hugo_Symbol,Variant_Classification,Protein_Change,COSMIC_overlapping_mutations,ONCOGENIC,MUTATION_EFFECT,Pathogenic,Pathogenic_notes,MUTATION_EFFECT_DESCRIPTION)%>%
  group_by(Hugo_Symbol,Protein_Change)%>%
  mutate(count = n())%>%distinct(Hugo_Symbol, Protein_Change, .keep_all = TRUE)

#Add Q-value from cancer hotspots database:
hotspots_db = read.xlsx("T:/Barbara/WES_analyses/Methods/Pathogenic_mutations/cancerhotspotsdb/hotspots_v2.xlsx")
hotspots_db$Protein_Change= paste0("p.",gsub("\\:.*","",hotspots_db$Reference_Amino_Acid),hotspots_db$Amino_Acid_Position,gsub("\\:.*","",hotspots_db$Variant_Amino_Acid))

oncokb_oncogenic=merge(oncokb_oncogenic,hotspots_db[,c("Hugo_Symbol","Protein_Change","qvalue")],by=c("Hugo_Symbol","Protein_Change"),all.x=TRUE)
oncokb_oncogenic$qvalue = as.numeric(oncokb_oncogenic$qvalue )

# Exclude mutations that are not found by intogen:
oncokb_oncogenic=oncokb_oncogenic %>%
  mutate(Intogen_driver= if_else(Hugo_Symbol %in%intogen_drivers,"Yes","No"))%>%
  as.data.frame()


write.xlsx(oncokb_oncogenic ,file.path(dir_results_script,"04_02_01_oncokb.xlsx"))
