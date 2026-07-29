#---------------------------------------------------
# Aim: Mutational signature analysis
# Author: B. Rentroia Pacheco
# Input: MM lists
# Output: Muattional signatures
#---------------------------------------------------

#-----------------------
# 0. Setup
#-----------------------

if(IS_classification=="IS.at.cSCC"){
  #-----------------------
  # 1. Setup
  #-----------------------
  #mega_MML_TiN_analysis = read.xlsx(paste0("T:/Barbara/WES_analyses/Results/MML/Tumor_only_analyses/b_02_01_10_Mega_MML_tonly_nonly_comparison_MasterMutationList.xlsx"))
  #gr_analysis ="Somatic_Rescued_Clonal_AFTER_rescuing"
  #sample_summary = read.xlsx("T:/Barbara/WES_analyses/Results/Meta_data/DA_04_Sample_summary_extended_version_incl_exclusion.xlsx")
  print("Only clonal mutations are used for mutational signatures")
  mega_MML_selected = mega_MML %>%filter(Clonal=="Yes")%>%as.data.frame()
  
  #-----------------------
  # 2. Prepare input matrix to SigProfiler
  #-----------------------
  dir.create( file.path(dir_results_script,"Data"))
  # Convert into MAF format:
  laml <- read_maf(maf = mega_MML_selected)
  
  # Tally components
  # BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
  
  # Build matrix fo SBS most common 96 components
  mt_tally <- sig_tally(laml, 
                        ref_genome = "BSgenome.Hsapiens.UCSC.hg19", 
                        genome_build = "hg19",
                        useSyn = TRUE)
  str(mt_tally$all_matrices, max.level = 1)
  
  out <- mt_tally$nmf_matrix
  out <-as.data.frame(t(out))
  colnames(out) = paste0("S",colnames(out))
  out[, "MutationType"]= rownames(out)
  #Need to reshuffle the columns:
  out = out[,c("MutationType",setdiff(colnames(out),"MutationType"))]
  rownames(out)=NULL
  sig_matrix_file_name = file.path(dir_results_script,"Data",paste0("S04_01_02_mutational_sig_matrix.txt"))
  write.table(out, file=sig_matrix_file_name, sep='\t', quote=F, row.names=F)
  
  #-----------------------
  # 3. Run SigProfilerAssignment in the terminal, use python 3.9!
  #-----------------------
  #import SigProfilerAssignment as spa
  #from SigProfilerAssignment import Analyzer as Analyze
  sig_matrix_file_name_command = gsub("T:/","",sig_matrix_file_name)
  sig_profiler_command='Analyze.cosmic_fit(samples="/Volumes/shainlab/sig_matrix_file_name_command", 
                     output="/Volumes/shainlab/sig_matrix_file_name_commandSigProfilerAssignment/ExomeRN_exome_parameter/gr_analysis",
                     exome = exome_parameter,
                     input_type="matrix",
                     genome_build="GRCh37",
                     cosmic_version=3.4)'
  # Create directories:
  dir.create(file.path(dir_results_script,"Data","SigProfilerAssignment"))
  dir.create(file.path(dir_results_script,"Data","SigProfilerAssignment","ExomeRN_True"))
  dir.create(file.path(dir_results_script,"Data","SigProfilerAssignment","ExomeRN_True",gr_analysis))
  for (ex_par in c("False","True")){
    cat(gsub("S04_01_02_mutational_sig_matrix.txtSig","Sig",
             gsub("sig_matrix_file_name_command",sig_matrix_file_name_command,
                  gsub("exome_parameter",ex_par,
                       gsub("gr_analysis",gr_analysis,sig_profiler_command)))))
    cat("\n")
    
  }
}
#-----------------------
# 4. Analyze results
#-----------------------
# Find sample order according to mutational tumor burden:
# Create group and arrange
sample_summary_plot <- sample_summary_plot%>%
  mutate(Sample_Id_S = paste0("S",Sample.ID))%>%as.data.frame()

for(exome_par in c("TRUE","FALSE")){
  for(gr_analysis in c(gr_analysis)){
    # Delete signatures that are not present in the dataset:
    sig_results = read.csv(file.path(dir_results_script,"Data","SigProfilerAssignment",paste0("ExomeRN_",exome_par),gr_analysis,"Assignment_Solution","Activities","Assignment_Solution_Activities.txt"),sep="\t")
    rownames(sig_results) = sig_results$Samples
    sig_results$Samples = NULL
    
    # Remove signatures that are not present in the dataset
    sig_results = sig_results[,colSums(sig_results)!=0]
    
    signature_names = c("Other","SBS58","SBS14", "SBS20","SBS15","SBS6", "SBS5","SBS1","SBS13","SBS2","SBS32","SBS87","SBS7d","SBS7c","SBS7a","SBS7b","SBS40b","SBS11","SBS31","SBS30","SBS10b","SBS23","SBS54","SBS19")
    signature_palette = c("lightgray","#696969", "#2E8B57","#708238", "#32CD32","#3F704D","#5EA1C7","#3F5AE2","#008080","#40BFBF","#C8A2C8","#66D9CC","#F6D45B", "#FCF0C2","#D99C39","#A64C4C","#DB7093","#A398CC","#581845","#F4A6A6","#C2B280","#2C3E50","#A9A9A9","lightgray")
    names(signature_palette) = signature_names
    
  
    # Calculate proportions:
    sig_results_props = sig_results %>%
      mutate(across(everything())/rowSums(across(everything())))%>%
      t()%>%
      as.data.frame()%>%
      mutate(Signature = rownames(.))%>%
      select(Signature, everything())
    
    sbs_data_long <- reshape2::melt(sig_results_props)%>% filter(value!=0)%>%
      mutate(variable = factor(as.character(variable),levels=paste0("S",sample_summary_plot$Sample.ID)),
             gr_analysis = gr_analysis)%>%
      as.data.frame()
  
    
    # Merge less frequent signatures:
    #subset_signatures = names(sort(rowMeans(sig_results_props[,-1]),decreasing=TRUE)[1:10])
    sbs_data_long =  sbs_data_long  %>% 
      mutate(Signature_factor = fct_lump(Signature, n = 9,other_level = "Other")) %>%
      #mutate(Signature_factor =ifelse(Signature%in%subset_signatures,Signature, "Other")) %>%
      as.data.frame()
    # Make sure that color scheme is consistent:
    sign_present = as.character(unique(sbs_data_long$Signature_factor))
    sbs_data_long = sbs_data_long %>%
      mutate(Signature_factor=factor(Signature_factor,levels=signature_names[signature_names%in%sign_present]))%>%
      as.data.frame()
    
    # Find patients that are classified as having signature SBS32:
    pts_with_SBS_32=sbs_data_long%>%filter(Signature=="SBS32")%>%pull(variable)%>%as.character()
    write.csv(pts_with_SBS_32,file.path(dir_results_script,"Data","SigProfilerAssignment",paste0("ExomeRN_",exome_par),gr_analysis,"/S04_01_02_SBS32_patients.csv"))
    sample_summary_plot$IS.at.cSCC.updatedS32 = ifelse(sample_summary_plot$Sample_Id_S%in%pts_with_SBS_32,"Yes",sample_summary_plot$IS.at.cSCC)
    
    # Now with panels - Immunosuppressed according to clinical records:
    sbs_data_long_wclin = merge(sbs_data_long,sample_summary_plot[,c("Sample_Id_S","Metastasis","IS.at.cSCC","OTR.at.cSCC","IS.at.cSCC.updatedS32")],by.x="variable",by.y="Sample_Id_S",all.x=TRUE)
    
    # We might want to specify signatures of interest
    sbs_plot <- ggplot(sbs_data_long, aes(fill = Signature_factor, y = value, x = variable))+ 
      geom_bar(position = "fill", stat = "identity", width = 0.95)+
      scale_fill_manual(values=signature_palette[signature_names%in%sign_present])+
      theme_minimal() +
      theme_bw() +ylab("Signature proportion")+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      xlab("")+ theme(axis.ticks.x = element_blank(),  axis.text.x = element_blank())+ggtitle(paste0(gr_analysis," mutations"))
    ggsave(file.path(dir_results_script,"Data","SigProfilerAssignment",paste0("ExomeRN_",exome_par),gr_analysis,paste0("/S04_01_02_Mut_Sig_proportion",IS_classification,".pdf")),sbs_plot,width=16,height=8)
    
    write.xlsx(sbs_data_long_wclin,file.path(dir_results_script,"Data","SigProfilerAssignment",paste0("ExomeRN_",exome_par),gr_analysis,paste0("/S04_01_02_Mut_Sig_proportion_",IS_classification,"_source_data.xlsx")))
    for(IS_pts in c("IS.at.cSCC","IS.at.cSCC.updatedS32")){
      sbs_plot_facet = ggplot(sbs_data_long_wclin, aes(fill = Signature_factor, y = value, x = variable))+ 
        geom_bar(position = "fill", stat = "identity", width = 0.95)+
        scale_fill_manual(values=signature_palette[signature_names%in%sign_present])+
        theme_minimal() +
        theme_bw() +ylab("Signature Proportion")+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        xlab("")+ theme(axis.ticks.x = element_blank(),  axis.text.x = element_blank())+ggtitle(paste0(gr_analysis," mutations"))+facet_wrap(~Metastasis+.data[[IS_pts]],scales = "free")
      ggsave(file.path(dir_results_script,"Data","SigProfilerAssignment",paste0("ExomeRN_",exome_par),gr_analysis,paste0("/S04_01_02_Mut_Sig_proportion_clinical_facets_",IS_pts,".pdf")),sbs_plot_facet,width=16,height=10)
      
    }
      
  }
}
