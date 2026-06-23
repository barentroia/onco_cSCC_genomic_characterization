#---------------------------------------------------
# Aim: Create CNV folders for easy analysis
# Author: B. Rentroia Pacheco
# Input: CNV seg files
# Output: New organization folders
#---------------------------------------------------

#-----------------------
# 1. Create separate folders with all cns files:
#-----------------------
dir.create(file.path(dir_results_script,"cns"))
for(incl_samples_description in c("All_CC_but_tcel_0","All_CC_but_tcel_0.135","All_CC_but_tcel_0.268")){
  cnr_outputs = ".*deduplicated.cns$"
  dir_cnv_res = file.path(dir_results_script,"cns",incl_samples_description,"/")
  dir.create(dir_cnv_res)
  dir.create(paste0(dir_cnv_res,"Tumor/"))
  dir.create(paste0(dir_cnv_res,"Normal/"))
  
  # Sample exclusion
  accepted_samples = samples_cutoff[[incl_samples_description]]
  
  for (material_type in c("Tumor","Normal")){
    for(smp_id in unique(sample_summary$Sample.ID)){
      if(smp_id%in%accepted_samples){
        dir_final_batch = paste0("T:/Barbara/EMC_cSCC/Seq_Batches_ALL/",smp_id,"/CopyNumber/",material_type,"/")
        files_folder=list.files(dir_final_batch)
        
        # rename file: 
        #coverage_file = files_folder[grepl("fastp_coverage",files_folder)]
        #correct_prefix = gsub("_MasterMutationList.*","",sub("_1","",gsub("_fastp","_1_fastp",mml_file))) 
        
        cnv_file = files_folder[grepl(cnr_outputs,files_folder)]
        
        if(length(cnv_file)==0){print(smp_id)}
        #print(cnv_file)
        #print(paste0(paste0(dir_cnv_res,"/",material_type,"/",cnv_file_new)))
        #print(paste0(dir_final_batch,cnv_file))
        #print(paste0(paste0(dir_cnv_res,material_type,"/",cnv_file)))
        file.copy(paste0(dir_final_batch,cnv_file),paste0(paste0(dir_cnv_res,material_type,"/",cnv_file)))
      }
    }
  }
}

# Generate separate cns for cases and controls:
immuno_sup = sample_summary$Sample.ID[which(sample_summary[,IS_classification]=="Yes")]
cnr_outputs = ".*deduplicated.cns$"
for(incl_samples_description in c("All_CC_but_tcel_0","All_CC_but_tcel_0.135","All_CC_but_tcel_0.268")){
  dir_cnv_res = file.path(dir_results_script,"cns",incl_samples_description,"/")
  dir.create(paste0(dir_cnv_res,"/Controls/"))
  dir.create(paste0(dir_cnv_res,"/Cases/"))
  dir.create(paste0(dir_cnv_res,"/Immunosuppressed/"))
  dir.create(paste0(dir_cnv_res,"/Immunocompetent/"))
  
  # Sample exclusion
  accepted_samples = samples_cutoff[[incl_samples_description]]
  
  for(fld in c("Cases","Controls","Immunosuppressed","Immunocompetent")){
    dir.create(paste0(dir_cnv_res,"/",fld,"/CNS_files"))
    
    if(fld %in%c("Immunosuppressed","Immunocompetent")){
      dir.create(paste0(dir_cnv_res,"/",fld,"/Cases"))
      dir.create(paste0(dir_cnv_res,"/",fld,"/Controls"))
    }
  } 
  
  for(smp_id in sample_summary$Sample.ID){
    if(smp_id%in%accepted_samples){
      # Get the file:
      dir_final_batch = paste0("T:/Barbara/EMC_cSCC/Seq_Batches_ALL/",smp_id,"/CopyNumber/Tumor/")
      files_folder=list.files(dir_final_batch)
      cnv_file = files_folder[grepl(cnr_outputs,files_folder)]
      
      if(!file.exists(paste0(dir_cnv_res,"Tumor/",cnv_file))){print(smp_id)
      }else{
        # Save it in the appropriate group:
        # Case vs control
        if(grepl("_0",smp_id)){
          file.copy(paste0(dir_cnv_res,"Tumor/",cnv_file),paste0(paste0(dir_cnv_res,"/Controls/CNS_files/",cnv_file)))
        }else{
          file.copy(paste0(dir_cnv_res,"Tumor/",cnv_file),paste0(paste0(dir_cnv_res,"/Cases/CNS_files/",cnv_file)))
        }
        
        # Immunosuppressed vs immunocompetent:
        if(smp_id %in% immuno_sup){
          file.copy(paste0(dir_cnv_res,"Tumor/",cnv_file),paste0(paste0(dir_cnv_res,"/Immunosuppressed/CNS_files/",cnv_file)))
          
          if(grepl("_0",smp_id)){
            file.copy(paste0(dir_cnv_res,"Tumor/",cnv_file),paste0(paste0(dir_cnv_res,"/Immunosuppressed/Controls/CNS_files/",cnv_file)))
            
            
          }else{
            file.copy(paste0(dir_cnv_res,"Tumor/",cnv_file),paste0(paste0(dir_cnv_res,"/Immunosuppressed/Cases/CNS_files/",cnv_file)))
            
          }
        }else{
          file.copy(paste0(dir_cnv_res,"Tumor/",cnv_file),paste0(paste0(dir_cnv_res,"/Immunocompetent/CNS_files/",cnv_file)))
          if(grepl("_0",smp_id)){
            file.copy(paste0(dir_cnv_res,"Tumor/",cnv_file),paste0(paste0(dir_cnv_res,"/Immunocompetent/Controls/CNS_files/",cnv_file)))
            
          }else{
            file.copy(paste0(dir_cnv_res,"Tumor/",cnv_file),paste0(paste0(dir_cnv_res,"/Immunocompetent/Cases/CNS_files/",cnv_file)))
            
          }
          
        }
      } 
    }
  }
}

#-----------------------
# 2. Create .seg files
#-----------------------
# In the command line:
for(cns_folder in c("cns")){
  for(incl_samples_description in c("All_CC_but_tcel_0","All_CC_but_tcel_0.135","All_CC_but_tcel_0.268")){
    for(gr in c("Tumor","Normal","Cases","Controls","Immunocompetent","Immunosuppressed")){
      if(gr%in%c("Cases","Controls","Immunocompetent","Immunosuppressed")){
        cat(gsub("incl_samples_description",incl_samples_description,paste0('cd "/Volumes/shainlab/Barbara/WES_analyses/Publication/Tot_147_with_rescued_muts/04_CNV/cns/incl_samples_description/',gr,'/CNS_files/"')))
      }else{
        cat(gsub("incl_samples_description",incl_samples_description,paste0('cd "/Volumes/shainlab/Barbara/WES_analyses/Publication/Tot_147_with_rescued_muts/04_CNV/cns/incl_samples_description/',gr,'/"')))
      }
      cat("\n")
      cat(gsub("incl_samples_description",incl_samples_description,paste0('cnvkit.py export seg *.cns -o Samples_',gr,'_',cns_folder,'_incl_samples_description.seg')))
      cat("\n")
    }
  }
}

#-----------------------
# 3. Merge all cns files together into one xlsx:
#-----------------------
for(incl_samples_description in c("All_CC_but_tcel_0","All_CC_but_tcel_0.135","All_CC_but_tcel_0.268")){
  dir_cns_files = file.path(dir_results_script,"cns",incl_samples_description,"/")
  cns_files =list.files(paste0(dir_cns_files,"Tumor"))
  cns_files = cns_files[grepl("cns$",cns_files)]
  # Save all segment data together
  cns_files_total=data.frame()
    for (smp_name in cns_files){
    
      cns_sample_file = read.table(file.path(dir_cns_files,"Tumor",smp_name),header=TRUE)
      cns_sample_file$Sample_id = gsub("-","_",sub("^S([0-9]+-[0-9]+)-.*", "\\1", smp_name))
      cns_files_total = rbind(cns_files_total,cns_sample_file)
    }
  
  #Add Tumor cellularity to the cns files:
  cns_files_total = merge(cns_files_total,sample_summary[,c("Sample.ID","Tumor.cellularity.avg.pct")],by.x="Sample_id",by.y="Sample.ID",all.x=TRUE)
  write.xlsx(cns_files_total, file.path(dir_cns_files,"Tumor",paste0("04_02_02_cns_files_total_",incl_samples_description,".xlsx")))
}
