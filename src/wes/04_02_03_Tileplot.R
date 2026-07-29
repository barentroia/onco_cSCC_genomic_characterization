#---------------------------------------------------
# Aim: Tileplot with cohort alterations
# Author: B. Rentroia Pacheco, based on Harsh's script
# Input: Gene list of interest + results from CGA
# Output: Output plots
#---------------------------------------------------


#---------------------------------------------------
# 1. Small changes and checks to the global variables in the masterscript:
#---------------------------------------------------

# We focus only on MM lists on clonal mutations, and exclude all mutations outside TERT that are outside exons:
mega_MML_clonal_visualization = mega_MML %>%dplyr::filter(Clonal=="Yes")%>%
  mutate(Exclude = ifelse(Hugo_Symbol!="TERT"& Variant_Classification%in%setdiff(synonymous,"Silent"),"Yes","No"))%>%
  filter(Exclude =="No")%>%
  as.data.frame()
mega_MML_clonal_visualization $Exclude= NULL

# Sanity check that we have annotations for all oncogenes and tsps
gene_path_dictionary = list("genes"=genes_and_pathways_correspondence$Gene,"pathways"=genes_and_pathways_correspondence$Pathway_publication)
stopifnot(length(setdiff(gene_path_dictionary$genes,c(oncogenes,tsps)))==0)

# Extract genes from intogen:
intogen_list  = driver_genes_oficial_list %>%filter(Driver=="Yes")%>%pull(SYMBOL)
stopifnot(length(setdiff(intogen_list,c(oncogenes,tsps)))==0)

# Check that for intogen driver genes, nonsynonymous mutations are shown as pathogenic:
tsps_need_hotspot = pathogenic_mutations_shain_lab$Gene[which(pathogenic_mutations_shain_lab$Type=="TSP" &pathogenic_mutations_shain_lab$Need_to_be_in_hotspot=="Yes")]
print(paste0("Pathogenic mutations: ",(mega_MML_clonal_visualization  %>%pull(Pathogenic)%>%table())["Probably"]))
mega_MML_clonal_visualization = mega_MML_clonal_visualization  %>%
  mutate(add_probably = ifelse(Hugo_Symbol%in%intersect(intogen_list,setdiff(tsps,tsps_need_hotspot)) &Variant_Classification%in%nonsynonymous ,"Probably","No"))%>%
  mutate(Pathogenic=ifelse(add_probably=="Probably"&!Pathogenic%in%"Probably","Probably",Pathogenic))%>%as.data.frame() #%>%pull(Pathogenic)%>%table()
print(paste0("TODO:Pathogenic mutations after updating with intogen this number should be ideally the same as the one before: ",(mega_MML_clonal_visualization  %>%pull(Pathogenic)%>%table())["Probably"]))
#write.csv(mega_MML_clonal_visualization,file.path(dir_results_script,"04_02_03_mega_MML_clonal_visualization.csv"))

# We also investigate the genes that were not considered as drivers:
intogen_list_not_drivers  = driver_genes_oficial_list %>%filter(Driver=="Possibly")%>%pull(SYMBOL)
oncogenes_extra=c("TPM3","KDR","KNSTRN","SDHC","CNTRL")
tsps_extra = c("MLH1","LAMC1", "IGF2BP2")
stopifnot(length(setdiff(intogen_list_not_drivers ,c(oncogenes_extra,tsps_extra)))==0)

#---------------------------------------------------
# 2. Check if there are arm-level changes affecting the chromosomes
#---------------------------------------------------
# Add chromosomal arms to each pathogenic mutation:
genes_pathogenic = unique(c(mega_MML_clonal_visualization$Hugo_Symbol[which(mega_MML_clonal_visualization$Pathogenic=="Probably")],focal_copy_nr$Gene))
genes_pathogenic_plus_removed = c(genes_pathogenic, intogen_list_not_drivers)
coords <- getBM(
  attributes = c(
    "hgnc_symbol",
    "chromosome_name",
    "start_position",
    "end_position",
    "band"          # cytoband (gives p/q arm info)
  ),
  filters = "hgnc_symbol",
  values = genes_pathogenic_plus_removed ,
  mart = mart
)%>%
  dplyr::filter(chromosome_name %in% c(as.character(1:22), "X", "Y"))
coords$arm <- paste0(coords$chromosome_name,substr(coords$band, 1, 1))

# Create a copy number table that is useful:
copy_number_drivers_alt = coords %>%
  select(hgnc_symbol,band,arm)%>%
  tidyr::crossing(
    Sample_ID = mega_MML_clonal_visualization$Tumor_Sample_Barcode
  )%>% dplyr::select(Sample_ID, everything())%>%
  mutate(Del_Amp =NA, 
         Chr_alt_type = NA)%>%
  dplyr::rename(Gene = hgnc_symbol)

# Add arm-level copy number alterations:
for(smp in c(copy_number_data_tileplot$Tumor_Sample_Barcode)){
  for(chr in setdiff(colnames(copy_number_data_tileplot),"Tumor_Sample_Barcode")){
    if(!chr%in%c("chrXp","chrXq","chrYp","chrYq")){
      index = which(copy_number_drivers_alt$Sample_ID==smp & copy_number_drivers_alt$arm== gsub("chr","",chr)) 
      alt_direction = copy_number_data_tileplot[which(copy_number_data_tileplot$Tumor_Sample_Barcode==smp),chr]
    }else if(chr%in%c("chrXp","chrXq")){
      female_samples = sample_summary$Sample.ID[which(sample_summary$Sex=="Female")]
      if(smp %in%female_samples){
        index = which(copy_number_drivers_alt$Sample_ID==smp & copy_number_drivers_alt$arm== gsub("chr","",chr)) 
        alt_direction = copy_number_data_tileplot[which(copy_number_data_tileplot$Tumor_Sample_Barcode==smp),chr]
      }else{
        alt_direction=0
      }
    }
    if(alt_direction==-1){
      copy_number_drivers_alt$Del_Amp[index] <- "Del"
    }else if(alt_direction==1){
      copy_number_drivers_alt$Del_Amp[index] <- "Amp"
    }
    
    if(alt_direction!=0){
      copy_number_drivers_alt$Chr_alt_type[index] <- "Arm_level"
    }
    }
}

# Add focal copy number alterations:
for(i in c(1:nrow(focal_copy_nr))){
  index = which(copy_number_drivers_alt$Sample_ID==focal_copy_nr$Sample_ID[i] & copy_number_drivers_alt$Gene==focal_copy_nr$Gene[i] )
  copy_number_drivers_alt$Del_Amp[index] <- focal_copy_nr$Del_Amp[i]
  copy_number_drivers_alt$Chr_alt_type[index] <- "Focal"
}

copy_number_drivers_and_others_alt = copy_number_drivers_alt
copy_number_drivers_alt = copy_number_drivers_alt %>%filter(Gene%in%genes_pathogenic)%>%as.data.frame()
#---------------------------------------------------
# 2. Tileplot functions
#---------------------------------------------------
# This function annotates the mutations that will be used in the tileplot
# mml: Data frame with all mutations found in the samples of interest (each row corresponds to one mutation). Columns should at least contain:Tumor_Sample_Barcode,Hugo_Symbol,Variant_Classification,Pathogenic,MAF,Adjusted_MAF
# oncogenes: character vector with the Hugo_symbols of the genes that should be considered as oncogenes
# tsps: character vector with the Hugo_symbols of the genes that should be considered as tumor suppressor genes
# path_mutation_labels: string vector with the levels that should be considered to identify pathogenic mutations: "Probably", "Possibly"
annotate_mutations_for_tileplot<-function(mml,oncogenes,tsps,path_mutations_labels,sample_summary,file_name=NULL){
  
  driver_chr_df = mml %>%
    filter(Hugo_Symbol %in% c(oncogenes,tsps)) %>%
    group_by(Hugo_Symbol) %>%
    dplyr::slice(1) %>%
    ungroup()%>%
    select(Hugo_Symbol, Chromosome)%>%as.data.frame()
  
  # Remove MML columns that are not required:
  mml = mml %>% dplyr::select(Tumor_Sample_Barcode,Hugo_Symbol,Variant_Classification,Pathogenic,MAF,Adjusted_MAF)%>%as.data.frame()
  all_samples <- unique(mml$Tumor_Sample_Barcode)
  
  # Label pathogenic mutations of interest:
  mml$Pathogenic[is.na(mml$Pathogenic)]="Unlikely"
  mml$Driver = "No"
  mml$Driver[which(mml[,"Pathogenic"]%in%path_mutations_labels)]="Driver"
  
  # Only retain genes of interest:
  goi = c(oncogenes,tsps)
  mml = mml %>% dplyr::filter(Hugo_Symbol%in%goi)%>%as.data.frame()
  all_genes   <- goi  # genes of interest (oncogenes + tsps)
  
  # Extra annotation on the genes:
  # Mutation subtype: variant classification 	
  Mutation_type <- onc_tsg <- driver_gain_of_function_onc <- driver_missense_tsg <- driver_others_tsg <- passenger_synonymous <- passenger_nonsynonymous <- NA
  
  data <- data.frame( mml, Mutation_type, onc_tsg, driver_gain_of_function_onc, driver_missense_tsg, driver_others_tsg, passenger_synonymous, passenger_nonsynonymous)
  
  # Synonymous or non synonymous mutation classification:
  data$Variant_Classification[is.na(data$Variant_Classification)]="NA"
  data$Mutation_type = ifelse(data$Variant_Classification%in%nonsynonymous,"Nonsynonymous","")
  data$Mutation_type = ifelse(data$Variant_Classification%in%synonymous,"Synonymous",data$Mutation_type)
  
  # Note that if TERT is present, pathogenic mutations should not be considered synonymous!
  if("TERT" %in%data$Hugo_Symbol){
    data$Mutation_type[intersect(which(data$Hugo_Symbol=="TERT"),which(data$Pathogenic=="Probably"))]="Nonsynonymous"
  }
  
  # Oncogene or TSP:
  data$onc_tsg = ifelse(data$Hugo_Symbol%in%oncogenes,"Oncogene","")
  data$onc_tsg = ifelse(data$Hugo_Symbol%in%tsps,"Tumor_suppressor_gene",data$onc_tsg)
  
  # Label type of alteration in nonsynonymous mutations:
  data$driver_gain_of_function_onc[which(data$onc_tsg == "Oncogene" & data$Driver == "Driver")] <- "gain_of_function"
  
  data$driver_missense_tsg[which(data$onc_tsg == "Tumor_suppressor_gene" & data$Driver == "Driver" & data$Variant_Classification == "Missense_Mutation")] <- "missense_mutation"
  
  data$driver_others_tsg[which(data$onc_tsg == "Tumor_suppressor_gene" & data$Driver == "Driver" & data$Variant_Classification != "Missense_Mutation" & data$Mutation_type == "Nonsynonymous")] <- "others_non_synonymous"
  
  # Label type of alteration in passenger mutations:
  data$passenger_synonymous[which(data$Mutation_type == "Synonymous")] <- "passenger_synonymous"
  data$passenger_nonsynonymous[which(data$Driver == "No" & 	data$onc_tsg == "Oncogene" &	data$Mutation_type == "Nonsynonymous")] <- "passenger_nonsynonymous"
  
  # Note: we are assuming that driver+onco = nonsynonymous
  #stopifnot(length(which(data$Driver=="Driver"&data$Mutation_type=="Synonymous"))==0) # No synonymous mutation is considered pathogenic
  print("Note: this TSPs have nonsynonymous mutations that were not considered as drivers")
  print(unique(data$Hugo_Symbol[which(data$Driver=="No"&data$Mutation_type=="Nonsynonymous"&data$onc_tsg == "Tumor_suppressor_gene")])) # All nonsynonymous mutation in TSP is considered oncogenic
  

  multiple_hits_df = data %>%dplyr::filter(Driver=="Driver")%>%add_count(Tumor_Sample_Barcode,Hugo_Symbol, name = "total")%>%distinct(Tumor_Sample_Barcode,Hugo_Symbol,total)%>%
    mutate(multiple_hits_in_gene = ifelse(total>1,"Yes","No"))%>%as.data.frame()
  data = merge(data,multiple_hits_df,by=c("Tumor_Sample_Barcode","Hugo_Symbol"),all.x=TRUE)
  
  # Note that a mutation in a X chromosome in males it is equivalent to a double hit. 
  # Get the chromosomes of all genes:
  genes_in_x_chr = driver_chr_df$Hugo_Symbol[which(driver_chr_df$Chromosome=="X")]
  for(gn in genes_in_x_chr){
    # Male samples:
    male_samples = sample_summary$Sample.ID[which(sample_summary$Sex=="Male")]
    i.males.mut.x = intersect(which(data$Hugo_Symbol==gn&data$Tumor_Sample_Barcode%in%male_samples),which(!is.na(data$total)))
    if(length(i.males.mut.x)>0){
      data[i.males.mut.x,"multiple_hits_in_gene"] ="Yes"
    }
  }
  
  
  # If a mutation has an adjusted MAF >1.5X, it probably underwent LOH of some sort:
  data$Adjusted_MAF = ifelse(data$Adjusted_MAF>1,1,data$Adjusted_MAF)
  data$LOH_imbalance = data$Adjusted_MAF/0.5
  data$LOH_presence = ifelse(data$LOH_imbalance>1.5,"Yes","No")
  
  # add multiple hits
  data$multiple_hits_in_gene = ifelse(data$LOH_presence=="Yes"& !(is.na(data$multiple_hits_in_gene)),"Yes",data$multiple_hits_in_gene) 
  
  # We give a rank for the color as follows:
  data$ranking = NA 
  # Driver on Oncogenes 
  data$ranking[which(data$driver_gain_of_function_onc=="gain_of_function")]=1
  # Driver on TSP
  data$ranking[which(data$driver_others_tsg=="others_non_synonymous")]=1
  data$ranking[which(data$driver_missense_tsg=="missense_mutation")]=2
  stopifnot(length(intersect(c(which(data$driver_missense_tsg=="missense_mutation"),which(data$driver_others_tsg=="others_non_synonymous")),which(data$passenger_nonsynonymous=="passenger_nonsynonymous")))==0)
  #stopifnot(length(intersect(c(which(data$driver_missense_tsg=="missense_mutation"),which(data$driver_others_tsg=="others_non_synonymous")),which(data$driver_others_tsg=="others_non_synonymous")))==0)
  # Driver on Passenger:
  data$ranking[which(data$passenger_synonymous=="passenger_synonymous")]=3
  data$ranking[which(data$passenger_nonsynonymous=="passenger_nonsynonymous")]=2
  
  # Only keep the mutation of highest ranking for each gene in each sample
  data = data %>%arrange(ranking)%>%
    group_by(Tumor_Sample_Barcode,Hugo_Symbol)%>%
    slice_head(n=1)%>%
    as.data.frame()
  
  # Make that the level for the color:
  # final_mutation_subtype
  data$Mutation_subtype <- NA
  for(i in c("driver_gain_of_function_onc","driver_missense_tsg","driver_others_tsg","passenger_synonymous","passenger_nonsynonymous")){
    for(j in 1:nrow(data)){
      if(!is.na(data[j , i])){
        data$Mutation_subtype[j] <- data[j , i]
      }
    } 
  }
  
  # Ensure one row per sample × gene (even if all NA)

  
  data <- data %>%
    tidyr::complete(
      Tumor_Sample_Barcode = all_samples,
      Hugo_Symbol = all_genes
    )
  
  if(!is.null(file_name)){
    write.xlsx(data,file_name)
  }
  return(data)
  
}

# This function converts the annotated master mutation list into a dataframe that can be plotted in a tileplot
# data: dataframe annotated using the annotate_mutations_for_tileplot function
# gene_path_dictionary: list with 1 element "genes containing all genes, and $ pathways with the corresponding pathways
# focal_copy_nr:  dataframe with the focal copy number alterations
# oncogenes: character vector with the Hugo_symbols of the genes that should be considered as oncogenes
# tsps: character vector with the Hugo_symbols of the genes that should be considered as tumor suppressor genes
# file_name: string indicating where to save the final table
generate_tileplot_dataframe <- function(data,gene_path_dictionary,focal_copy_nr=NULL,oncogenes=NULL,tsps=NULL,file_name=NULL){
     
  data$Mutation_subtype <- as.vector(data$Mutation_subtype)
  data$multiple_hits_in_gene <- as.vector(data$multiple_hits_in_gene)
  
  final_data <- data %>%
    dplyr::rename(
      Sample = Tumor_Sample_Barcode,
      Genes  = Hugo_Symbol
    )%>%
    select(Sample,Genes,Mutation_subtype,multiple_hits_in_gene)
  
  # Generate the final data with results of ALL driver genes in ALL samples
  uniq_genes <- unique(data$Hugo_Symbol)
  
  # Expand with the Genes from the focal copy number:
  if(!is.null(focal_copy_nr)){
    uniq_genes = unique(c(uniq_genes,unique(focal_copy_nr$Gene)))
  }
  
  samplenames <- unique(final_data$Sample)
  full_grid <- expand.grid(Sample = samplenames,
                           Genes = uniq_genes,
                           stringsAsFactors = FALSE)
  final_data <- full_grid %>%
    left_join(final_data, by = c("Sample", "Genes"))
  # add focal copy number:
  if(!is.null(focal_copy_nr)){
    # remove all rows without any copy number alteration:
    focal_copy_nr = focal_copy_nr[!is.na(focal_copy_nr$Del_Amp),]
    
    for(i in 1:nrow(final_data)){
      index <- which(final_data[ i, 1] == focal_copy_nr$Sample_ID &  final_data[ i, 2] == focal_copy_nr$Gene)
      if(length(index) != 0){
        #print(focal_copy_nr$Gene[index])
        
        # If there are no previous mutations on that gene in that sample:
        if(is.na(final_data$Mutation_subtype[i])){
          # We will only consider a hit in a gene without any mutation, if it is a focal copy number alteration, not an arm level alteration, as it is difficult to attribute those to selection in that specific gene:
          if( focal_copy_nr$Chr_alt_type[index]=="Focal"){
            final_data$Mutation_subtype[i] = focal_copy_nr$Del_Amp[index]
            if( final_data$Mutation_subtype[i] =="Del"){
              final_data$Mutation_subtype[i] = "others_non_synonymous"
            }else if (final_data$Mutation_subtype[i] =="Amp"){
              final_data$Mutation_subtype[i] = "gain_of_function"
            }
            final_data$multiple_hits_in_gene[i] = "No"
          }
        
        }else{
          two_alterations = paste0(final_data$Mutation_subtype[i]," + ",focal_copy_nr$Del_Amp[index])
          # This means that there is a mutation and a deletion/amplification
          # we consider two hits if we have a non synonymous mutation + focal deletion on a TSP:
          
          if(focal_copy_nr$Del_Amp[index]=="Del"){
            if(final_data$Genes[i]%in%tsps & two_alterations=="others_non_synonymous + Del"){
              final_data$multiple_hits_in_gene[i] = "Yes"
              final_data$Mutation_subtype[i] = "others_non_synonymous"
            }else if(final_data$Genes[i]%in%tsps & two_alterations=="missense_mutation + Del"){
              final_data$multiple_hits_in_gene[i] = "Yes"
              final_data$Mutation_subtype[i] = "others_non_synonymous"
            }else{
              if(focal_copy_nr$Chr_alt_type[index]=="Focal"){
                final_data$Mutation_subtype[i] = "others_non_synonymous"
                final_data$multiple_hits_in_gene[i] = "No"
              }
            }
          }
          
          if(focal_copy_nr$Chr_alt_type[index]=="Focal"){
          if(focal_copy_nr$Del_Amp[index]=="Amp"){
            if (final_data$Genes[i]%in%oncogenes & two_alterations=="gain_of_function + Amp"){
              final_data$multiple_hits_in_gene[i] = "Yes"
              final_data$Mutation_subtype[i] = "gain_of_function"
            }else{
             
                final_data$Mutation_subtype[i] = "gain_of_function"
                final_data$multiple_hits_in_gene[i] = "No"
              }

            }
          }
          #final_data$Mutation_subtype[i] = focal_copy_nr$Del_Amp[index]
          
        }
      }
    }
  }
  # Reorder factors:
  final_data$Mutation_subtype <- factor(final_data$Mutation_subtype, levels = c("gain_of_function", "others_non_synonymous", "missense_mutation", "passenger_nonsynonymous", "passenger_synonymous", "Del","Amp",NA))
  final_data$multiple_hits_in_gene <- factor(final_data$multiple_hits_in_gene,levels=c("No","Yes"))
  
  # Add pathway information:
  final_data$Pathway <- NA
  genes = gene_path_dictionary[["genes"]]
  pathways = gene_path_dictionary[["pathways"]]
  for(i in seq_along(genes)){
    index <- which(final_data$Genes == genes[i])
    final_data$Pathway[index] <- pathways[i]
  }
  
  # Factors to decide order on the plot:
  # final_data$Genes <- factor(final_data$Genes , levels=rev(genes))
  
  #print(table(final_data$Sample))
  stopifnot(sum(is.na(final_data$Sample))==0)
  
  if(!is.null(file_name)){
    final_data %>%
      dplyr::mutate(dplyr::across(where(is.character), ~ na_if(.x, ""))) %>%
      arrange(Genes, Sample) %>%
      write.xlsx(file_name)
  }
  return(final_data)
}

generate_tileplot <-function(data,gene_path_dictionary,sample_order,gene_order=NULL,focal_copy_nr=NULL,oncogenes=NULL,tsps=NULL,remove_only_passenger=TRUE,file_path=NULL,sample_summary=NULL){
  
  final_data <-generate_tileplot_dataframe(data,gene_path_dictionary,focal_copy_nr,oncogenes,tsps,file.path(file_path,"04_02_03_source_data_tileplot.xlsx"))


  if(remove_only_passenger){
    # Identify genes to keep: those that have at least one non-passenger mutation
    genes_to_keep <- final_data %>%
      group_by(Genes) %>%
      dplyr::summarise(has_non_passenger = any(!is.na(Mutation_subtype) & 
                                          !Mutation_subtype %in% c("passenger_nonsynonymous", "passenger_synonymous"))) %>%
      dplyr::filter(has_non_passenger) %>%
      pull(Genes)
    
    final_data=final_data %>%dplyr::filter(Genes%in%genes_to_keep)%>%as.data.frame()
  }
  
  if(is.null(gene_order)){
    final_data$Genes = factor(final_data$Genes,levels=unique(final_data$Genes))
    final_data = final_data %>%dplyr::filter(Genes%in%gene_order)%>%as.data.frame()
  }else if (gene_order == "Pathway"){
    gene_path_df = data.frame(gene_path_dictionary)
    count_path_alt = final_data %>%
      dplyr::filter(!is.na(Mutation_subtype) & 
               !Mutation_subtype %in% c("passenger_nonsynonymous", "passenger_synonymous"))%>%
      left_join(gene_path_df, by = c("Genes" = "genes"))%>%
      group_by(pathways) %>%
      dplyr::summarise(Pathway_frequency = n_distinct(Sample)) %>%
      ungroup()
    
    count_gene_alt = final_data %>%
      dplyr::filter(!is.na(Mutation_subtype) & 
               !Mutation_subtype %in% c("passenger_nonsynonymous", "passenger_synonymous"))%>%
      dplyr::count(Genes, name = "Count")%>%as.data.frame()
    
    merged_gn_path_counts = merge(gene_path_df,count_path_alt,by="pathways",all.x=TRUE)
    
    ordered_genes = merge(gene_path_df, count_gene_alt, by.x="genes", by.y="Genes", all.x=TRUE) %>%
      left_join(count_path_alt, by = "pathways") %>%
      arrange(desc(Pathway_frequency), pathways, desc(Count))%>%as.data.frame()
    
    gene_order = ordered_genes$genes
    
    final_data = final_data %>%dplyr::filter(Genes%in%gene_order)%>%as.data.frame()
    
    final_data$Genes = factor(final_data$Genes,levels=rev(gene_order))
    
    print(paste0("The following genes were removed from the plot: ",setdiff(data$Hugo_Symbol,gene_order)))
  }
  
  final_data$Genes <- droplevels(final_data$Genes)
  if(is.null(sample_order)){
    is_altered <- function(x) {
      !is.na(x) & !x %in% c("passenger_nonsynonymous", "passenger_synonymous")
    }
    
    gene_levels <- rev(levels(final_data$Genes))
    
    sample_order <- final_data %>%
      mutate(altered = as.integer(
        !is.na(Mutation_subtype) &
          !Mutation_subtype %in% c("passenger_nonsynonymous", "passenger_synonymous")
      )) %>%
      select(Sample, Genes, altered) %>%
      pivot_wider(
        names_from = Genes,
        values_from = altered,
        values_fill = 0
      ) %>%
      # force column order to match gene order
      select(Sample, all_of(gene_levels)) %>%
      mutate(total_alterations = rowSums(dplyr::across(-Sample))) %>%
      dplyr::arrange(
        across(all_of(gene_levels), desc),
        desc(total_alterations)
      ) %>%
      pull(Sample)
    
    final_data$Sample <- factor(final_data$Sample, levels = sample_order)
    
    
  }else{
    final_data$Sample = factor(final_data$Sample,levels=sample_order)
  }
  
  final_data = final_data %>%dplyr::filter(Sample%in%sample_order)%>%as.data.frame()
  
  n_sampl = length(unique(levels(final_data$Sample)))
  print(n_sampl)
  
  # Add percentage of driver mutations in each gene
  gene_freq <- final_data %>%
    group_by(Genes,Pathway) %>%
    dplyr::summarise(
      N_alterations =sum(
        !is.na(Mutation_subtype) &
          !Mutation_subtype %in% c("passenger_nonsynonymous", "passenger_synonymous")      ) 
    )%>%
    mutate(Pct_alterations = 100 *N_alterations/n_sampl)%>%as.data.frame()
 
  # Re-apply factor with same levels but new labels
  current_levels <- levels(final_data$Genes)
  new_labels <- gene_freq %>%
    mutate(
      label = paste0(
        Genes,
        " (",
        round(N_alterations, 1),
        "%)"
      )
    )%>%
    filter(Genes %in% current_levels) %>%
    arrange(match(Genes, current_levels))
  final_data$Genes <- factor(
    final_data$Genes,
    levels = current_levels,
    labels = paste0(
      current_levels,
      " (",
      round(new_labels$Pct_alterations, 1),
      "%)"
    )
  )
  
  # Compute the percentage of mutations within a pathway, and save it to add them to the plot in illustrator
  path_freq <- final_data %>%
    mutate(Mutation_subtype = as.character(Mutation_subtype)) %>%
    
    # Step 1: one row per sample per pathway
    group_by(Pathway, Sample) %>%
    dplyr::summarise(
      Has_alteration = any(
        !is.na(Mutation_subtype) &
          !Mutation_subtype %in% c("passenger_nonsynonymous",
                                   "passenger_synonymous")
      ),
      .groups = "drop"
    ) %>%
    
    # Step 2: count samples per pathway
    group_by(Pathway) %>%
    dplyr::summarise(
      N_samples_altered = sum(Has_alteration),
      .groups = "drop"
    ) %>%
    
    # Step 3: compute percentage
    mutate(
      Pct_alterations = 100 * N_samples_altered / n_sampl
    )
  
  sheet_source_data <- list("Genes" = gene_freq, "Pathways" = path_freq)
  
  write.xlsx( sheet_source_data,file.path(file_path,"04_02_03_pct_driver_muts_within_gen_pathways.xlsx"))
  
  theme1 <-theme(
    # LABELS APPEARANCE
    plot.title = element_text(colour= "black" ),
    axis.title.x = element_text(colour = "black"),    
    axis.title.y = element_text(colour = "black"),    
    axis.text.x = element_text(colour = "black"), 
    axis.text.y = element_text(colour = "black"),
    strip.text.x = element_text(colour = "black" ),
    strip.text.y = element_text(colour = "black"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank()
  ) 
  
  p = ggplot(final_data)+
    geom_tile(mapping = aes(x = Sample, y = Genes, fill = Mutation_subtype), height=1, width=1, color = "black",
              lwd = 0.25, linetype = 1, na.rm = T) +
    # scale_fill_manual( values = c("deeppink2", "dodgerblue2", "darkseagreen3", "burlywood2", "cornsilk2"), na.value= "white") +
    scale_fill_manual(
      values = c(
        "gain_of_function" = "#E41A1C",
        "others_non_synonymous" = "#6A51A3",
        "missense_mutation" = "dodgerblue3",
        "passenger_nonsynonymous" = "#FDB678",
        "passenger_synonymous" = "#FFFF99",
        "Del" = "orchid1",
        "Amp" = "tomato",
        "Case" = group_colors_Mets["Case"],
        "Control"=group_colors_Mets["Control"] 
      ),
      na.value = "white"
    )+
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme1 +
    geom_point(
      data = final_data %>% dplyr::filter(multiple_hits_in_gene == "Yes"),
      aes(x = Sample, y = Genes),
      size = 2.5,
      color = "white",
      na.rm = TRUE
    )+ 
    ylab("")+
    xlab("")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
  
  if(!is.null(sample_summary)){

    clinical_df = sample_summary[,c("Sample.ID","Metastasis",IS_classification)]
    clinical_df = clinical_df %>%filter(Sample.ID%in%unique(final_data$Sample))%>%as.data.frame()
    clinical_df$Sample.ID = factor(clinical_df$Sample.ID,levels=levels(final_data$Sample))
    clinical_long <- clinical_df %>%
      pivot_longer(cols = -Sample.ID, names_to = "Variable", values_to = "Value")
    clinical_long$y <- as.numeric(factor(clinical_long$Variable, levels = rev(unique(clinical_long$Variable))))
    
    # To have the same
    offset_y =  length(unique(clinical_long$Variable))
    # Clinical tiles below
    p=ggplot() +
      # Main mutation tiles
      geom_tile(
        data = final_data,
        mapping = aes(x = Sample, y = as.numeric(Genes) + offset_y , fill = Mutation_subtype),
        height = 1, width = 1, color = "black"
      ) +
      
      # Multiple hits overlay (white dot)
      geom_point(
        data = final_data %>% filter(multiple_hits_in_gene == "Yes"),
        mapping = aes(x = Sample, y = as.numeric(Genes) + offset_y ),
        size = 2.5,
        color = "white",
        na.rm = TRUE
      ) +
      
      # Clinical track below
      geom_tile(
        data = clinical_long,
        mapping = aes(x = Sample.ID, y = y, fill = Value),
        height = 1, width = 1, color = "black"
      ) +
      
      # Axis & themes
      scale_y_continuous(
        breaks = c(seq_along(levels(final_data$Genes)) + offset_y , unique(clinical_long$y)),
        labels = c(levels(final_data$Genes), unique(clinical_long$Variable)),
        expand = c(0,0)
      ) +
      scale_fill_manual(
        values = c(
          "gain_of_function" = "#E41A1C",
          "others_non_synonymous" = "#6A51A3",
          "missense_mutation" = "dodgerblue3",
          "passenger_nonsynonymous" = "#FDB678",
          "passenger_synonymous" = "#FFFF99",
          "Del" = "orchid1",
          "Amp" = "tomato",
          "Case" = group_colors_Mets[["Case"]],
          "Control" = group_colors_Mets[["Control"]], 
          "Yes"="darkred",
          "No" = "lightgray"
        ),
        na.value = "white"
      ) +
      theme_bw() +
      theme1 +
      theme(
        axis.text.y=element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = 18),
        legend.text  = element_text(size = 18)
      ) +
      ylab("") +
      xlab("")
  }
  
  ggsave(file.path(file_path,"04_02_03_tileplot.pdf"),p,width=18,height=10)
  return(p)
}

#-----------------------
# 3. Generate the plot and save associated data
#-----------------------
genes_to_show = unique(c(intogen_list,setdiff(oncokb_muts$Hugo.Symbol,NA)))
write.csv(genes_to_show,file.path(dir_results_script,"04_02_03_genes_to_show.csv"))
source_data_mutations_fpath = file.path(dir_results_script,"04_02_03_source_data_mutations_tileplot.xlsx")
data =annotate_mutations_for_tileplot(mega_MML_clonal_visualization ,intersect(oncogenes,genes_to_show),intersect(tsps,genes_to_show),c("Probably"),sample_summary = sample_summary,source_data_mutations_fpath )
sample_order = sample_summary %>%
  arrange(Metastasis, desc(`Burden.(Mutations/megabase).(MB)`))%>%pull(Sample.ID)
  #arrange(.data[[IS_classification]],Metastasis, desc(`Burden.(Mutations/megabase).(MB)`))%>%pull(Sample.ID)#unique(c(data$Tumor_Sample_Barcode[grepl("_1",data$Tumor_Sample_Barcode)],data$Tumor_Sample_Barcode[grepl("_0",data$Tumor_Sample_Barcode)]))
#gene_order = unique(c(oncogenes,tsps,unique(focal_copy_nr$Gene)))
generate_tileplot(data,gene_path_dictionary, NULL,gene_order="Pathway",copy_number_drivers_alt,oncogenes,tsps,file_path=dir_results_script,sample_summary=sample_summary)
# Add clinical information to the source data:
final_data=read.xlsx(file.path(dir_results_script,"04_02_03_source_data_tileplot.xlsx"))
write.xlsx(merge(final_data,sample_summary[,c("Sample.ID","Metastasis",IS_classification)],by.x="Sample",by.y="Sample.ID",all.x=TRUE),file.path(dir_results_script,"04_02_03_source_data_tileplot.xlsx"))

#-----------------------
# 4. Compare number of driver genes between cases and controls
#-----------------------
# Add follow-up information:
final_data=merge(final_data,sample_summary[,c("Sample.ID","Metastasis","FU.metastasis.years","Percent_Bait_with_callable_coverage",IS_classification)],by.x="Sample",by.y="Sample.ID")

# Prepare the data
mut_data_drivers = final_data %>%
  mutate(
    Met_binary = ifelse(Metastasis == "Control", 0, 1),
    Mutation_subtype = as.character(Mutation_subtype),
    
    # Create weight per mutation
    gene_weight = case_when(
      is.na(Mutation_subtype) ~ 0,
      Mutation_subtype %in% c("passenger_nonsynonymous",
                              "passenger_synonymous") ~ 0,
      multiple_hits_in_gene == "Yes" ~ 2,
      TRUE ~ 1
    ),
    Any_mutation = factor(ifelse(gene_weight>0,"Yes","No"),levels=c("No","Yes"))
  )

onc_score_df <- mut_data_drivers%>%
  group_by(Sample)%>%
  dplyr::summarise(
    Nr_hits= sum(gene_weight),
    Met_factor = dplyr::first(Metastasis),
    Percent_Bait_with_callable_coverage=dplyr::first(Percent_Bait_with_callable_coverage),
    IS_classification = dplyr::first(.data[[IS_classification]]),
    .groups = "drop"
  )

# Visual comparison of the the number of hits:
# Fit linear model adjusting for callable coverage
lm_model <- lm(Nr_hits ~Met_factor+ Percent_Bait_with_callable_coverage, data = onc_score_df )

# Extract p-value for Metastasis
p_val <- summary(lm_model)$coefficients["Met_factorControl", "Pr(>|t|)"]  # adjust name if needed
# Format for display
p_label <- ifelse(p_val < 0.001, "p < 0.001", paste0("p adj call.cov= ", signif(p_val, 2)))

stat.test <- compare_means(
  Nr_hits ~ Met_factor,
  data = onc_score_df ,
  method = "wilcox.test"
)%>%
  mutate(label = paste0("p-val: ", p.format))


p_dhits_Mets = ggplot(onc_score_df,
    aes(y = Nr_hits,
        x = Met_factor)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0,alpha = 0.8, size = 1)  +
  labs(y = "Oncogenic score", fill = "Group") +
  theme_bw(base_size = 18) +
  labs(
    x = "Metastatic Outcome",
    y = "Number of driver hits")+
  ylim(c(0, max(onc_score_df$Nr_hits) * 1.13))+
  stat_pvalue_manual(
    stat.test,
    label = "label",
    y.position = max(onc_score_df$Nr_hits) * 1.1,size=5)+
  annotate("text", x = 1.5, y = max(onc_score_df$Nr_hits) * 1.05,
           label = p_label, size = 5)

ggsave(paste0(dir_results_script,"/S04_02_03_dhits_Metastasis.pdf"),p_dhits_Mets,width=5,height=6)
write.xlsx(onc_score_df[,c("Sample","Nr_hits","Met_factor")],paste0(dir_results_script,"/S04_02_03_dhits_Metastasis_sourcedata.xlsx"))
plot_driver_hits <- function(df, group_col, output_file, x_label = NULL) {
  # df: dataframe with Nr_hits and other variables
  # group_col: string, column name to group by (e.g., "Met_factor" or "IS.at.cSCC.updatedS32")
  # output_file: path to save PDF
  # x_label: optional x-axis label
  
  # Summarize per Sample
  onc_score_df <- df %>%
    group_by(Sample) %>%
    summarise(
      Nr_hits = sum(gene_weight),
      Group = first(.data[[group_col]]),
      Percent_Bait_with_callable_coverage = first(Percent_Bait_with_callable_coverage),
      .groups = "drop"
    )%>%
    mutate(Group = factor(Group,levels=levels(df[,group_col])))
  
  # Linear model adjusting for callable coverage
  lm_model <- lm(Nr_hits ~ Group + Percent_Bait_with_callable_coverage, data = onc_score_df)
  
  # Extract p-value for the first level of the group
  coef_name <- paste0("Group", levels(onc_score_df$Group)[2])
  p_val <- summary(lm_model)$coefficients[coef_name, "Pr(>|t|)"]
  p_label <- ifelse(p_val < 0.001, "p < 0.001", paste0("p adj call.cov= ", signif(p_val, 2)))
  
  # Wilcoxon test
  stat.test <- compare_means(
    Nr_hits ~ Group,
    data = onc_score_df,
    method = "wilcox.test"
  ) %>%
    mutate(label = paste0("p-val: ", p.format))
  
  # Plot
  p <- ggplot(onc_score_df, aes(x = Group, y = Nr_hits)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, height = 0, alpha = 0.8, size = 1) +
    theme_bw(base_size = 18) +
    labs(
      x = ifelse(is.null(x_label), group_col, x_label),
      y = "Number of driver hits"
    ) +
    coord_cartesian(ylim = c(0, max(onc_score_df$Nr_hits) * 1.13)) +
    stat_pvalue_manual(
      stat.test,
      label = "label",
      y.position = max(onc_score_df$Nr_hits) * 1.1,
      size = 5
    ) +
    annotate(
      "text",
      x = 1.5,
      y = max(onc_score_df$Nr_hits) * 1.05,
      label = p_label,
      size = 5
    )
  
  # Save plot
  ggsave(output_file, p, width = 5, height = 6)
  
  return(p)
}

# Metastasis plot
mut_data_drivers = mut_data_drivers %>%
  mutate(Metastasis=factor(Metastasis,levels=c("Case","Control")),
         IS_classification=factor(.data[[IS_classification]],levels=c("No","Yes")))
plot_driver_hits(
  df = mut_data_drivers,
  group_col = "Metastasis",
  output_file = paste0(dir_results_script, "/S04_02_03_dhits_Metastasis.pdf"),
  x_label = "Metastatic Outcome"
)

# IS classification plot
plot_driver_hits(
  df = mut_data_drivers,
  group_col = "IS_classification",
  output_file = paste0(dir_results_script, "/S04_02_03_dhits_IS_classification.pdf"),
  x_label = "IS Classification"
)

write.xlsx(onc_score_df[,c("Sample","Nr_hits","IS_classification")],paste0(dir_results_script,"/S04_02_03_dhits_IS_sourcedata.xlsx"))

# Focus on genes altered instead of driver genes:
plot_driver_hits_and_genes <- function(df, group_col, output_file, x_label = NULL) {
  # df: dataframe with Nr_hits and other variables
  # group_col: string, column name to group by (e.g., "Met_factor" or "IS.at.cSCC.updatedS32")
  # output_file: path to save PDF
  # x_label: optional x-axis label
  
  # Summarize per Sample
  onc_score_df <- df %>%
    mutate(driver_alteration_present = ifelse(Any_mutation == "Yes", 1, 0)) %>%
    group_by(Sample) %>%
    summarise(
      Nr_hits = sum(gene_weight),
      Nr_alt_genes = sum(driver_alteration_present),
      Group = first(.data[[group_col]]),
      Percent_Bait_with_callable_coverage = first(Percent_Bait_with_callable_coverage),
      .groups = "drop"
    )%>%
    mutate(Group = factor(Group,levels=levels(df[,group_col])))
  
  write.xlsx(onc_score_df,gsub("pdf","xlsx",output_file))
  # Number of hits:
  # Linear model adjusting for callable coverage, for the number of hits:
  lm_model <- lm(Nr_hits ~ Group + Percent_Bait_with_callable_coverage, data = onc_score_df)
  # Extract p-value for the first level of the group
  coef_name <- paste0("Group", levels(onc_score_df$Group)[2])
  p_val <- summary(lm_model)$coefficients[coef_name, "Pr(>|t|)"]
  p_label_nr_hits<- ifelse(p_val < 0.001, "p < 0.001", paste0("p adj call.cov= ", signif(p_val, 2)))
  
  # Same but for Number of altered genes:
  # Linear model adjusting for callable coverage, for the number of altered genes:
  lm_model <- lm(Nr_alt_genes ~ Group + Percent_Bait_with_callable_coverage, data = onc_score_df)
  # Extract p-value for the first level of the group
  coef_name <- paste0("Group", levels(onc_score_df$Group)[2])
  p_val <- summary(lm_model)$coefficients[coef_name, "Pr(>|t|)"]
  p_label_nr_alt_genes <- ifelse(p_val < 0.001, "p < 0.001", paste0("p adj call.cov= ", signif(p_val, 2)))
  
  lm_labels <- tibble(
    Metric  = c("Nr_hits", "Nr_alt_genes"),
    p.lm    = c(p_label_nr_hits, p_label_nr_alt_genes)
  )
  
  # Wilcoxon test
  plot_data <- onc_score_df %>%
    select(Nr_hits, Nr_alt_genes, Group) %>%
    pivot_longer(cols = c(Nr_hits, Nr_alt_genes),
                 names_to  = "Metric",
                 values_to = "Value")
  
  wilcox_results <- plot_data %>%
    group_by(Metric) %>%
    summarise(
      p.value = wilcox.test(Value ~ Group)$p.value,
      .groups = "drop"
    ) %>%
    mutate(
      p.wilcox = paste0("p = ", round(p.value, 2))
      )
  
  # Compute bracket y positions per Metric (just above the max value)
  bracket_pos <- plot_data %>%
    group_by(Metric) %>%
    summarise(y.pos = max(Value) * 1.15, .groups = "drop")
  
  # Merge annotation data
  annot <- wilcox_results %>%
    left_join(lm_labels,    by = "Metric") %>%
    left_join(bracket_pos,  by = "Metric") %>%
    mutate(
      full_label = paste0(p.wilcox, "\n", p.lm)  # stacked label
    )
  
  # Get group levels for bracket x positions
  groups <- unique(plot_data$Group)
  x1 <- 1   # position of first group on x-axis
  x2 <- 2   # position of second group
  
  # Plot
  p <- ggplot(plot_data, aes(x = Group, y = Value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_jitter(width = 0.15),
                alpha = 0.8, size = 1) +
    geom_segment(data = annot,
                 aes(x = x1, xend = x2, y = y.pos, yend = y.pos),
                 inherit.aes = FALSE) +
    geom_segment(data = annot,
                 aes(x = x1, xend = x1, y = y.pos * 0.98, yend = y.pos),
                 inherit.aes = FALSE) +
    geom_segment(data = annot,
                 aes(x = x2, xend = x2, y = y.pos * 0.98, yend = y.pos),
                 inherit.aes = FALSE) +
    geom_text(data = annot,
              aes(x = (x1 + x2) / 2, y = y.pos * 1.08, label = full_label),
              inherit.aes = FALSE, size = 3.5, lineheight = 0.9) +
    facet_wrap(~ Metric) +
    labs(x = x_label %||% "", y = "Frequency") +
    theme_bw(base_size = 18)
  
  # Save plot
  ggsave(output_file, p, width = 6.5, height = 5)
  
  return(p)
}

plot_driver_hits_and_genes(
  df = mut_data_drivers,
  group_col = "Metastasis",
  output_file = paste0(dir_results_script, "/S04_02_03_daltgenes_hits_Metastasis.pdf"),
  x_label = "Metastatic Outcome"
)

# IS classification plot
plot_driver_hits_and_genes(
  df = mut_data_drivers,
  group_col = "IS_classification",
  output_file = paste0(dir_results_script, "/S04_02_03_daltgenes_hits_IS_classification.pdf"),
  x_label = "IS Classification"
)

# combine the source data into one table:
onco_1 = read.xlsx(paste0(dir_results_script, "/S04_02_03_daltgenes_hits_Metastasis.xlsx"))
onco_2 = read.xlsx(paste0(dir_results_script, "/S04_02_03_daltgenes_hits_IS_classification.xlsx"))

onco_2 %>%
  rename(Immunossuppression = Group)%>%
  left_join(onco_1)%>%
  select(Sample,Nr_alt_genes,Nr_hits,Group,Immunossuppression,Percent_Bait_with_callable_coverage )%>%
  rename(Metastasis = Group)%>%
  write.xlsx(paste0(dir_results_script, "/S04_02_03_daltgenes_hits_combined.xlsx"))
#-----------------------
# 5. Supplementary Figure, focusing on the genes that we excluded:
#-----------------------
mega_MML_clonal_visualization_not_drivers = mega_MML_clonal_visualization  %>%
  mutate(add_probably = ifelse(Hugo_Symbol%in%intersect(intogen_list_not_drivers ,setdiff(tsps_extra,tsps_need_hotspot)) &Variant_Classification%in%nonsynonymous ,"Probably","No"))%>%
  mutate(Pathogenic=ifelse(add_probably=="Probably"&!Pathogenic%in%"Probably","Probably",Pathogenic))%>%
  as.data.frame() #%>%pull(Pathogenic)%>%table()

source_data_mutations_fpath_exc_genes = file.path(dir_results_script,"04_02_03_source_data_mutations_tileplot.xlsx")
data_excgenes =annotate_mutations_for_tileplot(mega_MML_clonal_visualization_not_drivers  ,oncogenes_extra,tsps_extra ,c("Probably"),sample_summary = sample_summary,source_data_mutations_fpath_exc_genes  )
gene_path_dictionary_extra = list(genes = c(gene_path_dictionary$genes,setdiff(c(oncogenes_extra,tsps_extra),gene_path_dictionary$genes)),"pathways"=c(gene_path_dictionary$pathways,rep("Other",length(setdiff(c(oncogenes_extra,tsps_extra),gene_path_dictionary$genes)))))
generate_tileplot(rbind(data,data_excgenes),gene_path_dictionary_extra, NULL,gene_order="Pathway",copy_number_drivers_and_others_alt,c(oncogenes,oncogenes_extra),c(tsps,tsps_extra),file_path=file.path(dir_results_script,"Excluded_genes"),sample_summary=sample_summary)
final_data_all = read.xlsx(file.path(dir_results_script,"04_02_03_source_data_tileplot.xlsx"))

# Now join the two:
pathway_order = c("Notch","P53","Hippo","Rb","SWI_SNF_PRC2","Ras_MAPK","TGF_BETA")


# 1. Prepare plot data
df_plot <- final_data_all %>%
  filter(!is.na(multiple_hits_in_gene)) %>%
  count(Genes, multiple_hits_in_gene) %>%
  complete(Genes, multiple_hits_in_gene = c("Yes", "No"), fill = list(n = 0)) %>%
  group_by(Genes) %>%
  filter(sum(n) >= 5) %>% # remove genes with very few total observations
  mutate(frac = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    onc_tsg = ifelse(Genes %in% c(tsps, tsps_extra), "TSP", "Oncogene"),
    Nominated = ifelse(Genes %in% genes_to_show, "Yes", "No")
  )

# 2. Compute n_yes / n_total and join chromosome arm
df_labels <- df_plot %>%
  group_by(Genes, Nominated) %>%
  summarise(
    n_yes = sum(n[multiple_hits_in_gene == "Yes"]),
    n_total = sum(n),
    order_frac = ifelse(n_yes > 0, n_yes / n_total, 0),
    .groups = "drop"
  ) %>%
  # join chromosome arm from coords
  left_join(coords %>% select(hgnc_symbol, arm),
            by = c("Genes" = "hgnc_symbol")) %>%
  # create label: Gene (arm) + n_yes/n_total
  mutate(
    Genes_label = paste0(Genes, "\n (", arm, ")\n", n_yes, "/", n_total)
  )

# 3. Join labels back and reorder
df_plot <- df_plot %>%
  left_join(df_labels %>% select(Genes, Nominated, Genes_label, order_frac),
            by = c("Genes", "Nominated")) %>%
  mutate(
    Genes_label = reorder_within(Genes_label, -order_frac, Nominated),
    multiple_hits_in_gene = factor(multiple_hits_in_gene, levels = c("No", "Yes")),
    Nominated = factor(ifelse(Nominated=="Yes","Nominated","Not nominated"),
                       levels=c("Nominated","Not nominated"))
  )

# 4. Plot
p_tsp = ggplot(df_plot %>% filter(onc_tsg != "Oncogene", !is.na(frac)),
       aes(x = Genes_label, y = frac, fill = multiple_hits_in_gene)) +
  geom_col(width = 0.8) +
  scale_x_reordered() +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("No" = "grey80", "Yes" = "darkred")) +
  labs(
    x = "Gene (chromosome arm) + number of multiple hits / total",
    y = "Mutations with more than 1 hit (%)",
    fill = "Multiple hits"
  ) +
  theme_minimal(base_size = 16) +                       # <-- set base_size here
  theme(
    axis.text.x = element_text(size = 14),  # bigger labels
    axis.text.y = element_text(size = 14),
    panel.spacing = unit(2, "lines"), strip.text = element_text(size = 16, face = "bold") 
  )+
  facet_grid(~Nominated, scales = "free_x", space = "free_x")

ggsave(file.path(dir_results_script,"Excluded_genes","04_02_03_TSP_stacked_barplot_multiplehits.pdf"),p_tsp,width=20,height=5)
df_plot %>% filter(onc_tsg != "Oncogene", !is.na(frac)) %>%
  write.xlsx(file.path(dir_results_script,"Excluded_genes","04_02_03_TSP_stacked_barplot_multiplehits_sourcedata.xlsx"))

p_onco = ggplot(df_plot %>% filter(onc_tsg == "Oncogene", !is.na(frac)),
               aes(x = Genes_label, y = frac, fill = multiple_hits_in_gene)) +
  geom_col(width = 0.8) +
  scale_x_reordered() +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("No" = "grey80", "Yes" = "darkred")) +
  labs(
    x = "Gene (chromosome arm) + number of multiple hits / total",
    y = "Mutations with more than 1 hit (%)",
    fill = "Multiple hits"
  ) +
  theme_minimal(base_size = 16) +                       # <-- set base_size here
  theme(
    axis.text.x = element_text(size = 14),  # bigger labels
    axis.text.y = element_text(size = 14),
    panel.spacing = unit(2, "lines"), strip.text = element_text(size = 16, face = "bold") 
  )
ggsave(file.path(dir_results_script,"Excluded_genes","04_02_03_ONCO_stacked_barplot_multiplehits.pdf"),p_onco,width=6.5,height=4.5)
df_plot %>% filter(onc_tsg == "Oncogene", !is.na(frac)) %>%
  write.xlsx(file.path(dir_results_script,"Excluded_genes","04_02_03_ONCO_stacked_barplot_multiplehits_sourcedata.xlsx"))

#-----------------------
# 6. Supplementary Figure: lollipop plot for the excluded genes:
#-----------------------

# Loliplots for each candidate gene:

#Lolliplot function:
make_lolliplot = function(gene_symbol,mega_MML,file_name=NULL){
  gene_id <- get(gene_symbol, org.Hs.egSYMBOL2EG)
  # Get gene ensemblE:
  gene_ensb = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                    filters = "hgnc_symbol",
                    values = gene_symbol,
                    mart = mart)[,1]
  
  # Extract all mutations in the coding region:
  mml_gene = mega_MML %>%dplyr::filter(Hugo_Symbol==gene_symbol,!Variant_Classification%in%c("Intron","3'UTR","5'UTR","5'Flank"))%>%arrange(Start_Position)%>%as.data.frame()
  mml_gene$aa_pos = as.numeric(str_extract(mml_gene$Protein_Change, "\\d+"))
  mml_gene = mml_gene[!is.na(mml_gene$aa_pos),]
  
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
    if(sample.gr$type[z] == "DE_NOVO_START_OUT_FRAME") return("#9400d3")
    if(sample.gr$type[z] == "Nonsense_Mutation") return("#9400d3")
    if(sample.gr$type[z] == "Splice_Site") return("#9400d3")
    if(sample.gr$type[z] == "START_CODON_SNP") return("#9400d3")
    if(sample.gr$type[z] == "Frame_Shift_Ins") return("#9400d3")
    if(sample.gr$type[z] == "Frame_Shift_Del") return("#9400d3")
  }))
  
  
  
  sample.gr$score = mutation_counts$score[match(sample.gr$mutation, mutation_counts$Protein_Change)]

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
      domains$fill <- 1+seq_along(domains)
      domains$height <- 0.05
      domains$type <- "domain"
      
      if(!is.null(file_name)){
        pdf(file_name,width=10,height=4)
      }
      lolliplot(sample.gr, domains, ranges = gr,score = sample.gr$score,legend = NULL,domainCol =  c("domain" = "darkgray"))
      grid.text(gene_symbol, x=.5, y=.6, just="top", 
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
      
      if(!is.null(file_name)){
        pdf(file_name,width=10,height=4)
      }
      lolliplot(sample.gr, features = features,legend = NULL)
      grid.text(gene_symbol, x=.5, y=.9, just="top", 
                gp=gpar(cex=1.5, fontface="bold"))
    }
    if(!is.null(file_name)){
    dev.off()}
  }

}

for(gn_name in c("LAMC1","IGF2BP2")){
  make_lolliplot(gn_name,mega_MML_clonal_visualization,file.path(dir_results_script,paste0("S04_02_03_",gn_name,"_lollipop.pdf")))
}

