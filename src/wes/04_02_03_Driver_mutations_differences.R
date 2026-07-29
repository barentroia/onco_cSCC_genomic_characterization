#---------------------------------------------------
# Aim: Show alterations within pathway
# Author: B. Rentroia Pacheco
# Input: Mutations within each pathway
# Output: Output plots
#---------------------------------------------------

#---------------------------------------------------
# 1. Use the tileplot data to show associations with prognosis
#---------------------------------------------------
final_data=read.xlsx(file.path(gsub("06_Forest_plots",paste0("05_Tileplot/",IS_classification),dir_results_script),"04_02_03_source_data_tileplot.xlsx"))

# Add follow-up information:
final_data=merge(final_data,sample_summary[,c("Sample.ID","FU.metastasis.years","Percent_Bait_with_callable_coverage")],by.x="Sample",by.y="Sample.ID")

#if (add_HLA){
#  final_data = rbind(final_data,sample_summary %>%select(Sample.ID,Metastasis,!!IS_classification,FU.metastasis.years,Percent_Bait_with_callable_coverage)%>%
#    mutate(Genes = "HLA",
#           Mutation_subtype = ifelse(Sample.ID%in%HLA_info$Point.Mutation[!is.na(HLA_info$Point.Mutation)],"missense_mutation",NA),
          # Mutation_subtype = ifelse(Sample.ID%in%HLA_info$LOH[!is.na(HLA_info$LOH)],"others_non_synonymous",Mutation_subtype),
#           multiple_hits_in_gene = ifelse(Sample.ID%in%intersect(HLA_info$LOH[!is.na(HLA_info$LOH)],HLA_info$Point.Mutation[!is.na(HLA_info$Point.Mutation)]),"Yes",NA))%>%
#    dplyr::rename(Sample = Sample.ID)%>%
#    select(Sample,Genes,Pathway, Mutation_subtype,multiple_hits_in_gene,Metastasis,!!IS_classification,FU.metastasis.years,Percent_Bait_with_callable_coverage))
#}

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

path_sample = mut_data_drivers %>%  dplyr::filter(!is.na(Pathway)) %>%
  group_by(Sample, Pathway) %>%
  dplyr::summarise(
    Weighted_burden = sum(gene_weight),
    Met_binary = dplyr::first(Met_binary),
    FU.metastasis.years = dplyr::first(FU.metastasis.years),
    Percent_Bait_with_callable_coverage = dplyr::first(Percent_Bait_with_callable_coverage),
    IS_classification = dplyr::first(.data[[IS_classification]]),
    .groups = "drop"
  )%>%
  mutate(
    Met_factor = factor(Met_binary,
                        levels = c(0,1),
                        labels = c("Control","Case")),
    Any_mutation = factor(ifelse(Weighted_burden>0,"Yes","No"),levels=c("No","Yes"))
  )

# Save a dataframe that can be used to compare with RNA-seq data:
# This dataframe will contain mutation, pathway and copy number data:
# Include mutation data:
genomic_df_summary = mut_data_drivers %>%
  pivot_wider(
    id_cols = Sample,
    names_from = Genes,
    values_from = gene_weight,
    values_fill = 0
  )%>%as.data.frame()

# Include pathway information:
pathway_info = path_sample %>%
  pivot_wider(
    id_cols = Sample,
    names_from = Pathway,
    values_from = Weighted_burden,
    values_fill = 0
  )%>%as.data.frame()

# now remove overlapping columns
cols_to_remove <- setdiff(intersect(colnames(pathway_info), colnames(genomic_df_summary)),"Sample")

pathway_info <- pathway_info %>%
  select(-c(any_of(cols_to_remove),"Metabolism","Jak_Stat","Chromatin_remodeling","Spliceosome","NF_kb")) %>%
  as.data.frame()

# Add pathway data:
genomic_df_summary = merge(genomic_df_summary,pathway_info,by="Sample")

# Merge with copy number data:
genomic_df_summary = merge(genomic_df_summary,copy_number_data_tileplot,by.x="Sample",by.y="Tumor_Sample_Barcode")

# Add basic information:
# Metastasis, Tumor cellularity, Mutation burden, pct of cnv altered
genomic_df_summary = merge(genomic_df_summary,sample_summary[,c("Sample.ID","IS.at.cSCC.updatedS32","Tumor.cellularity.avg.pct","Dominant.clone.UV.mutations.(%)","Burden.(Mutations/megabase).(MB)","pct_altered","Percent_Bait_with_callable_coverage")],by.x="Sample",by.y="Sample.ID")
genomic_df_summary = merge(genomic_df_summary,df_map[,c("UCSF_short_ID","universal_patient_ID","SkylineDx.ID")],by.x="Sample",by.y="UCSF_short_ID",all.x=TRUE)
write.xlsx(genomic_df_summary,file.path(dir_results_script,paste0("S04_02_03_","genomic_df_summary.xlsx")))
#---------------------------------------------------
# 2. Build a forest plot with the Hazard ratios for each pathway
#---------------------------------------------------
# Plot results in a forest plot:
pathways_with_enough_mutations = path_sample %>%
  group_by(Pathway,Sample)%>%dplyr::summarise(count = sum(Any_mutation=="Yes"),.groups = "drop")%>%
  group_by(Pathway)%>%dplyr::summarise(count = sum(count))%>%filter(count>2)%>%pull(Pathway)
genes_with_enough_mutations = mut_data_drivers %>%group_by(Genes,Sample)%>%dplyr::summarise(count = sum(gene_weight>0),.groups = "drop")%>%
  group_by(Genes)%>%dplyr::summarise(count = sum(count))%>%filter(count>2)%>%pull(Genes)

#Function to perform Cox models + plot them in forest plot:
apply_coxmodels = function(df, var_of_interest,group_of_interest,adjustment=NULL,group_levels_show=NULL,cox_file_path){
  if(!is.null(adjustment)){
    formula_txt <- paste0(
      "Surv(FU.metastasis.years, Met_binary) ~ ",
      var_of_interest, " + ", adjustment
    )
    
  }else{
    formula_txt <- paste0(
      "Surv(FU.metastasis.years, Met_binary) ~ ",
      var_of_interest
    )
  }
  formula_obj <- as.formula(formula_txt)
  
  if(var_of_interest=="Any_mutation"){
    var_of_interest_cox_r = paste0(var_of_interest,"Yes")
  }else{
    var_of_interest_cox_r = var_of_interest
  }
  
  
  cox_results <- df%>%
    group_by(.data[[group_of_interest]]) %>%
    nest() %>%
    mutate(
      model = map(data, ~coxph(formula_obj, data = .x)),
      tidy = map(model, ~tidy(.x, exponentiate = TRUE, conf.int = TRUE))
    ) %>%
    unnest(tidy) %>%
    filter(term == var_of_interest_cox_r) %>%   # keep mutation effect only
    select(.data[[group_of_interest]], estimate, conf.low, conf.high, p.value) %>%
    ungroup()%>%as.data.frame()
  group_levels_show = intersect(group_levels_show,cox_results[!is.infinite(cox_results$conf.high),1])
  x_max <- cox_results %>%filter(.data[[group_of_interest]]%in%group_levels_show)%>%pull(conf.high) %>%max(na.rm = TRUE)
  adj_txt = ifelse(is.null(adjustment),", no adjustment",paste0(",with adjustment:",adjustment))
  p=cox_results %>%filter(.data[[group_of_interest]]%in%group_levels_show)%>% # Fewer than 2 mutations%>%
    mutate(Significant = ifelse(p.value < 0.1,"signif","ns"),  
           label = sprintf("%.2f (%.2f–%.2f), p = %.3f",estimate, conf.low, conf.high, p.value))%>%
    mutate(Significant= ifelse(p.value<0.05,"ssignif",Significant))%>%
    
    ggplot(aes(x = estimate,
               y = reorder(.data[[group_of_interest]], estimate),
               color = Significant)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                   height = 0.2,
                   linewidth = 0.8) +
    
    # Add text to the right of the CI
    geom_text(aes(x = x_max * 1.4, label = label),
              hjust = 0,
              size = 3.5,
              color = "black") +
    
    geom_vline(xintercept = 1, linetype = "dashed") +
    
    scale_x_log10(
      limits = c(min(cox_results$conf.low[cox_results[,group_of_interest]%in%group_levels_show], na.rm = TRUE) * 0.8,
                 x_max * 1.6)
    ) +
    
    scale_color_manual(values = c("ssignif" = "darkred",
                                  "signif" = "goldenrod2",
                                  "ns" = "black")) +
    
    labs(
      x = "Hazard Ratio (log scale)",
      y = "",
      title = paste("Cox Models per", gsub("s$","",group_of_interest),adj_txt)
    ) +
    
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(size = 12), 
      legend.position = "none",
      plot.margin = margin(5.5, 120, 5.5, 5.5)
    ) +
    coord_cartesian(clip = "off")
  
  if(!is.null(cox_file_path)){
    write.xlsx(cox_results,cox_file_path)
  }
  
  return(list("Cox_results"=cox_results,forest_plot = p))
}

# Weighted mutations:
cox_pathways = apply_coxmodels(path_sample,"Weighted_burden","Pathway","Percent_Bait_with_callable_coverage",setdiff(pathways_with_enough_mutations,c("MYC","CASP8","ZNF750","TERT","KMT2D","CHUK")),cox_file_path =file.path(dir_results_script,"04_02_03_forest_plots_pathways_adjusted_source_data.xlsx") )
cox_pathways_unadj = apply_coxmodels(path_sample,"Weighted_burden","Pathway",NULL,pathways_with_enough_mutations,cox_file_path =file.path(dir_results_script,"04_02_03_forest_plots_pathways_unadjusted_source_data.xlsx") )
cox_pathways_anymut = apply_coxmodels(path_sample,"Any_mutation","Pathway","Percent_Bait_with_callable_coverage",setdiff(pathways_with_enough_mutations,c("MYC","CASP8","ZNF750","TERT","KMT2D","CHUK")),cox_file_path =file.path(dir_results_script,"04_02_03_forest_plots_pathways_adjusted_ANY_MUTATION_source_data.xlsx") )

# Cox results for all pathways:
cox_genes = apply_coxmodels(mut_data_drivers,"gene_weight","Genes","Percent_Bait_with_callable_coverage",genes_with_enough_mutations,cox_file_path =file.path(dir_results_script,"04_02_03_forest_plots_genes_adjusted_source_data.xlsx") )
cox_genes_unadj = apply_coxmodels(mut_data_drivers,"gene_weight","Genes",NULL,genes_with_enough_mutations,cox_file_path =file.path(dir_results_script,"04_02_03_forest_plots_genes_unadjusted_source_data.xlsx") )
cox_genes_anymut = apply_coxmodels(mut_data_drivers,"Any_mutation","Genes","Percent_Bait_with_callable_coverage",genes_with_enough_mutations,cox_file_path =file.path(dir_results_script,"04_02_03_forest_plots_genes_adjusted_ANY_MUTATION_source_data.xlsx") )


# Save plots:
# Adjusted plot:
combined_plot <- cox_genes$forest_plot+ cox_pathways$forest_plot
  

ggsave(file.path(dir_results_script,"04_02_03_forest_plots_genes_adjusted.pdf"),
       combined_plot,
       width = 12,
       height = 5,
       dpi = 300)

# Unadjusted plot:
combined_plot <- cox_pathways_unadj$forest_plot +
  cox_genes_unadj$forest_plot

ggsave(file.path(dir_results_script,"04_02_03_forest_plots_genes_not_adjusted.pdf"),
       combined_plot,
       width = 14,
       height = 6,
       dpi = 300)

 
pathways_of_interest=c("CASP8","TERT","KMT2D","ZNF750","CHUK",
                      "Hippo","Notch","P53","Ras_MAPK","Rb","SWI_SNF_PRC2","TGF_BETA")

# Visualize each pathway individually:
results <- path_sample %>%
  filter(Pathway%in%pathways_of_interest)%>%
  group_by(Pathway) %>%
  do(tidy(glm(Met_binary ~ Weighted_burden,
              data = .,
              family = binomial))) %>%
  filter(term == "Weighted_burden")%>%
  mutate(
    label = paste0("p = ", signif(p.value, 3))
  )

path_sample%>%
  filter(Pathway%in%pathways_with_enough_mutations)%>%
  ggplot(
       aes(x = Weighted_burden,
           y = Met_binary)) +
  geom_jitter(aes(col=Met_factor),height = 0.05, alpha = 0.3) +
  geom_smooth(method = "glm",
              method.args = list(family = "binomial"),
              se = TRUE) +
  geom_text(data = results ,
            aes(x = Inf, y = 0.45, label = label),
            inherit.aes = FALSE,
            hjust = 1.5, vjust = 0, size = 3) +
  facet_wrap(~ Pathway, scales = "free_x", nrow = 3, ncol = 5) +
  theme_classic()+xlab("Number of driver mutations \n(weighted for double hits)")



results_any_mut <- path_sample %>%
  filter(Pathway%in%pathways_of_interest)%>%
  group_by(Pathway) %>%
  do(tidy(glm(Met_binary ~ Any_mutation,
              data = .,
              family = binomial))) %>%
  filter(term == "Any_mutationYes")

#---------------------------------------------------
# 3. Compare the number of driver mutations between groups
#---------------------------------------------------

prognostic_pathways = cox_pathways$Cox_results %>%filter(Pathway%in%pathways_with_enough_mutations,p.value<0.1)%>%pull(Pathway)
negative_estimate_pathways = cox_pathways$Cox_results %>%filter(Pathway%in%pathways_with_enough_mutations,estimate<1)%>%pull(Pathway)
plot_data = path_sample %>%
  mutate(
    prognostic_group = ifelse(Pathway %in% prognostic_pathways,
                           "Prognostic",
                           "Other")
  )%>%
  mutate(Oncogenic_contribution = ifelse(Pathway%in%negative_estimate_pathways,-Weighted_burden,Weighted_burden))%>%
  mutate(Oncogenic_contribution = ifelse(Pathway=="KMT2D"&Oncogenic_contribution ==(-2),-1,Oncogenic_contribution ))%>%
  group_by(Sample, prognostic_group) %>%
  dplyr::summarise(
    Weighted_burden = sum(Weighted_burden),
    Oncogenic_contribution = sum(Oncogenic_contribution),
    Unweighted_burden = sum(Any_mutation=="Yes"),
    Met_factor = dplyr::first(Met_factor),
    .groups = "drop"
  )

# Oncogenic score: considering double hits in pathways that are prognostic vs not, and decreasing score if there is a negative association
plot_data%>%filter(prognostic_group=="Prognostic")%>%
  ggplot(aes(x = prognostic_group, y = Oncogenic_contribution, fill = Met_factor)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              alpha = 0.5, size = 1.5)  +
  labs(y = "Oncogenic score", fill = "Group") +
  stat_compare_means(aes(group = Met_factor),
                     method = "wilcox.test",
                     label = "p.format",
                     label.y = max(plot_data$Oncogenic_contribution) + 0.5) +
  theme_classic()+scale_fill_manual(values=group_colors_Mets)+xlab("Pathway")+ylab("Oncogenic score")

# Summing all alterations, accounting for double hit:
plot_data%>%group_by(Sample)%>%
  dplyr::summarise(Weighted_burden=sum(Weighted_burden),
                   Unweighted_burden = sum(Unweighted_burden),
            Met_factor = dplyr::first(Met_factor),
            .groups = "drop"
  )%>%
  ggplot(aes(x=Met_factor, y = Weighted_burden, fill = Met_factor)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
                alpha = 0.5, size = 1.5)  +
    labs(y = "Driver mutation burden (weighted)", fill = "Group") +
    stat_compare_means(
                       method = "wilcox.test",
                       label = "p.format",
                       label.y = max(plot_data$Weighted_burden) + 0.5,  aes(x = 1.5)) +
  theme_classic()+scale_fill_manual(values=group_colors_Mets)+xlab("")

# Summing all alterations, not accounting for double hits:
plot_data%>%group_by(Sample)%>%
  dplyr::summarise(Weighted_burden=sum(Weighted_burden),
                   Unweighted_burden = sum(Unweighted_burden),
                   Met_factor = dplyr::first(Met_factor),
                   .groups = "drop"
  )%>%
  ggplot(aes(x=Met_factor, y = Unweighted_burden, fill = Met_factor)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              alpha = 0.5, size = 1.5)  +
  labs(y = "Driver mutation burden (unweighted)", fill = "Group") +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    label.y = max(plot_data$Weighted_burden) -2,  aes(x = 1.5)) +
  theme_classic()+scale_fill_manual(values=group_colors_Mets)+xlab("")

#---------------------------------------------------
# 3. Forest plot with copy number alterations
#---------------------------------------------------
copy_number_data_fplot <- merge(copy_number_data[,setdiff(colnames(copy_number_data),c("chrXp","chrXq","chrYp","chrYq"))],sample_summary[,c("Sample.ID","Metastasis","FU.metastasis.years","Percent_Bait_with_callable_coverage",IS_classification)],by.x="Tumor_Sample_Barcode",by.y="Sample.ID")
copy_number_data_fplot<-copy_number_data_fplot %>%
  mutate(
    event = ifelse(Metastasis == "Case", 1, 0)
  )

# Identify chromosome arm columns:
chr_cols <- grep("^chr", names(copy_number_data_fplot), value = TRUE)
# Identify chromosome arm with enough alterations:
chr_with_enough_alterations_gain= apply(copy_number_data_fplot %>%select(all_of(chr_cols))%>%as.data.frame(),2,function(x)(sum(x>0)))
chr_with_enough_alterations_loss= apply(copy_number_data_fplot %>%select(all_of(chr_cols))%>%as.data.frame(),2,function(x)(sum(x<0)))
chr_with_enough_alterations_gain = paste(names(chr_with_enough_alterations_gain)[chr_with_enough_alterations_gain>9],"Gain", sep = " - ")
chr_with_enough_alterations_loss = paste(names(chr_with_enough_alterations_loss)[chr_with_enough_alterations_loss>9],"Loss", sep = " - ")
chr_with_enough_alterations_any= apply(copy_number_data_fplot %>%select(all_of(chr_cols))%>%as.data.frame(),2,function(x)(sum(abs(x)>0)))
chr_with_enough_alterations_any=names(chr_with_enough_alterations_any)[chr_with_enough_alterations_any>9]

apply_coxmodels_cnv = function(copy_number_data_fplot,adjustment=NULL,split_cnv=TRUE,cox_file_path=NULL){
  chr_cols <- grep("^chr", names(copy_number_data_fplot), value = TRUE)
  cox_results <- map_df(chr_cols, function(chr) {
    
    if(split_cnv){
      # make CN categorical with 0 as reference
      copy_number_data_fplot[[chr]] <- factor(copy_number_data_fplot[[chr]], levels = c(0, -1, 1))
    }
    
      if(!is.null(adjustment)){
        formula_txt <- paste0(
          "Surv(FU.metastasis.years, event) ~ ",
          chr, " + ", adjustment
        )
    
      }else{
        formula_txt <- paste0(
          "Surv(FU.metastasis.years, event) ~ ",
          chr
        )
    }

    model <- coxph(
        as.formula(formula_txt),
        data = copy_number_data_fplot
    )
      

    
    tidy(model, conf.int = TRUE,exponentiate=TRUE) %>%
      filter(!term %in%c("(Intercept)",adjustment)) %>%
      mutate(
        Chromosome = chr
      )
  })

  if(split_cnv){
    cox_results_clean <- cox_results %>%
      mutate(
        CN_status = case_when(
          grepl("-1", term) ~ "Loss",
          grepl("1", term)  ~ "Gain"
        ),
        Arm = Chromosome,
        
        q_value = p.adjust(p.value, method = "BH"),
        Chr_arm_status = paste(Arm, CN_status, sep = " - ")
      ) %>%
      select(Arm, CN_status, Chr_arm_status,estimate, conf.low, conf.high, 
             p.value, q_value)
    x_max <- cox_results_clean %>%filter(Chr_arm_status %in% c(chr_with_enough_alterations_gain,chr_with_enough_alterations_loss))%>%pull(conf.high) %>%max(na.rm = TRUE)
    x_min = cox_results_clean %>%filter(Chr_arm_status %in% c(chr_with_enough_alterations_gain,chr_with_enough_alterations_loss))%>%pull(conf.high) %>%min(na.rm = TRUE)
    
    p_data = cox_results_clean %>% filter(Chr_arm_status %in% c(chr_with_enough_alterations_gain,chr_with_enough_alterations_loss))%>%
      mutate(Significant = ifelse(p.value < 0.1,"signif","ns"),  
             label = sprintf("%.2f (%.2f–%.2f), p = %.3f",estimate, conf.low, conf.high, p.value),
             Ylab = paste0(Arm,"-",CN_status))%>%
      mutate(Significant= ifelse(p.value<0.05,"ssignif",Significant))
  }else{
    cox_results_clean <- cox_results%>%
      mutate( q_value = p.adjust(p.value, method = "BH"))%>%
      select(Chromosome,estimate,conf.low, conf.high, 
             p.value, q_value)
    x_max <- cox_results_clean %>%filter(Chromosome %in% chr_with_enough_alterations_any)%>%pull(conf.high) %>%max(na.rm = TRUE)
    x_min = cox_results_clean %>%filter(Chromosome %in% chr_with_enough_alterations_any)%>%pull(conf.high) %>%min(na.rm = TRUE)
    
    p_data = cox_results_clean %>% filter(Chromosome %in% chr_with_enough_alterations_any)%>%
      mutate(Significant = ifelse(p.value < 0.1,"signif","ns"),  
             label = sprintf("%.2f (%.2f–%.2f), p = %.3f",estimate, conf.low, conf.high, p.value),
             Ylab = paste0(Chromosome))%>%
      mutate(Significant= ifelse(p.value<0.05,"ssignif",Significant))
  }

 p =p_data%>%
  
  ggplot(aes(x = estimate,
             y = reorder(Ylab, estimate),
             color = Significant)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.2,
                 linewidth = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  # Add text to the right of the CI
  geom_text(aes(x = x_max * 1.4, label = label),
            hjust = 0,
            size = 3.5,
            color = "black") +
  
  geom_vline(xintercept = 1, linetype = "dashed") +
  
  scale_x_log10() +
  
   scale_color_manual(values = c("ssignif" = "darkred",
                                 "signif" = "goldenrod2",
                                 "ns" = "black")) +
  
  labs(
    x = "Hazard Ratio (log scale)",
    y = "",
    title = paste("Cox Models per chr arm")
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(size = 12), 
    legend.position = "none",
    plot.margin = margin(5.5, 120, 5.5, 5.5)
  ) +
  coord_cartesian(clip = "off")
 
 if(!is.null(cox_file_path)){
   write.xlsx(cox_results_clean,cox_file_path)
 }
 
  return(list("Cox_results"=cox_results_clean,forest_plot = p))
}

cox_cnv = apply_coxmodels_cnv(copy_number_data_fplot,adjustment=NULL,cox_file_path=file.path(dir_results_script,"04_02_03_forest_plots_cnv_notadjusted_source_data.xlsx") )
cox_cnv_adjusted = apply_coxmodels_cnv(copy_number_data_fplot,adjustment="Percent_Bait_with_callable_coverage",cox_file_path=file.path(dir_results_script,"04_02_03_forest_plots_cnv_adjusted_source_data.xlsx") )

cox_cnv_nosplit = apply_coxmodels_cnv(copy_number_data_fplot,adjustment=NULL,split_cnv=FALSE,cox_file_path=NULL)
cox_cnv_adjusted_nosplit = apply_coxmodels_cnv(copy_number_data_fplot,adjustment="Percent_Bait_with_callable_coverage",split_cnv=FALSE,cox_file_path=NULL)

# Save plot: 
ggsave(file.path(dir_results_script,"04_02_03_forest_plots_cnv_not_adjusted.pdf"),cox_cnv$forest_plot,width=5,height=7.5)
ggsave(file.path(dir_results_script,"04_02_03_forest_plots_cnv_adjusted.pdf"),cox_cnv_adjusted$forest_plot,width=5,height=7.5)

# Combine both:
prognostic_cnv_alt = cox_cnv_adjusted$Cox_results %>%filter(Chr_arm_status %in%c(chr_with_enough_alterations_gain,chr_with_enough_alterations_loss),p.value<0.1)%>%pull(Chr_arm_status)

copy_number_split <- copy_number_data %>%
  
  # Keep barcode separate
  mutate(across(-Tumor_Sample_Barcode,
                
                # For each chr column, create Gain and Loss
                .fns = list(
                  Gain = ~ as.numeric(. > 0),
                  Loss = ~ as.numeric(. < 0)
                ),
                
                # Rename columns to "chr1p - Gain"
                .names = "{.col} - {.fn}"
  )) %>%
  
  # Remove original chr columns
  select(Tumor_Sample_Barcode, contains(" - ")) %>%
  dplyr::rename(Sample = Tumor_Sample_Barcode)
copy_number_split=merge(copy_number_split,sample_summary[,c("Sample.ID","Metastasis")],by.x="Sample",by.y="Sample.ID")

plot_data_cnv = copy_number_split %>%
  mutate(
    Metastasis = factor(Metastasis,
                        levels = c("Control","Case")))%>%
  pivot_longer(
    cols = -c(Sample,Metastasis),
    names_to = "Chromosome_Event",
    values_to = "Status"
  )%>%
  mutate(
    prognostic_group = ifelse(Chromosome_Event %in% prognostic_cnv_alt ,
                           "Prognostic",
                           "Other"))%>%
  #mutate(Weighted_burden = ifelse(Pathway%in%negative_estimate_pathways,-Weighted_burden,Weighted_burden))%>%
  group_by(Sample, prognostic_group) %>%
  dplyr::summarise(
    Weighted_burden = sum(Status),
    Oncogenic_contribution = sum(Status),
    Unweighted_burden = sum(Status),
    Met_factor = dplyr::first(Metastasis),
    .groups = "drop"
  )

ggplot(plot_data_cnv,aes(x = Met_factor, y = Weighted_burden, fill = Met_factor)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              alpha = 0.5, size = 1.5)  +
  labs(y = "Weighted mutation burden", fill = "Group") +
  stat_compare_means(aes(group = Met_factor, x = 1.5),
                     method = "wilcox.test",
                     label = "p.format",
                     label.y = max(plot_data_cnv$Weighted_burden) + 0.5) +
  theme_classic()+scale_fill_manual(values=group_colors_Mets)+xlab("Pathway")+ylab("Sum of chromosome alterations")


ggplot(plot_data_cnv,aes(x = prognostic_group, y = Weighted_burden, fill = Met_factor)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              alpha = 0.5, size = 1.5)  +
  labs(y = "Weighted mutation burden", fill = "Group") +
  stat_compare_means(aes(group = Met_factor),
                     method = "wilcox.test",
                     label = "p.format",
                     label.y = max(plot_data_cnv$Weighted_burden) + 0.5) +
  theme_classic()+scale_fill_manual(values=group_colors_Mets)+xlab("Pathway")+ylab("Sum of chromosome alterations")

# Combine both:
common_samples = intersect(plot_data$Sample,plot_data_cnv$Sample)
combined_plot_df = bind_rows(plot_data, plot_data_cnv) %>%
  group_by(Sample, prognostic_group) %>%
  dplyr::summarise(
    Weighted_burden = sum(Weighted_burden, na.rm = TRUE),
    Oncogenic_contribution = sum(Oncogenic_contribution, na.rm = TRUE),
    Met_factor = dplyr::first(Met_factor),
    .groups = "drop"
  )%>%
  filter(Sample%in%common_samples)

combined_plot_df %>%
  filter(prognostic_group=="Prognostic")%>%
  ggplot(aes(x = prognostic_group, y = Oncogenic_contribution, fill = Met_factor)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              alpha = 0.5, size = 1.5)  +
  labs(y = "Weighted mutation burden", fill = "Group") +
  stat_compare_means(aes(group = Met_factor),
                     method = "wilcox.test",
                     label = "p.format",
                     label.y = max(combined_plot_df$Oncogenic_contribution) + 0.5) +
  theme_classic()+scale_fill_manual(values=group_colors_Mets)+xlab("Pathway")+ylab("Sum of significant alterations")

#---------------------------------------------------
# 4. Tileplot
#---------------------------------------------------

path_sample_combination = path_sample %>%
  select(Sample,Pathway,Weighted_burden)%>%
  dplyr::rename(Alteration = Pathway,
                Value = Weighted_burden)

copy_nr_data_combination <- copy_number_data %>%
  pivot_longer(
    cols = -Tumor_Sample_Barcode,
    names_to = "arm",
    values_to = "value"
  ) %>%
  mutate(
    chr = str_extract(arm, "\\d+"),    # extract the number
    arm_type = str_extract(arm, "[pq]") # p or q
  )

chr_order = copy_nr_data_combination %>%distinct(chr, arm_type,arm) %>%
  arrange(as.integer(chr), arm_type) %>%pull(arm)

copy_nr_data_combination = copy_nr_data_combination%>%
  mutate(
    chr = as.integer(chr),
    arm = factor(arm, levels = chr_order))%>%
  mutate(
    value = as.factor(value)   # <— REQUIRED for scale_fill_manual()
  )%>%
  dplyr::rename(Alteration = arm,
                Value = value,
                Sample = Tumor_Sample_Barcode)%>%
  select(Sample,Alteration,Value)

combined_df_tileplot = rbind(path_sample_combination,copy_nr_data_combination)

# We only want samples for which we trust both driver mutations and copy number alterations:
combined_df_tileplot = combined_df_tileplot %>%
  filter(Sample %in%common_samples)

# Order according to Hazard ratios:
order_by_hr = cox_cnv_adjusted$Cox_results %>%
  select(-c(Arm,CN_status))%>%
  dplyr::rename(Pathway = Chr_arm_status)%>%
  bind_rows(cox_pathways$Cox_results)%>%
  dplyr::rename(Alteration = Pathway)%>%
  filter(p.value<0.1,Alteration %in% c(pathways_with_enough_mutations,chr_with_enough_alterations_gain,chr_with_enough_alterations_loss))%>%
  arrange(estimate)%>%
  pull(Alteration)

combined_df_tileplot %>%
  filter(Alteration%in%c(prognostic_pathways,gsub(" -.*","",prognostic_cnv_alt)))%>%
  mutate(Alteration = factor(Alteration,levels = unique(gsub(" -.*","",order_by_hr))))%>%
  ggplot() +
  geom_tile(
    aes(y = Alteration, x = Sample, fill = Value),
    color = "black",
    lwd = 0.25,
    width = 1,
    height = 1,
    na.rm = TRUE
  ) +
  scale_fill_manual(
    values = c(
      "-1" = "dodgerblue3",   # Del
      "1"  = "indianred",# Amp or 1 hit
      "2"  = "firebrick1", 
      "3"  = "firebrick3",
      "4"  = "firebrick4",
      "5"  = "coral4",
      "0"  = "white"
    ),
    na.value = "white"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  xlab("") +
  ylab("") 
