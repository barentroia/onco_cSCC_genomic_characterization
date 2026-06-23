#---------------------------------------------------
# Aim: Summarize MML generated
# Author: B. Rentroia Pacheco
# Input: Folders with MML lists
# Output: Summary files of MMLs
#---------------------------------------------------

#-----------------------
# 0. Setup
#-----------------------
# source(file.path(dir_scripts,"Pub_00_required_libraries_WES.R"))
#batch_nr = 1
#inspection_folders = c("2024_08_06_MML_for_inspection","2024_08_07_MML_for_inspection","2024_08_26_MML_for_inspection","2024_08_30_MML_for_inspection","2024_09_06_MML_for_inspection")
#batch_info = read.xlsx("T:/Barbara/WES_analyses/Methods/WES_cSCC_Analysis_progress.xlsx")
#tert_mutations_summary = read.xlsx("T:/Barbara/WES_analyses/Results/Manual_inspections/Indels_pathogenic/2024_08_14_TERT/2024_08_14_TERT_mutations.xlsx")
#sample_summary = read.xlsx("T:/Barbara/WES_analyses/Results/Meta_data/DA_04_Sample_summary_extended_version_incl_exclusion.xlsx")
dir.create(dir_results_script)


#-----------------------
# 1. Summarize Tumor mutational burden and fraction of UV mutations:
#-----------------------

# Create IS group and arrange
sample_summary_plot <- sample_summary %>%
  mutate(Group = paste(Metastasis, paste0("IS_",.data[[IS_classification]]), sep = "_")) %>%
  mutate(Group = factor(Group,levels=c("Case_IS_No","Control_IS_No","Case_IS_Yes","Control_IS_Yes")))%>%
  #arrange(Group, desc(`Burden.(Mutations/megabase).(MB)`)) %>%
  arrange(.data[[IS_classification]], desc(`Burden.(Mutations/megabase).(MB)`)) %>%
  mutate(Sample_ID = row_number())


# Make Sample_ID a factor with fixed levels to preserve order
sample_summary_plot$Sample_ID <- factor(sample_summary_plot$Sample_ID, levels = sample_summary_plot$Sample_ID)


# Top plot: Mutation burden
plot_burden <- ggplot(sample_summary_plot, aes(x = Sample_ID, y = `Burden.(Mutations/megabase).(MB)`,fill=Metastasis)) +
  geom_col() +
  scale_fill_manual(values = group_colors_Mets) +
  labs(y = "Mutation Burden (Mut/Mb)", x = NULL) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top"
  )
plot_burden

ggsave(file.path(dir_results_script,"S04_01_01_01_Mutation_burden.pdf"),plot_burden,width=15,height=6)
ggsave(file.path(dir_results_script,"S04_01_01_01_Mutation_burden.png"),plot_burden,width=15,height=6)


# Prepare long-format annotation data
plot_burden <- ggplot(sample_summary_plot, aes(x = Sample_ID, y = `Burden.(Mutations/megabase).(MB)`)) +
  geom_col() +
  #scale_fill_manual(values = group_colors_Mets) +
  labs(y = "Mutation Burden (Mut/Mb)", x = NULL) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top"
  )
annotation_df <- sample_summary_plot %>%
  select(Sample_ID, Metastasis) %>%
  pivot_longer(cols = c(Metastasis),
               names_to = "Annotation",
               values_to = "Group") %>%
  mutate(Annotation = factor(Annotation, levels = c( "Metastasis")))

# Plot tiles
p_tile <- ggplot(annotation_df, aes(x = Sample_ID, y = Annotation, fill = Group)) +
  geom_tile(width = 0.9, height = 0.9) +
  scale_fill_manual(values = c(group_colors_Mets)) +
  theme_void() +                  # remove axes
  theme(
    legend.position = "none"      # remove legend if you prefer
  )

plot_burden_monochromatic = plot_burden / p_tile + plot_layout(heights = c(7, 1))

ggsave(file.path(dir_results_script,"S04_01_01_01_Mutation_burden_monochr.pdf"),plot_burden_monochromatic,width=15,height=6)
ggsave(file.path(dir_results_script,"S04_01_01_01_Mutation_burden_monochr.png"),plot_burden_monochromatic,width=15,height=6)

# Bottom plot: % UV mutations
#plot_uv <- ggplot(sample_summary_plot, aes(x = Sample_ID, y = `Dominant.clone.UV.mutations.(%)`)) +
#  geom_col() +
#  scale_fill_manual(values = group_colors) +
#  labs(y = "% UV Mutations", x = "Samples") +
#  theme_bw() +
#  theme(
#    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
#    axis.ticks.x = element_blank(),
#    legend.position = "none"
#  )+ylim(c(0,100))
#ggarrange(plot_burden,plot_uv,nrow=2)

# Same but with boxplots, and statistical comparisons:
# Comparison of Mutation burden:
my_comparisons_2gr_mets <- list(
  c("Control","Case")
)
p_boxplot_mets_groups = ggplot(sample_summary_plot,
                             aes(x = Metastasis,
                                 y = `Burden.(Mutations/megabase).(MB)`))+
                                 #fill = Metastasis)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 1)  +
  labs(y = "Mutation Burden (Mut/Mb)", x = NULL) +
  theme_bw(base_size = 18)+
  stat_compare_means(comparisons = my_comparisons_2gr_mets,
                     method = "wilcox.test",
                     label = "p.format") +ylim(c(0,270))
ggsave(file.path(dir_results_script,"S04_01_01_01_A1_Boxplot_Mutation_burden_2gr_Mets.pdf"),p_boxplot_mets_groups,width=5,height=6)
ggsave(file.path(dir_results_script,"S04_01_01_01_A1_Boxplot_Mutation_burden_2gr_Mets.png"),p_boxplot_mets_groups,width=5,height=6)

# Comparison IS vs IC
my_comparisons_2gr <- list(
  c("No","Yes")
)

p_boxplot_IS_groups = ggplot(sample_summary_plot,
       aes(x = .data[[IS_classification]],
           y = `Burden.(Mutations/megabase).(MB)`))+
           #fill = .data[[IS_classification]])) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 1) +
  #scale_fill_manual(values = group_colors_IS,name = "Immunocompromised" ) +
  labs(y = "Mutation Burden (Mut/Mb)", x = NULL) +
  theme_bw(base_size = 18) +ylim(c(0,270))


# 4 groups:
my_comparisons_4gr <- list(
  c("Case_IS_No", "Control_IS_No"),
  c("Case_IS_Yes", "Control_IS_Yes")
)

p_boxplot_4_groups = ggplot(sample_summary_plot,
                            aes(x = Group,
                                y = `Burden.(Mutations/megabase).(MB)`,
                                fill = Group)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 1)  +
  scale_fill_manual(values = group_colors_Mets_IS) +
  labs(y = "Mutation Burden (Mut/Mb)", x = NULL) +
  theme_bw(base_size = 18) 

for(type_test in c("wilcox","adjusted")){
  if(type_test =="wilcox"){
    p_boxplot_IS_groups_t = p_boxplot_IS_groups +
      stat_compare_means(comparisons = my_comparisons_2gr,
                         method = "wilcox.test",
                         label = "p.format")
    
    p_boxplot_4_groups_t = p_boxplot_4_groups +
      stat_compare_means(comparisons = my_comparisons_4gr,
                         method = "wilcox.test",
                         label = "p.format")
    
    
  }else if(type_test =="adjusted"){
    y_max <- max(sample_summary_plot$`Burden.(Mutations/megabase).(MB)`)
    
    formula_IS <- as.formula(
      paste0("`Burden.(Mutations/megabase).(MB)` ~ ",
             IS_classification,
             " + Percent_Bait_with_callable_coverage")
    )
    
    model_IS <- lm(formula_IS, data = sample_summary_plot)
    p_group_IS <- summary(model_IS)$coefficients[paste0(IS_classification,"Yes"), "Pr(>|t|)"]
    p_group_formatted_IS <- format.pval(p_group_IS, digits = 3)
    
    p_boxplot_IS_groups_t = p_boxplot_IS_groups +
      annotate("text",
               x = mean(c(1,2)),   
               y = y_max * 1.05,
               label = paste0("Adj p = ", p_group_formatted_IS),
               size = 3.5)
    #p_boxplot_IS_groups_t
    
    # 4 groups: 
    model_IS_no <- lm(`Burden.(Mutations/megabase).(MB)` ~ Metastasis +Percent_Bait_with_callable_coverage,
                data = sample_summary_plot[which(sample_summary_plot[,IS_classification]=="No"),])
    p_group_IS_no <- summary(model_IS_no)$coefficients["MetastasisControl", "Pr(>|t|)"]
    p_group_formatted_IS_no <- format.pval(p_group_IS_no, digits = 3)
    
    model_IS_yes <- lm(`Burden.(Mutations/megabase).(MB)` ~ Metastasis +Percent_Bait_with_callable_coverage,
                      data = sample_summary_plot[which(sample_summary_plot[,IS_classification]=="Yes"),])
    p_group_IS_yes <- summary(model_IS_yes)$coefficients["MetastasisControl", "Pr(>|t|)"]
    p_group_formatted_IS_yes <- format.pval(p_group_IS_yes, digits = 3)
    
   
    
    p_boxplot_4_groups_t = p_boxplot_4_groups +
      annotate("text",
               x = mean(c(1,2)),   # middle of first two groups (IS = No)
               y = y_max * 1.05,
               label = paste0("Adj p = ", p_group_formatted_IS_no),
               size = 3.5) +
      annotate("text",
               x = mean(c(3,4)),   # middle of last two groups (IS = Yes)
               y = y_max * 1.05,
               label = paste0("Adj p = ", p_group_formatted_IS_yes),
               size = 3.5)
  }
  # 2 groups:
  ggsave(file.path(dir_results_script,paste0("S04_01_01_01_A1_Boxplot_Mutation_burden_2gr_",IS_classification,"_",type_test,".pdf")),p_boxplot_IS_groups_t ,width=5,height=6)
  ggsave(file.path(dir_results_script,paste0("S04_01_01_01_A1_Boxplot_Mutation_burden_2gr_",IS_classification,"_",type_test,".png")),p_boxplot_IS_groups_t ,width=5,height=6)
  # 4 groups:
  ggsave(file.path(dir_results_script,paste0("S04_01_01_01_A1_Boxplot_Mutation_burden_4gr_",IS_classification,"_",type_test,".pdf")),p_boxplot_4_groups_t ,width=8,height=5)
  ggsave(file.path(dir_results_script,paste0("S04_01_01_01_A1_Boxplot_Mutation_burden_4gr_",IS_classification,"_",type_test,".png")),p_boxplot_4_groups_t,width=8,height=5)
}


# % of UV mutations:

p_boxplot_4_groups_UV_muts = ggplot(sample_summary_plot, aes(y = `Dominant.clone.UV.mutations.(%)` , x = Group,fill = Group)) +
  geom_boxplot(outliers = FALSE) +geom_jitter(position = position_jitterdodge(jitter.width = 0.1))+
  scale_fill_manual(values = group_colors_Mets_IS) +
  labs(y = "% UV Mutations", x = NULL) +
  theme_bw()+
  stat_compare_means(comparisons = my_comparisons_4gr,
                     method = "wilcox.test",
                     label = "p.format")

ggsave(file.path(dir_results_script,paste0("S04_01_01_01_A2_Boxplot_UV_muts_4gr_",IS_classification,".pdf")),p_boxplot_4_groups_UV_muts,width=15,height=6)
ggsave(file.path(dir_results_script,paste0("S04_01_01_01_A2_Boxplot_UV_muts_4gr_",IS_classification,".png")),p_boxplot_4_groups_UV_muts,width=15,height=6)

# Plot relationship with tumor cellularity
p_scatter_tcel_MB = sample_summary%>%
  ggplot(aes(x=Tumor.cellularity.avg.pct,y=`Burden.(Mutations/megabase).(MB)`))+
  geom_point()+xlab("Tumor cellularity (%)")+ylab("Mutation burden (Mut/MB)") +
  stat_cor(method = "spearman",
           label.x = 0.1 * max(sample_summary$Tumor.cellularity.avg.pct),
           label.y = 0.9 * max(sample_summary$`Burden.(Mutations/megabase).(MB)`))+geom_smooth(method = "loess", se = FALSE)+theme_bw()
ggsave(file.path(dir_results_script,paste0("S04_01_01_01_A3_Scatterplot_Tcel_MB.pdf")),p_scatter_tcel_MB,width=10,height=6)
ggsave(file.path(dir_results_script,paste0("S04_01_01_01_A3_Scatterplot_Tcel_MB.png")),p_scatter_tcel_MB,width=10,height=6)

# Plot relationship with percentage of callable coverage:
p_scatter_CCoverage_MB = sample_summary%>%
  ggplot(aes(x=Percent_Bait_with_callable_coverage,y=`Burden.(Mutations/megabase).(MB)`))+
  geom_point()+xlab("Callable coverage(%)")+ylab("Mutation burden (Mut/MB)") +
  stat_cor(method = "spearman",
           label.x = 0.9 * max(sample_summary$Percent_Bait_with_callable_coverage),
           label.y = 0.9 * max(sample_summary$`Burden.(Mutations/megabase).(MB)`))+geom_smooth(method = "loess", se = FALSE)+theme_bw()
ggsave(file.path(dir_results_script,paste0("S04_01_01_A2_Scatterplot_CallableCov_MB.pdf")),p_scatter_CCoverage_MB ,width=10,height=6)
ggsave(file.path(dir_results_script,paste0("S04_01_01_A2_Scatterplot_CallableCov_MB.png")),p_scatter_CCoverage_MB ,width=10,height=6)

# Stratify per group:
group_vars <- c("Metastasis", IS_classification)
for (grp in group_vars) {
  if(grp=="Metastasis"){
    color_scale_plot = group_colors_Mets
  }else{
    color_scale_plot = group_colors_IS
  }
  
  # Tumor cellularity:
  p <- sample_summary %>%
    ggplot(aes(x = Tumor.cellularity.avg.pct,
               y = `Burden.(Mutations/megabase).(MB)`)) +
    geom_point(aes(color = .data[[grp]])) +
    geom_smooth(aes(color = .data[[grp]]),method = "loess", se = FALSE) +
    scale_color_manual(values =  color_scale_plot) +
    stat_cor(method = "spearman",
             label.x = 0.1 * max(sample_summary$Tumor.cellularity.avg.pct),
             label.y = 0.9 * max(sample_summary$`Burden.(Mutations/megabase).(MB)`)) +
    xlab("Tumor cellularity (%)") +
    ylab("Mutation burden (Mut/MB)") +
    theme_bw()
  
  # Save PDF
  ggsave(
    file.path(dir_results_script,
              paste0("S04_01_01_01_A3_Scatterplot_Tcel_MB_by_", grp, ".pdf")),
    p, width = 10, height = 6
  )
  
  # Save PNG
  ggsave(
    file.path(dir_results_script,
              paste0("S04_01_01_01_A3_Scatterplot_Tcel_MB_by_", grp, ".png")),
    p, width = 10, height = 6
  )
  
  # Percent Bait:
  p_bait = sample_summary%>%
    ggplot(aes(x=Percent_Bait_with_callable_coverage,y=`Burden.(Mutations/megabase).(MB)`))+
    geom_point(aes(color = .data[[grp]])) +
    scale_color_manual(values =  color_scale_plot) +
    xlab("Callable coverage(%)")+ylab("Mutation burden (Mut/MB)") +
    geom_vline(xintercept = 97.5, linetype = "dashed", linewidth = 0.5, color = "black")+
    stat_cor(method = "spearman",
             label.x = 0.9 * max(sample_summary$Percent_Bait_with_callable_coverage),
             label.y = 0.9 * max(sample_summary$`Burden.(Mutations/megabase).(MB)`))+geom_smooth(aes(color = .data[[grp]]),method = "loess", se = FALSE)+theme_bw()
  
  ggsave(
    file.path(dir_results_script,
              paste0("S04_01_01_A2_Scatterplot_CallableCov_MB_by_", grp, ".pdf")),
    p_bait, width = 10, height = 6
  )
  
  # Save PNG
  ggsave(
    file.path(dir_results_script,
              paste0("S04_01_01_A2_Scatterplot_CallableCov_MB_by_", grp, ".png")),
    p_bait, width = 10, height = 6
  )

}
  
my_comparisons_4gr_pct_mets <- list(
  c("Case.<=95%", "Control.<=95%"),
  c("Case.>95%", "Control.>95%")
)

p_boxplot_CCoverage_MB_Mets = sample_summary%>%
  mutate(Percent_Bait_binary = ifelse(Percent_Bait_with_callable_coverage>95,">95%","<=95%"))%>%
  ggplot(aes(x=interaction(Metastasis,Percent_Bait_binary), fill=Metastasis,y=`Burden.(Mutations/megabase).(MB)`))+
  geom_boxplot(outliers = FALSE)+geom_jitter(position = position_jitterdodge(jitter.width = 0.1))+xlab("Callable coverage(%)")+ylab("Mutation burden (Mut/MB)") + theme_bw()+stat_compare_means(comparisons = my_comparisons_4gr_pct_mets,
                                                                                                                                                                                method = "wilcox.test",
                                                                                                                                                                                label = "p.format")
  
ggsave(file.path(dir_results_script,paste0("S04_01_01_A2_Boxplot_CallableCov_MB_Mets.pdf")),p_boxplot_CCoverage_MB_Mets,width=10,height=6)
ggsave(file.path(dir_results_script,paste0("S04_01_01_A2_Boxplot_CallableCov_MB_Mets.png")),p_boxplot_CCoverage_MB_Mets,width=10,height=6)

my_comparisons_4gr_pct <- list(
  c("No.<=95%", "Yes.<=95%"),
  c("No.>95%", "Yes.>95%")
)

p_boxplot_CCoverage_MB_IS = sample_summary%>%
  mutate(Percent_Bait_binary = ifelse(Percent_Bait_with_callable_coverage>95,">95%","<=95%"))%>%
  ggplot(aes(x=interaction(.data[[IS_classification]],Percent_Bait_binary), fill=.data[[IS_classification]],y=`Burden.(Mutations/megabase).(MB)`))+
  geom_boxplot(outliers = FALSE)+geom_jitter(position = position_jitterdodge(jitter.width = 0.1)) +xlab("Callable coverage(%)")+ylab("Mutation burden (Mut/MB)") + theme_bw() +stat_compare_means(comparisons = my_comparisons_4gr_pct,
                                                                                                               method = "wilcox.test",
                                                                                                               label = "p.format")

ggsave(file.path(dir_results_script,paste0("S04_01_01_A2_Boxplot_CallableCov_MB_IS.pdf")),p_boxplot_CCoverage_MB_IS,width=10,height=6)
ggsave(file.path(dir_results_script,paste0("S04_01_01_A2_Boxplot_CallableCov_MB_IS.png")),p_boxplot_CCoverage_MB_IS,width=10,height=6)
