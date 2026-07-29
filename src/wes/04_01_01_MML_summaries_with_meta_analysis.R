#---------------------------------------------------
# Aim: Summarize mutations in Chang and shain meta-analysis + our cohort
# Author: B. Rentroia Pacheco
# Input: Meta-analysis and EMC's MM lists
# Output: Plots summarizing mutation burdens
#---------------------------------------------------

#-----------------------
# 1. Summarize Tumor mutational burden and fraction of UV mutations in meta-analysis and EMC cohort
#-----------------------
# First we need to harmonize the datasets:
sample_summary_plot_reduced = sample_summary_plot %>%
  select(Sample.ID,Tumor.cellularity.avg.pct,`Burden.(Mutations/megabase).(MB)`,Percent_Bait_with_callable_coverage,.data[[IS_classification]])%>%
  rename(Immune.Status = !!sym(IS_classification))%>%
  mutate(Cohort="EMC")%>%
  as.data.frame()

chang_plot_reduces_sample_summary = original_data_chang%>%
  filter(RDEB!="Yes")%>%
  select(Sample,`Tumor.Cellularity.(%)`,`Burden.(Mutations/.megabase)`,`Fraction.of.Footprint.to.Target.Territory.for.Reference/Tumor.Pair.(%)`,`Burden.(Mutations/.megabase)`,Immune.Status)%>%
  mutate(Cohort="Meta_analysis",
         Immune.Status = ifelse(Immune.Status=="Suppressed","Yes","No"))%>%
  rename(Percent_Bait_with_callable_coverage = `Fraction.of.Footprint.to.Target.Territory.for.Reference/Tumor.Pair.(%)`,
         Sample.ID = Sample,
         Tumor.cellularity.avg.pct = `Tumor.Cellularity.(%)`,
         `Burden.(Mutations/megabase).(MB)`=`Burden.(Mutations/.megabase)`)%>%
  arrange(Immune.Status, desc(`Burden.(Mutations/megabase).(MB)`)) %>%
  as.data.frame()
  
sample_summary_chang_EMC = rbind(sample_summary_plot_reduced,chang_plot_reduces_sample_summary)%>%
  arrange(Immune.Status, desc(`Burden.(Mutations/megabase).(MB)`)) %>%
  mutate(Sample_ID_plot = row_number())%>%
  mutate(Sample_ID_plot=factor(Sample_ID_plot, levels = Sample_ID_plot))%>%
  as.data.frame()


# Top plot: Mutation burden
plot_burden_2cohorts <- ggplot(sample_summary_chang_EMC, aes(x = Sample_ID_plot, y = `Burden.(Mutations/megabase).(MB)`,fill=Cohort)) +
  geom_col() +
  #scale_fill_manual(values = group_colors_Mets) +
  labs(y = "Mutation Burden (Mut/Mb)", x = NULL) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top"
  )
plot_burden_2cohorts


my_comparisons_2gr <- list(
  c("No","Yes")
)

p_boxplot_IS_groups_2cohorts = ggplot(sample_summary_chang_EMC,
                             aes(x = Immune.Status,
                                 y = `Burden.(Mutations/megabase).(MB)`,
                                 fill = Immune.Status)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1)) +
  scale_fill_manual(values = group_colors_IS) +
  labs(y = "Mutation Burden (Mut/Mb)", x = NULL) +
  theme_bw() +
  stat_compare_means(comparisons = my_comparisons_2gr,
                     method = "wilcox.test",
                     label = "p.format")
p_boxplot_IS_groups_2cohorts + facet_grid(~Cohort)

sample_summary_chang_EMC%>%
  mutate(Percent_Bait_binary = ifelse(Percent_Bait_with_callable_coverage>95,">95%","<=95%"))%>%
  ggplot(aes(x=interaction(Immune.Status,Percent_Bait_binary), fill=Immune.Status,y=`Burden.(Mutations/megabase).(MB)`))+
  geom_boxplot(outliers = FALSE)+geom_jitter(position = position_jitterdodge(jitter.width = 0.1)) +xlab("Callable coverage(%)")+ylab("Mutation burden (Mut/MB)") + theme_bw() + facet_grid(~Cohort)+stat_compare_means(comparisons = my_comparisons_4gr_pct,
                                                                                                                                                                                                 method = "wilcox.test",
                                                                                                                                                                                                  label = "p.format")
