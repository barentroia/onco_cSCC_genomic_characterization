#---------------------------------------------------
# Aim: Integrate WES and RNA-seq results
# Author: B. Rentroia Pacheco
# Input: Results from WES
# Output: Integration
#---------------------------------------------------

#---------------------------------------------------
# Load genomic clusters:
#---------------------------------------------------
genomic_clusters = read.xlsx(file.path(dir_results,"08_RNAseq_integration","SfigXA_PCA.xlsx"))
# Add to genomic_df_summary:
genomic_df_summary_wRNAseq = merge(genomic_df_summary,genomic_clusters[,c("Sample_id","cluster")],by.x="SkylineDx.ID",by.y="Sample_id")

#---------------------------------------------------
# 1. Explore associations:
#---------------------------------------------------
# Veriables of interest: All

variables_of_int = pathways_with_enough_mutations #c("MYC","CASP8","ZNF750","KMT2D","Hippo","Notch","P53","Ras_MAPK","Rb","SWI_SNF_PRC2","TGF_BETA")
# 1. Set order of cluster
genomic_df_summary_wRNAseq <- genomic_df_summary_wRNAseq%>%
  mutate(cluster = ifelse(cluster =="Well-differentiated","Differentiated",cluster))%>%
  mutate(cluster = ifelse(cluster =="Moderately differentiated","Basal-like",cluster))%>%
  mutate(cluster = ifelse(cluster =="Poorly differentiated","Mesenchymal-like",cluster))%>%
  mutate(cluster = factor(cluster,
                          levels = c("Differentiated",
                                     "Basal-like",
                                     "Mesenchymal-like")))
# 2. Pivot to long format
df_long <- genomic_df_summary_wRNAseq %>%
  select(cluster, Percent_Bait_with_callable_coverage,all_of(variables_of_int)) %>%
  pivot_longer(
    cols = -c(cluster,Percent_Bait_with_callable_coverage),
    names_to = "variable",
    values_to = "value"
  )

# P-values computed with Poisson regression, adjusted by callable coverage:
pvals_df_poiss <- df_long %>%
  group_by(variable) %>%
  do({
    model <- glm(value ~ Percent_Bait_with_callable_coverage+cluster, data = ., family = "poisson")
    model_null<- glm(value ~ Percent_Bait_with_callable_coverage, data = ., family = "poisson")
    # global p-value for cluster effect
    pval <- anova(model, test = "Chisq")$`Pr(>Chi)`[3]
    #pval_null<-anova(model_null,model)$`Pr(>Chi)`[2]
    data.frame(p_value = pval)#,pvalue_null =pval_null)
  }) %>%
  ungroup() %>%
  mutate(label = paste0("p = ", signif(p_value, 3)))

# Pathways
# Compute max per variable
label_positions <- df_long %>%
  group_by(variable) %>%
  summarise(y_pos = max(value, na.rm = TRUE) * 1.2, .groups = "drop")

# Merge positions with p-values
pvals_df <- left_join(pvals_df_poiss , label_positions, by = "variable")

p1=ggplot(df_long, aes(x = cluster, y = value, fill = cluster)) +
  geom_violin(trim = FALSE, scale = "width")+
  geom_jitter(width = 0.15,height = 0, alpha = 0.4, size = 0.8) +
  facet_wrap(~ variable, scales = "free_y") +
  geom_text(
    data = pvals_df,
    aes(x = 0.8, y = y_pos, label = label),  # x=2 centers above middle cluster
    inherit.aes = FALSE,
    size = 3
  ) +
  scale_fill_manual(values = c(
    "Differentiated" = "wheat2",
    "Basal-like" = "tan1",
    "Mesenchymal-like" = "tan4"
  )) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave(paste0(dir_results_script,"/04_04_alterations_per_gcluster_violin.pdf"),p1,width=10,height=10)

# Now check summary variables, using linear regressions:
variables_of_int = c("pct_altered","Dominant.clone.UV.mutations.(%)","Tumor.cellularity.avg.pct","Burden.(Mutations/megabase).(MB)")
df_long_linreg <- genomic_df_summary_wRNAseq %>%
  select(cluster, Percent_Bait_with_callable_coverage,all_of(variables_of_int)) %>%
  pivot_longer(
    cols = -c(cluster,Percent_Bait_with_callable_coverage),
    names_to = "variable",
    values_to = "value"
  )

pvals_df_linreg <- df_long_linreg %>%
  group_by(variable) %>%
  do({
    model <- glm(value ~ Percent_Bait_with_callable_coverage+cluster, data = .)
    
    # global p-value for cluster effect
    pval <- anova(model, test = "Chisq")$`Pr(>Chi)`[3]
    
    data.frame(p_value = pval)
  }) %>%
  ungroup() %>%
  mutate(label = paste0("p = ", signif(p_value, 3)))

label_positions <- df_long_linreg %>%
  group_by(variable) %>%
  summarise(y_pos = max(value, na.rm = TRUE) * 1.2, .groups = "drop")

# Merge positions with p-values
pvals_df <- left_join(pvals_df_linreg, label_positions, by = "variable")


p2=ggplot(df_long_linreg, aes(x = cluster, y = value, fill = cluster)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_jitter(width = 0.15, alpha = 0.4, size = 0.8) +
  facet_wrap(~ variable, scales = "free_y") +
  geom_text(
    data = pvals_df,
    aes(x = 0.7, y = y_pos, label = label),  # x=2 centers above middle cluster
    inherit.aes = FALSE,
    size = 3
  )+
  scale_fill_manual(values = c(
    "Differentiated" = "wheat2",
    "Basal-like" = "tan1",
    "Mesenchymal-like" = "tan4"
  )) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
ggsave(paste0(dir_results_script,"/04_04_summary_stats_per_gcluster_violin.pdf"),p2,width=10,height=10)

# Save p-values:
rbind(pvals_df_poiss,pvals_df_linreg)%>%
  write.xlsx(paste0(dir_results_script,"/04_04_gassociation_tests_gclusters.xlsx"))

# Copy number alterations:
# Example: one gene/alteration column
chr_cols <- grep("^chr", colnames(genomic_df_summary), value = TRUE)
chr_order <- chr_cols %>%
  # create a tibble with chromosome number and arm
  tibble::tibble(chr = .) %>%
  mutate(
    chr_num = as.numeric(str_extract(chr, "(?<=chr)\\d+")),
    arm = str_extract(chr, "[pq]$")
  ) %>%
  arrange(chr_num, arm) %>%
  pull(chr)%>%setdiff(c("chrXp","chrXq","chrYp","chrYq"))

df_long <- genomic_df_summary_wRNAseq %>%
  select(cluster, Percent_Bait_with_callable_coverage,pct_altered,all_of(chr_order)) %>%
  pivot_longer(
    cols = -c(cluster,Percent_Bait_with_callable_coverage,pct_altered),
    names_to = "gene",
    values_to = "alteration"
  ) %>%
  mutate(alteration = factor(alteration, levels = c(-1, 0, 1),
                             labels = c("Loss", "Neutral", "Gain")),
         gene = factor(gene, levels = chr_order)) 

# Code gains and losses separately:
df_long_bin <- df_long %>%
  mutate(gain = ifelse(alteration == "Gain", 1, 0),
         loss =  ifelse(alteration == "Loss", 1, 0))

pvals_logreg = df_long_bin%>% 
  group_by(gene) %>%summarise(
  pvalue_gain_adjCallableCov = {
    fit <- glm(gain ~ Percent_Bait_with_callable_coverage+cluster, family = binomial)
    # overall effect of cluster (not individual coefficients)
    anova(fit, test = "Chisq")["cluster", "Pr(>Chi)"]
  },
  pvalue_loss_adjCallableCov = {
    fit <- glm(loss ~ Percent_Bait_with_callable_coverage+cluster, family = binomial)
    # overall effect of cluster (not individual coefficients)
    anova(fit, test = "Chisq")["cluster", "Pr(>Chi)"]
  },
  pvalue_gain_PctAltered = {
    fit <- glm(gain ~ pct_altered+cluster, family = binomial)
    # overall effect of cluster (not individual coefficients)
    anova(fit, test = "Chisq")["cluster", "Pr(>Chi)"]
  },
  pvalue_loss_PctAltered = {
    fit <- glm(loss ~ pct_altered+cluster, family = binomial)
    # overall effect of cluster (not individual coefficients)
    anova(fit, test = "Chisq")["cluster", "Pr(>Chi)"]
  },
  .groups = "drop"
)

pvals_logreg%>%
  write.xlsx(paste0(dir_results_script,"/04_04_gassociation_CNV_tests_gclusters.xlsx"))

p3 = ggplot(df_long, aes(x = cluster, fill = alteration)) +
  geom_bar(position = "fill") +   # proportions per cluster
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("Loss" = "red", "Neutral" = "gray80", "Gain" = "blue")) +
  facet_wrap(~ gene, scales = "free_y") +
  theme_bw() +
  labs(y = "% of samples", fill = "Alteration") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(dir_results_script,"/04_04_CNV_per_gcluster_barplots.pdf"),p3,width=10,height=10)