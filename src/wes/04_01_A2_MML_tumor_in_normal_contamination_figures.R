#---------------------------------------------------
# Aim: Report tumor in normal contamination:
# Author: B. Rentroia Pacheco
# Input: Results from WES
# Output: Statistics on tumor in normal contamination
#---------------------------------------------------

# Check distribution of rescued mutations per sample:
mega_MML_clonal%>%
  group_by(Tumor_Sample_Barcode) %>%
  mutate(
    n_mutations = n()
  )%>%
  ungroup() %>%
  ggplot(
    aes(
      x = fct_reorder(Tumor_Sample_Barcode,
                      n_mutations,
                      .desc = TRUE),
      fill = Rescued_SNP
    )
  ) +
  geom_bar()+xlab("Sample")+theme_bw()

# Identify % of samples with contamination and what is the average contamination:
rescued_summary <- mega_MML_clonal %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(
    n_mutations = n(),
    n_rescued = sum(Rescued_SNP=="Yes", na.rm = TRUE),
    pct_rescued = 100 * n_rescued / n_mutations,
    n_pathogenic_rescued = sum(Rescued_SNP=="Yes"&Pathogenic=="Probably", na.rm = TRUE),
    n_pathogenic = sum(Pathogenic=="Probably", na.rm = TRUE),
    pct_rescued_path = 100 * n_pathogenic_rescued / n_pathogenic,
    .groups = "drop"
  ) %>% as.data.frame()

# Estimate contamination based on normal MAF (all clonal mutations):
contamination_all <- mega_MML_clonal %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(
    mean_MAF_normal = mean(MAF_normal, na.rm = TRUE),
    estimated_contamination_mean = 2 * mean_MAF_normal,
    max_MAF_normal = max(MAF_normal, na.rm = TRUE),
    n_variants = n(),
    .groups = "drop"
  )

# Compare contamination estimate with % rescued mutations per sample:
rescued_contamination <- merge(rescued_summary, contamination_all)

ggplot(rescued_contamination, aes(x = estimated_contamination_mean, y = pct_rescued)) +
  geom_point() +
  xlab("Contamination using mean MAF (%)") +
  ylab("Rescued mutations (%)") +
  geom_vline(xintercept = 0.05, linetype = 2, color = "red")

# Number of samples with estimated contamination (mean-based) > 5%:
sum(rescued_contamination$estimated_contamination_mean > 0.05)/147

# Maximum average contamination (~16%):
max(rescued_contamination$estimated_contamination_mean)
