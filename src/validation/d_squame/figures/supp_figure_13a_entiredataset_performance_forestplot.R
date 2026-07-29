#-------------------------------------------------------------------------------
# Aim: Plot of SCCore-GEP performance in entire D-SQUAME validation
# Author: L.Pozza
# Input: SCCore-GEP performances
# Output: Supplementary figure 13A: forest plot of SCCore-GEP entire D-SQUAME validation performances
#-------------------------------------------------------------------------------

# Load libraries
#-------------------------------------------------------------------------------
library(conflicted)
library(tidyverse)
library(forestplot)
library(patchwork)
conflicted::conflicts_prefer(dplyr::filter)
#-------------------------------------------------------------------------------

# Filenames and settings
#-------------------------------------------------------------------------------
# File with results
res_fn <- file.path(results_dir, "intermediate", "validation", "d_squame", "SCCoreGEP_BWH_AJCC8_EMCmodel_performances_entiredataset_with_imputation.csv")
# Output folder
output_dir <- file.path(results_dir, "publication")
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}
#-------------------------------------------------------------------------------

# Source functions
#-------------------------------------------------------------------------------
source(file.path(code_dir, "functions", "colors.R"))
#-------------------------------------------------------------------------------

# Read data
#-------------------------------------------------------------------------------
res <- read.csv(res_fn)
#-------------------------------------------------------------------------------

# Prepare data for forest plot
#-------------------------------------------------------------------------------
df <- res %>%
    filter(model %in% c("AJCC8_model", "BWH_model", "CP_risk", "sccore_gep") &
            variable %in% c("All", "BWH_staging", "AJCC_staging", "Type_of_material_bin") &
            !((model == "AJCC8_model" & variable == "BWH_staging") | (model == "BWH_model" & variable == "AJCC_staging"))) %>%
    select(model, variable, group, wcindex) %>%
    mutate(wcindex_pe = as.numeric(sapply(strsplit(wcindex, "[(|)|-]"), "[[", 1)),
           wcindex_95ci_low = as.numeric(sapply(strsplit(wcindex, "[(|)|-]"), "[[", 2)),
           wcindex_95ci_up = as.numeric(sapply(strsplit(wcindex, "[(|)|-]"), "[[", 3)),
           model = recode(model, "CP_risk" = "EMC model", "sccore_gep" = "SCCore-GEP", "AJCC8_model" = "AJCC8", "BWH_model" = "BWH"),
           group = recode(group, "All" = "Entire cohort", "Excision/Reexcision/Resectie/Mohs" = "Excision"),
           variable = recode(variable, "All" = "Entire cohort", "BWH_staging" = "BWH", "AJCC_staging" = "AJCC8", "Type_of_material_bin" = "Sample type"),
           est = round(wcindex_pe, 2),
           lower = round(wcindex_95ci_low, 2),
           upper = round(wcindex_95ci_up, 2),
           wcindex_95ci = paste0(format(round(wcindex_pe, digits=2), nsmall=2), " (", format(round(wcindex_95ci_low, digits=2), nsmall=2), "-", format(round(wcindex_95ci_up, digits=2), nsmall=2), ")"))
fp_df <- df %>%
    mutate(group_model = paste0(variable, ":", group, ":", model),
           model = factor(model, levels = c("SCCore-GEP", "EMC model", "BWH", "AJCC8"))) %>%
    mutate(category = model) %>%
    select(group_model, group, category, est, lower, upper) %>%
    rbind(data.frame(group_model = c("Entire cohort",
                                     "Sample type", "Sample type:Biopsy", "Sample type:Excision",
                                     "BWH", "BWH:T1", "BWH:T2a", "BWH:T2b",
                                     "AJCC8", "AJCC8:T1", "AJCC8:T2", "AJCC8:T3"),
                     group = c("Entire cohort",
                               "Sample type", "Biopsy", "Excision",
                               "BWH", "T1", "T2a", "T2b",
                               "AJCC8", "T1", "T2", "T3"),
                     category = c("header", 
                                  "header", "subgroup", "subgroup",
                                  "header", "subgroup", "subgroup", "subgroup",
                                  "header", "subgroup", "subgroup", "subgroup"),
                     est = rep(NA, 12),
                     lower = rep(NA, 12),
                     upper = rep(NA, 12))) %>%
    mutate(group_model = factor(group_model,
                                levels = c("AJCC8:T3:EMC model", "AJCC8:T3:SCCore-GEP",
                                           "AJCC8:T3",
                                           "AJCC8:T2:EMC model", "AJCC8:T2:SCCore-GEP",
                                           "AJCC8:T2",
                                           "AJCC8:T1:EMC model", "AJCC8:T1:SCCore-GEP",
                                           "AJCC8:T1",
                                           "AJCC8",
                                           "BWH:T2b:EMC model", "BWH:T2b:SCCore-GEP",
                                           "BWH:T2b",
                                           "BWH:T2a:EMC model", "BWH:T2a:SCCore-GEP",
                                           "BWH:T2a",
                                           "BWH:T1:EMC model", "BWH:T1:SCCore-GEP",
                                           "BWH:T1",
                                           "BWH",
                                           "Sample type:Excision:AJCC8", "Sample type:Excision:BWH",
                                           "Sample type:Excision:EMC model", "Sample type:Excision:SCCore-GEP", 
                                           "Sample type:Excision",
                                           "Sample type:Biopsy:AJCC8", "Sample type:Biopsy:BWH",
                                           "Sample type:Biopsy:EMC model", "Sample type:Biopsy:SCCore-GEP",
                                           "Sample type:Biopsy",
                                           "Sample type",
                                           "Entire cohort:Entire cohort:AJCC8", "Entire cohort:Entire cohort:BWH",
                                           "Entire cohort:Entire cohort:EMC model", "Entire cohort:Entire cohort:SCCore-GEP",
                                           "Entire cohort")),
           group = ifelse(category == "header", as.character(group),
                          ifelse(category == "subgroup", paste0("   ", as.character(group)),
                          ifelse(grepl("Entire cohort", group_model), paste0("   ", as.character(category)), paste0("      ", as.character(category)))))) %>%
    rename(unique_label = group_model,
           label = group) %>%
    arrange(unique_label)
# Split dataset in 2 to make 2 separate forestplot
fp_df1 <- fp_df %>% filter(grepl("Entire cohort|Sample type", unique_label))
fp_df2 <- fp_df %>% filter(!grepl("Entire cohort|Sample type", unique_label))
#-------------------------------------------------------------------------------

# Forest plot #1
#-------------------------------------------------------------------------------
# Identify header rows for separator lines
header_positions1 <- which(unique(fp_df1$unique_label) %in% fp_df1$unique_label[fp_df1$category %in% c("header")])
subgroup_positions1 <- which(unique(fp_df1$unique_label) %in% fp_df1$unique_label[fp_df1$category %in% c("subgroup")])
# Left panel (Subset + weighted C-index with CI)
p_table1 <- ggplot(fp_df1, aes(y = unique_label)) +
  geom_hline(yintercept = header_positions1 - 0.4, color = "grey40", linewidth = 0.7) +
  geom_hline(yintercept = subgroup_positions1 - 0.4, color = "grey80", linewidth = 0.5) +
  geom_text(aes(x = 0, label = label, fontface = ifelse(category %in% c("header", "subgroup"), "bold", "plain")), hjust = 0, size = 5) +
  geom_text(aes(x = 12.6, label = ifelse(!grepl("header|subgroup", category), paste0(est, " (", lower, "-", upper, ")"),"")), hjust = 0.5, size = 5) +
  scale_x_continuous(limits = c(0, 18)) +
  theme_void() +
  theme(plot.margin = margin(r = 10))
# Right panel (Forest plot)
p_forest1 <- ggplot(fp_df1, aes(y = unique_label)) +
  geom_hline(yintercept = header_positions1 - 0.4, color = "grey40", linewidth = 0.7) +
  geom_hline(yintercept = subgroup_positions1 - 0.4, color = "grey80", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = category), linewidth = 0.8, na.rm = TRUE) +
  geom_point(aes(x = est, color = category), size = 5, shape = 15, na.rm = TRUE) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = cols_strat_systems) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(x = "Weighted C-index (95% CI)", y = NULL, color = "") +
  theme_classic()+
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        plot.margin = margin(l = 5),
        legend.position = "bottom",
        text = element_text(size = 15)
  )
#-------------------------------------------------------------------------------

# Forest plot #2
#-------------------------------------------------------------------------------
# Identify header rows for separator lines
header_positions2 <- which(unique(fp_df2$unique_label) %in% fp_df2$unique_label[fp_df2$category %in% c("header")])
subgroup_positions2 <- which(unique(fp_df2$unique_label) %in% fp_df2$unique_label[fp_df2$category %in% c("subgroup")])
# Left panel (Subset + weighted C-index with CI)
p_table2 <- ggplot(fp_df2, aes(y = unique_label)) +
  geom_hline(yintercept = header_positions2 - 0.4, color = "grey40", linewidth = 0.7) +
  geom_hline(yintercept = subgroup_positions2 - 0.4, color = "grey80", linewidth = 0.5) +
  geom_text(aes(x = 0, label = label, fontface = ifelse(category %in% c("header", "subgroup"), "bold", "plain")), hjust = 0, size = 5) +
  geom_text(aes(x = 12.6, label = ifelse(!grepl("header|subgroup", category), paste0(est, " (", lower, "-", upper, ")"),"")), hjust = 0.5, size = 5) +
  scale_x_continuous(limits = c(0, 18)) +
  theme_void() +
  theme(plot.margin = margin(r = 10))
# Right panel (Forest plot)
p_forest2 <- ggplot(fp_df2, aes(y = unique_label)) +
  geom_hline(yintercept = header_positions2 - 0.4, color = "grey40", linewidth = 0.7) +
  geom_hline(yintercept = subgroup_positions2 - 0.4, color = "grey80", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = category), linewidth = 0.8, na.rm = TRUE) +
  geom_point(aes(x = est, color = category), size = 5, shape = 15, na.rm = TRUE) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = cols_strat_systems) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(x = "Weighted C-index (95% CI)", y = NULL, color = "") +
  theme_classic()+
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        plot.margin = margin(l = 5),
        legend.position = "bottom",
        text = element_text(size = 15)
  )
# Combine panels with patchwork
fp <- p_table1 + p_forest1 + p_table2 + p_forest2 +
  plot_layout(widths = c(2, 2.5, 2, 2.5))
#-------------------------------------------------------------------------------

# Source data that needs to be saved
#-------------------------------------------------------------------------------
source_data_df <- df %>%
  select(variable, group, model, est, lower, upper) %>%
  rename(`Stratification variable` = variable,
         Group = group,
         Model = model,
         Weighted_C_index = est,
         CI_95_lower = lower,
         CI_95_upper = upper) %>% 
  mutate(`Stratification variable` = factor(`Stratification variable`,
                                            levels = c("Entire cohort", "Sample type", "Immunosuppression status", "BWH", "AJCC8")),
         Model = factor(Model, levels = c("SCCore-GEP", "EMC model", "BWH", "AJCC8"))) %>%
  arrange(`Stratification variable`, Group, Model)
#-------------------------------------------------------------------------------

# Save results
#-------------------------------------------------------------------------------
write.csv(source_data_df, file.path(output_dir, "supp_figure_13a_entiredataset_performance_forestplot_source_data.csv"), row.names = F)
pdf(file.path(output_dir, "supp_figure_13a_entiredataset_performance_forestplot.pdf"), width = 17, height = 7.5)
print(fp)
dev.off()
#-------------------------------------------------------------------------------