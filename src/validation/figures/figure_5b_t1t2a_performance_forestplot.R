#-------------------------------------------------------------------------------
# Aim: Plot of SCCore-GEP performance in T1-T2a subset of D-SQUAME validation and Nassir et al. datasets
# Author: L.Pozza
# Input: SCCore-GEP performances
# Output: Figure 5B: forest plot of SCCore-GEP performances in T1-T2a D-SQUAME validation and Nassir et al. datasets
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
res_dsquame_fn <- file.path(results_dir, "intermediate", "validation", "d_squame", "SCCoreGEP_BWH_AJCC8_EMCmodel_performances_t1t2a_with_imputation.csv")

res_nassir_fn <- file.path(results_dir, "intermediate", "validation", "nassir", "SCCoreGEP_BWH_performances_max_samples.csv")
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
res_dsquame <- read.csv(res_dsquame_fn)
res_nassir <- read.csv(res_nassir_fn)
#-------------------------------------------------------------------------------

# Prepare data for forest plot
#-------------------------------------------------------------------------------
# D-SQUAME
df_dsquame <- res_dsquame %>%
    filter(model %in% c("AJCC8_model", "BWH_model", "sccore_gep") & variable == "All") %>%
    select(model, variable, group, wcindex) %>%
    mutate(wcindex_pe = as.numeric(sapply(strsplit(wcindex, "[(|)|-]"), "[[", 1)),
           wcindex_95ci_low = as.numeric(sapply(strsplit(wcindex, "[(|)|-]"), "[[", 2)),
           wcindex_95ci_up = as.numeric(sapply(strsplit(wcindex, "[(|)|-]"), "[[", 3)),
           model = recode(model, "sccore_gep" = "SCCore-GEP", "BWH_model" = "BWH", "AJCC8_model" = "AJCC8"),
           group = recode(group, "All" = "D-SQUAME"),
           variable = recode(variable, "All" = "D-SQUAME"),
           est = round(wcindex_pe, 2),
           lower = round(wcindex_95ci_low, 2),
           upper = round(wcindex_95ci_up, 2),
           wcindex_95ci = paste0(format(round(wcindex_pe, digits=2), nsmall=2), " (", format(round(wcindex_95ci_low, digits=2), nsmall=2), "-", format(round(wcindex_95ci_up, digits=2), nsmall=2), ")"))

fp_df_dsquame <- df_dsquame %>%
    mutate(group_model = paste0(group, ":", model),
           model = factor(model, levels = c("SCCore-GEP", "BWH", "AJCC8"))) %>%
    rename(category = model) %>%
    select(group_model, category, est, lower, upper) %>%
    rbind(data.frame(group_model = c("D-SQUAME"),
                     category = c("header"),
                     est = NA,
                     lower = NA,
                     upper = NA)) %>%
    mutate(group_model = factor(group_model,
                                levels = c("D-SQUAME:AJCC8",
                                           "D-SQUAME:BWH",
                                           "D-SQUAME:SCCore-GEP",
                                           "D-SQUAME")),
           group = ifelse(category == "header", as.character(group_model),
                          ifelse(category == "subgroup", paste0("   ", as.character(group_model)), paste0("      ", as.character(category))))) %>%
    rename(unique_label = group_model,
           label = group) %>%
    arrange(unique_label)
# Nassir et al.
df_nassir <- res_nassir %>%
    filter(model %in% c("bwh", "sccore_gep") & variable == "BWH_staging_T1T2a_T2bT3" & group == "T1-T2a") %>%
    mutate(group = "Nassir et al.",
           auc_pe = as.numeric(sapply(strsplit(auc, "[(|)|-]"), "[[", 1)),
           auc_95ci_low = as.numeric(sapply(strsplit(auc, "[(|)|-]"), "[[", 2)),
           auc_95ci_up = as.numeric(sapply(strsplit(auc, "[(|)|-]"), "[[", 3)),
           model = recode(model, "sccore_gep" = "SCCore-GEP", "bwh" = "BWH"),
           est = round(auc_pe, 2),
           lower = round(auc_95ci_low, 2),
           upper = round(auc_95ci_up, 2),
           auc_95ci = paste0(format(round(auc_pe, digits=2), nsmall=2), " (", format(round(auc_95ci_low, digits=2), nsmall=2), "-", format(round(auc_95ci_up, digits=2), nsmall=2), ")")) %>%
    select(-variable)
fp_df_nassir <- df_nassir %>%
    mutate(group_model = paste0(group, ":", model),
           category = factor(model, levels = c("SCCore-GEP", "BWH", "AJCC8"))) %>%
    select(group_model, group, category, est, lower, upper) %>%
    rbind(data.frame(group_model = c("Nassir et al."),
                     group = c("Nassir et al."),
                     category = c("header"),
                     est = NA,
                     lower = NA,
                     upper = NA)) %>%
    mutate(group_model = factor(group_model,
                                levels = c("Nassir et al.:BWH", "Nassir et al.:SCCore-GEP", "Nassir et al.")),
           group = ifelse(category == "header", as.character(group),
                          ifelse(category == "subgroup", paste0("   ", as.character(group)), paste0("      ", as.character(category))))) %>%
    droplevels() %>%
    rename(unique_label = group_model,
           label = group) %>%
    arrange(unique_label)
# Combined
fp_df <- fp_df_dsquame %>%
    rbind(fp_df_nassir) %>%
    mutate(unique_label = factor(unique_label, levels = c("Nassir et al.:BWH", "Nassir et al.:SCCore-GEP", "Nassir et al.", 
                                                          "D-SQUAME:AJCC8", "D-SQUAME:BWH", "D-SQUAME:SCCore-GEP", "D-SQUAME"))) %>%
    arrange(unique_label)
#-------------------------------------------------------------------------------

# Forest plot
#-------------------------------------------------------------------------------
# Identify header rows for separator lines
header_positions <- which(fp_df$label %in% fp_df$label[fp_df$category %in% c("header")])
# Left panel (Subset + weighted C-index with CI)
p_table <- ggplot(fp_df, aes(y = unique_label)) +
  geom_hline(yintercept = header_positions - 0.4, color = "grey40", linewidth = 0.7) +
  geom_text(aes(x = 0, label = label, fontface = ifelse(category %in% c("header", "subgroup"), "bold", "plain")), hjust = 0, size = 5) +
  geom_text(aes(x = 12.6, label = ifelse(!grepl("header|subgroup", category), paste0(est, " (", lower, "-", upper, ")"),"")), hjust = 0.5, size = 5) +
  scale_x_continuous(limits = c(0, 18)) +
  theme_void() +
  theme(plot.margin = margin(r = 10))
# Right panel (Forest plot)
p_forest <- ggplot(fp_df, aes(y = unique_label)) +
  geom_hline(yintercept = header_positions - 0.4, color = "grey40", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = category), linewidth = 0.8, na.rm = TRUE) +
  geom_point(aes(x = est, color = category), size = 5, shape = 15, na.rm = TRUE) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = cols_strat_systems) +
  scale_x_continuous(limits = c(0.3, 1), breaks = seq(0.3, 1, 0.1)) +
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
fp <- p_table + p_forest +
  plot_layout(widths = c(2, 3))
#-------------------------------------------------------------------------------

# Source data that needs to be saved
#-------------------------------------------------------------------------------
source_data_df <- fp_df %>%
  mutate(Dataset = ifelse(grepl("D-SQUAME", unique_label), "D-SQUAME", "Nassir et al.")) %>%
  filter(category != "header") %>%
  select(Dataset, category, est, lower, upper) %>%
  rename(Model = category,
         Weighted_C_index = est,
         CI_95_lower = lower,
         CI_95_upper = upper) %>% 
  mutate(Dataset = factor(Dataset, levels = c("D-SQUAME", "Nassir et al.")),
         Model = factor(Model, levels = c("SCCore-GEP", "BWH", "AJCC8"))) %>%
  arrange(Dataset, Model)
#-------------------------------------------------------------------------------

# Save results
#-------------------------------------------------------------------------------
write.csv(source_data_df, file.path(output_dir, "figure_5b_t1t2a_performance_forestplot_source_data.csv"), row.names = F)
pdf(file.path(output_dir, "figure_5b_t1t2a_performance_forestplot.pdf"),
    width = 10, height = 3)
print(fp)
dev.off()
#-------------------------------------------------------------------------------