#===============================================================================
# Aim:      Generate Supplementary Figure 13C: forest plot of the multivariable
#           analyses of SCCore-GEP combined with BWH, AJCC8, the EMC model or
#           the DvP score, in the BWH T1-T2a subset of D-SQUAME validation dataset.
# Author:   Adapted from code of L. Pozza
# Input:    - Multivariable analysis results
#             (multivariable_analysis_t1t2a_SCCoreGEP_vs_<model>.csv)
# Output:   - supp_figure_13C_t1t2a_multivariable_analysis_ext.pdf and _source_data.csv in results_dir/publication
# Requires: code_dir, results_dir (main.R)
# Usage:    Sourced from main.R
#===============================================================================

# 1. Load libraries
#-------------------------------------------------------------------------------
library(conflicted)
library(openxlsx)
library(tidyverse)
library(forestplot)
library(patchwork)
conflicted::conflicts_prefer(dplyr::filter)
#-------------------------------------------------------------------------------

# 2. Filenames and settings
#-------------------------------------------------------------------------------
# Multivariable analysis results
mv_dir <- file.path(results_dir, "intermediate", "validation", "d_squame")
sccoregep_bwh_fn      <- file.path(mv_dir, "multivariable_analysis_t1t2a_SCCoreGEP_vs_BWH.csv")
sccoregep_ajcc8_fn    <- file.path(mv_dir, "multivariable_analysis_t1t2a_SCCoreGEP_vs_AJCC8.csv")
sccoregep_emcmodel_fn <- file.path(mv_dir, "multivariable_analysis_t1t2a_SCCoreGEP_vs_EMCmodel.csv")
sccoregep_dvp_fn      <- file.path(mv_dir, "multivariable_analysis_t1t2a_SCCoreGEP_vs_DvPscore.csv")
# Output folder
output_dir <- file.path(results_dir, "publication")
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}
#-------------------------------------------------------------------------------

# 3. Source functions
#-------------------------------------------------------------------------------
source(file.path(code_dir, "functions", "colors.R"))
#-------------------------------------------------------------------------------

# 4. Read data
#-------------------------------------------------------------------------------
sccoregep_bwh      <- read.csv(sccoregep_bwh_fn)
sccoregep_ajcc8    <- read.csv(sccoregep_ajcc8_fn)
sccoregep_emcmodel <- read.csv(sccoregep_emcmodel_fn)
sccoregep_dvp      <- read.csv(sccoregep_dvp_fn)
#-------------------------------------------------------------------------------

# 5. Settings for the plot
#-------------------------------------------------------------------------------
# Order of the analyses from top to bottom in the plot
analysis_order <- c("SCCore-GEP vs EMC model", "SCCore-GEP vs BWH",
                    "SCCore-GEP vs AJCC8", "SCCore-GEP vs DvP score")
model_levels <- c("SCCore-GEP", "EMC model", "BWH", "AJCC8", "DvP score")

#-------------------------------------------------------------------------------

# 6. Pre-process data
#-------------------------------------------------------------------------------
# Combine multivariable analysis results
df <- bind_rows(
  sccoregep_bwh      %>% mutate(analysis = "SCCore-GEP vs BWH"),
  sccoregep_ajcc8    %>% mutate(analysis = "SCCore-GEP vs AJCC8"),
  sccoregep_emcmodel %>% mutate(analysis = "SCCore-GEP vs EMC model"),
  sccoregep_dvp      %>% mutate(analysis = "SCCore-GEP vs DvP score")) %>%
  mutate(mean = round(hr, 2),
         lower = round(hr_ci95_low, 2),
         upper = round(hr_ci95_up, 2),
         pvalue = ifelse(pvalue < 0.01, pixiedust::pvalString(pvalue, format = "scientific", digits = 3), round(pvalue, 2)),
         hr_ci95 = paste0(format(round(hr, digits = 2), nsmall = 2), " (", format(round(hr_ci95_low, digits = 2), nsmall=2), "-", format(round(hr_ci95_up, digits=2), nsmall=2), ")"),
         model = case_when(grepl("BWH", term)        ~ "BWH",
                           grepl("AJCC8", term)      ~ "AJCC8",
                           grepl("SCCore-GEP", term) ~ "SCCore-GEP",
                           grepl("DvP", term)        ~ "DvP score",
                           grepl("EMC", term)        ~ "EMC model"),
         term = sub("^.*:", "", term))
stopifnot(!anyNA(df$model))

# Rows per analysis: header, SCCore-GEP, comparator term(s)
fp_df <- df %>%
  mutate(analysis_term = paste0(analysis, ":", term),
         category = as.character(model)) %>%
  select(analysis, analysis_term, term, category, mean, lower, upper, hr_ci95, pvalue) %>%
  bind_rows(data.frame(analysis = analysis_order,
                       analysis_term = analysis_order,
                       category = "header")) %>%
  mutate(analysis = factor(analysis, levels = analysis_order),
         is_gep = term %in% "SCCore-GEP",
         category = factor(category, levels = c(model_levels, "header"))) %>%
  # within each analysis: header first, then SCCore-GEP, then the comparator
  arrange(analysis, category != "header", !is_gep) %>%
  # ggplot draws the first factor level at the bottom -> reverse
  mutate(analysis_term = factor(analysis_term, levels = rev(unique(analysis_term)))) %>%
  arrange(analysis_term) %>%
  select(-is_gep)

#-------------------------------------------------------------------------------

# 7. Forest plot
#-------------------------------------------------------------------------------
# Identify header rows for separator lines
header_positions <- which(levels(fp_df$analysis_term) %in% fp_df$analysis_term[fp_df$category %in% "header"])
fp_right <- fp_df %>%
  ggplot(aes(y = analysis_term))+
  geom_hline(yintercept = header_positions - 0.45, color = "grey40", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = category), linewidth = 0.7, na.rm = TRUE)+
  geom_point(aes(x = mean, color = category), size = 3, shape = 15, na.rm = TRUE) +
  geom_vline(xintercept = 1, linetype = 2) +
  labs(x = "Hazard ratio (95% CI)", color = "", shape = "", linetype = "")+    
  scale_x_continuous(breaks = seq(0, 11, 1), limits = c(0, 11))+
  scale_color_manual(values = cols_strat_systems, breaks = model_levels)+
  theme_classic()+
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        plot.margin = margin(l = 0))
fp_left <- fp_df %>%
  ggplot(aes(y = analysis_term)) + 
  geom_hline(yintercept = header_positions - 0.45, color = "grey40", linewidth = 0.4) +
  geom_text(aes(x = 0, label = ifelse(category == "header", as.character(analysis_term), ""), fontface = ifelse(category == "header", "bold", "plain")), hjust = 0, size = 3) +
  geom_text(aes(x = 3.5, label = ifelse(category != "header", term, "")), hjust = 0, size = 3) +
  geom_text(aes(x = 6.5, label = ifelse(category != "header", paste0(mean, " (", lower, "-", upper, ")"), "")), hjust = 0.5, size = 3) +
  geom_text(aes(x = 8.5, label = ifelse(category != "header", pvalue, "")), hjust = 0.5, size = 3) +
  theme_void() +
  theme(plot.margin = margin(r = 0))+
  coord_cartesian(xlim = c(0, 9))
fp <- fp_left +
  fp_right + 
  plot_layout(widths = c(2.5, 4))
#-------------------------------------------------------------------------------

# 8. Source data that needs to be saved
#-------------------------------------------------------------------------------
source_data_df <- df %>%
  select(analysis, term, mean, lower, upper, pvalue) %>%
  rename(Analysis = analysis,
         Term = term,
         HR = mean,
         CI_95_lower = lower,
         CI_95_upper = upper,
         p_value = pvalue) %>%
  mutate(Analysis = factor(Analysis, levels = analysis_order)) %>%
  arrange(Analysis)
#-------------------------------------------------------------------------------

# 9. Save results
#-------------------------------------------------------------------------------
write.csv(source_data_df, file.path(output_dir, "supp_figure_13C_t1t2a_multivariable_analysis_ext_source_data.csv"), row.names = F)
pdf(file.path(output_dir, "supp_figure_13C_t1t2a_multivariable_analysis_ext.pdf"), height = 3.3, width = 16)
print(fp)
dev.off()
#-------------------------------------------------------------------------------