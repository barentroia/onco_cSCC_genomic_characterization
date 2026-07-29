#-------------------------------------------------------------------------------
# Aim: Forest plot of multivariable analysis in T1-T2a subset of D-SQUAME validation
# Author: L.Pozza
# Input: Multivariable analysis results of SCCore-GEP vs BWH, AJCC8 and EMC model
# Output: Supplememntary Figure 12C: forest plot of multivariable analysis in T1-T2a subset of D-SQUAME validation
#-------------------------------------------------------------------------------

# Load libraries
#-------------------------------------------------------------------------------
library(conflicted)
library(openxlsx)
library(tidyverse)
library(forestplot)
library(patchwork)
conflicted::conflicts_prefer(dplyr::filter)
#-------------------------------------------------------------------------------

# Filenames and settings
#-------------------------------------------------------------------------------
# Multivaraite analysis results
sccoregep_bwh_fn <- file.path(results_dir, "intermediate", "validation", "d_squame", "multivariable_analysis_t1t2a_SCCoreGEP_vs_BWH.csv")
sccoregep_ajcc8_fn <- file.path(results_dir, "intermediate", "validation", "d_squame", "multivariable_analysis_t1t2a_SCCoreGEP_vs_AJCC8.csv")
sccoregep_emcmodel_fn <- file.path(results_dir, "intermediate", "validation", "d_squame", "multivariable_analysis_t1t2a_SCCoreGEP_vs_EMCmodel.csv")
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
sccoregep_bwh <- read.csv(sccoregep_bwh_fn)
sccoregep_ajcc8 <- read.csv(sccoregep_ajcc8_fn)
sccoregep_emcmodel <- read.csv(sccoregep_emcmodel_fn)
#-------------------------------------------------------------------------------

# Pre-process data
#-------------------------------------------------------------------------------
# Combine multivariable analysis results
df <- sccoregep_bwh %>%
    mutate(analysis = "SCCore-GEP vs BWH") %>%
    rbind(sccoregep_ajcc8 %>%
            mutate(analysis = "SCCore-GEP vs AJCC8")) %>%
    rbind(sccoregep_emcmodel %>%
            mutate(analysis = "SCCore-GEP vs EMC model")) %>%
    mutate(mean = round(hr, 2),
           lower = round(hr_ci95_low, 2),
           upper = round(hr_ci95_up, 2),
           pvalue = ifelse(pvalue < 0.01, pixiedust::pvalString(pvalue, format = "scientific", digits = 3), round(pvalue, 2)),
           hr_ci95 = paste0(format(round(hr, digits = 2), nsmall = 2), " (", format(round(hr_ci95_low, digits = 2), nsmall=2), "-", format(round(hr_ci95_up, digits=2), nsmall=2), ")"),
           model = ifelse(grepl("BWH", term), "BWH", ifelse(grepl("AJCC8", term), "AJCC8", ifelse(grepl("SCCore-GEP", term), "SCCore-GEP", "EMC model"))),
           term = sub("^.*:", "", term))
fp_df <- df %>%
    mutate(analysis_term = paste0(analysis, ":", term),
           model = factor(model, levels = c("SCCore-GEP", "EMC model", "BWH", "AJCC8")),
           category = model) %>%
    select(analysis_term, term, category, mean, lower, upper, hr_ci95, pvalue) %>%
    rbind(data.frame(analysis_term = c("SCCore-GEP vs EMC model", "SCCore-GEP vs BWH", "SCCore-GEP vs AJCC8"),
                     term = rep(NA, 3),
                     category = c("header", "header", "header"),
                     mean = rep(NA, 3),
                     lower = rep(NA, 3),
                     upper = rep(NA, 3),
                     hr_ci95 = rep(NA, 3),
                     pvalue = rep(NA, 3))) %>%
    mutate(analysis_term = factor(analysis_term,
                                  levels = c("SCCore-GEP vs AJCC8:T2-T3", "SCCore-GEP vs AJCC8:SCCore-GEP", "SCCore-GEP vs AJCC8",
                                             "SCCore-GEP vs BWH:T2a", "SCCore-GEP vs BWH:SCCore-GEP", "SCCore-GEP vs BWH",
                                             "SCCore-GEP vs EMC model:EMC model", "SCCore-GEP vs EMC model:SCCore-GEP", "SCCore-GEP vs EMC model")),
           group = ifelse(category == "header", as.character(term),
                          ifelse(category == "subgroup", paste0("   ", as.character(group)), paste0("      ", as.character(category))))) %>%
    arrange(analysis_term)
#-------------------------------------------------------------------------------

# Forest plot
#-------------------------------------------------------------------------------
# Identify header rows for separator lines
header_positions <- which(unique(fp_df$analysis_term) %in% fp_df$analysis_term[fp_df$category %in% c("header")])
fp_right <- fp_df %>%
    ggplot(aes(y = analysis_term))+
    geom_hline(yintercept = header_positions - 0.45, color = "grey40", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = lower, xmax = upper, color = category), linewidth = 0.7)+
    geom_point(aes(x = mean, color = category), size = 3, shape = 15, na.rm = TRUE) +
    geom_vline(xintercept = 1, linetype = 2) +
    labs(x = "Hazard ratio (95% CI)", color = "", shape = "", linetype = "")+    
    scale_x_continuous(breaks = seq(0, 11, 1), limits = c(0, 11))+
    scale_fill_manual(values = c("white", "#C8C8C8"))+
    scale_color_manual(values = cols_strat_systems)+
    theme_classic()+
    theme(axis.line.y = element_blank(),
          axis.ticks.y= element_blank(),
          axis.text.y= element_blank(),
          axis.title.y= element_blank(),
          plot.margin = margin(l = 0))+
    guides(fill = "none")
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

# Source data that needs to be saved
#-------------------------------------------------------------------------------
source_data_df <- df %>%
  select(analysis, term, mean, lower, upper, pvalue) %>%
  rename(Analysis = analysis,
         Term = term,
         HR = mean,
         CI_95_lower = lower,
         CI_95_upper = upper,
         p_value = pvalue) %>%
  mutate(Analysis = factor(Analysis, levels = c("SCCore-GEP vs EMC model", "SCCore-GEP vs BWH", "SCCore-GEP vs AJCC8"))) %>%
  arrange(Analysis)
#-------------------------------------------------------------------------------

# Save results
#-------------------------------------------------------------------------------
write.csv(source_data_df, file.path(output_dir, "supp_figure_12c_t1t2a_multivariable_analysis_ext_source_data.csv"), row.names = F)
pdf(file.path(output_dir, "supp_figure_12c_t1t2a_multivariable_analysis_ext.pdf"), height = 2.5, width = 16)
print(fp)
dev.off()
#-------------------------------------------------------------------------------