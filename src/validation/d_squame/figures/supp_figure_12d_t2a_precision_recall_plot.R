#-------------------------------------------------------------------------------
# Aim: Plot precision-recall plot for T2a in D-SQUAME validation
# Author: L.Pozza
# Input: SCCore-GEP threshold-based metrics in T2a
# Output: Supplementary Figure 12D: T2a precision recall plot
#-------------------------------------------------------------------------------

# Load libraries
#-------------------------------------------------------------------------------
library(conflicted)
library(tidyverse)
library(patchwork)
library(scales)
conflicted::conflicts_prefer(dplyr::filter)
#-------------------------------------------------------------------------------

# Filenames and settings
#-------------------------------------------------------------------------------
# Path to threshold-based metrics
res_fn <- file.path(results_dir, "intermediate", "validation", "d_squame", "threshold_based_metrics.csv")
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

# Prevalences from literature (Zakhem GA, Qiblawi S, Shelton E, Xu YG. Prevalence of poor outcomes
# in cutaneous squamous cell carcinoma by AJCC and BWH tumor stages: A systematic review and meta-analysis.
# J Am Acad Dermatol. 2025 May;92(5):1064-1071. doi: 10.1016/j.jaad.2024.11.082. Epub 2025 Jan 29. PMID: 39889854.​)
#-------------------------------------------------------------------------------
literature_based_prev_bwh <- data.frame(BWH = c("T3", "T2b", "T2a", "T1"),
                                        Stratum = rep("T2a", 4),
                                        prev_min = c(22.4, 12.6, 1.7, 0.6),
                                        prev_max = c(53.6, 19.5, 8.9, 1.9))
#-------------------------------------------------------------------------------

# Subset data, risk by stages and regions of interest
#-------------------------------------------------------------------------------
t2a_res <- res %>%
       filter(stratification_system == "BWH_staging" & stratum == "T2a") %>%
       mutate(Threshold = threshold, 
              Stratum = stratum,
              `Post-test risk` = 100 * PPV_w,
              Sensitivity = 100 * SE_w,
              `% patients high-risk` = 100 * high_risk_w,
              `% patients low-risk` = 100 - `% patients high-risk`) %>%
       select(Stratum, Threshold, `Post-test risk`, Sensitivity, `% patients high-risk`, `% patients low-risk`)
stages_risk_bwh <- data.frame(Stratum = c("T2a"),
                              risk = c(2.63),
                              prevalence = c(2.62))
regions_of_interest <- t2a_res %>%
       filter(`Post-test risk` > 8.9) %>%
       head(1) %>%
       select(Stratum, Sensitivity)
#-------------------------------------------------------------------------------

# Make and save plot
#-------------------------------------------------------------------------------
ppv_vs_se <- ggplot(data = t2a_res %>% arrange(Threshold),
                          aes(x = Sensitivity,
                              y = `Post-test risk`)) +
                        geom_rect(data = literature_based_prev_bwh,
                                  aes(ymin = prev_min, ymax = prev_max, fill = BWH, xmin = -Inf, xmax = Inf), inherit.aes = F, alpha = .2)+
                        geom_vline(data = regions_of_interest, aes(xintercept = Sensitivity), linetype = 2, linewidth = 0.4)+
                        scale_fill_manual(values = cols_clinical_data$BWH)+
                        geom_text(data =literature_based_prev_bwh %>% mutate(mean_prev = (prev_min+prev_max)/2) ,
                                  aes(y = mean_prev, x = 5, label = BWH, color = BWH), size = 3, fontface = "bold") +
                        scale_color_manual(values = cols_clinical_data$BWH)+
                        geom_hline(data = stages_risk_bwh, aes(yintercept = prevalence), linetype = 3, linewidth = 0.4)+
                        geom_path()+
                        theme_bw() +
                        theme(axis.text.x = element_blank(),
                              axis.title.x = element_blank(),
                              axis.ticks.x = element_blank(),
                              panel.grid.minor = element_blank(),
                              axis.line = element_line(colour = "black"),
                              legend.position = "none",
                              strip.background = element_rect(fill = "white")) +
                        labs(y = "Post-test risk (%)", fill = "BWH")+
                        scale_y_continuous(trans = log2_trans(),
                                           breaks = c(0, 1, 2, 4, 8, 16, 32, 64, 100))+
                        facet_wrap("Stratum", nrow = 1, scales = "free")
perc_hr_vs_se <- ggplot(data = t2a_res %>%
                              group_by(`% patients high-risk`, Sensitivity) %>%
                              unique() %>%
                              arrange(`% patients high-risk`),
                          aes(x = Sensitivity,
                              y = `% patients high-risk`)) +
                        geom_ribbon(aes(ymin = 0, ymax = `% patients high-risk`, fill = "High-risk"), alpha = 0.6)+
                        geom_ribbon(aes(ymax = 100, ymin = `% patients high-risk`, fill = "Low-risk"), alpha = 0.6)+
                        geom_vline(data = regions_of_interest, aes(xintercept = Sensitivity), linetype = 2, linewidth = 0.4)+
                        scale_fill_manual(values = c("High-risk" = "#CC79A7", "Low-risk" = "#009E73"))+
                        theme_bw() +
                        theme(strip.text = element_blank(),
                              panel.grid.minor = element_blank(),
                              axis.line = element_line(colour = "black"),
                              legend.position = "bottom",
                              strip.background = element_rect(fill = "white")) +
                        labs(x = "% Metastases captured",
                             y = "% patients",
                             fill = "")+
                        facet_wrap("Stratum", nrow = 1, scales = "free")
p <- ppv_vs_se + perc_hr_vs_se + plot_layout(nrow = 2, heights = c(3, 1))
#-------------------------------------------------------------------------------


# Source data that needs to be saved
#-------------------------------------------------------------------------------
source_data_df <- t2a_res %>% select(-Threshold)
#-------------------------------------------------------------------------------

# Save results
#-------------------------------------------------------------------------------
write.csv(source_data_df, file.path(output_dir, "supp_figure_12d_t2a_precision_recall_plot_source_data.csv"), row.names = F)
pdf(file.path(output_dir, "supp_figure_12d_t2a_precision_recall_plot.pdf"), height = 5, width = 3.5)
print(p)
dev.off()
#-------------------------------------------------------------------------------