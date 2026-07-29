#-------------------------------------------------------------------------------
# Aim: Plot precision-recall plot for T1-T2 in D-SQUAME validation
# Author: L.Pozza
# Input: SCCore-GEP threshold-based metrics in T2a
# Output: Supplementary Figure 14B: T1-T2 precision recall plot
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
literature_based_prev_ajcc8 <- data.frame(AJCC8 = c("T4", "T4", "T3", "T3", "T2", "T2", "T1", "T1"),
                                          Stratum = factor(c("T2", "T1", "T2", "T1", "T2", "T1", "T2", "T1"), levels = c("T1", "T2")),
                                          prev_min = c(15.9, 15.9, 8.6, 8.6, 5.5, 5.5, 0.5, 0.5),
                                          prev_max = c(81.1, 81.1, 38.5, 38.5, 17, 17, 4, 4))
#-------------------------------------------------------------------------------

# Subset data, risk by stages and regions of interest
#-------------------------------------------------------------------------------
t1t2_res <- res %>%
       filter(stratification_system == "AJCC_staging" & stratum %in% c("T1", "T2")) %>%
       mutate(Threshold = threshold, 
              Stratum = stratum,
              `Post-test risk` = 100 * PPV_w,
              Sensitivity = 100 * SE_w,
              `% patients high-risk` = 100 * high_risk_w,
              `% patients low-risk` = 100 - `% patients high-risk`) %>%
       select(Stratum, Threshold, `Post-test risk`, Sensitivity, `% patients high-risk`, `% patients low-risk`)
stages_risk_ajcc8 <- data.frame(Stratum = c("T1", "T2"),
                                prevalence = c(0.83, 1.73))
regions_of_interest <- t1t2_res %>%
       filter(Stratum == "T1" & `Post-test risk` > (literature_based_prev_ajcc8 %>% filter(AJCC8 == "T2") %>% pull(prev_min))) %>%
       head(1) %>%
       select(Stratum, Sensitivity) %>%
       rbind(t1t2_res %>%
              filter(Stratum == "T2" & `Post-test risk` > (literature_based_prev_ajcc8 %>% filter(AJCC8 == "T3") %>% pull(prev_min))) %>%
              head(1) %>%
              select(Stratum, Sensitivity))
#-------------------------------------------------------------------------------

# Make and save plot
#-------------------------------------------------------------------------------
# T1
ppv_vs_se_t1 <- ggplot(data = t1t2_res %>% filter(Stratum == "T1") %>% arrange(Threshold),
                          aes(x = Sensitivity,
                              y = `Post-test risk`)) +
                        geom_rect(data = literature_based_prev_ajcc8 %>% filter(Stratum == "T1") %>% droplevels(),
                                  aes(ymin = prev_min, ymax = prev_max, fill = AJCC8, xmin = -Inf, xmax = Inf), inherit.aes = F, alpha = .2)+
                        geom_vline(data = regions_of_interest %>% filter(Stratum == "T1"), aes(xintercept = Sensitivity), linetype = 2, linewidth = 0.4)+
                        scale_fill_manual(values = cols_clinical_data$AJCC8)+
                        geom_text(data =literature_based_prev_ajcc8 %>% filter(Stratum == "T1") %>% droplevels() %>% mutate(mean_prev = (prev_min+prev_max)/2) ,
                                  aes(y = mean_prev, x = 5, label = AJCC8, color = AJCC8), size = 3, fontface = "bold") +
                        scale_color_manual(values = cols_clinical_data$AJCC8)+
                        geom_hline(data = stages_risk_ajcc8 %>% filter(Stratum == "T1"), aes(yintercept = prevalence), linetype = 3, linewidth = 0.4)+
                        geom_path()+
                        theme_bw() +
                        theme(axis.text.x = element_blank(),
                              axis.title.x = element_blank(),
                              axis.ticks.x = element_blank(),
                              panel.grid.minor = element_blank(),
                              axis.line = element_line(colour = "black"),
                              legend.position = "none",
                              strip.background = element_rect(fill = "white")) +
                        labs(y = "Post-test risk (%)", fill = "AJCC8")+
                        scale_y_continuous(trans = log2_trans(),
                                           breaks = c(0, 1, 2, 4, 8, 16, 32, 64, 100))+
                        facet_wrap("Stratum", nrow = 1)
perc_hr_vs_se_t1 <- ggplot(data = t1t2_res %>%
                              filter(Stratum == "T1") %>%
                              group_by(`% patients high-risk`, Sensitivity) %>%
                              unique() %>%
                              arrange(`% patients high-risk`),
                          aes(x = Sensitivity,
                              y = `% patients high-risk`)) +
                        geom_ribbon(aes(ymin = 0, ymax = `% patients high-risk`, fill = "High-risk"), alpha = 0.6)+
                        geom_ribbon(aes(ymax = 100, ymin = `% patients high-risk`, fill = "Low-risk"), alpha = 0.6)+
                        geom_vline(data = regions_of_interest %>% filter(Stratum == "T1"), aes(xintercept = Sensitivity), linetype = 2, linewidth = 0.4)+
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
                        facet_wrap("Stratum", nrow = 1)
# T2
ppv_vs_se_t2 <- ggplot(data = t1t2_res %>% filter(Stratum == "T2") %>% arrange(Threshold),
                          aes(x = Sensitivity,
                              y = `Post-test risk`)) +
                        geom_rect(data = literature_based_prev_ajcc8 %>% filter(Stratum == "T2") %>% droplevels(),
                                  aes(ymin = prev_min, ymax = prev_max, fill = AJCC8, xmin = -Inf, xmax = Inf), inherit.aes = F, alpha = .2)+
                        geom_vline(data = regions_of_interest %>% filter(Stratum == "T2"), aes(xintercept = Sensitivity), linetype = 2, linewidth = 0.4)+
                        scale_fill_manual(values = cols_clinical_data$AJCC8)+
                        geom_text(data =literature_based_prev_ajcc8 %>% filter(Stratum == "T2") %>% droplevels() %>% mutate(mean_prev = (prev_min+prev_max)/2) ,
                                  aes(y = mean_prev, x = 5, label = AJCC8, color = AJCC8), size = 3, fontface = "bold") +
                        scale_color_manual(values = cols_clinical_data$AJCC8)+
                        geom_hline(data = stages_risk_ajcc8 %>% filter(Stratum == "T2"), aes(yintercept = prevalence), linetype = 3, linewidth = 0.4)+
                        geom_path()+
                        theme_bw() +
                        theme(axis.text.x = element_blank(),
                              axis.title.x = element_blank(),
                              axis.ticks.x = element_blank(),
                              panel.grid.minor = element_blank(),
                              axis.line = element_line(colour = "black"),
                              legend.position = "none",
                              strip.background = element_rect(fill = "white")) +
                        labs(y = "Post-test risk (%)", fill = "AJCC8")+
                        scale_y_continuous(trans = log2_trans(),
                                           breaks = c(0, 1, 2, 4, 8, 16, 32, 64, 100))+
                        facet_wrap("Stratum", nrow = 1)
perc_hr_vs_se_t2 <- ggplot(data = t1t2_res %>%
                              filter(Stratum == "T2") %>%
                              group_by(`% patients high-risk`, Sensitivity) %>%
                              unique() %>%
                              arrange(`% patients high-risk`),
                          aes(x = Sensitivity,
                              y = `% patients high-risk`)) +
                        geom_ribbon(aes(ymin = 0, ymax = `% patients high-risk`, fill = "High-risk"), alpha = 0.6)+
                        geom_ribbon(aes(ymax = 100, ymin = `% patients high-risk`, fill = "Low-risk"), alpha = 0.6)+
                        geom_vline(data = regions_of_interest %>% filter(Stratum == "T2"), aes(xintercept = Sensitivity), linetype = 2, linewidth = 0.4)+
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
                        facet_wrap("Stratum", nrow = 1)
p <- ppv_vs_se_t1 + perc_hr_vs_se_t1 + ppv_vs_se_t2 + perc_hr_vs_se_t2 + plot_layout(nrow = 4, heights = c(3, 1, 3, 1))
#-------------------------------------------------------------------------------


# Source data that needs to be saved
#-------------------------------------------------------------------------------
source_data_df <- t1t2_res %>% select(-Threshold)
#-------------------------------------------------------------------------------

# Save results
#-------------------------------------------------------------------------------
write.csv(source_data_df, file.path(output_dir, "supp_figure_14b_t1t2_precision_recall_plot_source_data.csv"), row.names = F)
pdf(file.path(output_dir, "supp_figure_14b_t1t2_precision_recall_plot.pdf"), height = 10, width = 4)
print(p)
dev.off()
#-------------------------------------------------------------------------------