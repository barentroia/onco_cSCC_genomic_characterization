#===============================================================================
# Aim:      Generate Figure 5D: precision-recall plot of SCCore-GEP in the BWH
#           T1-T2a subset of the D-SQUAME validation dataset, with clinical
#           decision bands (IQRs) from the SCOUT expert survey and the
#           sensitivity at the selected operating point.
# Author:   L. Pozza, adjusted by R. Ruiter
# Input:    - SCCore-GEP threshold-based metrics (threshold_based_metrics.csv)
#           - Source data of the PPV vs sensitivity scatter plot
#             (output of supp_figure_13D_t1t2a_scatterplot.R)
# Output:   - figure_5d_t1t2a_precision_recall_plot.pdf and _source_data.csv
#             in results_dir/publication
# Requires: code_dir, results_dir (main.R); supp_figure_13D run first
# Usage:    Sourced from main.R
#===============================================================================

# 1. Load libraries
#-------------------------------------------------------------------------------
library(conflicted)
library(tidyverse)
library(patchwork)
library(scales)
conflicted::conflicts_prefer(dplyr::filter)
#-------------------------------------------------------------------------------

# 2. Filenames and settings
#-------------------------------------------------------------------------------
# Path to threshold-based metrics
res_fn <- file.path(results_dir, "intermediate", "validation", "d_squame", "threshold_based_metrics.csv")

# Sensitivity in validation at the operating point chosen on discovery (0.362745098039216).
src <- read.csv(file.path(results_dir, "publication", "supp_figure_13D_t1t2a_scatterplot_PPV_vs_SE_scatter_sourcedata.csv"))
se_op1  <- src$SE_w[src$Rule == "OP1"]    # x position of the SCCore-GEP point
se_op1_pct <- 100 * se_op1   # this figure uses %, the scatter source data uses a proportion
op_line <- data.frame(Stratum = "T1-T2a", Sensitivity = se_op1_pct)
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
res <- read.csv(res_fn)
#-------------------------------------------------------------------------------

# 5. Clinical decision thresholds
#-------------------------------------------------------------------------------
# SCOUT consortium survey (n = 53 experts), 5-year metastatic risk (%):
# Rentroia-Pacheco B, Tokez S, Bramer EM, et al. Personalised decision making to
# predict absolute metastatic risk in cutaneous squamous cell carcinoma: development
# and validation of a clinico-pathological model. EClinicalMedicine. 2023;63:102150.
# Supplementary Table S6.
scout <- data.frame(decision = c("Follow-up", "Adjuvant RT", "Adjuvant systemic"),
                    median   = c(10, 20, 30),
                    q1       = c(5, 10, 20),
                    q3       = c(20, 20, 40))

# Standard of care: below Q1 of the follow-up threshold:
# at least 75% of experts would not include the patient in a follow-up schedule
soc_cutoff <- scout$q1[scout$decision == "Follow-up"]   # 5
#-------------------------------------------------------------------------------

# 6. Subset data
#-------------------------------------------------------------------------------
t1t2a_res <- res %>%
  filter(stratification_system == "BWH_staging_T1T2a_T2bT3" & stratum == "T1-T2a") %>%
  mutate(Threshold = threshold,
         Stratum = stratum,
         `Post-test risk` = 100 * PPV_w,
         Sensitivity = 100 * SE_w,
         `% patients high-risk` = 100 * high_risk_w,
         `% patients low-risk` = 100 - `% patients high-risk`) %>%
  select(Stratum, Threshold, `Post-test risk`, Sensitivity,
         `% patients high-risk`, `% patients low-risk`)
#-------------------------------------------------------------------------------

# 7. Stratum prevalence (from the data) and region of interest
#-------------------------------------------------------------------------------
# prevalence = PPV x % high-risk / sensitivity (valid at any threshold)
stages_risk_bwh <- t1t2a_res %>%
  filter(Sensitivity > 0, !is.na(`Post-test risk`)) %>%
  mutate(prevalence = `Post-test risk` * (`% patients high-risk` / 100) / (Sensitivity / 100)) %>%
  slice_max(Sensitivity, n = 1, with_ties = FALSE) %>%
  transmute(Stratum, prevalence)
message("T1-T2a prevalence from data: ", round(stages_risk_bwh$prevalence, 2), "%")

# Vertical line: highest sensitivity at which the PPV is still above 5% (lower limit of follow-up band)
regions_of_interest <- t1t2a_res %>%
  filter(`Post-test risk` > scout$q1[scout$decision == "Follow-up"]) %>%   # 5%
  slice_max(Sensitivity, n = 1, with_ties = FALSE) %>%
  select(Stratum, Sensitivity)
#-------------------------------------------------------------------------------

# 8. Decision bands (IQRs of the SCOUT thresholds)
#-------------------------------------------------------------------------------
# Bottom of the log axis: nearest power of 2 below the lowest plotted value
risk_vals <- c(t1t2a_res$`Post-test risk`, stages_risk_bwh$prevalence, soc_cutoff)
y_min <- 2^floor(log2(min(risk_vals[risk_vals > 0], na.rm = TRUE)))

# Rows are drawn in this order: follow-up first, RT on top (the two IQRs overlap at 10-20%)
decision_bands <- data.frame(
  band     = c("Modulation of follow-up", "Adjuvant radiotherapy",
               "Adjuvant systemic treatment", "Standard of care"),
  prev_min = c(scout$q1[scout$decision == "Follow-up"],           # 5
               scout$q1[scout$decision == "Adjuvant RT"],         # 10
               scout$q1[scout$decision == "Adjuvant systemic"],   # 20
               0),                                                # 0 -> bottom of panel (log axis)
  prev_max = c(scout$q3[scout$decision == "Follow-up"],           # 20
               scout$q3[scout$decision == "Adjuvant RT"],         # 20
               100,                                               # top of panel
               soc_cutoff)                                        # 5
) %>%
  mutate(band = factor(band, levels = band))

cols_bands <- c("Modulation of follow-up"     = unname(cols_clinical_data$BWH["T2a"]),
                "Adjuvant radiotherapy"       = unname(cols_clinical_data$BWH["T3"]),
                "Adjuvant systemic treatment" = "#6E5A93",
                "Standard of care"            = unname(cols_clinical_data$BWH["T1"]))

# Legend order: top to bottom, matching the vertical order of the bands
legend_order <- c("Adjuvant systemic treatment", "Adjuvant radiotherapy",
                  "Modulation of follow-up", "Standard of care")
#-------------------------------------------------------------------------------

# 9. Make and save plot
#-------------------------------------------------------------------------------
ppv_vs_se <- ggplot(data = t1t2a_res %>% arrange(Threshold),
                    aes(x = Sensitivity, y = `Post-test risk`)) +
  geom_rect(data = decision_bands,
            aes(ymin = prev_min, ymax = prev_max, fill = band, xmin = -Inf, xmax = Inf),
            inherit.aes = F, alpha = .3) +
  scale_fill_manual(values = cols_bands, breaks = legend_order,
                    name = "Considered by SCOUT:") +
  geom_hline(data = stages_risk_bwh, aes(yintercept = prevalence), linetype = 3, linewidth = 0.4) +
  geom_vline(data = op_line, aes(xintercept = Sensitivity),
             linetype = 2, linewidth = 0.4, colour = "black") +
  geom_path() +
  coord_cartesian(ylim = c(y_min, 100)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        strip.background = element_rect(fill = "white")) +
  labs(y = "PPV (%)") +
  scale_y_continuous(trans = log2_trans(),
                     breaks = c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 100)) +
  facet_wrap("Stratum", nrow = 1, scales = "free")

perc_hr_vs_se <- ggplot(data = t1t2a_res %>%
                          group_by(`% patients high-risk`, Sensitivity) %>%
                          unique() %>%
                          arrange(`% patients high-risk`),
                        aes(x = Sensitivity,
                            y = `% patients high-risk`)) +
  geom_ribbon(aes(ymin = 0, ymax = `% patients high-risk`,
                  fill = "High-risk"), alpha = 0.6) +
  geom_ribbon(aes(ymax = 100, ymin = `% patients high-risk`,
                  fill = "Low-risk"), alpha = 0.6) +
  geom_vline(data = op_line, aes(xintercept = Sensitivity),
             linetype = 2, linewidth = 0.4, colour = "black") +
  
  scale_fill_manual(values = c("High-risk" = "#CC79A7",
                               "Low-risk"  = "#009E73"),
                    breaks = c("Low-risk", "High-risk"),   # top-to-bottom, matching the stack
                    name = "SCCore-GEP:") +
  theme_bw() +
  theme(strip.text = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right",
        legend.justification = "left",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        strip.background = element_rect(fill = "white")) +
  labs(x = "Sensitivity (%)",
       y = "% patients") +
  facet_wrap("Stratum", nrow = 1, scales = "free")

p <- ppv_vs_se + perc_hr_vs_se + plot_layout(nrow = 2, heights = c(3, 1))
#-------------------------------------------------------------------------------

# 10. Source data that needs to be saved
#-------------------------------------------------------------------------------
source_data_df <- t1t2a_res %>% select(-Threshold) %>%
  mutate(op1_sensitivity = se_op1_pct)
#-------------------------------------------------------------------------------

# 11. Save results
#-------------------------------------------------------------------------------
write.csv(source_data_df, file.path(output_dir, "figure_5d_t1t2a_precision_recall_plot_source_data.csv"), row.names = F)
pdf(file.path(output_dir, "figure_5d_t1t2a_precision_recall_plot.pdf"), height = 4.2, width = 4.7)
print(p)
dev.off()
#-------------------------------------------------------------------------------