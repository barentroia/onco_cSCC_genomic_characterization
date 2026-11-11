#===============================================================================
# Aim:      Generate Supplementary Figure 15B: precision-recall plot of
#           SCCore-GEP in the AJCC8 T1-T2 subset of the D-SQUAME validation
#           dataset, with clinical decision bands (IQRs) from the SCOUT expert
#           survey and the operating point (from discovery BWH T1-T2a) shown
#           as a dashed vertical line.
# Author:   L. Pozza, adjusted by R. Ruiter
# Input:    - SCCore-GEP threshold-based metrics in validation
#             (intermediate/validation/d_squame/threshold_based_metrics.csv)
#           - Discovery T1-T2a threshold metrics, for the operating point
#             (intermediate/sccore_gep/discovery_T1_T2a_threshold_metrics_discoveryCal_OOB.csv)
#           - Imputed validation clinical data (mids object) and merged
#             imputed validation clinical data (intermediate/clinical/)
#           - Validation NCC weights (e7_weights)
#           - Recalibrated SCCore-GEP validation predictions
#             (intermediate/validation/d_squame/recalibrated_SCCoreGEP_predictions.csv)
#           - Colour definitions (src/functions/colors.R)
# Output:   - publication/supp_figure_15B_t1t2_precision_recall_plot.pdf
#           - publication/supp_figure_15B_t1t2_precision_recall_plot_source_data.csv
# Requires: code_dir, results_dir, dir_results_intermediate_sccore_gep (main.R);
#           e7_weights (config.R)
# Usage:    Sourced from main.R
# Notes:    - The weighted sensitivity at the operating point is pooled across
#             imputations in the same way as Supplementary Figure 13D.
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
res_fn <- file.path(results_dir ,"intermediate", "validation", "d_squame", "threshold_based_metrics.csv")
# Output folder
output_dir <- file.path(results_dir, "publication")
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}

# Stratification system and stratum to plot
strat_system_sel <- "AJCC_staging_T1T2_T3T4"
stratum_sel      <- "T1-T2"

# Variables and identifiers used for the operating point
pred_var    <- "sccore_gep"
outcome_var <- "Metastasis_num"
weights_var <- "Weight_rescaled"
AJCC_LEVELS <- c("T1", "T2")
ID_CANDIDATES <- c("Patient_ID_SKY", "Skyline_ID", "Sample_id", "Skyline_Sample_ID", "Set_id")

# Validation (val102) input files
p102 <- list(
  clinical = file.path(results_dir, "intermediate", "clinical", "merged_clinical_data_imputed_d_squame_validation.csv"),
  mimp     = file.path(results_dir, "intermediate", "clinical", "clinical_data_imputed_d_squame_validation.Rdata"),
  weights  = e7_weights,
  preds    = file.path(results_dir, "intermediate", "validation", "d_squame", "recalibrated_SCCoreGEP_predictions.csv")
)

#-------------------------------------------------------------------------------

# 3. Operating point from discovery
#-------------------------------------------------------------------------------
# The operating point is taken from the discovery BWH T1-T2a threshold analysis
# (maximum finite weighted LR+) and applied to AJCC8 T1-T2 in val102
thr_tab <- readr::read_csv(file.path(dir_results_intermediate_sccore_gep,
                                     "discovery_T1_T2a_threshold_metrics_discoveryCal_OOB.csv"))
op_row <- thr_tab[is.finite(thr_tab$LR_pos_w), ]
op_row <- op_row[which.max(op_row$LR_pos_w), ]
OPERATING_POINT <- op_row$threshold
# Guard against silent changes of the operating point
stopifnot(abs(OPERATING_POINT - 0.069823237) < 1e-8)

#-------------------------------------------------------------------------------

# 4. Helper functions
#-------------------------------------------------------------------------------

#' Division that returns NA for missing or non-positive denominators
#'
#' @param a numeric. Numerator.
#' @param b numeric. Denominator.
#'
#' @return a / b, or NA_real_ if a or b is NA or b <= 0.
safeDiv <- function(a, b) if (is.na(a) || is.na(b) || b <= 0) NA_real_ else a / b

#' First of a set of candidate column names present in a data frame
#'
#' @param df data.frame. Data to search.
#' @param candidates character vector. Column names, in order of preference.
#'
#' @return The first matching column name, or NULL if none is present.
pickCol <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) NULL else hit[1]
}

#' First ID column shared by two data frames
#'
#' @param a data.frame. First data frame.
#' @param b data.frame. Second data frame.
#' @param candidates character vector. ID column names, in order of preference.
#' @param what character. Description used in the error message.
#'
#' @return The first shared column name; stops if there is none.
shared_id <- function(a, b, candidates, what) {
  hit <- candidates[candidates %in% colnames(a) & candidates %in% colnames(b)]
  if (length(hit) == 0) stop("No shared ", what)
  hit[1]
}

#' Load the first mids object from an .Rdata file
#'
#' @param path character. Path to the .Rdata file.
#'
#' @return The mids object (multiple imputations).
load_mids <- function(path) {
  e <- new.env(); load(path, envir = e)
  objs <- ls(e)
  is_mids <- vapply(objs, function(o) inherits(get(o, envir = e), "mids"), logical(1))
  if (!any(is_mids)) stop("No mids object in ", path)
  get(objs[is_mids][1], envir = e)
}

#' Read SCCore-GEP predictions on the percent scale
#'
#' The scale is detected from the data: values <= 1 are multiplied by 100.
#'
#' @param path character. Path to the predictions file.
#'
#' @return data.frame with an ID column and sccore_gep (percent), one row per ID.
readPreds <- function(path) {
  raw <- read.csv(path)
  pred_col <- pickCol(raw, c("sccore_gep_cox_calibrated", "cox", "sccore_gep"))
  vals <- suppressWarnings(as.numeric(raw[[pred_col]]))
  scale_factor <- if (max(vals, na.rm = TRUE) <= 1) 100 else 1   # to percent scale
  id_col <- pickCol(raw, ID_CANDIDATES)
  out <- data.frame(raw[[id_col]], vals * scale_factor)
  names(out) <- c(id_col, "sccore_gep")
  out[!duplicated(out[[id_col]]), ]
}

#' Weighted confusion matrix values for a binary rule in one completed dataset
#'
#' @param d data.frame. One completed dataset with the rule, outcome and weights.
#' @param rule character. Name of the binary (0/1) rule column.
#'
#' @return Named numeric vector with weighted TP, FP, FN and TN.
cells_w <- function(d, rule) {
  ok <- !is.na(d[[rule]]) & !is.na(d[[outcome_var]]) & !is.na(d[[weights_var]])
  r <- d[[rule]][ok]; y <- d[[outcome_var]][ok]; w <- d[[weights_var]][ok]
  c(TP = sum(w[r == 1 & y == 1]), FP = sum(w[r == 1 & y == 0]),
    FN = sum(w[r == 0 & y == 1]), TN = sum(w[r == 0 & y == 0]))
}

#-------------------------------------------------------------------------------

# 5. Sensitivity at the operating point in val102 AJCC8 T1-T2
#-------------------------------------------------------------------------------
# 5.1 Load completed datasets with staging, weights and predictions
mimp     <- load_mids(p102$mimp)
clinical <- read.csv(p102$clinical)
weights  <- read.csv(p102$weights)
preds    <- readPreds(p102$preds)

datasets <- lapply(seq_len(mimp$m), function(i) {
  d <- as.data.frame(complete(mimp, i))
  if (!("AJCC_staging" %in% colnames(d))) {                       # put back if mice dropped it
    id <- shared_id(d, clinical, ID_CANDIDATES, "ID for clinical")
    raw_small <- clinical[!duplicated(clinical[[id]]), c(id, "AJCC_staging")]
    d <- merge(d, raw_small, by = id, all.x = TRUE, sort = FALSE)
  }
  if (!("Metastasis_num" %in% colnames(d))) {
    d$Metastasis_num <- ifelse(as.character(d$Metastasis) == "Case", 1, 0)
  }
  id_w <- shared_id(d, weights, ID_CANDIDATES, "ID for weights")
  w_small <- weights[!is.na(weights$Weight_rescaled), c(id_w, "Weight_rescaled")]
  w_small <- w_small[!duplicated(w_small[[id_w]]), ]
  d <- merge(d, w_small, by = id_w, all.x = TRUE, sort = FALSE)
  id_p <- shared_id(d, preds, ID_CANDIDATES, "ID for predictions")
  merge(d, preds, by = id_p, all.x = TRUE, sort = FALSE)
})
print(table(datasets[[1]]$AJCC_staging, useNA = "ifany"))         # check level names

# 5.2 Restrict to AJCC8 T1-T2, apply the cut, pool the 2x2 cells across imputations
op_percent <- 100 * OPERATING_POINT                               # sccore_gep is in %
cells <- sapply(datasets, function(d) {
  s <- d[!is.na(d$AJCC_staging) & d$AJCC_staging %in% AJCC_LEVELS, , drop = FALSE]
  s$OP1 <- as.integer(suppressWarnings(as.numeric(s[[pred_var]])) > op_percent)
  cells_w(s, "OP1")
})
pooled <- rowMeans(cells)

se_op1  <- safeDiv(pooled["TP"], pooled["TP"] + pooled["FN"])
ppv_op1 <- safeDiv(pooled["TP"], pooled["TP"] + pooled["FP"])
message("AJCC8 T1-T2 at OP: weighted SE = ", signif(se_op1, 6),
        " | weighted PPV = ", signif(100 * ppv_op1, 4), "%")

# 5.3 Position of the operating point line (sensitivity, %)
se_op1_pct <- 100 * unname(se_op1)
op_line <- data.frame(Stratum = stratum_sel, Sensitivity = se_op1_pct)
#-------------------------------------------------------------------------------

# 6. Source functions and read data
#-------------------------------------------------------------------------------
source(file.path(code_dir, "functions", "colors.R"))
#-------------------------------------------------------------------------------

# Read data
#-------------------------------------------------------------------------------
res <- read.csv(res_fn)
#-------------------------------------------------------------------------------

# 7. Clinical decision thresholds (SCOUT survey)
#-------------------------------------------------------------------------------
# SCOUT consortium survey (n = 53 experts), 5-year metastatic risk (%):
# Rentroia-Pacheco B, Tokez S, Bramer EM, et al. Personalised decision making to
# predict absolute metastatic risk in cutaneous squamous cell carcinoma: development
# and validation of a clinico-pathological model. EClinicalMedicine. 2023;63:102150.
scout <- data.frame(decision = c("Follow-up", "Adjuvant RT", "Adjuvant systemic"),
                    median   = c(10, 20, 30),
                    q1       = c(5, 10, 20),
                    q3       = c(20, 20, 40))

# Standard of care: below Q1 of the follow-up threshold:
# at least 75% of experts would not include the patient in a follow-up schedule
soc_cutoff <- scout$q1[scout$decision == "Follow-up"]   # 5
#-------------------------------------------------------------------------------

# 8. Prepare the T1-T2 data
#-------------------------------------------------------------------------------
# 8.1 Subset the threshold metrics to AJCC8 T1-T2
t1t2_raw <- res %>%
  filter(stratification_system == strat_system_sel & stratum == stratum_sel)

t1t2_res <- t1t2_raw %>%
  mutate(Threshold = threshold,
         Stratum = stratum,
         `Post-test risk` = 100 * PPV_w,
         Sensitivity = 100 * SE_w,
         `% patients high-risk` = 100 * high_risk_w,
         `% patients low-risk` = 100 - `% patients high-risk`) %>%
  select(Stratum, Threshold, `Post-test risk`, Sensitivity,
         `% patients high-risk`, `% patients low-risk`)

# 8.2 Weighted prevalence of the combined T1-T2 stratum, derived from the metrics:
# TP fraction = PPV * P(high-risk) = SE * prevalence  ->  prevalence = PPV * P(high-risk) / SE
# (constant across thresholds; median taken to absorb rounding)
prev_t1t2 <- t1t2_raw %>%
  filter(SE_w > 0) %>%
  summarise(prevalence = median(100 * PPV_w * high_risk_w / SE_w, na.rm = TRUE)) %>%
  pull(prevalence)

stages_risk_ajcc8 <- data.frame(Stratum = stratum_sel,
                                prevalence = prev_t1t2)
message("Weighted T1-T2 prevalence (%): ", round(prev_t1t2, 2))

# 8.3 Decision bands (IQRs of the SCOUT thresholds)
# Bottom of the log axis: nearest power of 2 below the lowest plotted value
risk_vals <- c(t1t2_res$`Post-test risk`, stages_risk_ajcc8$prevalence, soc_cutoff)
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

# 9. Plot functions
#-------------------------------------------------------------------------------

#' PPV (post-test risk) versus sensitivity, with SCOUT decision bands
#'
#' Draws the precision-recall curve on a log2 PPV axis, the decision bands,
#' the stratum prevalence (dotted) and the operating point (dashed).
#'
#' @param stratum_sel character. Stratum to plot (e.g. "T1-T2").
#'
#' @return ggplot object (top panel of the figure).
make_ppv_plot <- function(stratum_sel) {
  ggplot(data = t1t2_res %>% filter(Stratum == stratum_sel) %>% arrange(Threshold),
         aes(x = Sensitivity, y = `Post-test risk`)) +
    geom_rect(data = decision_bands,
              aes(ymin = prev_min, ymax = prev_max, fill = band, xmin = -Inf, xmax = Inf),
              inherit.aes = F, alpha = .3) +
    scale_fill_manual(values = cols_bands, breaks = legend_order,
                      name = "Considered by SCOUT:") +
    geom_hline(data = stages_risk_ajcc8 %>% filter(Stratum == stratum_sel),
               aes(yintercept = prevalence), linetype = 3, linewidth = 0.4) +
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
    facet_wrap("Stratum", nrow = 1)
}

#' Percentage of patients classified high-risk versus sensitivity
#'
#' Stacked areas of high-risk and low-risk patients, with the operating point
#' (dashed).
#'
#' @param stratum_sel character. Stratum to plot (e.g. "T1-T2").
#'
#' @return ggplot object (bottom panel of the figure).
make_perc_hr_plot <- function(stratum_sel) {
  ggplot(data = t1t2_res %>%
           filter(Stratum == stratum_sel) %>%
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
    facet_wrap("Stratum", nrow = 1)
}
#-------------------------------------------------------------------------------

# 10. Make and save the figure
#-------------------------------------------------------------------------------
# 10.1 Combine the two panels
ppv_vs_se_t1t2     <- make_ppv_plot(stratum_sel)
perc_hr_vs_se_t1t2 <- make_perc_hr_plot(stratum_sel)

p <- ppv_vs_se_t1t2 + perc_hr_vs_se_t1t2 +
  plot_layout(nrow = 2, heights = c(3, 1))
# 10.2 Source data: the values plotted
source_data_df <- t1t2_res %>% select(-Threshold)
# 10.3 Save
write.csv(source_data_df, file.path(output_dir, "supp_figure_15B_t1t2_precision_recall_plot_source_data.csv"), row.names = F)
pdf(file.path(output_dir, "supp_figure_15B_t1t2_precision_recall_plot.pdf"), height = 3.8, width = 4.5)
print(p)
dev.off()
#-------------------------------------------------------------------------------
