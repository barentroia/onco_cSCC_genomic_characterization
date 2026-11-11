#===============================================================================
# Aim:      Generate Supplementary Figure 12: performance (C-index) of the GEP,
#           WES and combined GEP + WES models, and bootstrap inclusion
#           frequency and direction of the selected features.
# Author:   R. Ruiter
# Input:    - GEP_model, WES_model and Combined_model .rds files
#             (dir_results_intermediate_integration)
#           - GTF annotation (p9_gtf), to map Ensembl IDs to gene names
# Output:   - supp_figure_12_WESGEP_integration.pdf in dir_results_publication
# Requires: dir_results_intermediate_integration, dir_results_publication
#           (main.R); packages loaded by Model_input_processing.R
# Usage:    Sourced from main.R
#===============================================================================

# 1. Load annotation
#-------------------------------------------------------------------------------
gtf <- import(p9_gtf)
#-------------------------------------------------------------------------------

# 2. Gene name mapping, helper functions and models
#-------------------------------------------------------------------------------
# Create Ensembl gene ID -> external gene name mapping
gene_mappings <- as.data.frame(mcols(gtf)) %>%
  filter(type == "gene") %>%
  select(
    ensembl_gene_id = gene_id,
    external_gene_name = gene_name
  ) %>%
  distinct()


conflicts_prefer(dplyr::count)

# Named vector to map unversioned Ensembl IDs to gene names
gene_mappings <- gene_mappings %>%
  mutate(gene_id = str_split(ensembl_gene_id, "\\.", simplify = TRUE)[,1]) %>%
  select(gene_id, external_gene_name)

gene_map_vec <- setNames(
  gene_mappings$external_gene_name,
  gene_mappings$gene_id
)

#' Rename WES feature labels for plotting
#'
#' @param x character vector. Feature names.
#'
#' @return character vector with readable labels (e.g. "_del_bin" -> "Loss",
#'   "_amp_bin" -> "Gain").
clean_features <- function(x) {
  x <- str_replace(x, "Ras_MAPK", "RAS_MAPK/PI3K")
  x <- str_replace(x, "_del_bin$", "Loss")
  x <- str_replace(x, "_amp_bin$", "Gain")
  x
}


#' Bootstrap selection frequency of each feature
#'
#' @param obj list. Saved model object (output of bootstrap_validation()).
#' @param gene_map_vec named character vector. Ensembl ID -> gene name.
#'
#' @return tibble with feature (readable name) and n (number of bootstrap
#'   models in which the feature was selected).
get_freq <- function(obj, gene_map_vec) {
  
  if (!is.null(obj$Selected_Features)) {
    
    df <- tibble(
      feature = rownames(obj$Selected_Features),
      n = rowSums(obj$Selected_Features)
    )
    
  } else {
    
    feats <- unlist(lapply(obj$stacked_bootstrap_models, function(x) {
      names(x$coefficients)
    }))
    
    df <- tibble(feature = feats) %>%
      count(feature, name = "n")
  }
  
  df %>%
    mutate(
      feature_clean = str_replace(feature, "\\..*", ""),
      feature_clean = clean_features(feature_clean),
      feature_clean = ifelse(feature_clean == "Tumor", "Cellularity", feature_clean),
      
      mapped_name = ifelse(
        str_detect(feature_clean, "^ENSG"),
        gene_map_vec[feature_clean],
        feature_clean
      ),
      
      mapped_name = ifelse(is.na(mapped_name), feature_clean, mapped_name)
    ) %>%
    select(feature = mapped_name, n) %>%
    group_by(feature) %>%
    summarise(n = sum(n), .groups = "drop")
}

#' Coefficients of the final model with readable feature names
#'
#' @param model_obj list. Saved model object with Final_model2.
#' @param gene_map_vec named character vector. Ensembl ID -> gene name.
#'
#' @return data.frame with feature_name and coefficient.
extract_coefficients <- function(model_obj, gene_map_vec) {
  
  coefs <- model_obj$Final_model2$coefficients
  
  df <- data.frame(
    feature = names(coefs),
    coefficient = as.numeric(coefs)
  )
  
  df$feature_clean <- str_replace(df$feature, "\\..*", "")
  df$feature_clean <- clean_features(df$feature_clean)
  
  df <- df %>%
    mutate(
      is_ensg = str_detect(feature_clean, "^ENSG"),
      feature_name = ifelse(is_ensg, gene_map_vec[feature_clean], feature_clean)
    )
  
  df$feature_name[is.na(df$feature_name)] <- df$feature_clean[is.na(df$feature_name)]
  
  df %>%
    select(feature_name, coefficient)
}

# Load the models and extract final coefficients and selection frequencies
RNAseq <- readRDS(file.path(dir_results_intermediate_integration, "GEP_model_robj.rds"))
WES <-readRDS(file.path(dir_results_intermediate_integration, "WES_model_robj.rds"))
COMBINED <- readRDS(file.path(dir_results_intermediate_integration, "Combined_model_robj.rds"))

df_gep <- extract_coefficients(RNAseq, gene_map_vec)
df_wes <- extract_coefficients(WES, gene_map_vec)
df_combined <- extract_coefficients(COMBINED, gene_map_vec)

freq_gep <- get_freq(RNAseq, gene_map_vec)
freq_wes <- get_freq(WES, gene_map_vec)
freq_combined <- get_freq(COMBINED, gene_map_vec)

#-------------------------------------------------------------------------------

# 3. Performance (C-index) per model
#-------------------------------------------------------------------------------
#' C-index of a model: .632+ estimate and out-of-bag values
#'
#' @param obj list. Saved model object.
#'
#' @return list with boot632 (.632+ C-index) and oob (C-index per bootstrap,
#'   out-of-bag samples).
extract_metrics <- function(obj) {
  
  boot632 <- obj$`Optimism-Corrected0632plus` %>%
    filter(variable == "All", metric == "cind_model") %>%
    pull(value)
  
  oob <- obj$`All out of bag` %>%
    filter(variable == "All", metric == "cind_model") %>%
    pull(value)
  
  list(boot632 = boot632, oob = oob)
}

models <- list(
  GEP = RNAseq,
  WES = WES,
  Combined = COMBINED
)

oob_plot_df <- imap_dfr(models, function(obj, nm) {
  extract_metrics(obj)$oob %>%
    tibble(oob_value = .,
           group = nm)
})

metrics_df <- imap_dfr(models, function(obj, nm) {
  tibble(
    group = nm,
    boot632 = extract_metrics(obj)$boot632
  )
})

ref_line <- metrics_df %>%
  filter(group == "Combined") %>%
  pull(boot632)

#-------------------------------------------------------------------------------

# 4. Bootstrap inclusion frequency plots
#-------------------------------------------------------------------------------
#' Bar plot of bootstrap inclusion frequencies (top 75 features)
#'
#' Bars are coloured by the sign of the coefficient in the final model;
#' features in the final model are shown in bold.
#'
#' @param df tibble. Output of get_freq().
#' @param coef_df data.frame. Output of extract_coefficients().
#' @param model_name character. Plot title.
#' @param n_boot integer. Number of bootstrap repetitions.
#'
#' @return ggplot object.
plot_model <- function(df, coef_df, model_name, n_boot = 200) {
  
  df <- df %>%
    mutate(n_pct = 100 * (n / n_boot)) %>%
    left_join(coef_df, by = c("feature" = "feature_name")) %>%
    mutate(
      coef = ifelse(is.na(coefficient), 0, coefficient),
      fill_group = case_when(
        coef > 0 ~ "Positive selected",
        coef < 0 ~ "Negative selected",
        TRUE ~ "Not selected"
      ),
      label = ifelse(coef != 0, paste0("<b>", feature, "</b>"), feature)
    ) %>%
    arrange(desc(n_pct)) %>%
    slice_head(n = 75)
  
  ggplot(df, aes(x = reorder(label, n_pct), y = n_pct, fill = fill_group)) +
    geom_col() +
    coord_flip() +
    scale_y_continuous(limits = c(0, 100)) +
    scale_fill_manual(values = c(
      "Positive selected" = "red",
      "Negative selected" = "blue",
      "Not selected" = "grey80"
    )) +
    labs(title = model_name, x = NULL, y = "Bootstrap inclusion frequency (%)") +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.y = ggtext::element_markdown()
    )
}

p_freq_gep <- plot_model(freq_gep, df_gep, "GEP")
p_freq_wes <- plot_model(freq_wes, df_wes, "WES")
p_freq_combined <- plot_model(freq_combined, df_combined, "Combined")

#-------------------------------------------------------------------------------

# 5. Performance plots
#-------------------------------------------------------------------------------
ymin <- min(oob_plot_df$oob_value, metrics_df$boot632, na.rm = TRUE)
ymax <- max(oob_plot_df$oob_value, metrics_df$boot632, na.rm = TRUE)
ylim_low <- max(0, ymin - 0.02)
ylim_high <- min(1, ymax + 0.02)
#' Box plot of out-of-bag C-indices with the .632+ estimate
#'
#' The dashed line marks the .632+ C-index of the combined model.
#'
#' @param group_name character. Model to plot ("GEP", "WES" or "Combined").
#' @param ylim_low numeric. Lower y-axis limit.
#' @param ylim_high numeric. Upper y-axis limit.
#'
#' @return ggplot object.
plot_cindex <- function(group_name, ylim_low, ylim_high) {
  
  ggplot() +
    geom_boxplot(
      data = oob_plot_df %>% filter(group == group_name),
      aes(x = group, y = oob_value),
      fill = "white",
      colour = "black",
      width = 0.6,
      outlier.shape = NA
    ) +
    geom_jitter(
      data = oob_plot_df %>% filter(group == group_name),
      aes(x = group, y = oob_value),
      width = 0.25,
      alpha = 0.35,
      size = 2
    ) +
    geom_point(
      data = metrics_df %>% filter(group == group_name),
      aes(x = group, y = boot632),
      colour = "#8a181a",
      size = 4
    ) +
    geom_hline(yintercept = ref_line, linetype = "dashed", colour = "#8a181a") +
    scale_y_continuous(limits = c(ylim_low, ylim_high)) +
    labs(title = group_name, x = NULL, y = "C-index") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 20, hjust = 1),
      panel.border = element_rect(colour = "black", fill = NA),
      
      # add background y-axis grid lines
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor.y = element_line(color = "grey95", linewidth = 0.3),
      
      # optional: keep plot clean horizontally
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}

p_box_gep <- plot_cindex("GEP", ylim_low, ylim_high)
p_box_wes <- plot_cindex("WES", ylim_low, ylim_high)
p_box_combined <- plot_cindex("Combined", ylim_low, ylim_high)
#-------------------------------------------------------------------------------

# 6. Final aligned figure
#-------------------------------------------------------------------------------
final_plot <-
  (p_box_wes | p_box_gep | p_box_combined) /
  (p_freq_wes | p_freq_gep | p_freq_combined) +
  plot_layout(heights = c(1, 3))

#-------------------------------------------------------------------------------

# 7. Save figure
#-------------------------------------------------------------------------------
ggsave(
  file.path( dir_results_publication ,paste0("supp_figure_12_WESGEP_integration.pdf") ),
  final_plot,
  width = 14,
  height = 15
)