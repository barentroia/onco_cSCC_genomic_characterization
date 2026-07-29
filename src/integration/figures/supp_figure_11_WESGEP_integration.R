#---------------------------------------------------
# Aim: Make supplementary figure 11
# Author: R.Ruiter
# Input: gtf file
#        integration models
#        
# Output: supplementary figure 11
#---------------------------------------------------
gtf <- import(p9_gtf)
#-------------------------------------------------------------------------------

# 1. Preprocess DvP signature
#-------------------------------------------------------------------------------
# Create Ensembl gene ID -> external gene name mapping
gene_mappings <- as.data.frame(mcols(gtf)) %>%
  filter(type == "gene") %>%
  select(
    ensembl_gene_id = gene_id,
    external_gene_name = gene_name
  ) %>%
  distinct()

# View result
head(gene_mappings)

conflicts_prefer(dplyr::count)
# gene/mutation mapping and cleaning

gene_mappings <- gene_mappings %>%
  mutate(gene_id = str_split(ensembl_gene_id, "\\.", simplify = TRUE)[,1]) %>%
  select(gene_id, external_gene_name)

gene_map_vec <- setNames(
  gene_mappings$external_gene_name,
  gene_mappings$gene_id
)

clean_features <- function(x) {
  x <- str_replace(x, "Ras_MAPK", "RAS_MAPK/PI3K")
  x <- str_replace(x, "_del_bin$", "Loss")
  x <- str_replace(x, "_amp_bin$", "Gain")
  x
}


# frequency function
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

#Coefficients extraction
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

# Models
RNAseq <- readRDS(file.path(dir_results_intermediate_integration, "GEP_model_robj.rds"))
WES <-readRDS(file.path(dir_results_intermediate_integration, "WES_model_robj.rds"))
COMBINED <- readRDS(file.path(dir_results_intermediate_integration, "Combined_model_robj.rds"))

df_gep <- extract_coefficients(RNAseq, gene_map_vec)
df_wes <- extract_coefficients(WES, gene_map_vec)
df_combined <- extract_coefficients(COMBINED, gene_map_vec)

freq_gep <- get_freq(RNAseq, gene_map_vec)
freq_wes <- get_freq(WES, gene_map_vec)
freq_combined <- get_freq(COMBINED, gene_map_vec)

# =========================
# performance
# =========================
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

# =========================
# frequency
# =========================
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

# =========================
# plot performance
# =========================
ymin <- min(oob_plot_df$oob_value, metrics_df$boot632, na.rm = TRUE)
ymax <- max(oob_plot_df$oob_value, metrics_df$boot632, na.rm = TRUE)
ylim_low <- max(0, ymin - 0.02)
ylim_high <- min(1, ymax + 0.02)
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
# =========================
# final aligned figure
# =========================
final_plot <-
  (p_box_wes | p_box_gep | p_box_combined) /
  (p_freq_wes | p_freq_gep | p_freq_combined) +
  plot_layout(heights = c(1, 3))

# =========================
# SAVE
# =========================
ggsave(
  file.path( dir_results_publication ,paste0("supp_figure_11_WESGEP_integration.pdf") ),
  final_plot,
  width = 14,
  height = 15
)