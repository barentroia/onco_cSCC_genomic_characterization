#---------------------------------------------------
# Aim: Show associations of gene alterations and immunosuppression
# Author: B. Rentroia Pacheco
# Input: Master mutation list
# Output: Forest plots
#---------------------------------------------------

#---------------------------------------------------
# 2. Build a forest plot with the Odds ratio for IS
#---------------------------------------------------

apply_models <- function(df,
                         var_of_interest,
                         group_of_interest,
                         model_type = "cox",
                         outcome = NULL,
                         time = NULL,
                         adjustment = NULL,
                         group_levels_show = NULL,
                         output_file = NULL){
  
  # Build formula
  rhs <- if(is.null(adjustment)){
    var_of_interest
  } else {
    paste(var_of_interest, adjustment, sep = " + ")
  }
  
  if(model_type == "cox"){
    formula_txt <- paste0("Surv(", time, ", ", outcome, ") ~ ", rhs)
  } else if(model_type == "logistic"){
    formula_txt <- paste0(outcome, " ~ ", rhs)
  }
  
  formula_obj <- as.formula(formula_txt)
  
  # mutation naming fix
  if(var_of_interest == "Any_mutation"){
    var_of_interest_r <- paste0(var_of_interest,"Yes")
  } else{
    var_of_interest_r <- var_of_interest
  }
  
  # model fitting
  model_results <- df %>%
    group_by(.data[[group_of_interest]]) %>%
    nest() %>%
    mutate(
      model = map(data, ~{
        if(model_type == "cox"){
          survival::coxph(formula_obj, data = .x)
        } else{
          glm(formula_obj, data = .x, family = binomial)
        }
      }),
      tidy = map(model, ~broom::tidy(.x, exponentiate = TRUE, conf.int = TRUE))
    ) %>%
    unnest(tidy) %>%
    filter(term == var_of_interest_r) %>%
    select(.data[[group_of_interest]], estimate, conf.low, conf.high, p.value) %>%
    ungroup() %>%
    as.data.frame()
 
  # restrict groups shown: no infinite coefficients
  group_levels_show <- intersect(
    group_levels_show,
    model_results[!is.infinite(model_results$conf.high)&!is.na(model_results$conf.high),1]
  )
  
  
  x_max <- model_results %>%
    filter(.data[[group_of_interest]] %in% group_levels_show) %>%
    pull(conf.high) %>%
    max(na.rm = TRUE)
  
  adj_txt <- ifelse(
    is.null(adjustment),
    ", no adjustment",
    paste0(", with adjustment: ", adjustment)
  )
  
  # plotting dataframe
  plot_df <- model_results %>%
    filter(.data[[group_of_interest]] %in% group_levels_show) %>%
    mutate(
      Significant = ifelse(p.value < 0.1,"signif","ns"),
      Significant = ifelse(p.value < 0.05,"ssignif", Significant),
      label = sprintf("%.2f (%.2f–%.2f), p = %.3f",
                      estimate, conf.low, conf.high, p.value)
    )
  
  xlab <- ifelse(model_type == "cox",
                 "Hazard Ratio (log scale)",
                 "Odds Ratio (log scale)")
  
  p <- ggplot(plot_df,
              aes(x = estimate,
                  y = reorder(.data[[group_of_interest]], estimate),
                  color = Significant)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                   height = 0.2,
                   linewidth = 0.8) +
    geom_text(aes(x = x_max * 1.4, label = label),
              hjust = 0,
              size = 3.5,
              color = "black") +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_x_log10(
      limits = c(min(plot_df$conf.low, na.rm = TRUE) * 0.8,
                 x_max * 1.6)
    ) +
    scale_color_manual(values = c(
      "ssignif" = "darkred",
      "signif" = "goldenrod2",
      "ns" = "black"
    )) +
    labs(
      x = xlab,
      y = "",
      title = paste(toupper(model_type), "models per",
                    gsub("s$","",group_of_interest),
                    adj_txt)
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(size = 12),
      legend.position = "none",
      plot.margin = margin(5.5, 120, 5.5, 5.5)
    ) +
    coord_cartesian(clip = "off")
  
  if(!is.null(output_file)){
    openxlsx::write.xlsx(model_results, output_file)
  }
  
  return(list(results = model_results,
              forest_plot = p))
}

# Adjust the variables so that we can run the logistic regression models:
path_sample$IS_classification_bin = ifelse(path_sample$IS_classification=="Yes",1,0)
mut_data_drivers$IS_classification_bin = ifelse(mut_data_drivers[,IS_classification]=="Yes",1,0)

# Find the pathways with enough mutations:
IS_pathways_adj <- apply_models(
  df = path_sample,
  var_of_interest = "Weighted_burden",
  group_of_interest = "Pathway",
  model_type = "logistic",
  outcome = "IS_classification_bin",
  adjustment = "Percent_Bait_with_callable_coverage",
  group_levels_show = setdiff(pathways_with_enough_mutations,c("MYC","CASP8","ZNF750","TERT","KMT2D","CHUK")),
  output_file = file.path(dir_results_script,"04_02_03_forest_plots_IS_pathways_adjusted_source_data.xlsx") 
)

IS_pathways_unadj <- apply_models(
  df = path_sample,
  var_of_interest = "Weighted_burden",
  group_of_interest = "Pathway",
  model_type = "logistic",
  outcome = "IS_classification_bin",
  adjustment = NULL,
  group_levels_show = setdiff(pathways_with_enough_mutations,c("MYC","CASP8","ZNF750","TERT","KMT2D","CHUK")),
  output_file = file.path(dir_results_script,"04_02_03_forest_plots_IS_pathways_unadjusted_source_data.xlsx") 
)

IS_genes_adj <- apply_models(
  df = mut_data_drivers,
  var_of_interest = "gene_weight",
  group_of_interest = "Genes",
  model_type = "logistic",
  outcome = "IS_classification_bin",
  adjustment = "Percent_Bait_with_callable_coverage",
  group_levels_show = genes_with_enough_mutations,
  output_file = file.path(dir_results_script,"04_02_03_forest_plots_IS_genes_adjusted_source_data.xlsx") 
)

IS_genes_unadj <- apply_models(
  df = mut_data_drivers,
  var_of_interest = "gene_weight",
  group_of_interest = "Genes",
  model_type = "logistic",
  outcome = "IS_classification_bin",
  adjustment = NULL,
  group_levels_show = genes_with_enough_mutations,
  output_file = file.path(dir_results_script,"04_02_03_forest_plots_IS_genes_unadjusted_source_data.xlsx") 
)



# Save plots:
# Adjusted plot:
combined_plot <- IS_genes_adj$forest_plot+ IS_pathways_adj$forest_plot
ggsave(file.path(dir_results_script,"04_02_04_forest_plots_IS_genes_adjusted.pdf"),
       combined_plot,
       width = 12,
       height = 5,
       dpi = 300)

# Now the same for CN alterations:
copy_number_data_fplot$IS_bin = ifelse(copy_number_data_fplot[,IS_classification]=="Yes",1,0)
apply_logisticmodels_cnv = function(copy_number_data_fplot,
                                    outcome = "event",
                                    adjustment = NULL,
                                    split_cnv = TRUE,
                                    glm_file_path = NULL){
  
  chr_cols <- grep("^chr", names(copy_number_data_fplot), value = TRUE)
  
  glm_results <- map_df(chr_cols, function(chr) {
    
    df <- copy_number_data_fplot
    
    if(split_cnv){
      df[[chr]] <- factor(df[[chr]], levels = c(0, -1, 1))
    }
    
    if(!is.null(adjustment)){
      formula_txt <- paste0(outcome, " ~ ", chr, " + ", adjustment)
    } else {
      formula_txt <- paste0(outcome, " ~ ", chr)
    }
    
    model <- glm(
      as.formula(formula_txt),
      data = df,
      family = "binomial"
    )
    
    tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
      filter(!term %in% c("(Intercept)", adjustment)) %>%
      mutate(Chromosome = chr)
  })
  
  if(split_cnv){
    
    glm_results_clean <- glm_results %>%
      mutate(
        CN_status = case_when(
          grepl("-1", term) ~ "Loss",
          grepl("1", term) ~ "Gain"
        ),
        Arm = Chromosome,
        q_value = p.adjust(p.value, method = "BH"),
        Chr_arm_status = paste(Arm, CN_status, sep = " - ")
      ) %>%
      select(Arm, CN_status, Chr_arm_status,
             estimate, conf.low, conf.high,
             p.value, q_value)
    
    x_max <- glm_results_clean %>%
      filter(Chr_arm_status %in% c(chr_with_enough_alterations_gain,
                                   chr_with_enough_alterations_loss)) %>%
      pull(conf.high) %>% max(na.rm = TRUE)
    
    p_data <- glm_results_clean %>%
      filter(Chr_arm_status %in% c(chr_with_enough_alterations_gain,
                                   chr_with_enough_alterations_loss)) %>%
      mutate(
        Significant = ifelse(p.value < 0.1,"signif","ns"),
        label = sprintf("%.2f (%.2f–%.2f), p = %.3f",
                        estimate, conf.low, conf.high, p.value),
        Ylab = paste0(Arm,"-",CN_status)
      ) %>%
      mutate(Significant=ifelse(p.value<0.05,"ssignif",Significant))
    
  } else {
    
    glm_results_clean <- glm_results %>%
      mutate(q_value = p.adjust(p.value, method="BH")) %>%
      select(Chromosome,estimate,conf.low,conf.high,p.value,q_value)
    
    x_max <- glm_results_clean %>%
      filter(Chromosome %in% chr_with_enough_alterations_any) %>%
      pull(conf.high) %>% max(na.rm = TRUE)
    
    p_data <- glm_results_clean %>%
      filter(Chromosome %in% chr_with_enough_alterations_any) %>%
      mutate(
        Significant = ifelse(p.value < 0.1,"signif","ns"),
        label = sprintf("%.2f (%.2f–%.2f), p = %.3f",
                        estimate, conf.low, conf.high, p.value),
        Ylab = Chromosome
      ) %>%
      mutate(Significant=ifelse(p.value<0.05,"ssignif",Significant))
  }
  
  p = p_data %>%
    
    ggplot(aes(x = estimate,
               y = reorder(Ylab, estimate),
               color = Significant)) +
    
    geom_point(size = 3) +
    
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                   height = 0.2,
                   linewidth = 0.8) +
    
    geom_vline(xintercept = 1, linetype = "dashed") +
    
    geom_text(aes(x = x_max * 1.4, label = label),
              hjust = 0,
              size = 3.5,
              color = "black") +
    
    scale_x_log10() +
    
    scale_color_manual(values = c(
      "ssignif" = "darkred",
      "signif" = "goldenrod2",
      "ns" = "black"
    )) +
    
    labs(
      x = "Odds Ratio (log scale)",
      y = "",
      title = "Logistic Models per chr arm"
    ) +
    
    theme_minimal(base_size = 13) +
    
    theme(
      plot.title = element_text(size = 12),
      legend.position = "none",
      plot.margin = margin(5.5,120,5.5,5.5)
    ) +
    
    coord_cartesian(clip="off")
  
  if(!is.null(glm_file_path)){
    write.xlsx(glm_results_clean, glm_file_path)
  }
  
  return(list(
    "GLM_results" = glm_results_clean,
    "forest_plot" = p
  ))
}

IS_cnv_adjusted = apply_logisticmodels_cnv(copy_number_data_fplot,
                                           outcome = "IS_bin",
                                           adjustment="Percent_Bait_with_callable_coverage",
                                           glm_file_path=file.path(dir_results_script,"04_02_03_forest_plots_IS_cnv_adjusted_source_data.xlsx"))
IS_cnv_unadjusted = apply_logisticmodels_cnv(copy_number_data_fplot,
                                           outcome = "IS_bin",
                                           adjustment=NULL,
                                           glm_file_path=file.path(dir_results_script,"04_02_03_forest_plots_IS_cnv_unadjusted_source_data.xlsx"))

ggsave(file.path(dir_results_script,"04_02_03_forest_plots_IS_cnv_unadjusted.pdf"),IS_cnv_unadjusted$forest_plot,width=5,height=7.5)
ggsave(file.path(dir_results_script,"04_02_03_forest_plots_IS_cnv_adjusted.pdf"),IS_cnv_adjusted$forest_plot,width=5,height=7.5)

