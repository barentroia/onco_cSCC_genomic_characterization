#---------------------------------------------------
# Aim: Make supplementary figure 10A 
# Author: R.Ruiter
# Input: Processed performance dataset of models of DvP, non-DvP, 
#        Early integration,Late feature integration,
#        Late prediction integration
# Output: supplementary figure 10A 
#---------------------------------------------------


# -------------------------------
# Colors
# -------------------------------

base_cols <- c("#A03A6FFF", "#266986FF")
fill_cols <- alpha(base_cols, 0.4)


# -------------------------------
# Load data
# -------------------------------
perf_ds <- read_csv(
  file.path(dir_results_intermediate_sccore_gep , "perf_OOB_apparent_bootstrap.csv"))

perf_632plus <- read_csv(
  file.path(dir_results_intermediate_sccore_gep , "perf_632plus.csv"))

perf_apparent <- read_csv(
  file.path(dir_results_intermediate_sccore_gep , "perf_apparent.csv"))

gene_count <- read_csv(
  file.path(dir_results_intermediate_sccore_gep , "model_gene_counts.csv"))


# -------------------------------
# Feature levels
# -------------------------------

feature_levels <- c(
  "DvP",
  "non-DvP",
  "Early integration",
  "Late feature integration",
  "Late prediction integration"
)


# -------------------------------
# Model parser
# -------------------------------

parse_model <- function(df){
  df %>%
    mutate(
      # -----------------------------
      # Integration strategy
      # -----------------------------
      integration = case_when(
        str_detect( model,"DvP" ) &
          !str_detect( model, "non_DvP|Early|Late") ~ "DvP",
        str_detect( model, "non_DvP" ) ~"non-DvP",
        str_detect( model, "Early_integration" ) ~ "Early integration",
        str_detect( model, "Late_feature_integration") ~
          "Late feature integration",
        str_detect(model,"Late_prediction_integration") ~
          "Late prediction integration",
        TRUE ~ NA_character_
      ),
      integration = factor( integration,levels = feature_levels ),
      # -----------------------------
      # Model type
      # -----------------------------
      model_names = case_when(
        str_detect( model,"^RSF_" ) ~"RSF",
        str_detect(model, "^coxnet_" ) ~ "Regularized Cox",
        TRUE ~ NA_character_),
      facet_group = factor(model_names, levels = c( "Regularized Cox","RSF" )))
}

# -------------------------------
# Variables to loop
# -------------------------------
variables_to_plot <- c(
  "All"
)

# -------------------------------
# Output directory
# -------------------------------



dir.create(dir_results_publication , recursive = TRUE, showWarnings = FALSE)

# -------------------------------
# Plot function
# -------------------------------
make_plot <- function(var_name){
  message("Processing: ", var_name)
  # -------------------------------
  # Determine groups dynamically
  # -------------------------------
  groups_var <- unique( perf_ds[ perf_ds$variable == var_name,]$group)
  
  # -------------------------------
  # Special handling by variable
  # -------------------------------
  if(var_name == "All"){
    groups_use <- "All"
    facet_by_group <- FALSE
  } 
  # -------------------------------
  # Filter OOB dataset
  # -------------------------------
  perf_ds_f <- perf_ds %>%
    filter(
      variable == var_name,
      group %in% groups_use,
      metric == "wcind_model",
      type == "OOB"
    ) %>%
    parse_model() %>%
    filter(
      model_names %in% c(
        "Regularized Cox",
        "RSF"
      )
    ) %>%
    mutate(
      group = factor(group)
    )
  
  # -------------------------------
  # Filter 0.632+ dataset
  # -------------------------------
  
  perf_632plus_f <- perf_632plus %>%
    filter(
      variable == var_name,
      group %in% groups_use
    ) %>%
    parse_model() %>%
    filter(
      model_names %in% c(
        "Regularized Cox",
        "RSF"
      )
    ) %>%
    mutate(
      legend_var = ".632+ estimate",
      group = factor(group)
    )
  
  
  # -------------------------------
  # Filter apparent dataset
  # -------------------------------
  
  perf_apparent_f <- perf_apparent %>%
    filter(
      variable == var_name,
      group %in% groups_use
    ) %>%
    parse_model() %>%
    filter(
      model_names %in% c(
        "Regularized Cox",
        "RSF"
      )
    ) %>%
    mutate(
      legend_var = "Apparent",
      group = factor(group)
    )
  
  
  # -------------------------------
  # Gene counts
  # -------------------------------
  
  gene_counts_f <- gene_count %>%
    parse_model() %>%
    filter(
      model_names %in% c(
        "Regularized Cox",
        "RSF"
      )
    ) %>%
    group_by(
      model_names,
      integration,
      facet_group
    ) %>%
    summarise(
      n_genes = dplyr::first(gene_count),
      .groups = "drop"
    )
  
  
  # -------------------------------
  # Build facet formula dynamically
  # -------------------------------
  
  if(facet_by_group){
    facet_formula <- as.formula(
      "group ~ facet_group"
    )
  } else {
    facet_formula <- as.formula(
      "~ facet_group" )
  }
  
  
  # -------------------------------
  # Reference line
  #
  # Reference:
  # Late feature integration
  # Regularized Cox
  # .632+ estimate
  # -------------------------------
  if(facet_by_group){
    ref_line_df <- perf_632plus_f %>%
      filter(
        integration ==
          "Late feature integration",
        model_names ==
          "Regularized Cox"
      ) %>%
      dplyr::select(
        group,
        ref_value = value
      )
  } else {
    ref_line_df <- perf_632plus_f %>%
      filter(
        integration ==
          "Late feature integration",
        model_names ==
          "Regularized Cox"
      ) %>%
      summarise(
        ref_value = mean(
          value,
          na.rm = TRUE
        ))
  }
  
  # -------------------------------
  # Attach reference value
  # -------------------------------
  ref_value <- ref_line_df$ref_value[1]

  # -------------------------------
  # Plot source data
  # -------------------------------
  
  plot_source_df <- bind_rows(
    perf_ds_f %>%
      mutate(dataset = "OOB" ),
    perf_632plus_f %>%
      mutate( dataset = ".632plus" ),
    perf_apparent_f %>%
      mutate(dataset = "Apparent" )
  ) %>%
    mutate( variable = var_name,
      ref_632plus_final = ref_value
    ) %>%
    dplyr::select(
      variable,
      group,
      integration,
      model,
      model_names,
      value,
      dataset,
      ref_632plus_final
    )
  
  # -------------------------------
  # Plot 1:
  # Performance
  # -------------------------------
  
  p1 <- ggplot(
    perf_ds_f,
    aes( x = integration, y = value)) +
    # OOB boxplot
    geom_boxplot(
      aes( colour = "OOB" ),
      fill = alpha( "#266986FF",  0.4  ),
      outlier.shape = NA) +
    # OOB individual values
    geom_jitter(
      aes(colour = "OOB", fill = "OOB"),width = 0.15,size = 1, alpha = 0.8) +
    # .632+ estimate
    geom_point(
      data = perf_632plus_f,aes( x = integration, y = value,
                                 colour = ".632+ estimate",
                                 fill = ".632+ estimate"
                                 ),
      size = 4,shape = 23, stroke = 1) +
    # Apparent estimate
    geom_point(
      data = perf_apparent_f,
      aes( x = integration, y = value,colour = "Apparent",fill = "Apparent"),
      size = 4,shape = 23, stroke = 1) +
    # Reference line
    geom_hline(
      data = ref_line_df,
      aes(yintercept = ref_value, colour = ".632+ estimate final signature"),
      linetype = "dashed",linewidth = 1, inherit.aes = FALSE) +
    # Facets
    facet_grid(facet_formula ) +
    # Axis labels
    ylab("Weighted C-index") +
    xlab(  NULL) +
    # Theme
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text( size = 18),
      axis.title.y = element_text( size = 22 ),
      strip.text = element_text( size = 20),
      legend.text =element_text( size = 18 ),
      legend.title = element_text( size = 20 ),
      legend.position ="right" ) +
    
    # -----------------------------
  # Colour legend
  # -----------------------------
  
  scale_colour_manual(
    name = "Estimate type",
    values = c("OOB" ="#266986FF",
               ".632+ estimate" = "red",
               "Apparent" =  "#A03A6FFF",
               ".632+ estimate final signature" =  "blue" ) ) +
  # -----------------------------
  # Fill legend
  # -----------------------------
  scale_fill_manual( name = "Estimate type",
    values = c( "OOB" = "#266986FF", 
                ".632+ estimate" = "red", 
                "Apparent" = "#A03A6FFF") ) +
  # -----------------------------
  # Guides
  # -----------------------------
  guides(
    fill = "none",
    shape = "none",
    linetype = "none",
    colour = guide_legend(
      title = "Estimate Type",
      override.aes = list(
        fill = c(  "red", NA, base_cols[1],fill_cols[2] ),
        linetype = c( "solid", "dashed", "solid",  "solid" ),
        shape = c(23, NA, 23, NA ))))
  
  
  # -------------------------------
  # Plot 2:
  # Gene counts
  # -------------------------------
  
  p2 <- ggplot(
    gene_counts_f,
    aes(x = integration, y = n_genes)) +
    geom_col(
      fill = "grey",
      width = 0.6) +
    facet_grid(
      ~facet_group) +
    ylab( "Number of genes") +
    xlab("Feature integration strategy" ) +
    theme_bw() +
    theme(axis.title.x = element_text( size = 22),
      axis.title.y = element_text( size = 22 ),
      axis.text.x =element_text( angle = 45,hjust = 1, size = 16 ),
      axis.text.y =element_text( size = 18 ),
      strip.text =element_blank(),
      legend.position ="none" )
  
  # -------------------------------
  # Save PDF
  # -------------------------------
  pdf(file.path( dir_results_publication ,paste0("supp_figure_10a_", var_name,
                                                       "_performance_boxplot.pdf") ),width = 20, height =10 +(length(unique(perf_ds_f$group)) * 2))
  print(p1 /p2 +plot_layout( heights = c( 3,1)) )
  dev.off()

}

# -------------------------------
# Run loop
# -------------------------------

for(v in variables_to_plot){
  print(v)
  make_plot(v)
  
}
