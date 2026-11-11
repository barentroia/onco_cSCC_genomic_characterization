#===============================================================================
# Aim:      Generate the spatial transcriptomics panel of Figure 5A: per
#           SCCore-GEP gene, the fraction of expression contributed by each
#           cell type in the spatial sequencing data (14 patients).
# Input:    - Final SCCore-GEP model (coxnet_Late_feature_integration_robj.rds)
#           - Spatial sequencing Seurat object (seurat_object_fn)
#           - GTF annotation (p9_gtf), to map Ensembl IDs to gene symbols
# Output:   - figure_5a_ST_barplot_contributions_N14.pdf and .xlsx (source data)
#             in results_dir/publication
# Requires: code_dir, results_dir (main.R) and config.R paths
# Usage:    Sourced from main.R
# Notes:    Model genes absent from the spatial panel are shown as "NA" bars so
#           the gene order matches the Figure 5A heatmap.
#===============================================================================

# 1. Load libraries
#-------------------------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(reshape2)
library(conflicted)
library(rtracklayer)

source(file.path(code_dir,"functions","colors.R"))
source(file.path(code_dir,"functions","single_cell_inputFuncs.R"))

conflicts_prefer(base::intersect)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(base::as.factor)
conflicts_prefer(base::setdiff)

#-------------------------------------------------------------------------------

# 2. General setup
#-------------------------------------------------------------------------------
# N = 14 patients (SCC09 and SCC10 not included)
patients <- sprintf("SCC%02d", c(1:8, 11:16))

# Merge fine annotations into broader cell types (anything not listed keeps its name)
celltype_merge <- c(
  # immune cells
  large_pre_B_cell  = "immune_cell",
  macrophage        = "immune_cell",
  monocyte          = "immune_cell",
  neutrophils       = "immune_cell",
  plasma_cell       = "immune_cell",
  T_cell            = "immune_cell",
  others            = "immune_cell",
  # pilosebaceous unit
  sebaceous_gland   = "pilosebaceous_cells",
  hair_follicle     = "pilosebaceous_cells",
  # undefined stromal cells + myofibroblasts
  stormal_cell      = "fibroblast",
  myofibroblast     = "fibroblast",
  # tumor melanocytes counted as tumor keratinocytes
  tumor_melanocytes = "tumor_keratinocytes"
)

# Plot order + colours; "NA" = model gene not present in the ST panel
celltype_colors <- c(
  immune_cell         = "#A53C73FF",
  sweat_gland         = "#823575FF",
  pilosebaceous_cells = "#B46727FF",
  keratinocytes       = "#C2CBCDFF",
  endothelial_cell    = "#A7A7A7FF",
  fibroblast          = "#7E7568FF",
  smooth_muscle       = "#A4D8E5FF",
  tumor_keratinocytes = "#2FBCE3FF",
  tumor_fibroblast    = "#1F6582FF",
  "NA"                = "white"
)
celltype_levels <- names(celltype_colors)

# output folder
output_dir <- file.path(results_dir, "publication")
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}
#-------------------------------------------------------------------------------

# 3. Final model
#-------------------------------------------------------------------------------
sccoregep_fn <- file.path(results_dir, "intermediate", "sccore_gep", "Discovery_models", "coxnet_Late_feature_integration_robj.rds") 

GEP_final_model <- readRDS(file.path(sccoregep_fn))
coeff_genes <- GEP_final_model$Final_model2$coefficients
stopifnot(length(coeff_genes) == 23)

# Ensembl ID -> gene symbol 
gene_map <- import(p9_gtf) %>%
  as.data.frame() %>%
  distinct(gene_id, gene_name) %>%
  mutate(gene_id_nover = sub("\\..*$", "", gene_id))

names(coeff_genes) <- gene_map$gene_name[match(names(coeff_genes), gene_map$gene_id_nover)]
stopifnot(!anyNA(names(coeff_genes)))
coeff_genes <- sort(coeff_genes)   # gene order follows coefficients (aligned with heatmap)

#-------------------------------------------------------------------------------

# 4. Load and prepare ST data
#-------------------------------------------------------------------------------
# Object has no spatial coordinates and only contains the genes of interest,
# so the usual QC checks are not possible.
ST_CSCC_seurat <- LoadSeuratRds(seurat_object_fn)
all_RNA <- LayerData(ST_CSCC_seurat, assay = "RNA", layer = "data")

genes_in_ST   <- intersect(names(coeff_genes), rownames(all_RNA))
genes_missing <- setdiff(names(coeff_genes), genes_in_ST)
message("Model genes not in ST data: ", paste(genes_missing, collapse = ", "))  # LCE1F, DEFB103A

meta <- ST_CSCC_seurat@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(cell,
            patient_id = sample_id,
            custom_id  = as.character(annotate),
            custom_id  = coalesce(unname(celltype_merge[custom_id]), custom_id),
            custom_id  = factor(custom_id, levels = setdiff(celltype_levels, "NA"))) %>%
  filter(patient_id %in% patients)

stopifnot(all(patients %in% meta$patient_id))

#-------------------------------------------------------------------------------

# 5. Generate cell type plot
#-------------------------------------------------------------------------------
expr_long <- as.matrix(all_RNA[genes_in_ST, meta$cell, drop = FALSE]) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  pivot_longer(-cell, names_to = "gene", values_to = "value") %>%
  inner_join(meta, by = "cell")

# Summed expression per patient x cell type x gene
counts_per_patient <- expr_long %>%
  group_by(patient_id, custom_id, gene) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(patient_id, gene) %>%
  mutate(frac_within_patient = value / sum(value)) %>%
  group_by(gene) %>%
  mutate(frac_of_total = value / sum(value)) %>%
  ungroup() %>%
  mutate(gene = factor(gene, levels = rev(names(coeff_genes))))

# Pooled distribution
to_plot <- counts_per_patient %>%
  group_by(custom_id, gene) %>%
  summarise(frac_counts_total = sum(frac_of_total), .groups = "drop") %>%
  mutate(across(c(custom_id, gene), as.character)) %>%
  # placeholder bars for model genes absent from ST, to align with heatmap
  bind_rows(tibble(custom_id = "NA", gene = genes_missing, frac_counts_total = 1)) %>%
  mutate(custom_id = factor(custom_id, levels = celltype_levels),
         gene      = factor(gene,      levels = rev(names(coeff_genes))))

# Plot barplots
p_pooled <- ggplot(to_plot, aes(x = gene, y = frac_counts_total, fill = custom_id)) +
  geom_col() +
  scale_fill_manual(values = celltype_colors) +
  coord_flip() +
  theme_bw()

pdf(file.path(output_dir, "figure_5a_ST_barplot_contributions_N14.pdf"), width = 7, height = 6)
print(p_pooled)
dev.off()

#-------------------------------------------------------------------------------

# 6. Save source data
#-------------------------------------------------------------------------------
writexl::write_xlsx(to_plot,file.path(output_dir, "figure_5a_ST_barplot_contributions_N14.xlsx"))


