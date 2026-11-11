#===============================================================================
# Aim:      Generate Supplementary Figure 11A: per SCCore-GEP gene, the
#           fraction of expression contributed by each cell type in the
#           single-cell data of Ji et al. 2020 (GSE144236).
# Input:    - Final SCCore-GEP model (coxnet_Late_feature_integration_robj.rds)
#           - Single-cell counts and metadata (SC_counts_fn, SC_metadata_fn)
#           - GTF annotation (p9_gtf), to map Ensembl IDs to gene symbols
# Output:   - AndrewLJi_seurat_processed.rds (cached processed Seurat object)
#           - supp_figure_11A_SC_Andrew_barplot_contributions_counts.pdf and
#             .xlsx (source data)
# Requires: results_dir (main.R) and config.R paths
# Usage:    Sourced from main.R
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
conflicts_prefer(base::union)

#-------------------------------------------------------------------------------

# 2. General setup
#-------------------------------------------------------------------------------
# Gene symbols that differ between the bulk GTF and the single-cell data (bulk = sc)
sc_gene_rename <- c(CCN5 = "WISP2", CYSRT1 = "C9orf169")

# Map annotations (level2_celltype:tum.norm) to broader cell types.
# Smooth muscle cells are not defined in this data set; Multiplet and Melanocyte are excluded.
immune_types <- c("LC", "CLEC9A", "PDC", "CD1C", "B Cell", "Mac", "ASDC", "Tcell", "MDSC", "NK")
celltype_groups <- list(
  immune_cell         = c(paste0(immune_types, ":Normal"), paste0(immune_types, ":Tumor")),
  keratinocytes       = c("Keratinocyte:Normal", "Normal_KC_Basal:Normal", "Normal_KC_Diff:Normal",
                          "Normal_KC_Cyc:Normal"),
  endothelial_cell    = c("Endothelial Cell:Normal", "Endothelial Cell:Tumor"),
  fibroblast          = "Fibroblast:Normal",
  tumor_keratinocytes = c("Keratinocyte:Tumor", "Tumor_KC_Diff:Tumor", "Tumor_KC_Basal:Tumor",
                          "Tumor_KC_Cyc:Tumor", "TSK:Tumor"),
  tumor_fibroblast    = "Fibroblast:Tumor",
  sweat_gland         = c("Eccrine:Normal", "Eccrine:Tumor"),
  pilosebaceous_cells = c("Pilosebaceous:Normal", "Pilosebaceous:Tumor")
)
celltype_lookup <- setNames(rep(names(celltype_groups), lengths(celltype_groups)),
                            unlist(celltype_groups))

# Plot order + colours
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

# Output folder
output_dir <- file.path(results_dir, "publication")
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}

SC_processed_fn <- file.path(output_dir, "AndrewLJi_seurat_processed.rds")

#-------------------------------------------------------------------------------

# 3. Final model
#-------------------------------------------------------------------------------
sccoregep_fn <- file.path(results_dir, "intermediate", "sccore_gep", "Discovery_models", "coxnet_Late_feature_integration_robj.rds") 

GEP_final_model <- readRDS(sccoregep_fn)
coeff_genes <- GEP_final_model$Final_model2$coefficients
stopifnot(length(coeff_genes) == 23)

# Ensembl ID -> gene symbol (strip version suffix from the GTF IDs)
gene_map <- import(p9_gtf) %>%
  as.data.frame() %>%
  distinct(gene_id, gene_name) %>%
  mutate(gene_id_nover = sub("\\..*$", "", gene_id))

names(coeff_genes) <- gene_map$gene_name[match(names(coeff_genes), gene_map$gene_id_nover)]
stopifnot(!anyNA(names(coeff_genes)))

# Rename to single-cell symbols (matched by name, not by position)
to_rename <- names(coeff_genes) %in% names(sc_gene_rename)
names(coeff_genes)[to_rename] <- sc_gene_rename[names(coeff_genes)[to_rename]]

coeff_genes <- sort(coeff_genes)   # gene order follows coefficients (aligned with heatmap)


#-------------------------------------------------------------------------------

# 4. Preprocess single cell data
#-------------------------------------------------------------------------------
if (!file.exists(SC_processed_fn)) {
  counts   <- read.table(SC_counts_fn,   header = TRUE, sep = "\t")
  metadata <- read.table(SC_metadata_fn, header = TRUE, sep = "\t")
  
  seurat_obj_Andrew <- CreateSeuratObject(counts = counts, meta.data = metadata)
  #seurat_obj_Andrew <- PercentageFeatureSet(seurat_obj_Andrew, pattern = "^RP[SL]", col.name = "percent_ribo")
  #seurat_obj_Andrew <- PercentageFeatureSet(seurat_obj_Andrew, pattern = "^MT.",    col.name = "percent_mito")
  
  seurat_obj_Andrew <- default_seurat_processing_PCA(seurat_obj_Andrew)
  seurat_obj_Andrew <- default_seurat_UMAP(seurat_obj_Andrew, dims = 15)  # based on PCA heatmaps and elbow plot
  
  saveRDS(seurat_obj_Andrew, SC_processed_fn)
  rm(counts, metadata)
} else {
  seurat_obj_Andrew <- readRDS(SC_processed_fn)
}

#-------------------------------------------------------------------------------

# 5. Annotation
#-------------------------------------------------------------------------------
seurat_obj_Andrew <- subset(seurat_obj_Andrew,
                            subset = !(level2_celltype %in% c("Multiplet", "Melanocyte")))

meta <- seurat_obj_Andrew@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(cell,
            patient,
            custom_id = unname(celltype_lookup[paste0(level2_celltype, ":", tum.norm)]),
            custom_id = factor(custom_id, levels = setdiff(celltype_levels, "NA")))

# every cell must have received a new annotation
stopifnot(!anyNA(meta$custom_id))

#-------------------------------------------------------------------------------

# 6. Per patient contribution
#-------------------------------------------------------------------------------
all_RNA <- LayerData(seurat_obj_Andrew, assay = "RNA", layer = "data")

genes_in_sc <- intersect(names(coeff_genes), rownames(all_RNA))

# subset to the model genes *before* densifying (avoids a full dense matrix)
expr_long <- as.matrix(all_RNA[genes_in_sc, meta$cell, drop = FALSE]) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  pivot_longer(-cell, names_to = "gene", values_to = "value") %>%
  inner_join(meta, by = "cell")

counts_per_patient <- expr_long %>%
  group_by(patient, custom_id, gene) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(patient, gene) %>%
  mutate(frac_within_patient = value / sum(value)) %>%
  group_by(gene) %>%
  mutate(frac_of_total = value / sum(value)) %>%
  ungroup() %>%
  mutate(gene = factor(gene, levels = rev(names(coeff_genes))))

# Not enough reads in LCE1F and DEFB103A (see plot at the end)
genes_na       <- c("LCE1F", "DEFB103A")
genes_missing  <- union(genes_na, setdiff(names(coeff_genes), genes_in_sc))
genes_detected <- setdiff(names(coeff_genes), genes_missing)
message("Model genes shown as NA: ", paste(genes_missing, collapse = ", "))

# Pooled contributions 
to_plot <- counts_per_patient %>%
  filter(gene %in% genes_detected) %>%
  group_by(custom_id, gene) %>%
  summarise(frac_counts_total = sum(frac_of_total), .groups = "drop") %>%
  mutate(across(c(custom_id, gene), as.character)) %>%
  bind_rows(tibble(custom_id = "NA", gene = genes_missing, frac_counts_total = 1)) %>%
  mutate(custom_id = factor(custom_id, levels = celltype_levels),
         gene      = factor(gene,      levels = rev(names(coeff_genes))))

# Plot
p_pooled <- ggplot(to_plot, aes(x = gene, y = frac_counts_total, fill = custom_id)) +
  geom_col() +
  scale_fill_manual(values = celltype_colors) +
  coord_flip() +
  theme_bw() +
  ylab("Fraction counts (%)")

pdf(file.path(output_dir, "supp_figure_11A_SC_Andrew_barplot_contributions_counts.pdf"), width = 6, height = 5)
print(p_pooled)
dev.off()

#-------------------------------------------------------------------------------

# 7. Save source data
#-------------------------------------------------------------------------------
writexl::write_xlsx(to_plot, paste0( output_dir, "/supp_figure_11A_SC_Andrew_barplot_contributions_counts.xlsx"))

# Not run: total expression per gene across all cells (used to check that
# LCE1F and DEFB103A have too few reads)
# gene_totals <- total_counts_all %>%
#   group_by(variable) %>%
#   summarise(
#     total_reads = sum(value, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   arrange(total_reads)
# ggplot(gene_totals, aes(x = reorder(variable, total_reads), y = total_reads)) +
#   geom_col() +
#   coord_flip() +
#   theme_bw() +
#   labs(
#     x = "Gene",
#     y = "Total expression",
#     title = "Total expression across all cells"
#   )
