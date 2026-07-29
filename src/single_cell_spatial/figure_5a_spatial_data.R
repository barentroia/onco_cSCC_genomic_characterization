#-------------------------------------------------------------------------------
# Aim: Generate Figure 5a with single cell data
# Input: Spatial sequencing data
# Output: Barplot with the fraction of reads expressed per cell type
#-------------------------------------------------------------------------------


# Load libraries
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

#-------------------------------------------------------------------------------


# General setup
#-------------------------------------------------------------------------------
# output folder
output_dir <- file.path(results_dir, "publication")
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}
#-------------------------------------------------------------------------------


# Final model
# TODO replace with final model in results section!
#-------------------------------------------------------------------------------
final_model <- file.path(project_dir, "rruiter", "combi", "step-ident-2023-genesig-development", "output", "model_dev","E00_combiGEP_coxnet_test",
                         "GEPBailey_remaining_record_based_removed_badQCsamples_coxnet_deseq2_Biopsy_excision_2f+Sex+Metastasis_num_50logfc_logfc0.5_robj.rds")
GEP_final_model <- readRDS(final_model)
final_model_genes <- names(GEP_final_model$Final_model2$coefficients)
stopifnot(length(final_model_genes)==23)
#-------------------------------------------------------------------------------


# Quick checks ST data
#-------------------------------------------------------------------------------
ST_CSCC_seurat <- LoadSeuratRds(seurat_object_fn)

# no spatial data in there
# only limited to the genes of interest, I can not do any regular checks

#-------------------------------------------------------------------------------


# Generate cell type plot
#-------------------------------------------------------------------------------
# mapping file bulk RNAseq
mapping <- import(p9_gtf)
out <- mapping %>%
  as.data.frame() %>%
  dplyr::select(gene_id, gene_name) %>%
  unique()

final_model_gene_names <- out[match(final_model_genes,do.call(rbind,strsplit(out$gene_id,"[.]"))[,1]),]$gene_name
stopifnot(length(unlist(final_model_gene_names))==23)
# LCE1F, DEFB103A not in ST data

# custom id
ST_CSCC_seurat$custom_id <- ST_CSCC_seurat$annotate
# immune cells merged
ST_CSCC_seurat$custom_id <- ifelse(ST_CSCC_seurat$custom_id %in% c("large_pre_B_cell","macrophage","neutrophils","plasma_cell","T_cell"),"immune_cell",as.character(ST_CSCC_seurat$custom_id))
# pilosebaceous cells merged
ST_CSCC_seurat$custom_id <- ifelse(ST_CSCC_seurat$custom_id %in% c("sebaceous_gland","hair_follicle"),"pilosebaceous_cells",as.character(ST_CSCC_seurat$custom_id))
# undefined stronmal cells merged
ST_CSCC_seurat$custom_id <- ifelse(ST_CSCC_seurat$custom_id %in% c("stormal_cell"),"fibroblast",as.character(ST_CSCC_seurat$custom_id))


# new barplot figure
ST_CSCC_seurat$custom_id <- as.factor(ST_CSCC_seurat$custom_id)
ST_CSCC_seurat$custom_id <- factor(ST_CSCC_seurat$custom_id, levels=rev(c("immune_cell", "keratinocytes","endothelial_cell",
                                                                          "fibroblast","smooth_muscle","tumor_keratinocytes",
                                                                          "tumor_endothelial cell","tumor_fibroblast",
                                                                          "sweat_gland","pilosebaceous_cells")))
Idents(ST_CSCC_seurat) <- "custom_id"
ST_CSCC_seurat$patient_id <- colnames(ST_CSCC_seurat)
ST_CSCC_seurat$patient_id <- substr(ST_CSCC_seurat$patient_id,1,5)

all_RNA <- LayerData(ST_CSCC_seurat,assay = "RNA", layer = "data")

# proportion of reads
# pat 1
meta_data_all <- ST_CSCC_seurat@meta.data %>% filter(patient_id == "SCC01")
RNA_data_all_final <- all_RNA[,rownames(meta_data_all)]
RNA_data_all_final <- RNA_data_all_final[rownames(RNA_data_all_final) %in% final_model_gene_names,]
scale_data_all_final <- t(RNA_data_all_final)
scale_data_all_final_meta <- merge(meta_data_all,scale_data_all_final,by=0)
scale_data_all_final_meta_sum <- scale_data_all_final_meta %>% group_by(custom_id) %>% summarise(across(c("PLA2G2A","PLOD2","CCN5","SPP1","KRT8",
                                                                                                          "KRT19","KRT18","TGM2","MMP13","IGF2BP3",
                                                                                                          "MMP11","INHBA","SCNN1B","WFDC5",
                                                                                                          "CYSRT1","EPGN","ALOX15B","USP2","RDH12",
                                                                                                          "TMEM86A","CD36"), ~ sum(.x, na.rm = TRUE)))
scale_data_all_final_meta_melt <- melt(scale_data_all_final_meta_sum)

total_counts_p1 <- scale_data_all_final_meta_melt %>% group_by(variable) %>% mutate(sum_gene = sum(value), frac_counts = value/sum_gene)

# pat 2
meta_data_all <- ST_CSCC_seurat@meta.data %>% filter(patient_id == "SCC02")
RNA_data_all_final <- all_RNA[,rownames(meta_data_all)]
RNA_data_all_final <- RNA_data_all_final[rownames(RNA_data_all_final) %in% final_model_gene_names,]
scale_data_all_final <- t(RNA_data_all_final)
scale_data_all_final_meta <- merge(meta_data_all,scale_data_all_final,by=0)
scale_data_all_final_meta_sum <- scale_data_all_final_meta %>% group_by(custom_id) %>% summarise(across(c("PLA2G2A","PLOD2","CCN5","SPP1","KRT8",
                                                                                                          "KRT19","KRT18","TGM2","MMP13","IGF2BP3",
                                                                                                          "MMP11","INHBA","SCNN1B","WFDC5",
                                                                                                          "CYSRT1","EPGN","ALOX15B","USP2","RDH12",
                                                                                                          "TMEM86A","CD36"), ~ sum(.x, na.rm = TRUE)))
scale_data_all_final_meta_melt <- melt(scale_data_all_final_meta_sum)

total_counts_p2 <- scale_data_all_final_meta_melt %>% group_by(variable) %>% mutate(sum_gene = sum(value), frac_counts = value/sum_gene)


# pat 3
meta_data_all <- ST_CSCC_seurat@meta.data %>% filter(patient_id == "SCC03")
RNA_data_all_final <- all_RNA[,rownames(meta_data_all)]
RNA_data_all_final <- RNA_data_all_final[rownames(RNA_data_all_final) %in% final_model_gene_names,]
scale_data_all_final <- t(RNA_data_all_final)
scale_data_all_final_meta <- merge(meta_data_all,scale_data_all_final,by=0)
scale_data_all_final_meta_sum <- scale_data_all_final_meta %>% group_by(custom_id) %>% summarise(across(c("PLA2G2A","PLOD2","CCN5","SPP1","KRT8",
                                                                                                          "KRT19","KRT18","TGM2","MMP13","IGF2BP3",
                                                                                                          "MMP11","INHBA","SCNN1B","WFDC5",
                                                                                                          "CYSRT1","EPGN","ALOX15B","USP2","RDH12",
                                                                                                          "TMEM86A","CD36"), ~ sum(.x, na.rm = TRUE)))
scale_data_all_final_meta_melt <- melt(scale_data_all_final_meta_sum)

total_counts_p3 <- scale_data_all_final_meta_melt %>% group_by(variable) %>% mutate(sum_gene = sum(value), frac_counts = value/sum_gene)


# pat 4
meta_data_all <- ST_CSCC_seurat@meta.data %>% filter(patient_id == "SCC04")
RNA_data_all_final <- all_RNA[,rownames(meta_data_all)]
RNA_data_all_final <- RNA_data_all_final[rownames(RNA_data_all_final) %in% final_model_gene_names,]
scale_data_all_final <- t(RNA_data_all_final)
scale_data_all_final_meta <- merge(meta_data_all,scale_data_all_final,by=0)
scale_data_all_final_meta_sum <- scale_data_all_final_meta %>% group_by(custom_id) %>% summarise(across(c("PLA2G2A","PLOD2","CCN5","SPP1","KRT8",
                                                                                                          "KRT19","KRT18","TGM2","MMP13","IGF2BP3",
                                                                                                          "MMP11","INHBA","SCNN1B","WFDC5",
                                                                                                          "CYSRT1","EPGN","ALOX15B","USP2","RDH12",
                                                                                                          "TMEM86A","CD36"), ~ sum(.x, na.rm = TRUE)))
scale_data_all_final_meta_melt <- melt(scale_data_all_final_meta_sum)

total_counts_p4 <- scale_data_all_final_meta_melt %>% group_by(variable) %>% mutate(sum_gene = sum(value), frac_counts = value/sum_gene)



# pat 5
meta_data_all <- ST_CSCC_seurat@meta.data %>% filter(patient_id == "SCC05")
RNA_data_all_final <- all_RNA[,rownames(meta_data_all)]
RNA_data_all_final <- RNA_data_all_final[rownames(RNA_data_all_final) %in% final_model_gene_names,]
scale_data_all_final <- t(RNA_data_all_final)
scale_data_all_final_meta <- merge(meta_data_all,scale_data_all_final,by=0)
scale_data_all_final_meta_sum <- scale_data_all_final_meta %>% group_by(custom_id) %>% summarise(across(c("PLA2G2A","PLOD2","CCN5","SPP1","KRT8",
                                                                                                          "KRT19","KRT18","TGM2","MMP13","IGF2BP3",
                                                                                                          "MMP11","INHBA","SCNN1B","WFDC5",
                                                                                                          "CYSRT1","EPGN","ALOX15B","USP2","RDH12",
                                                                                                          "TMEM86A","CD36"), ~ sum(.x, na.rm = TRUE)))
scale_data_all_final_meta_melt <- melt(scale_data_all_final_meta_sum)

total_counts_p5 <- scale_data_all_final_meta_melt %>% group_by(variable) %>% mutate(sum_gene = sum(value), frac_counts = value/sum_gene)

# pat 6
meta_data_all <- ST_CSCC_seurat@meta.data %>% filter(patient_id == "SCC06")
RNA_data_all_final <- all_RNA[,rownames(meta_data_all)]
RNA_data_all_final <- RNA_data_all_final[rownames(RNA_data_all_final) %in% final_model_gene_names,]
scale_data_all_final <- t(RNA_data_all_final)
scale_data_all_final_meta <- merge(meta_data_all,scale_data_all_final,by=0)
scale_data_all_final_meta_sum <- scale_data_all_final_meta %>% group_by(custom_id) %>% summarise(across(c("PLA2G2A","PLOD2","CCN5","SPP1","KRT8",
                                                                                                          "KRT19","KRT18","TGM2","MMP13","IGF2BP3",
                                                                                                          "MMP11","INHBA","SCNN1B","WFDC5",
                                                                                                          "CYSRT1","EPGN","ALOX15B","USP2","RDH12",
                                                                                                          "TMEM86A","CD36"), ~ sum(.x, na.rm = TRUE)))
scale_data_all_final_meta_melt <- melt(scale_data_all_final_meta_sum)

total_counts_p6 <- scale_data_all_final_meta_melt %>% group_by(variable) %>% mutate(sum_gene = sum(value), frac_counts = value/sum_gene)


# pat 7
meta_data_all <- ST_CSCC_seurat@meta.data %>% filter(patient_id == "SCC07")
RNA_data_all_final <- all_RNA[,rownames(meta_data_all)]
RNA_data_all_final <- RNA_data_all_final[rownames(RNA_data_all_final) %in% final_model_gene_names,]
scale_data_all_final <- t(RNA_data_all_final)
scale_data_all_final_meta <- merge(meta_data_all,scale_data_all_final,by=0)
scale_data_all_final_meta_sum <- scale_data_all_final_meta %>% group_by(custom_id) %>% summarise(across(c("PLA2G2A","PLOD2","CCN5","SPP1","KRT8",
                                                                                                          "KRT19","KRT18","TGM2","MMP13","IGF2BP3",
                                                                                                          "MMP11","INHBA","SCNN1B","WFDC5",
                                                                                                          "CYSRT1","EPGN","ALOX15B","USP2","RDH12",
                                                                                                          "TMEM86A","CD36"), ~ sum(.x, na.rm = TRUE)))
scale_data_all_final_meta_melt <- melt(scale_data_all_final_meta_sum)

total_counts_p7 <- scale_data_all_final_meta_melt %>% group_by(variable) %>% mutate(sum_gene = sum(value), frac_counts = value/sum_gene)


# pat 8
meta_data_all <- ST_CSCC_seurat@meta.data %>% filter(patient_id == "SCC08")
RNA_data_all_final <- all_RNA[,rownames(meta_data_all)]
RNA_data_all_final <- RNA_data_all_final[rownames(RNA_data_all_final) %in% final_model_gene_names,]
scale_data_all_final <- t(RNA_data_all_final)
scale_data_all_final_meta <- merge(meta_data_all,scale_data_all_final,by=0)
scale_data_all_final_meta_sum <- scale_data_all_final_meta %>% group_by(custom_id) %>% summarise(across(c("PLA2G2A","PLOD2","CCN5","SPP1","KRT8",
                                                                                                          "KRT19","KRT18","TGM2","MMP13","IGF2BP3",
                                                                                                          "MMP11","INHBA","SCNN1B","WFDC5",
                                                                                                          "CYSRT1","EPGN","ALOX15B","USP2","RDH12",
                                                                                                          "TMEM86A","CD36"), ~ sum(.x, na.rm = TRUE)))
scale_data_all_final_meta_melt <- melt(scale_data_all_final_meta_sum)

total_counts_p8 <- scale_data_all_final_meta_melt %>% group_by(variable) %>% mutate(sum_gene = sum(value), frac_counts = value/sum_gene)


total_counts_all <- rbind(total_counts_p1,total_counts_p2,total_counts_p3,total_counts_p4,
                          total_counts_p5,total_counts_p6,total_counts_p7,total_counts_p8)
total_gene_sum <- total_counts_all %>% group_by(variable) %>% summarise(sum_gene = sum(value))
total_counts_all$sum_gene <- total_gene_sum[match(total_counts_all$variable,total_gene_sum$variable),]$sum_gene
total_counts_all <- total_counts_all %>% group_by(variable) %>% mutate(frac_counts_total = value/sum_gene)

total_counts_all$custom_id <- factor(total_counts_all$custom_id, levels=c("immune_cell", "keratinocytes","endothelial_cell",
                                                                          "fibroblast","smooth_muscle","tumor_keratinocytes",
                                                                          "tumor_fibroblast",
                                                                          "sweat_gland","pilosebaceous_cells"))

coeff_genes <- GEP_final_model$Final_model2$coefficients
names(coeff_genes) <- out[match(names(coeff_genes),do.call(rbind,strsplit(out$gene_id,"[.]"))[,1]),]$gene_name
coeff_genes <- coeff_genes[order(coeff_genes)]

total_counts_all$variable <- factor(total_counts_all$variable, levels=names(coeff_genes))

to_plot <- total_counts_all %>% group_by(custom_id, variable) %>% summarise(frac_counts_total = sum(frac_counts_total))

# align with heatmap
setdiff(names(coeff_genes),as.character(unique(to_plot$variable)))
to_plot_extra  <- rbind(to_plot,as.data.frame(list(custom_id="NA",variable = "LCE1F", frac_counts_total = 1)))
to_plot_extra  <- rbind(to_plot_extra, as.data.frame(list(custom_id="NA",variable = "DEFB103A", frac_counts_total = 1)))
to_plot_extra$variable <- factor(to_plot_extra$variable,levels=rev(names(coeff_genes)))
to_plot_extra$custom_id <- factor(to_plot_extra$custom_id, levels= c(as.character(unique(to_plot$custom_id)),"NA"))

to_plot_extra$custom_id <- factor(to_plot_extra$custom_id, 
                                  levels= c("immune_cell", "sweat_gland","pilosebaceous_cells",
                                            "keratinocytes","endothelial_cell","fibroblast","smooth_muscle","tumor_keratinocytes",
                                            "tumor_fibroblast","NA"))

pdf(file=paste0(output_data_dir,"/figure_5a_ST_barplot_contributions_counts_immune_N8_flip_highres_v2.pdf"), width=7, height=6)
ggplot(to_plot_extra) + geom_bar(aes(x=variable , y=frac_counts_total, fill= custom_id),stat="identity") +
  scale_fill_manual(values=c("#A53C73FF", "#823575FF", "#B46727FF", "#C2CBCDFF", "#A7A7A7FF", "#7E7568FF", "#A4D8E5FF" ,"#2FBCE3FF" ,"#1F6582FF" ,"white")) + theme_bw() +
  coord_flip()
dev.off()

writexl::write_xlsx(to_plot_extra,paste0(output_data_dir,"/figure_5a_barplot_ST_N8_source_data_highres_v2.xlsx"))



