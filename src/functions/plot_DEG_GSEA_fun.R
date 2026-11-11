#===============================================================================
# Aim:      Plot differentially expressed genes (volcano plots) and run/plot
#           gene set enrichment analysis (GSEA).
# Contents: - plotsSignatures(): volcano plot, optionally with signature genes
#           - runGSEA(): fgsea on DESeq2 statistics + bar plot of top pathways
# Requires: output_dir (and color_fig when a signature is plotted) in the
#           calling environment
# Usage:    Sourced by the DEG and GSEA figure scripts in src/rna_seq/figures
#===============================================================================

#' Volcano plot of DEGs, optionally highlighting signature genes
#'
#' Without a signature, significant genes are labelled; with a signature,
#' signature genes are coloured using color_fig (global object).
#'
#' @param dgea_df data.frame. DGEA results from dgea_deseq2() (padj,
#'   log2FoldChange, padj_string, gene_name and, if used, signature).
#' @param design character. Design/contrast label used in title and file names.
#' @param signature character vector or NA. Signature gene names; NA = none.
#' @param n_samples numeric. Number of samples in the DGEA (for the title).
#' @param figure_name character. Prefix of the output file names.
#'
#' @return Called for its side effects: saves the volcano plot (pdf) and the
#'   plotted data (csv) in output_dir (global object).
plotsSignatures <- function(dgea_df,
                            design,
                            signature=NA,
                            n_samples=378,
                            figure_name){

  
  if (any(is.na(signature))){
    dgea_df <- dgea_df[order(dgea_df$padj,decreasing = F),]
    dgea_df <- dgea_df[!is.na(dgea_df$padj_string),]
    pdf(file=paste0(output_dir,"/",figure_name,"_DEG_",design,".pdf"), width=6, height=4.2)
    p1=dgea_df %>% ggplot() + 
      geom_point(aes(y=-log10(padj),x=log2FoldChange,color=padj_string),alpha=0.5) +
      scale_color_manual(values=list("Adj p-value < 0.05"="black","Not significant"="grey")) +
      ggtitle(paste0(design," (N=",n_samples,", #DEG=",nrow(dgea_df),")")) +
      theme_bw() + 
      xlab("fold change (log2)") +
      geom_vline(xintercept = 0) +   geom_hline(yintercept = -log10(0.05),linetype=2) +
      geom_label_repel(data=dgea_df[1:100,], box.padding = 0.3,min.segment.length = 0,
                       aes(y=-log10(padj),x=log2FoldChange,label=gene_name),size=2,
                       max.overlaps = getOption("ggrepel.max.overlaps", default = 15))
    print(p1)
    dev.off()
    
    write.csv(x=dgea_df,file=paste0(output_dir,"/",figure_name,"_DEG_",design,".csv"),
              quote=F, row.names = F)
  } else {
    
    if (is.list(signature)){
      if(length(signature)==2){
        color_fig <- c("#BE6E1BFF","#7C396B")
      }    
      if(length(signature)>2){
        color_fig <- skyColors("#A03A6FFF", "#D58000FF", "#C7C7C7FF", "#14B5E2FF", "#266986FF","darkgreen","#7C396B") 
      }
      
    } else {
      color_fig <- c("#BE6E1BFF")
    }
    
    dgea_df$signature <- ""
    for (sig in names(signature)){
      dgea_df[dgea_df$gene_name %in% signature[[sig]],]$signature <- sig
      print(paste0(sig,", # genes present in data: "))
      print(nrow(dgea_df[dgea_df$signature == sig,]))
    }
    
    dgea_df <- dgea_df[order(dgea_df$padj,decreasing = F),]
    dgea_df <- dgea_df[!is.na(dgea_df$padj_string),]
    pdf(file=paste0(output_dir,"/",figure_name,"_DEG_signature",design,".pdf"), width=6, height=4.2)
    p1=ggplot() + 
      geom_point(data=dgea_df[dgea_df$signature == "",],aes(y=-log10(padj),x=log2FoldChange,color=padj_string),alpha=0.5) +
      scale_color_manual(values=c("black","black")) +
      ggtitle(paste0(design," (N=",n_samples,", #DEG=",nrow(dgea_df),")")) +
      theme_bw() + 
      xlab("fold change (log2)") +
      new_scale_color() +
      geom_point(data=dgea_df[dgea_df$signature != "",],aes(y=-log10(padj),x=log2FoldChange,color=signature),alpha=1) +
      scale_color_manual(values=color_fig) +
      geom_vline(xintercept = 0) +   geom_hline(yintercept = -log10(0.05),linetype=2) +
      ggrepel::geom_label_repel(data=dgea_df[dgea_df$padj_string != "Not significant",][1:500,], box.padding = 0.3,min.segment.length = 0,
                                aes(y=-log10(padj),x=log2FoldChange,label=gene_name),size=2,
                                max.overlaps = getOption("ggrepel.max.overlaps", default = 15))
    print(p1)
    dev.off()
    
    write.csv(x=dgea_df,file=paste0(output_dir,"/",figure_name,"_DEG_signature",design,".csv"),
              quote=F, row.names = F)
  }
  
  
}

#' Run gene set enrichment analysis (fgsea) and plot the top pathways
#'
#' Genes are ranked by the DESeq2 Wald statistic. The n_top pathways with the
#' lowest and the highest normalized enrichment score (NES) are plotted.
#'
#' @param dgea_df data.frame. DGEA results from dgea_deseq2() (ensg_id, stat).
#' @param pathways_to_test list. Gene sets (Ensembl gene IDs) to test.
#' @param design character. Design/contrast label used in title and file names.
#' @param type_set character. Gene set collection label (e.g. "H" = Hallmarks).
#' @param n_top integer. Number of top pathways per direction to plot.
#' @param figure_name character. Prefix of the output file names.
#'
#' @return Called for its side effects: saves the GSEA bar plot (pdf) and all
#'   fgsea results (csv) in output_dir (global object).
runGSEA <- function(dgea_df,
                    pathways_to_test,
                    design,
                    type_set="H",
                    n_top=15,
                    figure_name){
  
  set.seed(123)
  
  dgea_df$ensembl_gene_id <- do.call(rbind,strsplit(dgea_df$ensg_id,"[.]"))[,1]
  
  res_cleaned <- dgea_df %>% 
    dplyr::select(ensembl_gene_id, stat) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(ensembl_gene_id) %>% 
    summarize(stat=mean(stat))
  
  ranks <- deframe(res_cleaned)
  
  fgseaRes <- fgseaMultilevel(pathways=pathways_to_test, stats=ranks)
  
  fgseaRes <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  fgseaRes_sig <- fgseaRes[order(fgseaRes$NES,decreasing = F),][1:n_top,]
  fgseaRes_sig <- rbind(fgseaRes_sig,fgseaRes[order(fgseaRes$NES,decreasing = T),][1:n_top,])
  width_fig <- 4 + max(unlist(lapply(fgseaRes_sig$pathway, nchar)))*0.07
  pdf(file=paste0(output_dir,"/",figure_name,"_fGSEA_",design,"_",type_set,".pdf"), width=width_fig, height=6.5+((n_top-15)*0.2))
  p1=ggplot(fgseaRes_sig, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=padj<0.05)) + 
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=design) + coord_flip() +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_fill_manual(values=list("FALSE"="grey","TRUE"="#A03A6FFF")) 
  print(p1)
  dev.off()
  
  write.csv(x=fgseaRes %>% select(-leadingEdge),file=paste0(output_dir,"/",figure_name,"_fGSEA_",design,"_",type_set,".csv"),
            quote=F, row.names = F)
}
