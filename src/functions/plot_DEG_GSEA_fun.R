
plotsSignatures <- function(dgea_df,
                            design,
                            signature=NA,
                            n_samples=378,
                            figure_name){
  
  #' Plot volcano plots with DEGs + signature genes
  #' 
  #' @param dgea_df data frame. Data frame with the DEG output from DESeq2
  #' @param design character. String with information on the design
  #' @param signature list. List of signatures (gene names)
  #' @param n_samples numeric. Number of samples used in DGEA (default is 378)
  #' @param figure_name character. Name of the figure
  #' 
  #' @return Volcano plots
  #'  
  
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


runGSEA <- function(dgea_df,
                    pathways_to_test,
                    design,
                    type_set="H",
                    n_top=15,
                    figure_name){
  
  #' Run and plot GSEA using fgsea
  #' 
  #' @param dgea_df data frame. Data frame with the DEG output from DESeq2
  #' @param pathways_to_test list. List(s) of pathways to test
  #' @param design character. Design model DGEA
  #' @param type_set character. Type of pathways
  #' 
  #' @return GSEA plot
  #'  
  
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
