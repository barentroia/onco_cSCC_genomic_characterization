# Load libraries
#-------------------------------------------------------------------------------
library(conflicted)
library(argparser)
library(tidyverse)
library(rtracklayer)
library(DESeq2)
#library(BiocParallel)
#-------------------------------------------------------------------------------

# Functions
#-------------------------------------------------------------------------------
dgea_deseq2 <- function(rnaseq_obj, design, output_folder,contrast_design=c("Metastasis","Case","Control"),subset="all"){
  #' DGEA with DESeq2
  #' 
  #' @param rnaseq_obj object. Object with RNA sequencing data
  #' @param design formula. Formula of the design for DESeq2
  #' @param output_folder character. Output folder
  #' @param contrast_design list. List with the contrasts
  #' @param subset character. Variable in the clinical data file on which to subset the data
  #'
  #' @return DEG list + corresponding DESeq2 objects
  #'  
  
  # Function
  #-----------------------------------------------------------------------------
  map_fun <- function(id, map_file, id_col = "gene_id", map_id_col = "gene_name"){
    map_file[map_file[, id_col] == id, map_id_col]
  }
  #-----------------------------------------------------------------------------
  
  
  # Prepare data
  #-----------------------------------------------------------------------------
  design_formula <- as.formula(design)
  
  rnaseq_obj@clinical <-  rnaseq_obj@clinical %>% mutate(Metastasis = factor(Metastasis, levels = c("Control", "Case")),
                                                         Set_id = as.factor(Set_id),
                                                         Biopsy_excision = factor(ifelse(Biopsy_excision %in% c("Curettage", "Excochleation"),
                                                                                         "Excision", Biopsy_excision),
                                                                                  levels = c("Excision", "Biopsy")),
                                                         Immunosupp_bin = factor(ifelse(OTR_at_cSCC == "Yes" | HM_at_cSCC == "Yes", "Yes","No"),
                                                                                 levels = c("Yes", "No")),
                                                         BWH = factor(BWH, levels = c("T1", "T2a", "T2b", "T3")),
                                                         AJCC_8 = factor(AJCC_8, levels = c("T1", "T2", "T3", "T4")))
  clinical <- rnaseq_obj@clinical
  
  #-----------------------------------------------------------------------------
  
  # Subset data, <colname>:<factor>
  #-----------------------------------------------------------------------------
  if (subset != "all"){
    subset_names <- strsplit(subset,":")
    if (length(grep("&",subset_names[[1]][2])) > 0){
      subset_names_vars <- strsplit(subset_names[[1]][2],"&")[[1]]
      clinical <- clinical[clinical[[subset_names[[1]][1]]] %in% subset_names_vars,]
    } else {
      clinical <- clinical[clinical[[subset_names[[1]][1]]] == subset_names[[1]][2],]
    }
    rnaseq_obj <- subset_samples(rnaseq_obj,
                                 samples = clinical$Skyline_ID,
                                 sample_id_column = "Skyline_ID")
    
    output_folder <- paste0(output_folder,"_",subset_names[[1]][1],"_",subset_names[[1]][2])
  }
  #-----------------------------------------------------------------------------
  
  
  # Process data
  #-----------------------------------------------------------------------------
  # Variables in design
  vars_in_design <- base::setdiff(unlist(strsplit(as.character(design), "~|[+]")), "")
  # Check if there are NAs for variables in design
  nas_in_vars <- sapply(vars_in_design,
                        function(x, df) any(is.na(df[[x]])),
                        df = clinical,
                        simplify = F) %>%
    any()
  print(paste0("There are NAs for variables in the design: ", nas_in_vars))
  # Remove samples with NAs
  clinical <- clinical %>%
    drop_na(any_of(vars_in_design))
  # Samples to keep
  samples_to_keep <- clinical$Skyline_ID
  if (grepl("Set_id", design)){
    samples_to_keep <- clinical %>%
      dplyr::filter(Skyline_ID %in% samples_to_keep) %>%
      arrange(Set_id) %>%
      group_by(Set_id) %>%
      mutate(n_samples_in_set = length(Skyline_ID)) %>%
      dplyr::filter(n_samples_in_set == 2) %>%
      pull(Skyline_ID)
  }
  print(paste0("Number of samples: ", length(samples_to_keep)))
  
  if (length(setdiff(samples_to_keep,rnaseq_obj@clinical$Skyline_ID))>0){
    rnaseq_obj <- subset_samples(rnaseq_obj,
                                 samples = samples_to_keep,
                                 sample_id_column = "Skyline_ID")
  }
  
  # Number of cases and controls
  print(paste0("Number of cases: ", clinical %>% dplyr::filter(Metastasis == "Case") %>% nrow()))
  print(paste0("Number of controls: ", clinical %>% dplyr::filter(Metastasis == "Control") %>% nrow()))
  
  # Design related
  print("Design related sample numbers:")
  print(table(clinical[[contrast_design[1]]]))
  
  # Sequencing data related
  print(paste0("Dimensions counts: ",paste0(dim(rnaseq_obj@sequencing_data$counts_filt),collapse = ", ")))
  
  #-----------------------------------------------------------------------------
  
  # DGEA
  #-----------------------------------------------------------------------------
  # Run DGEA
  dds <- DESeqDataSetFromMatrix(countData = round(rnaseq_obj@sequencing_data$counts_filt),
                                colData = rnaseq_obj@clinical,
                                design = design_formula)
  dds <- DESeq(dds)
  #-----------------------------------------------------------------------------
  
  # DGEA results
  #-----------------------------------------------------------------------------
  # Extract results with FDR thresh = alpha
  alpha <- 0.05
  res <- results(dds, alpha = alpha, contrast=contrast_design)
  res_df <- res %>% 
    as.data.frame() %>%
    mutate(ensg_id = rownames(.),
           gene_name = sapply(ensg_id, map_fun, out),
           padj_string = ifelse(padj < alpha & abs(log2FoldChange) > 2,
                                paste0("Adj p-value < ", alpha, " & |log2(FC)| > 2"),
                                ifelse(padj < alpha, paste0("Adj p-value < ", alpha),
                                       "Not significant")),
           neg_log10_padj = -log10(padj),
           neg_log10_pvalue = -log10(pvalue))
  #-----------------------------------------------------------------------------
  
  # Save results
  #-----------------------------------------------------------------------------
  if(!dir.exists(output_folder)) dir.create(paste0(output_folder,"/"))
  # Save DGEA results
  write.csv(res_df, paste0(output_folder,"/DGEA_results_DESeq2.csv"), row.names = F)
  # Save DGEA results object
  save(res, file = paste0(output_folder,"/DGEA_results_DESeq2.RData"))
  # Save dds object
  save(dds, file =  paste0(output_folder,"/dds_DESeq2.RData"))
  
  return(list(res,dds))
  #-----------------------------------------------------------------------------
}
#-------------------------------------------------------------------------------


