#===============================================================================
# Aim:      Load the RNA-seq datasets (TPMs, counts, clinical data, gene
#           annotation) into S4 data objects and save them as .rds files.
# Contents: - DiscoveryObject / ValidationObject: S4 classes and constructors
#           - subset_samples(): subset an object to a set of samples
#           - load_discovery(): discovery dataset, bad QC sets removed (N = 366)
#           - load_discovery_all(): discovery dataset, bad QC samples removed
#             (N = 378)
#           - load_validation_d_squame(): D-SQUAME validation dataset
#           - load_validation_nassir(): Nassir et al. dataset
# Usage:    Sourced from main.R; the .rds files are read by the RNA-seq,
#           heatmap and DGEA scripts.
#===============================================================================

# 1. Required libraries
#-------------------------------------------------------------------------------
library(rtracklayer)
library(tidyverse)
library(conflicted)

conflicts_prefer(base::intersect)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::rename)
#-------------------------------------------------------------------------------

# 2. Discovery data object
#-------------------------------------------------------------------------------
setClass(
  "DiscoveryObject",
  slots = list(
    sequencing_data = "list",
    clinical = "data.frame",
    out = "data.frame"
  )
)

#' Use the ID column as rownames (discovery data)
#'
#' @param df data.frame. Expression data with an ID column (gene IDs).
#'
#' @return data.frame with gene IDs as rownames and without the ID column.
.fix_gene_ids <- function(df) {
  
  if (!"ID" %in% colnames(df)) {
    stop("ID column is required.")
  }
  
  # Case 1: rownames are default
  default_rownames <- identical(rownames(df), as.character(seq_len(nrow(df))))
  
  # Case 2: rownames already equal ID
  already_set <- identical(rownames(df), as.character(df$ID))
  
  if (default_rownames) {
    rownames(df) <- df$ID
  } else if (!already_set) {
    warning("Rownames do not match ID; replacing with ID")
    rownames(df) <- df$ID
  }
  
  # Remove ID column
  df$ID <- NULL
  
  df
}

#' Filter lowly expressed and non protein-coding/lncRNA genes
#'
#' Keeps genes with non-zero variance, log2(TPM + 1) > 1 in more than half
#' of the samples, and biotype protein-coding or lncRNA.
#'
#' @param df data.frame. TPMs (genes x samples).
#'
#' @return Numeric matrix with the TPMs (not log-transformed) of the kept genes.
.filter_data <- function(df){
  # Prepare TPM data
  log2tpm <- log2(df + 1)
  
  # Filtered genes
  genes_to_keep <- readRDS(p12_genes_to_keep)
  gene_exp_data_for_filtering <- log2tpm
  feature_exclusion_df <- data.frame("Gene_name" = rownames(gene_exp_data_for_filtering),
                                     "Zero_variance" = 0,
                                     "Low_variance_0_01" = 0,
                                     "Filter_logtpm_p1" = 0,
                                     "Protein_lnc" = 0,
                                     "Low_mean_logtpm" = 0)
  ## Features with constant variance
  variance_GEfeats <- apply(gene_exp_data_for_filtering, 1, var)
  feature_exclusion_df$Zero_variance[which(variance_GEfeats == 0)] <- 1
  feature_exclusion_df$Low_variance_0_01[which(variance_GEfeats < 0.01)] <- 1
  
  # Filter genes where average logtpm should be 1 in most of the samples, instead of 0.5 in half of the samples
  filter_logtpm <- apply(gene_exp_data_for_filtering, 1, function(x) sum(x > 1)) > (ncol(gene_exp_data_for_filtering) * 0.5)
  feature_exclusion_df$Filter_logtpm_p1[!filter_logtpm] <- 1
  
  ## Mean log TPM filter: if average of gene in samples has log tpm of 1, it will be removed
  low_mean_logtpm <- apply(gene_exp_data_for_filtering, 1, function(x) mean(x) < 1)
  feature_exclusion_df$Low_mean_logtpm[low_mean_logtpm] <-1
  
  ## Identify protein coding genes
  feature_exclusion_df$Protein_lnc[!(gsub("\\..*","",rownames(gene_exp_data_for_filtering)) %in% genes_to_keep$ensembl_gene_id)] <- 1
  
  ## Combine filters
  feature_exclusion_df$filter_logp_and_protein_lnc <- ifelse(
    feature_exclusion_df$Filter_logtpm_p1 == 1 |
      feature_exclusion_df$Protein_lnc == 1,
    1,
    0
  )
  
  df_filt <- data.matrix(df)[intersect(which(feature_exclusion_df[, "filter_logp_and_protein_lnc"] == 0),which(feature_exclusion_df$Zero_variance == 0)),]
  print(paste0("Number of genes after filtering: ",nrow(df_filt)))
  
  df_filt
}

# Validity checks of DiscoveryObject
invisible(setValidity("DiscoveryObject", function(object) {
  
  seq_data <- object@sequencing_data
  
  # Check required elements exist
  if (!all(c("tpms", "counts") %in% names(seq_data))) {
    return("sequencing_data must contain 'tpms' and 'counts'")
  }
  
  tpms <- seq_data$tpms
  counts <- seq_data$counts
  
  # Check they are data.frames
  if (!is.data.frame(tpms) || !is.data.frame(counts)) {
    return("tpms and counts must be data.frames")
  }
    
  # Check number of columns match
  if (ncol(tpms) != ncol(counts)) {
    return("tpms and counts does not the same number of columns")
  }

  # Check against clinical rows
  if (ncol(tpms) != nrow(object@clinical)) {
    return("Number of sequencing samples does not match number of clinical samples")
  }
  
  # Check if IDs are identical
  if (!identical(rownames(tpms), rownames(counts))) {
    return("Gene IDs (rownames) of tpms and counts must match")
  }

  # Check against out
  if (nrow(tpms) == nrow(object@out)) {
    return("Number of gene ids is different between gtf and sequenced data")
  }
  
  TRUE
}))

#' Create a DiscoveryObject
#'
#' @param tpms data.frame. TPMs with an ID column (genes x samples).
#' @param counts data.frame. Raw counts with an ID column (genes x samples).
#' @param clinical data.frame. Clinical data, one row per sample.
#' @param out data.frame. Gene ID to gene name mapping (gene_id, gene_name).
#'
#' @return DiscoveryObject with sequencing_data (tpms, counts and the
#'   filtered tpms_filt, counts_filt), clinical and out.
DiscoveryObject <- function(tpms, counts, clinical, out) {
  
  if (!is.data.frame(tpms))
    stop("tpms must be a data.frame")
  if (!is.data.frame(counts))
    stop("counts must be a data.frame")
  if (!is.data.frame(clinical))
    stop("clinical must be a data.frame")
  
  tpms <- .fix_gene_ids(tpms)
  counts <- .fix_gene_ids(counts)
  
  print("Filtering sequencing data...")
  tpms_filt <- .filter_data(tpms)
  counts_filt <- counts[rownames(counts) %in% rownames(tpms_filt),]
    
  sequencing_data <- list(
    tpms = tpms,
    counts = counts,
    tpms_filt = tpms_filt,
    counts_filt = counts_filt
  )
  
  new("DiscoveryObject",
      sequencing_data = sequencing_data,
      clinical = clinical,
      out = out)
}
#-------------------------------------------------------------------------------

# 3. Validation data object
#-------------------------------------------------------------------------------
setClass(
  "ValidationObject",
  slots = list(
    sequencing_data = "list",
    clinical = "data.frame",
    out = "data.frame"
  )
)

#' Use the gene_id column as rownames (validation data, nf-core output)
#'
#' @param df data.frame. Expression data with gene_id and gene_name columns.
#'
#' @return data.frame with gene IDs as rownames, without gene_id/gene_name.
.fix_gene_ids_v2 <- function(df) {
  
  if (!"gene_id" %in% colnames(df)) {
    stop("gene_id column is required.")
  }
  
  # Case 1: rownames are default
  default_rownames <- identical(rownames(df), as.character(seq_len(nrow(df))))
  
  # Case 2: rownames already equal gene_id
  already_set <- identical(rownames(df), as.character(df$gene_id))
  
  if (default_rownames) {
    rownames(df) <- df$gene_id
  } else if (!already_set) {
    warning("Rownames do not match gene_id; replacing with gene_id")
    rownames(df) <- df$gene_id
  }
  
  # Remove gene_id and gene_name column
  df$gene_id <- NULL
  df$gene_name <- NULL
  
  df
}

# Validity checks of ValidationObject (same checks as DiscoveryObject)
invisible(setValidity("ValidationObject", function(object) {
  
  seq_data <- object@sequencing_data
  
  # Check required elements exist
  if (!all(c("tpms", "counts") %in% names(seq_data))) {
    return("sequencing_data must contain 'tpms' and 'counts'")
  }
  
  tpms <- seq_data$tpms
  counts <- seq_data$counts
  
  # Check they are data.frames
  if (!is.data.frame(tpms) || !is.data.frame(counts)) {
    return("tpms and counts must be data.frames")
  }
  
  # Check number of columns match
  if (ncol(tpms) != ncol(counts)) {
    return("tpms and counts does not the same number of columns")
  }
  
  # Check against clinical rows
  if (ncol(tpms) != nrow(object@clinical)) {
    return("Number of sequencing samples does not match number of clinical samples")
  }
  
  # Check if IDs are identical
  if (!identical(rownames(tpms), rownames(counts))) {
    return("Gene IDs (rownames) of tpms and counts must match")
  }
  
  # Check against out
  if (nrow(tpms) == nrow(object@out)) {
    return("Number of gene ids is different between gtf and sequenced data")
  }
  
  TRUE
}))

#' Create a ValidationObject
#'
#' @param tpms data.frame. TPMs with gene_id/gene_name columns.
#' @param counts data.frame. Raw counts with gene_id/gene_name columns.
#' @param clinical data.frame. Clinical data, one row per sample.
#' @param out data.frame. Gene ID to gene name mapping (gene_id, gene_name).
#'
#' @return ValidationObject with sequencing_data (tpms, counts), clinical and out.
ValidationObject <- function(tpms, counts, clinical, out) {
  
  if (!is.data.frame(tpms))
    stop("tpms must be a data.frame")
  if (!is.data.frame(counts))
    stop("counts must be a data.frame")
  if (!is.data.frame(clinical))
    stop("clinical must be a data.frame")
  
  tpms <- .fix_gene_ids_v2(tpms)
  counts <- .fix_gene_ids_v2(counts)
  
  sequencing_data <- list(
    tpms = tpms,
    counts = counts
  )
  
  new("ValidationObject",
      sequencing_data = sequencing_data,
      clinical = clinical,
      out = out)
}

#-------------------------------------------------------------------------------

# 4. General functions
#-------------------------------------------------------------------------------
#' Subset a data object to a set of samples
#'
#' @param obj DiscoveryObject or ValidationObject.
#' @param samples character vector. Sample IDs to keep.
#' @param sample_id_column character. Clinical column with the sample IDs.
#'
#' @return The object with clinical data, counts and tpms subset to samples.
#'   The filtered matrices (tpms_filt, counts_filt) are not subset.
subset_samples <- function(obj, samples, sample_id_column) {
  
  # Keep only requested samples that are present
  clinical <- obj@clinical[obj@clinical[[sample_id_column]] %in% samples, , drop = FALSE]
  samples <- clinical[[sample_id_column]]
  
  # Subset expression matrices
  obj@sequencing_data$counts <-
    obj@sequencing_data$counts[, samples, drop = FALSE]
  
  obj@sequencing_data$tpms <-
    obj@sequencing_data$tpms[, samples, drop = FALSE]
  
  # Update clinical data
  obj@clinical <- clinical
  
  obj
}


#-------------------------------------------------------------------------------

# 5. Load discovery dataset (N = 366)
#-------------------------------------------------------------------------------
#' Load the discovery dataset after removing bad QC sets (N = 366)
#'
#' Removes complete case-control sets that contain a bad QC sample.
#'
#' @param p_tpms character. Path to the TPM file.
#' @param p_counts character. Path to the count file.
#' @param p_clinical character. Path to the parsed clinical data.
#' @param p_badQC character. Path to the list of bad QC samples.
#' @param p_gtf character. Path to the GTF annotation file.
#' @param output_dir character. Folder where the object is saved.
#'
#' @return Called for its side effect: saves Discovery_data_set_N366.rds.
load_discovery <- function(p_tpms,p_counts,p_clinical,p_badQC,p_gtf,output_dir){

  print("Loading discovery data...")
  
  # Exclude patient sets with a bad QC
  excl_pats <- read.csv(p_badQC)
  excl_pats <- unique(excl_pats$Set_id)
  
  # Load/filter clinical data
  clinical_data <- read.csv(p_clinical)
  clinical_data_sel <- clinical_data[!(clinical_data$Set_id %in% excl_pats),]
  stopifnot(nrow(clinical_data_sel) == 366)
  
  # Load/filter counts and expression data
  tpm_data <- read.csv(p_tpms)
  count_data <- read.csv(p_counts)
  tpm_data_sel <- tpm_data[,c("ID",clinical_data_sel$Skyline_ID)]
  count_data_sel <- count_data[,c("ID",clinical_data_sel$Skyline_ID)]
  
  stopifnot(ncol(tpm_data_sel)==367)
  stopifnot(ncol(count_data_sel)==367)
  
  print("Loading gtf file...")
  # Mapping file
  mapping <- rtracklayer::import(p_gtf)
  out <- mapping %>%
    as.data.frame() %>%
    dplyr::select(gene_id, gene_name) %>%
    unique
  
  discovery_dataset <- DiscoveryObject(tpm_data_sel, count_data_sel, clinical_data_sel, out)
  
  write_rds(discovery_dataset, paste0(output_dir,"/Discovery_data_set_N366.rds"))
}
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------

# 6. Load discovery dataset including samples of incomplete sets (N = 378)
#-------------------------------------------------------------------------------
#' Load the discovery dataset after removing bad QC samples (N = 378)
#'
#' Removes only the bad QC samples themselves, so their matched samples are
#' kept (incomplete sets).
#'
#' @param p_tpms character. Path to the TPM file.
#' @param p_counts character. Path to the count file.
#' @param p_clinical character. Path to the parsed clinical data.
#' @param p_badQC character. Path to the list of bad QC samples.
#' @param p_gtf character. Path to the GTF annotation file.
#' @param output_dir character. Folder where the object is saved.
#'
#' @return Called for its side effect: saves Discovery_data_set_N378.rds.
load_discovery_all <- function(p_tpms,p_counts,p_clinical,p_badQC,p_gtf,output_dir){
  print("Loading discovery data...")
  
  # Exclude patient sets with a bad QC
  excl_pats <- read.csv(p_badQC)
  excl_pats <- unique(excl_pats$Skyline_ID)
  
  # Load/filter clinical data
  clinical_data <- read.csv(p_clinical)
  clinical_data_sel <- clinical_data[!(clinical_data$Skyline_ID %in% excl_pats),]
  stopifnot(nrow(clinical_data_sel) == 378)
  
  # Load/filter counts and expression data
  tpm_data <- read.csv(p_tpms)
  count_data <- read.csv(p_counts)
  tpm_data_sel <- tpm_data[,c("ID",clinical_data_sel$Skyline_ID)]
  count_data_sel <- count_data[,c("ID",clinical_data_sel$Skyline_ID)]
  
  stopifnot(ncol(tpm_data_sel)==379)
  stopifnot(ncol(count_data_sel)==379)
  
  print("Loading gtf file...")
  # Mapping file
  mapping <- rtracklayer::import(p_gtf)
  out <- mapping %>%
    as.data.frame() %>%
    dplyr::select(gene_id, gene_name) %>%
    unique
  
  discovery_dataset <- DiscoveryObject(tpm_data_sel, count_data_sel, clinical_data_sel, out)
  
  write_rds(discovery_dataset, paste0(output_dir,"/Discovery_data_set_N378.rds"))
}
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------

# 7. Load validation datasets
#-------------------------------------------------------------------------------
#' Load the D-SQUAME validation dataset (N = 102)
#'
#' For sets with both a biopsy and an excision, only the biopsy is kept.
#'
#' @param p_counts character. Path to the gene counts file (tsv).
#' @param p_tpms character. Path to the gene TPM file (tsv).
#' @param p_clinical character. Path to the parsed clinical data.
#' @param p_badQC character. Path to the list of bad QC samples (none in the
#'   final dataset).
#' @param p_gtf character. Path to the GTF annotation file.
#' @param output_dir character. Folder where the object is saved.
#'
#' @return Called for its side effect: saves D_SQUAME_validation_data_set.rds.
load_validation_d_squame <- function(p_counts,p_tpms,p_clinical,p_badQC,p_gtf,output_dir){

  
  print("Loading validation data...")
  
  # Exclude patient sets with a bad QC (none in the final validation set, N=102)
  excl_pats <- read.table(p_badQC)
  badQC_samples_id <- excl_pats$Skyline_ID
  
  # Load/filter/prepare clinical data
  clinical_seq_parsed <- read.csv(p_clinical)
  clinical_seq_parsed <- clinical_seq_parsed %>%
    filter(!(Set_info == "complete_set_double_sampletype" & Type_of_material_bin != "Biopsy"))
    
  # Load sequencing data
  tpm_data <- as.data.frame(read_tsv(p_tpms))
  count_data <- as.data.frame(read_tsv(p_counts))
  
  # Fix colnames
  clean_names <- function(x) {
    x |>
      gsub("\\.", "-", x = _) |>
      sub("^X(?=\\d)", "", x = _, perl = TRUE) |>
      sub("r\\d*$", "", x = _)
  }
  
  colnames(count_data) <- clean_names(colnames(count_data))
  colnames(tpm_data) <- clean_names(colnames(tpm_data))
  
  # Remove bad QC samples (none are excluded)
  clinical_data_val <- clinical_seq_parsed %>% filter(!clinical_seq_parsed$Skyline_ID %in% badQC_samples_id)
  
  # Rename sample type column
  clinical_data_val <- clinical_data_val %>% rename(Sample_type=Type_of_material)

  # Filter expression/count data
  tpm_data_val <- tpm_data[,c("gene_id","gene_name",clinical_data_val$Skyline_ID)]
  count_data_val <- count_data[,c("gene_id","gene_name",clinical_data_val$Skyline_ID)]
  
  stopifnot(ncol(tpm_data_val)==104)
  stopifnot(ncol(count_data_val)==104)
  
  print("Loading gtf file...")
  # Mapping file
  mapping <- rtracklayer::import(p_gtf)
  out <- mapping %>%
    as.data.frame() %>%
    dplyr::select(gene_id, gene_name) %>%
    unique
  
  validation_dataset <- ValidationObject(tpm_data_val, count_data_val, clinical_data_val, out)
  
  write_rds(validation_dataset, paste0(output_dir,"/D_SQUAME_validation_data_set.rds"))
}



#' Load the Nassir et al. validation dataset
#'
#' Two samples with bad QC (SRR31748002, SRR31748007) are removed.
#'
#' @param p_tpms character. Path to the gene TPM file (tsv).
#' @param p_counts character. Path to the gene counts file (tsv).
#' @param p_clinical character. Path to the clinical metadata.
#' @param p_gtf character. Path to the GTF annotation file.
#' @param output_dir character. Folder where the object is saved.
#'
#' @return Called for its side effect: saves Nassir_validation_data_set.rds.
load_validation_nassir <- function(p_tpms,p_counts,p_clinical,p_gtf,output_dir){
  print("Loading validation data...")
  
  # Exclude patient sets with a bad QC
  bad_samples <-  c("SRR31748002","SRR31748007") 
  
  # Load/filter clinical data
  clinical_data <- read.csv(p_clinical)
  clinical_data_nassir <- clinical_data[!(clinical_data$Run %in% bad_samples),]
  
  # Load/filter counts and expression data
  tpm_data <- as.data.frame(read_tsv(p_tpms))
  count_data <- as.data.frame(read_tsv(p_counts))
  
  # All samples, including controls
  tpm_data_nassir <- tpm_data[,c("gene_id","gene_name",clinical_data_nassir$Run)]
  count_data_nassir <- count_data[,c("gene_id","gene_name",clinical_data_nassir$Run)]
  
  print("Loading gtf file...")
  # Mapping file
  mapping <- rtracklayer::import(p_gtf)
  out <- mapping %>%
    as.data.frame() %>%
    dplyr::select(gene_id, gene_name) %>%
    unique
  
  validation_dataset <- ValidationObject(tpm_data_nassir, count_data_nassir, clinical_data_nassir, out)
  
  write_rds(validation_dataset, paste0(output_dir,"/Nassir_validation_data_set.rds"))
}

#-------------------------------------------------------------------------------