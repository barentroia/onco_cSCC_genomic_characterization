#---------------------------------------------------
# Aim: Analyze copy number data
# Author: B. Rentroia Pacheco
# Input: CNV estimates
# Output: Summary of CNV
#---------------------------------------------------

#-----------------------
# 0. Group all segment files together:
#-----------------------
df_n_cutoffs = data.frame(Cutoffs = c(0,13.5,26.8), N_total = NA,N_cases = NA,N_controls=NA)
samples_cutoff= list()
for(ct in df_n_cutoffs$Cutoffs){
  i.ct = which(df_n_cutoffs$Cutoffs==ct)
  df_n_cutoffs[i.ct,c("N_cases","N_controls")] = sample_summary%>%filter(Exclude_sample=="No")%>%
    filter(Tumor.cellularity.avg.pct>ct)%>%
    pull(Metastasis)%>%table()
  df_n_cutoffs[i.ct,"N_total"]=df_n_cutoffs[i.ct,"N_cases"]+df_n_cutoffs[i.ct,"N_controls"]
  samples_cutoff[[paste0("All_CC_but_tcel_",ct*0.01)]]=sample_summary%>%filter(Exclude_sample=="No")%>%
    filter(Tumor.cellularity.avg.pct>ct)%>%
    pull(Sample.ID)
}

write.xlsx(df_n_cutoffs,file.path(dir_results_script,"04_02_02_sample_size_exclusion_CNV.xlsx"))

#source(file.path(dir_scripts,"04_Publication_scripts","04_02_02_A1_CNV_folder_reorganization.R"))

# This contains all the cns files computed by cnvkit:
cns_files_total = read.xlsx(file.path(dir_results_script,"cns","All_CC_but_tcel_0","Tumor",paste0("04_02_02_cns_files_total_All_CC_but_tcel_0.xlsx")))

#-----------------------
# 1. Generate a genome-wide copy number plot:
#-----------------------
# This was generated in IGV, using the segment files generated in "04_02_02_A1_CNV_folder_reorganization.R"

#-----------------------
# 2. Merge chromosome band and cns data, to have log ratios per chromosome band
#-----------------------
# Convert dataframes to GRanges objects
cns_gr <- GRanges(
  seqnames = cns_files_total$chromosome,
  ranges = IRanges(start = cns_files_total$start, end = cns_files_total$end),
  Sample_id= cns_files_total$Sample_id,
  log2 = cns_files_total$log2,
  probes = cns_files_total$probes
)

# Convert chromosome bands into a GRange object:
chromosome_bands$chr_arm = gsub("[^a-zA-Z]", "", chromosome_bands$V4)
cytoband_gr <- GRanges(
  seqnames = chromosome_bands$V1,
  ranges = IRanges(start = chromosome_bands$V2, end = chromosome_bands$V3),
  band = chromosome_bands$V4,
  stain = chromosome_bands$V5,
  chr_arm = gsub("[^a-zA-Z]", "", chromosome_bands$V4),
  start_pos = chromosome_bands$V2,
  end_pos = chromosome_bands$V3
)

# Find overlaps, to summarize log ratios for each chr arm:
overlaps <- findOverlaps(cns_gr, cytoband_gr)
hits <- findOverlaps(cytoband_gr,cns_gr) # These overlaps are when the main data frame are the chromosomes

df <- data.frame(
  band_idx = queryHits(hits),
  chr = seqnames(cytoband_gr)[queryHits(hits)],
  chr_arm = mcols(cytoband_gr)$chr_arm[queryHits(hits)],
  start_band_pos = mcols(cytoband_gr)$start_pos[queryHits(hits)],
  end_band_pos = mcols(cytoband_gr)$end_pos[queryHits(hits)],
  sample = mcols(cns_gr)$Sample_id[subjectHits(hits)],
  log2 = mcols(cns_gr)$log2[subjectHits(hits)],
  probes = mcols(cns_gr)$probes[subjectHits(hits)],
  start_seg_pos = start(ranges(cns_gr))[subjectHits(hits)],
  end_seg_pos = end(ranges(cns_gr))[subjectHits(hits)]
 
)

df_minmax <- df %>%
  group_by(chr, chr_arm) %>%
  summarise(
    min_start = min(start_band_pos, na.rm = TRUE),
    max_end   = max(end_band_pos, na.rm = TRUE),
    .groups = "drop"
  )
#write.xlsx(df_minmax, paste0(dir_results_script,"S04_03_01_minmax_chr_bands.xlsx"))

result_arm_level = df%>%
  filter(probes>100)%>%
  mutate(
    overlap_length = pmin(end_seg_pos, end_band_pos) -
      pmax(start_seg_pos, start_band_pos)
  ) %>%
  group_by(band_idx, sample,chr,chr_arm,start_band_pos,end_band_pos) %>%
  summarise(log2_median_band = sum(log2 * overlap_length) / sum(overlap_length),
            Segment_length = sum(overlap_length)) %>%
  #mutate(band_length = end_band_pos-start_band_pos,
  #       log2_ratio_leng_weighted = band_length*log2_median_band)%>%
  mutate(log2_ratio_leng_weighted = Segment_length*log2_median_band)%>%
  group_by(chr,chr_arm,sample)%>%
  summarise(log2_median_arm = sum(log2_ratio_leng_weighted)/sum(Segment_length))%>%
  as.data.frame()

#-----------------------
# 3. Identify deletions and amplifications using a cutoff
#-----------------------
# Add some extra information on tumor cellularity:
result_arm_level$abs_log2_median = abs(result_arm_level$log2_median_arm)
result_arm_level=merge(result_arm_level,sample_summary[,c("Sample.ID","Tumor.cellularity.avg.pct")],by.x="sample",by.y="Sample.ID",all.x=TRUE)
cns_files_total =merge(cns_files_total,sample_summary[,c("Sample.ID","Metastasis","Percent_Bait_with_callable_coverage",IS_classification)],by.x="Sample_id",by.y="Sample.ID",all.x=TRUE)

# Now turn copy number results into a binary matrix:

# We save the original dataframes so that they can be used back
cns_files_total_all = cns_files_total 
result_arm_level_all = result_arm_level

for(tcel_cutoff in c("All_CC_but_tcel_0","All_CC_but_tcel_0.135","All_CC_but_tcel_0.268")){
  # Sample exclusion:
  keep_samples = samples_cutoff[[tcel_cutoff]]
  
  result_arm_level = result_arm_level_all %>%
    filter(sample %in%keep_samples)%>%
    as.data.frame()
  cns_files_total = cns_files_total_all %>%
    filter(Sample_id%in%keep_samples)%>%
    as.data.frame()
  
  
  dir_results_cutoff = file.path(dir_results_script,"cns",tcel_cutoff)
  
  # Cutoff determination:
  ct_nr = ct_cnv
  bin_matrix <-result_arm_level %>%
    mutate(chr_arm_combined = paste0(chr,chr_arm),
           altered = case_when(
             abs_log2_median > ct_nr & log2_median_arm > 0  ~  1,   # amplification
             abs_log2_median > ct_nr & log2_median_arm < 0  ~ -1,   # deletion
             TRUE                                     ~  0    # neutral
           ))%>%
    select(sample,chr_arm_combined,altered)%>%
    distinct()%>%
    pivot_wider(names_from = chr_arm_combined,
                values_from = altered, 
                values_fill = 0)%>%
    as.data.frame()
  
  cns_files_total$tcel_cutoff_del = -ct_nr
  cns_files_total$tcel_cutoff_amp = ct_nr
  
  
  # Same on the segments, for computation of the genome altered:
  cns_segments_summary = cns_files_total %>%
    filter(probes>100)%>%
    filter(!chromosome %in% c("chrX","chrY"))%>%
    mutate(seg_lenght = end-start,
           amp_det = ifelse(log2>abs(tcel_cutoff_amp),1,0),
           del_det =  ifelse(log2<tcel_cutoff_del,1,0))%>%
    group_by(Sample_id,Metastasis,Percent_Bait_with_callable_coverage,Tumor.cellularity.avg.pct, !!rlang::sym(IS_classification))%>%
    summarise(
      altered_genome = sum(seg_lenght*amp_det)+sum(seg_lenght*del_det),
      total_covered_genome = sum(seg_lenght))%>%
    ungroup()%>%
    mutate(pct_altered = altered_genome/total_covered_genome)
  
  # Paired test - all not significant
  #cns_segments_summary$Set_id = factor(gsub("_.*","",cns_segments_summary$Sample_id))
  #paired_df <- cns_segments_summary %>%
  #  select(Set_id, Metastasis, pct_altered) %>%
  # pivot_wider(names_from = Metastasis,
  #              values_from = pct_altered) %>%
  #  drop_na(Control, Case)
  #paired_wilcox = wilcox.test(
  #  paired_df$Case,
  #  paired_df$Control,
  #  paired = TRUE
  #)
  #print(paired_wilcox$p.value)
  #nrow(paired_df)
  
  # Fit linear model adjusting for callable coverage
  lm_model <- lm(pct_altered ~ Metastasis + Percent_Bait_with_callable_coverage, data = cns_segments_summary)
  
  # Extract p-value for Metastasis
  p_val <- summary(lm_model)$coefficients["MetastasisControl", "Pr(>|t|)"]  # adjust name if needed
  # Format for display
  p_label <- ifelse(p_val < 0.001, "p < 0.001", paste0("p adj call.cov= ", signif(p_val, 2)))
  
  # Visual comparison of the fraction of genome altered:
  stat.test <- compare_means(
    pct_altered ~ Metastasis,
    data = cns_segments_summary,
    method = "wilcox.test"
  )%>%
    mutate(label = paste0("p-val: ", p.format))
  
  
  p_box = ggplot(cns_segments_summary , aes(x = Metastasis,
                                            y = pct_altered)) +
    geom_boxplot( outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.8, size = 1) +
    theme_bw(base_size = 18) +
    labs(
      x = "",
      y = "Percentage of genome altered")+   stat_pvalue_manual(
      stat.test,
      label = "label",
      y.position = max(cns_segments_summary$pct_altered) * 1.1,size=5)+
    annotate("text", x = 1.5, y = max(cns_segments_summary$pct_altered) * 1.05,
             label = p_label, size = 5)
  
  ggsave(paste0(dir_results_cutoff,"/S04_02_02_alt_genome_boxplot",ct_cnv,".pdf"),p_box,width=5,height=6)
  p_dot = ggplot(cns_segments_summary,aes(x=Tumor.cellularity.avg.pct,y=pct_altered))+geom_point()+xlab("Tumor cellularity (%)")+ylab("Pct altered genome")+theme_bw()
  ggsave(paste0(dir_results_cutoff,"/S04_02_02_alt_genome_dotplot",ct_cnv,".pdf"),p_dot,width=5,height=6)
  
  write.xlsx(cns_segments_summary,paste0(dir_results_cutoff,"/S04_02_02_alt_genome_",ct_cnv,".xlsx"))
  
  # Same for Immunosuppresion vs immunocompetent patients:
  # Fit linear model adjusting for callable coverage
  formula_txt <- paste0( "pct_altered ~ ", IS_classification, " + Percent_Bait_with_callable_coverage")
  
  # Convert string to formula
  f <- as.formula(formula_txt)
  
  # Fit model (linear regression )
  lm_model <- lm(f, data = cns_segments_summary)
 
  # Extract p-value for Metastasis
  p_val <- summary(lm_model)$coefficients[paste0(IS_classification,"Yes"), "Pr(>|t|)"]  # adjust name if needed
  # Format for display
  p_label <- ifelse(p_val < 0.001, "p < 0.001", paste0("p adj call.cov= ", signif(p_val, 2)))
  
  # Visual comparison of the fraction of genome altered:
  stat.test <- compare_means(
    formula = as.formula(paste("pct_altered ~", IS_classification)),
    data = cns_segments_summary,
    method = "wilcox.test"
  )%>%
    mutate(label = paste0("p-val: ", p.format))
  
  p_box = ggplot(cns_segments_summary , aes(x =  !!rlang::sym(IS_classification),
                                            y = pct_altered)) +
    geom_boxplot( outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.8, size = 1) +
    theme_bw(base_size = 18) +
    labs(
      x = "",
      y = "Percentage of genome altered")+   stat_pvalue_manual(
        stat.test,
        label = "label",
        y.position = max(cns_segments_summary$pct_altered) * 1.1,size=5)+
    annotate("text", x = 1.5, y = max(cns_segments_summary$pct_altered) * 1.05,
             label = p_label, size = 5)
  
  ggsave(paste0(dir_results_cutoff,"/S04_02_02_alt_genome_boxplot",ct_cnv,"_IS_vs_IC.pdf"),p_box,width=5,height=6)
  
  # Extract column names (excluding 'sample')
  chr_cols <- colnames(bin_matrix)[-1]
  
  # Function to extract chromosome number and arm
  parse_chr <- function(x) {
    # Extract numeric part, but handle chrX / chrY
    chrom <- sub("chr", "", gsub("(p|q)$", "", x))
    chrom_num <- suppressWarnings(as.numeric(chrom))
    chrom_num[is.na(chrom_num) & chrom == "X"] <- 23
    chrom_num[is.na(chrom_num) & chrom == "Y"] <- 24
    data.frame(name = x, chrom_num, arm = ifelse(grepl("p$", x), "p", "q"))
  }
  
  # Create ordering dataframe
  order_df <- parse_chr(chr_cols)
  # Sort by chromosome number, then by arm (p before q)
  order_df <- order_df[order(order_df$chrom_num, order_df$arm), ]
  
  # Reorder columns in bin_matrix
  bin_matrix_ordered <- bin_matrix[, c("sample", order_df$name)]
  
  # Check new order
  #colnames(bin_matrix_ordered)
  colnames(bin_matrix_ordered)[which(colnames(bin_matrix_ordered)=="sample")]="Tumor_Sample_Barcode"
  write.xlsx(bin_matrix_ordered,paste0(dir_results_cutoff,"/S04_02_02_bin_matrix_cnv_armlevel_",ct_nr,".xlsx"))
  
  # Add information of fraction of genome altered to th sample summary:
  if(tcel_cutoff == include_CNV_samples){
    sample_summary = merge(sample_summary,cns_segments_summary[,c("Sample_id","pct_altered")],by.x="Sample.ID",by.y="Sample_id",all.x=TRUE)
  }
}

copy_number_data=read.xlsx(paste0(dir_results_script,"/cns/",include_CNV_samples,"/S04_02_02_bin_matrix_cnv_armlevel_",ct_nr,".xlsx"))

