#---------------------------------------------------
# Aim: Save the sample summary dataframe with all the information from samples (clinical + genomic)
# Author: B. Rentroia Pacheco
# Input: Sample summary
# Output: Excel and Upset plot
#---------------------------------------------------


#---------------------------------------------------
# 1. Tidy up the sample summary dataframe to only contain the information we need:
#---------------------------------------------------


sample_summary_publication = merge(sample_summary_total,sample_summary[,c("Sample.ID","pct_altered")],all.x=TRUE)%>%
  select(Sample.ID,Full_Sample_ID,Pt_ID,Starting.DNA.input.ng.Tumor,Starting.DNA.input.ng.Normal,Library.DNA.yield.ng.Tumor,Library.DNA.yield.ng.Normal,Mean.bait.coverage.Tumor,Mean.bait.coverage.Normal,Bait.territory,Deduplicated.Coverage.Tumor.gr5,Deduplicated.Eff.Coverage.Tumor.gr7,Deduplicated.Coverage.Normal.gr5,Deduplicated.Coverage,Percent_Bait_with_callable_coverage,Exclude_due_to_callable_coverage,RNA_Seq_mismatches_or_contamination,Exclude_sample,`RNA-seq_available`,Tumor.cellularity.avg.pct.with.RESCUED,Tumor.cellularity.Method.1,Value.1, Tumor.cellularity.Method.2,Value.2, Tumor.cellularity.Method.3,Value.3, Tumor.cellularity.Method.4, Value.4,Total.mutations.nr.with.RESCUED,Dominant.clone.mutations.nr.with.RESCUED,`Dominant.clone.UV.mutations.(%).with.RESCUED`,`Burden.(Mutations/megabase).(MB).with.RESCUED`,pct_altered,Sex,Age,IS.at.cSCC,Number.of.cSCC.before.culprit,Tumor.location,Tumor.location.morecats,AJCC.8,BWH,CP.score,Tumor.diameter,Tissue.involvement,Differentiation,PNI.or.LVI,Depth.of.Invasion,Breslow.thickness,Resection.margin.cat,Invasion.of.bones,Peritumoral.infiltration,Solar.elastosis,Tumor.budding,Mitotic.rate,Morphology.subtype,Metastasis,FU.metastasis.years,Vital.status,Vital.follow.up.years)%>%
  dplyr::rename(fraction_genome_altered_CN = pct_altered,
                WES_exclude_sample = Exclude_sample)%>%
  as.data.frame()

#---------------------------------------------------
# 2. Add necessary RNA-seq information
#---------------------------------------------------

# Add the samples that we are missing:
sample_summary_publication_all = sample_summary_publication %>%full_join(df_map %>% dplyr::rename(Pt_ID = universal_patient_ID)%>%select(Pt_ID,SkylineDx.ID), by="Pt_ID")%>%
  mutate(RNA_seq_available = ifelse(SkylineDx.ID%in%df_rnaseq_pairs$Skyline_ID,"Excluded_Pair",ifelse(SkylineDx.ID%in%df_map$SkylineDx.ID[!is.na(df_map$SkylineDx.ID)],"Yes","No")))%>%
  mutate(RNA_seq_available  = ifelse(SkylineDx.ID%in%df_rnaseq_exc$Skyline_ID,"Excluded",RNA_seq_available ))%>%
  relocate(RNA_seq_available , .after = `RNA-seq_available`)%>%
  select(-`RNA-seq_available`)%>%
  relocate(SkylineDx.ID, .after = Pt_ID)%>%
  mutate(RNA_Seq_Total_nr_reads = NA,
         RNA_Seq_Pct_reads_aligned_ref_genome=NA,
         RNA_Seq_Pct_coding_reads = NA)%>%
  relocate(RNA_Seq_Total_nr_reads,RNA_Seq_Pct_reads_aligned_ref_genome,RNA_Seq_Pct_coding_reads,RNA_seq_available, .after =SkylineDx.ID)%>%
  mutate(WES_available = ifelse(!is.na(WES_exclude_sample)&WES_exclude_sample=="Yes","Excluded",ifelse(is.na(WES_exclude_sample),"No","Yes")))%>%
  select(-WES_exclude_sample)%>%
  relocate(WES_available, .after = RNA_Seq_mismatches_or_contamination)

#---------------------------------------------------
# 3. Upset plot
#---------------------------------------------------

# Generate Upset plot
df_upset = sample_summary_publication_all%>%
  mutate(RNA_seq_available=ifelse(RNA_seq_available=="Excluded_Pair","Yes",RNA_seq_available))%>%
  mutate(
    RNAseq_PASS = RNA_seq_available == "Yes",
    RNAseq_NA  = RNA_seq_available == "No",
    RNAseq_FAIL = RNA_seq_available == "Excluded",
    WES_PASS = WES_available == "Yes",
    WES_NA  = WES_available == "No",
    WES_FAIL = WES_available == "Excluded"
  )%>%
  mutate(
    color_group = case_when(
      RNAseq_PASS & WES_PASS                          ~ "RNAseq + DNAseq",
      RNAseq_PASS & (WES_FAIL | WES_NA)               ~ "RNAseq only",
      WES_PASS    & (RNAseq_FAIL | RNAseq_NA)         ~ "WES only",
      (RNAseq_FAIL | RNAseq_NA) & (WES_FAIL | WES_NA) ~ "Neither"
    )
  )

set_order <- c("WES_FAIL", "WES_NA", "RNAseq_FAIL", "RNAseq_NA", "WES_PASS", "RNAseq_PASS")

upset_plot = ComplexUpset::upset(
  df_upset,
  intersect = set_order,
  sort_sets=FALSE,
  sort_intersections = FALSE,
  intersections = list(
    # RNAseq PASS group
    c("RNAseq_PASS", "WES_PASS"),
    c("RNAseq_PASS", "WES_FAIL"),
    c("RNAseq_PASS", "WES_NA"),
    c("RNAseq_FAIL", "WES_NA"),   # <-- moved here, after RNAseq_PASS + WES_NA
    
    # WES PASS group
    c("RNAseq_FAIL", "WES_PASS"),
    c("RNAseq_NA",   "WES_PASS"),
    
    # Neither
    c("RNAseq_FAIL", "WES_FAIL"),
    c("RNAseq_NA",   "WES_FAIL"),
    c("RNAseq_NA",   "WES_NA")
  ),set_sizes = upset_set_size() +
    geom_text(aes(label = after_stat(count)),
              hjust = 1.2,
              stat = 'count',
              size = 4)+ expand_limits(y = 450),
  base_annotations = list(
    'Intersection size' = intersection_size(
      mapping = aes(fill = color_group)
    ) +
      scale_fill_manual(values = c(
        "RNAseq + DNAseq" = "#2166ac",
        "RNAseq only" = "#92c5de",
        "WES only"    = "#f4a582",
        "Neither"     = "#d9d9d9"
      )) +
      geom_vline(xintercept = c(1.5,4.5),
                 linetype = "dashed",
                 color = "black",
                 linewidth = 0.7)
  ),
  min_size = 1,
  themes = upset_default_themes(text = element_text(size = 15))
) +
  ggtitle("Tumors by RNA-seq and WES availability")

ggsave(file.path(dir_results,"05_01_01_Upset_Samples.pdf"),upset_plot,width=10,height=6)

write.xlsx(sample_summary_publication_all,file.path(dir_results,"05_01_01_Sample_summary_v2.xlsx"))

#---------------------------------------------------
# 3. Add some additional information
#---------------------------------------------------
print("JJH Traets has added the information from the RNA-seq samples")
print("Reading: 05_01_01_Sample_summary_v3_updatedJT_unimputed")
sample_summary_publication_wRNAseq = read.xlsx(file.path(dir_results,"05_01_01_Sample_summary_v3_updatedJT_unimputed.xlsx"))

# Step 1: Check that previous information hasn't changed:
# Find positions where original dataframe was NOT NA
non_na_positions <- !is.na(sample_summary_publication_all)

# Compare only those positions
#differences <- sample_summary_publication_wRNAseq[non_na_positions] !=  sample_summary_publication_all[non_na_positions]
#diff_idx <- which(differences, arr.ind = TRUE)
# It is only rounding issue:
#all(round(as.numeric(sample_summary_publication_wRNAseq[non_na_positions][diff_idx]),4)==round(as.numeric(sample_summary_publication_all[non_na_positions][diff_idx]),4))
#differences_na <- is.na(sample_summary_publication_wRNAseq[non_na_positions]) 
#any(differences_na )

# Step 2: 
# Add information from IS SBS32 patients:
sample_summary_publication_wRNAseq$SBS32_presence = NA
sample_summary_publication_wRNAseq$SBS32_presence[which(sample_summary_publication_wRNAseq$WES_available=="Yes")] = "No"
sample_summary_publication_wRNAseq$SBS32_presence[sample_summary_publication_wRNAseq$Sample.ID%in%gsub("S","",pts_with_SBS_32)]="Yes"

# Add information of the genomic clusters: 
# Recoding:
genomic_clusters$cluster <- dplyr::recode(
  genomic_clusters$cluster,
  "Poorly differentiated" = "Mesenchymal-like",
  "Moderately differentiated" = "Basal-like",
  "Well-differentiated" = "Differentiated"
)
sample_summary_publication_wRNAseq = merge(sample_summary_publication_wRNAseq,genomic_clusters[,c("Sample_id","cluster","type")],by.x="SkylineDx.ID",by.y="Sample_id",all.x=TRUE)

# Change name of the Tumor cellularity methods:
for (tcel_column in c("Tumor.cellularity.Method.1","Tumor.cellularity.Method.2","Tumor.cellularity.Method.3","Tumor.cellularity.Method.4")){
  sample_summary_publication_wRNAseq[,tcel_column]=dplyr::recode(
    sample_summary_publication_wRNAseq[,tcel_column],
    "Allelic Imbalance of SNPs over Copy Neutral LOH" = "Allelic Imbalance of SNPs over copy number neutral LOH",
    "Chromosome XY MAF" = "Median Somatic MAF in sex chromosomes",
    "Deletion formula" = "Copy number deletion log ratios",
    "Allelic Imbalance of SNPs over Copy number deletions LOH"= "Allelic Imbalance of SNPs over copy number deletions")
}

# Reorder columns:
sample_summary_publication_wRNAseq_final = sample_summary_publication_wRNAseq %>%
  select(Sample.ID,Full_Sample_ID,Pt_ID,SkylineDx.ID,Sex,Age,IS.at.cSCC,SBS32_presence,Number.of.cSCC.before.culprit,Tumor.location,Tumor.location.morecats,AJCC.8,BWH,CP.score,Tumor.diameter,Tissue.involvement,Differentiation,PNI.or.LVI,Depth.of.Invasion,Breslow.thickness,Resection.margin.cat,Invasion.of.bones,Peritumoral.infiltration,Solar.elastosis,Tumor.budding,Mitotic.rate,Morphology.subtype,Metastasis,FU.metastasis.years,Vital.status,Vital.follow.up.years,
         RNA_Seq_Total_nr_reads,RNA_Seq_Pct_reads_aligned_ref_genome,RNA_Seq_Pct_coding_reads,RNA_seq_available,type,cluster,
         Starting.DNA.input.ng.Tumor,Starting.DNA.input.ng.Normal,Library.DNA.yield.ng.Tumor,Library.DNA.yield.ng.Normal,Mean.bait.coverage.Tumor,Mean.bait.coverage.Normal,Bait.territory,Deduplicated.Eff.Coverage.Tumor.gr7,Deduplicated.Coverage.Normal.gr5,Deduplicated.Coverage,Percent_Bait_with_callable_coverage,Exclude_due_to_callable_coverage,RNA_Seq_mismatches_or_contamination,WES_available,
         Tumor.cellularity.avg.pct.with.RESCUED,Tumor.cellularity.Method.1,Value.1,Tumor.cellularity.Method.2,Value.2,Tumor.cellularity.Method.3,Value.3,Tumor.cellularity.Method.4,Value.4,Total.mutations.nr.with.RESCUED,Dominant.clone.mutations.nr.with.RESCUED,'Dominant.clone.UV.mutations.(%).with.RESCUED','Burden.(Mutations/megabase).(MB).with.RESCUED',fraction_genome_altered_CN)%>%
  mutate(fraction_genome_altered_CN = 100*fraction_genome_altered_CN )%>%
  rename_with(~ gsub("\\.", " ", .x))%>%
  rename_with(~ gsub("_", " ", .x))%>%
  rename_with(~ gsub("years", "\\(years\\)", .x))%>%
  rename_with(~ gsub("(.+)Tumor", "\\1\\(Tumor\\)", .x))%>%
  rename_with(~ gsub("(.+)Normal", "\\1\\(Normal\\)", .x))%>%
  rename_with(~ gsub(" with RESCUED", "", .x))%>%
  rename("Immunossuppressed"= "IS at cSCC",
         "Tumor location (>3 categories)" ="Tumor location morecats",
         "AJCC8"="AJCC 8",
         "EMC model score"="CP score", 
         "Resection margin"="Resection margin cat",
         "RNA-Seq Total number of reads"="RNA Seq Total nr reads",
         "RNA-Seq reads aligned to reference genome (%)"="RNA Seq Pct reads aligned ref genome",
         "RNA-Seq coding reads (%)"="RNA Seq Pct coding reads",
         "Callable coverage (%)"="Percent Bait with callable coverage",
         "Tumor cellularity (%)" = "Tumor cellularity avg pct",
         "Number of somatic mutations" ="Total mutations nr",
         "Number of dominant clone mutations" = "Dominant clone mutations nr",
         "Genome altered (%)" = "fraction genome altered CN",
         "Transcriptomic cluster" =cluster,
         "RNA-Seq material"=type)

write.xlsx(sample_summary_publication_wRNAseq_final,file.path(dir_results,"05_01_01_Sample_summary_v4.xlsx"))
