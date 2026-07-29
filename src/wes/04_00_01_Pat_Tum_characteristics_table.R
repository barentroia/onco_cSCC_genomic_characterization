#---------------------------------------------------
# Aim: Prepare the patient and tumor characteristics tables
# Author: B. Rentroia Pacheco
# Input: Clinical dataset
# Output: Reduced dataset for publication
#---------------------------------------------------

#-----------------------
# Overview of scripts
#-----------------------
# Build first version of the sample summary:
# 02_06_01
# 02_06_01_a1: add information on rescued mutations
# DA_04: shows reasons for exclusion

#-----------------------
# 1. Produce patient and tumor characteristics table
#-----------------------
sample_summary_table = sample_summary_total %>% 
  mutate(Set.id = gsub("_0|_1","",Sample.ID))%>%
  mutate(Starting.DNA.input.ng.Tumor = as.numeric(Starting.DNA.input.ng.Tumor ),
         Starting.DNA.input.ng.Normal = as.numeric(Starting.DNA.input.ng.Normal),
         Library.DNA.yield.ng.Tumor = as.numeric(Library.DNA.yield.ng.Tumor),
         Library.DNA.yield.ng.Normal = as.numeric(Library.DNA.yield.ng.Normal),
         Mean.bait.coverage.Tumor = as.numeric(Mean.bait.coverage.Tumor),
         Mean.bait.coverage.Normal =as.numeric(Mean.bait.coverage.Normal ),
         Set.id = factor(Set.id),
         Sex = factor(Sex,levels=c("Female","Male")),
         Differentiation = factor(Differentiation, levels=c("Good/moderate","Poor/undifferentiated")),
         Tissue.involvement = factor(Tissue.involvement,levels=c("Dermis","Subcutaneous fat","Beyond subcutaneous fat")),
         HM.at.cSCC = factor(HM.at.cSCC,levels=c("No","Yes")),
         OTR.at.cSCC = factor(OTR.at.cSCC,levels=c("No","Yes")),
         IS.at.cSCC= factor(IS.at.cSCC,levels=c("No","Yes")),
         Metastasis = factor(Metastasis,levels=c("Control","Case")),
         Vital.status = factor(Vital.status,levels=c("Alive","Dead")),
         Tumor.location = factor(Tumor.location,levels=c("Face","Scalp/neck","Trunk/Extremities")),
         AJCC.8 = factor(AJCC.8,levels=c("T1","T2","T3","T4")),
         BWH = factor(BWH,levels=c("T1","T2a","T2b","T3")),
         PNI.bin = factor(PNI.bin,levels=c("No","Yes")),
         Lymphovascular.invasion.bin = factor(Lymphovascular.invasion.bin,levels=c("No","Yes")),
         PNI.or.LVI = factor(PNI.or.LVI ,levels=c("No","Yes")),
         Resection.margin.cat = factor(Resection.margin.cat,levels=c("R0 (Radical)","R1/R2 (irradicaal)")),
         Invasion.of.bones = factor(Invasion.of.bones,levels=c("No","Yes")),
         Peritumoral.infiltration = factor(Peritumoral.infiltration,levels=c("Absent/moderate","Abundant")),
         Solar.elastosis = factor(Solar.elastosis,levels=c("Absent/moderate","Extensive")),
         Morphology.subtype = factor(Morphology.subtype,levels=c("Acantholytic","Clear cell","dd Keratoacanthoma","Desmoplastic","Not classifiable","Not otherwise specified","Spindle cell")),
         Exclude_sample = factor(Exclude_sample,levels=c("No","Yes")))%>%
  mutate(Reason_exclusion = ifelse(Exclude_due_to_callable_coverage=="Pass"&Exclude_sample=="Yes","Other","None"))%>%
  mutate(Reason_exclusion = ifelse(Exclude_due_to_callable_coverage=="Fail","Low callable coverage",Reason_exclusion))%>%
  as.data.frame()

# Add reason for exclusion:
sample_summary_table%>%
  select(Reason_exclusion,Metastasis,FU.metastasis.years,Vital.status, Vital.follow.up.years,Sex,Age,Number.of.cSCC.before.culprit,Tumor.diameter,Differentiation,Breslow.thickness,PNI.or.LVI, Tissue.involvement,HM.at.cSCC,OTR.at.cSCC,IS.at.cSCC,Tumor.location,Tumor.location.morecats,CP.score,BWH,AJCC.8)%>%
  tbl_summary(by=Reason_exclusion,type = list(where(is.integer) ~ "continuous",where(is.numeric) ~ "continuous"))%>%
  add_n() %>% # add column with total number of non-missing observations
  add_p(list(all_continuous() ~ "wilcox.test")) %>%
  gtsummary::as_tibble() %>% 
  writexl::write_xlsx(., file.path(dir_results_script,"04_00_01_pat_characteristics_reduced_exc_vs_nexc.xlsx"))

# Excluded samples:
sample_summary_table%>%filter(Exclude_sample=="Yes")%>%
  select(Set.id,Metastasis,FU.metastasis.years,Vital.status, Vital.follow.up.years,Sex,Age,Number.of.cSCC.before.culprit,Tumor.diameter,Differentiation, Breslow.thickness,PNI.or.LVI,Tissue.involvement,HM.at.cSCC,OTR.at.cSCC,IS.at.cSCC,Tumor.location,Tumor.location.morecats,CP.score,BWH,AJCC.8)%>%
  tbl_summary(by=Metastasis,type = list(where(is.integer) ~ "continuous",where(is.numeric) ~ "continuous"),include = -Set.id)%>%
  add_n() %>% # add column with total number of non-missing observations
  add_p(list(all_continuous() ~ "paired.wilcox.test",all_categorical()~"mcnemar.test"),group=Set.id) %>% # test for a difference between groups
  gtsummary::as_tibble() %>% 
  writexl::write_xlsx(., file.path(dir_results_script,paste0("04_00_01_pat_characteristics_reduced_excluded.xlsx")))

sample_summary_table%>%filter(Exclude_sample=="No",`RNA-seq_available`!="No")%>%
  select(Set.id,Metastasis,FU.metastasis.years,Vital.status, Vital.follow.up.years,Sex,Age,Number.of.cSCC.before.culprit,Tumor.diameter,Differentiation, Breslow.thickness,PNI.or.LVI,Tissue.involvement,HM.at.cSCC,OTR.at.cSCC,IS.at.cSCC,Tumor.location,Tumor.location.morecats,CP.score,BWH,AJCC.8)%>%
  tbl_summary(by=Metastasis,type = list(where(is.integer) ~ "continuous",where(is.numeric) ~ "continuous"),include = -Set.id)%>%
  add_n() %>% # add column with total number of non-missing observations
  add_p(list(all_continuous() ~ "paired.wilcox.test",all_categorical()~"mcnemar.test"),group=Set.id) %>% # test for a difference between groups
  gtsummary::as_tibble() %>% 
  writexl::write_xlsx(., file.path(dir_results_script,paste0("04_00_01_pat_characteristics_reduced_afterexclusion_with_RNASEQ.xlsx")))

# Also extended set of characteristics:
n_total =  sample_summary_table%>%filter(Exclude_sample=="No")%>%nrow()
sample_summary_table[,setdiff(colnames(sample_summary_table),c("Sample.ID","Value.1","Value.2","Value.3","Value.4"))]%>%
  filter(Exclude_sample=="No")%>%
  mutate(Starting.DNA.input.ng.Tumor=as.numeric(Starting.DNA.input.ng.Tumor),
         Starting.DNA.input.ng.Normal=as.numeric(Starting.DNA.input.ng.Normal),
         Library.DNA.yield.ng.Tumor= as.numeric(Library.DNA.yield.ng.Tumor),
         Library.DNA.yield.ng.Normal = as.numeric(Library.DNA.yield.ng.Normal))%>%
  tbl_summary(by=Metastasis,type = list(where(is.integer) ~ "continuous",where(is.numeric) ~ "continuous"),include = -Set.id)%>%
  add_n() %>% # add column with total number of non-missing observations
  add_p(list(all_continuous() ~ "paired.wilcox.test",all_categorical()~"mcnemar.test"),group=Set.id) %>% # test for a difference between groups
  gtsummary::as_tibble() %>% 
  writexl::write_xlsx(., file.path(dir_results_script,paste0("04_00_01_pat_characteristics_reduced_afterexclusion_",n_total,"_ALL_charact.xlsx")))

#-----------------------
# 2. Plot showing CP-score
#-----------------------
# Imputation:
# Choose the variables you want to impute
vars_to_impute <- c("Tumor.diameter", "PNI.or.LVI", "Tissue.involvement")   
aux_vars <-c("Sex","Age","Number.of.cSCC.before.culprit","Metastasis","FU.metastasis.years","Vital.status","Vital.follow.up.years","Differentiation","Breslow.thickness","IS.at.cSCC","Tumor.location","AJCC.8","BWH","Depth.of.Invasion","Resection.margin.cat","Tumor.budding","Mitotic.rate")
ignore_cols <- c("Set.id", "Sample.ID","Exclude_sample","HM.at.cSCC","OTR.at.cSCC")  # replace with your column names

# Subset only those variables 
imp_data <- sample_summary_table[, c(ignore_cols,vars_to_impute,aux_vars)]%>%
  mutate(PNI.or.LVI = ifelse(is.na(PNI.or.LVI),"No",as.character(PNI.or.LVI)))

# Create a method vector for mice
meth <- make.method(imp_data )

# Set ignored columns to "" (no imputation)
meth[ignore_cols] <- ""

# Run multiple imputation
imp <- mice(
  data = imp_data,
  m = 25,          # number of imputed datasets
  maxit = 20,      # number of iterations per imputation
  method = meth,  
  seed = 123
)

# Check convergence
plot(imp)

# Summarize variables:
imputed_long <- complete(imp, "long")

# We now get the mode:
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# get the imputed:
final_imputed <- imputed_long %>%
  filter(.imp != 0) %>%     # keep only imputed datasets
  group_by(.id) %>%         # combine across imputations
  summarize(
    across(
      .cols = -c(.imp),     # summarize all non-imputation columns
      .fns = list(
        final = ~ if (is.numeric(.x)) median(.x) else Mode(.x)
      ),
      .names = "{.col}"     # keep original column names
    )
  ) %>% 
  ungroup()%>%
  as.data.frame()

# add CP-score:
emc_risk_pct = function(x){
  age_dec = as.numeric(x[["Age"]])/10
  sex_male = ifelse(x[["Sex"]]=="Male",1,0)
  number_cscc = ifelse(as.numeric(x[["Number.of.cSCC.before.culprit"]])>4,4,as.numeric(x[["Number.of.cSCC.before.culprit"]]))
  tdiam = ifelse(as.numeric(x[["Tumor.diameter"]])>4,4,as.numeric(x[["Tumor.diameter"]]))
  sc_neck = ifelse(x[["Tumor.location"]]=="Scalp/neck",1,0)
  face = ifelse(x[["Tumor.location"]]=="Face",1,0)
  sc = ifelse(x[["Tissue.involvement"]]=="Subcutaneous fat",1,0)
  bsc =ifelse(x[["Tissue.involvement"]]=="Beyond subcutaneous fat",1,0)
  diff =ifelse(x[["Differentiation"]]=="Poor/undifferentiated",1,0)
  pli=ifelse(x[["PNI.or.LVI"]]=="Yes",1,0)
  lp = -1.86+0.25*(age_dec-7.47)+0.51*sex_male+0.57*(number_cscc-0.29)+0.54*(tdiam-1.20)-0.72*sc_neck+0.52*face+0.33*sc+1.40*bsc+1.3*diff+0.80*pli
  s_met=1-0.973^(exp(lp))
  return(s_met*100)
}

final_imputed$CP.score = apply(final_imputed,1,emc_risk_pct)

final_imputed = merge(final_imputed,sample_summary_table[,c("Sample.ID","Invasion.of.bones","PNI.bin")])
#Add AJCC staging:
ajcc_calc = final_imputed %>% 
  mutate(Invasion.of.bones = as.character(Invasion.of.bones),
    Invasion.of.bones=ifelse(is.na(Invasion.of.bones),"No",Invasion.of.bones),
    PNI.bin=as.character(PNI.bin),
         PNI.bin=ifelse(is.na(PNI.bin),"No",PNI.bin))%>%
  mutate(T1 = ifelse(Tumor.diameter<2 &((PNI.bin=="No"|(PNI.bin=="Yes"&Tissue.involvement=="Dermis"))&Tissue.involvement!="Beyond subcutaneous fat"&Depth.of.Invasion<=6),1,0),
                 T2 = ifelse(Tumor.diameter>=2 & Tumor.diameter<4 & ((PNI.bin=="No"|(PNI.bin=="Yes"&Tissue.involvement=="Dermis"))&Tissue.involvement!="Beyond subcutaneous fat"&Depth.of.Invasion<=6),1,0),
                 T3 = ifelse((Tumor.diameter>=4|(PNI.bin=="Yes"&Tissue.involvement!="Dermis")|Tissue.involvement=="Beyond subcutaneous fat"|Depth.of.Invasion>6)&Invasion.of.bones=="No",1,0),
                 T4 = ifelse(Invasion.of.bones=="Yes",1,0))
final_imputed$AJCC.8 = ifelse(ajcc_calc$T4==1,"T4",ifelse(ajcc_calc$T3==1,"T3",ifelse(ajcc_calc$T2==1,"T2",ifelse(ajcc_calc$T1==1,"T1",NA)))) #This corrects mistakes done during the imputation

#Add BWH staging:
bwh_calc = final_imputed %>%
  mutate(Invasion.of.bones = as.character(Invasion.of.bones),
         Invasion.of.bones=ifelse(is.na(Invasion.of.bones),"No",Invasion.of.bones))%>%
  mutate(HR_Tumordiam = ifelse(Tumor.diameter>=2,1,0),
         HR_Diff = ifelse( Differentiation=="Poor/undifferentiated",1,0),
         PNI.bin =ifelse(is.na(PNI.bin),"No",as.character(PNI.bin)),
         HR_PNI = ifelse(PNI.bin=="Yes",1,0),
         HR_Tinv = ifelse(Tissue.involvement=="Beyond subcutaneous fat",1,0),
         Sum_HR = HR_Tumordiam+HR_Diff+HR_PNI+HR_Tinv,
         Sum_HR = ifelse(!(is.na(Invasion.of.bones))&Invasion.of.bones=="Yes",4,Sum_HR))


final_imputed$BWH=ifelse(bwh_calc$Sum_HR==0,"T1",ifelse(bwh_calc$Sum_HR==1,"T2a",ifelse(bwh_calc$Sum_HR%in%c(2,3),"T2b",ifelse(bwh_calc$Sum_HR>=4,"T3",NA)))) #This corrects mistakes done during the imputation

# Save final imputed dataset:
write.xlsx(final_imputed, file.path(dir_results_script, paste0("04_00_01_imputed_dataset.xlsx")))

# Check that cp scores for samples without any missing data are the same:
stopifnot(round(final_imputed$CP.score[!is.na(sample_summary_table$CP.score)],3)==round(sample_summary_table$CP.score[!is.na(sample_summary_table$CP.score)],3))
stopifnot(sum(!sample_summary_table$BWH[!is.na(sample_summary_table$BWH)]==final_imputed$BWH[!is.na(sample_summary_table$BWH)])==0)
# Make sure each set.id has exactly 2 observations (one per group)

compute_test = function(df,variable,test_type){
  df_wide <-df%>%
    select("Set.id", "Metastasis", variable)%>% 
    reshape(idvar = "Set.id",
            timevar = "Metastasis",
            direction = "wide")
  
  # Perform test
  if(test_type=="paired_wilcoxon"){
    test_result = wilcox.test(df_wide[,2], df_wide[,3], paired = TRUE)
    
  }else if(test_type=="paired_mcnemar"){
    test_result  = mcnemar.test(table(df_wide[,2], df_wide[,3]))
  }else if(test_type=="wilcoxon"){
    test_result = wilcox.test(df_wide[,2], df_wide[,3], paired = FALSE)
  }else if(test_type=="chisq"){
    test_result = chisq.test(table(df[,"Metastasis"],df[,variable]),correct=FALSE)}
  else if(test_type=="fisher"){
    test_result = fisher.test(table(df[,"Metastasis"],df[,variable]))}
  return(test_result$p.value)
}

df_wide <- final_imputed%>%
  filter(Exclude_sample=="No")%>%
  select("Set.id", "Metastasis", "CP.score")%>% 
  reshape(idvar = "Set.id",
                   timevar = "Metastasis",
                   direction = "wide")

wil_test_CP = wilcox.test(df_wide$CP.score.Control, df_wide$CP.score.Case, paired = TRUE)
boxplot(df_wide$CP.score.Case-df_wide$CP.score.Control,ylab="EMC CP risk Case-Control")
abline(h=10,lty=2)
abline(h=-10,lty=2)


#sample_summary_table %>% 
final_imputed%>%
  filter(Exclude_sample=="No")%>%
  ggplot(aes(x = Metastasis, y = log10(CP.score))) +geom_boxplot(aes(alpha=0.2))+
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.6, fatten = 0, color = "black") +
  theme_bw()+
  annotate("text",
           x = mean(c(1,2)),   # middle of last two groups (IS = Yes)
           y = 2.12,
           label = paste0("p = ", round(wil_test_CP$p.value,2)),
           size = 3.5)+
  annotate("segment", x = 1, xend = 2, y = 2.05, yend = 2.05) +      # horizontal line
  annotate("segment", x = 1, xend = 1, y = 2.00, yend = 2.05) +      # left vertical tick
  annotate("segment", x = 2, xend = 2, y = 2.00, yend = 2.05) +      # right vertical tick
  ylab("log10 (EMC model risk %)")



# Characterize cohort with and without imputation, using paired tests+added unpaired tests as sanity check:
for (paired_samples_only in c(TRUE,FALSE)){
  for (imputed in c(TRUE,FALSE)){
  if(imputed){
    base_table = final_imputed
    imp_info = "with_imp"
  }else{
    base_table = sample_summary_table
    imp_info = "without_imp"
  }
  
  base_table <-base_table %>%
    filter(Exclude_sample == "No")
  if(paired_samples_only){
    paired_sets = base_table%>%pull(Set.id)%>%table()
    paired_sets = names(paired_sets)[paired_sets==2]
    base_table <- base_table %>%filter(Set.id%in%paired_sets)
  }
  
  n_total =  base_table%>%nrow()
  
  base_table <- base_table %>%
    select(Set.id, Metastasis, FU.metastasis.years, Vital.status, 
           Vital.follow.up.years, Sex, Age, Number.of.cSCC.before.culprit,
           Tumor.diameter, Differentiation, Breslow.thickness, PNI.or.LVI,
           Tissue.involvement, HM.at.cSCC, OTR.at.cSCC, IS.at.cSCC,
           Tumor.location, CP.score, BWH, AJCC.8) %>%
    mutate(BWH_noT3 = factor(ifelse(BWH=="T3",NA,as.character(BWH)),levels=c("T1","T2a","T2b")),
           AJCC_noT4 =factor(ifelse(AJCC.8=="T4",NA,as.character(AJCC.8)),levels=c("T1","T2","T3")))%>%
    tbl_summary(
      by = Metastasis,
      type = list(
        where(is.integer) ~ "continuous",
        where(is.numeric) ~ "continuous"
      ),
      include = -Set.id
    ) %>%
    add_n()
  
  # Add paired tests
  paired_table <- base_table %>%
    add_p(
      list(
        all_continuous() ~ "paired.wilcox.test",
        all_categorical() ~ "mcnemar.test"
      ),pvalue_fun = ~style_pvalue(., digits = 2),
      group = Set.id
    ) %>%
    modify_header(p.value ~ "**Paired p-value**")
  
  # Add unpaired tests
  unpaired_table <- base_table %>%
    add_p(pvalue_fun = ~style_pvalue(., digits = 2)) %>%
    modify_header(p.value ~ "**Unpaired p-value**")
  
  # Merge the two tables
  merged_table <- tbl_merge(
    tbls = list(paired_table, unpaired_table),
    tab_spanner = c("**Paired Tests**", "**Unpaired Tests**")
  )
  
  # Export
  merged_table %>%
    gtsummary::as_tibble() %>%
    as.data.frame() %>%
    .[, -c(6:8)] %>%
    writexl::write_xlsx(
      .,
      file.path(dir_results_script, 
                paste0("04_00_01_pat_characteristics_reduced_afterexclusion_", n_total, "_",imp_info,".xlsx"))
    )
  
}
}

#Save types of tests:
# Extract paired test info
paired_test_info <- paired_table$table_body %>% 
  filter(row_type == "label") %>% 
  select(variable, label, test_name, p.value) %>%
  rename(paired_test = test_name, paired_p.value = p.value)

# Extract unpaired test info
unpaired_test_info <- unpaired_table$table_body %>% 
  filter(row_type == "label") %>% 
  select(variable, label, test_name, p.value) %>%
  rename(unpaired_test = test_name, unpaired_p.value = p.value)

# Combine both
combined_test_info <- paired_test_info %>%
  full_join(unpaired_test_info, by = c("variable", "label"))

# Save to Excel
writexl::write_xlsx(combined_test_info, 
                    file.path(dir_results_script, "04_00_01_statistical_tests_used.xlsx"))

# Quickly check if there is unbias:
final_imputed$imputed = ifelse(is.na(sample_summary_table$AJCC.8),"imputed","not")
final_imputed%>%ggplot(aes(x=interaction(AJCC.8,imputed),y=Breslow.thickness,fill=Metastasis))+geom_boxplot()
#-----------------------
# 3. Compare AJCC of Nationwide cohort with Cases and Controls
#-----------------------
# As reported in Steijlen et al 2026.
df_nation_wide_cohort <- data.frame(
  AJCC_8 = c( "T1", "T2", "T3", "T4", "Unknown"),
  Nationwide = c(14826, 1037, 95, 89, 727+2346)
  #Percent = c(12, 78, 5.4, 0.5, 0.5, 3.8)
)

# Define a blue gradient for T1 → T4 for the plot:
blue_palette <- c( "Unknown"="gray",
                   "T1" = "#add8e6",  # light blue
                  "T2" = "#6ca0dc",  # medium
                  "T3" = "#1f4e79",  # darker
                  "T4" = "#0b2f5a")  # darkest blue

fonts_plot = theme(
  
  axis.title.x = element_text(size = 20),  # X-axis label only
  axis.title.y = element_text(size = 20),   # Y-axis label only
  
  axis.text.x = element_text(size = 19),  # X-axis ticks only
  axis.text.y = element_text(size = 19)   # Y-axis ticks only
)

plot_list <- list()
plot_counter <- 1
for(remove_NA in c( TRUE, FALSE)){
  for(exclusion in c("no","yes")){
    for(inputation in c("not_imputed","imputed")){
      if(inputation =="imputed"){
        df_wes_cohort = final_imputed
      }else if(inputation =="not_imputed"){
        df_wes_cohort = sample_summary_table
      }
      
      if(exclusion =="yes"){
        df_wes_cohort = df_wes_cohort%>%filter(Exclude_sample=="No")%>%as.data.frame()
      }
      
      # Get the AJCC8 distribution on Cases and Controls:
      df_ncc_wes = as.data.frame(t(pivot_wider(data.frame(table(df_wes_cohort$AJCC.8,df_wes_cohort$Metastasis)),names_from= "Var1",values_from="Freq")))
      colnames(df_ncc_wes)=df_ncc_wes [1,]
      df_ncc_wes  =df_ncc_wes [-1,]
      df_ncc_wes["Unknown",] = c(sum(is.na(df_wes_cohort$AJCC.8[which(df_wes_cohort$Metastasis==colnames(df_ncc_wes)[1])])),
                                 sum(is.na(df_wes_cohort$AJCC.8[which(df_wes_cohort$Metastasis==colnames(df_ncc_wes)[2])])))
      
      # Join the two tables
      df=cbind(df_nation_wide_cohort,df_ncc_wes)%>%
        mutate(Control = as.numeric(Control),
               Case= as.numeric(Case))%>%
        as.data.frame()
      
      if(remove_NA){
        df = df%>%filter(AJCC_8!="Unknown")%>%as.data.frame()
        rem_NA="yes"
      }else{
        rem_NA="no"
      }
      
      # Convert to long format
      df_long <- df %>%
        pivot_longer(cols = c(Nationwide, Control, Case),
                     names_to = "Cohort",
                     values_to = "Count")%>%
        mutate(Cohort =factor(Cohort,levels=c("Nationwide","Control","Case")),
               AJCC_8 = factor(AJCC_8, levels = c("Unknown","T4", "T3", "T2", "T1")))
      
      # Compute percentage **within each group**
      df_long <- df_long %>%
        group_by(Cohort) %>%
        mutate(Percent = Count / sum(Count))
      
      
      # Plot
      plot_list[[plot_counter]] <- ggplot(df_long, aes(x = Cohort, y = Percent, fill = AJCC_8)) +
              geom_bar(stat = "identity") +
              scale_y_continuous(labels = percent_format()) +
              scale_fill_manual(values = blue_palette) +
              ylab("AJCC stage (%)") +
              xlab("Cohort") +
        theme_bw()+  theme(
          
            axis.title.x = element_text(size = 20),  # X-axis label only
            axis.title.y = element_text(size = 20),   # Y-axis label only
          
          axis.text.x = element_text(size = 19),  # X-axis label only
          axis.text.y = element_text(size = 19),
          legend.title = element_text(size = 18),  
          legend.text = element_text(size = 14) 
        ) +ggtitle(paste0("Settings: imp: ",inputation, " exc samples: ",exclusion," remove NA ",rem_NA))
      plot_counter <- plot_counter + 1
    }
  }
  
}
grid.arrange(grobs = plot_list, ncol = 2)
ggsave(filename = file.path(dir_results_script, "04_00_01_barplot_nationwide_AJCC_vs_NCC_all.pdf"), 
               plot = marrangeGrob(plot_list, nrow = 2, ncol = 4),
               width = 40, height = 15)
ggsave(filename = file.path(dir_results_script, "04_00_01_barplot_nationwide_AJCC_vs_NCC_chosen_subset.pdf"), 
       plot = marrangeGrob(plot_list[c(7,8)], nrow = 1, ncol = 2),
       width = 15, height = 5)

# Check individual variables:



#for(test_type in c("paired","unpaired")){
for (imputed in c(TRUE,FALSE)){
    if(imputed){
      base_table = final_imputed
      imp_info = "with_imp"
      
      
      
    }else{
      base_table = sample_summary_table
      imp_info = "without_imp"
    }
    
    # Save source data:
  
    
    print("Excluded samples are excluded")
    base_table <-base_table %>%
      filter(Exclude_sample == "No")
    
    # Source data:
    source_df = base_table%>%select(Metastasis, Sex,Age, Number.of.cSCC.before.culprit,Tumor.location,Tumor.diameter,Tissue.involvement,PNI.or.LVI,Breslow.thickness,CP.score,AJCC.8,BWH)
    write.xlsx(source_df, file.path(dir_results_script, paste0("04_00_01_data_source_",imp_info,"utation.xlsx")))
    
    
    #if(test_type=="paired"){
      p_value_tdiam <- compute_test(base_table, "Tumor.diameter", "paired_wilcoxon")
      p_value_bt <- compute_test(base_table, "Breslow.thickness", "paired_wilcoxon")
      p_value_doi <- compute_test(base_table, "Depth.of.Invasion", "paired_wilcoxon")
      p_value_age <- compute_test(base_table, "Age", "paired_wilcoxon")
      p_value_ncscc <- compute_test(base_table, "Number.of.cSCC.before.culprit", "paired_wilcoxon")
      p_value_diff<- compute_test(base_table, "Differentiation", "paired_mcnemar")
      p_value_tinv<- compute_test(base_table, "Tissue.involvement", "paired_mcnemar")
      p_value_sex<- compute_test(base_table, "Sex", "paired_mcnemar")
      p_value_pnilvi<- compute_test(base_table, "PNI.or.LVI", "paired_mcnemar")
      p_value_tloc<- compute_test(base_table, "Tumor.location", "paired_mcnemar")
      p_value_CP.score<-compute_test(base_table, "CP.score", "paired_wilcoxon")
    #}else if(test_type=="unpaired"){
      p_value_tdiam_unpaired <- compute_test(base_table, "Tumor.diameter", "wilcoxon")
      p_value_bt_unpaired <- compute_test(base_table, "Breslow.thickness", "wilcoxon")
      p_value_doi_unpaired <- compute_test(base_table, "Depth.of.Invasion", "wilcoxon")
      p_value_age_unpaired <- compute_test(base_table, "Age", "wilcoxon")
      p_value_ncscc_unpaired <- compute_test(base_table, "Number.of.cSCC.before.culprit", "wilcoxon")
      p_value_diff_unpaired<- compute_test(base_table, "Differentiation", "chisq")
      p_value_tinv_unpaired<- compute_test(base_table, "Tissue.involvement", "chisq")
      p_value_sex_unpaired<- compute_test(base_table, "Sex", "chisq")
      p_value_pnilvi_unpaired<- compute_test(base_table, "PNI.or.LVI", "chisq")
      p_value_tloc_unpaired<- compute_test(base_table, "Tumor.location", "chisq")
      p_value_CP.score_unpaired<-compute_test(base_table, "CP.score", "wilcoxon")
      
    #}
    
    # Tumor diameter: 
    max_tdiam <- max(base_table$Tumor.diameter, na.rm = TRUE)
    p_tdiam=base_table%>%
      ggplot(aes(x = Metastasis, y = Tumor.diameter)) +geom_boxplot(outliers = FALSE)+
      geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
      stat_summary(fun = median, geom = "crossbar",
                   width = 0.6, fatten = 0, color = "black") +
      theme_bw()+fonts_plot+
      geom_signif(
        comparisons = list(c("Control", "Case")),
        annotations = paste0(
          "Paired p-val = ", format.pval(p_value_tdiam, digits = 2), "\n",
          "Unpaired p-val = ", format.pval(p_value_tdiam_unpaired, digits = 2)
        ),
        textsize = 5,
        vjust = -0.1,
        y_position = max_tdiam * 1.0
      ) +
      coord_cartesian(ylim = c(NA, max_tdiam * 1.15))+
      ylab(gsub("\\.", " ", "Tumor.diameter"))+xlab("")
    
    # Breslow Thickness:
    max_bt = max(base_table$Breslow.thickness, na.rm = TRUE)
    p_bt=base_table%>%
      ggplot(aes(x = Metastasis, y = Breslow.thickness)) +geom_boxplot(outliers = FALSE)+
      geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
      stat_summary(fun = median, geom = "crossbar",
                   width = 0.6, fatten = 0, color = "black") +
      theme_bw()+fonts_plot+
      geom_signif(
        comparisons = list(c("Control", "Case")),
        annotations = paste0(
          "Paired p-val = ", format.pval(p_value_bt, digits = 2), "\n",
          "Unpaired p-val = ", format.pval(p_value_bt_unpaired, digits = 2)
        ),
        textsize = 5,
        vjust = -0.1,
        y_position = max_bt * 1
      ) +
      coord_cartesian(ylim = c(NA,  max_bt* 1.15))+
      ylab(gsub("\\.", " ", "Breslow.Thickness"))+xlab("")
    
    # Depth of invasion:
    #max_doi = max(base_table$Depth.of.Invasion, na.rm = TRUE)
    #p_doi=base_table%>%
    #  ggplot(aes(x = Metastasis, y = Depth.of.Invasion)) +geom_boxplot(outliers = FALSE)+
    #  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
    #  stat_summary(fun = median, geom = "crossbar",
    #               width = 0.6, fatten = 0, color = "black") +
    #  theme_bw()+fonts_plot+
    #  geom_signif(
    #    comparisons = list(c("Control", "Case")),
    #    annotations = paste0(
    #      "Paired p-val = ", format.pval(p_value_doi, digits = 2), "\n",
    #      "Unpaired p-val = ", format.pval(p_value_doi_unpaired, digits = 2)
    #    ),
    #    textsize = 5,
    #    vjust = -0.1,
    #    y_position = max_doi * 1
    #  ) +
    #  coord_cartesian(ylim = c(NA, max_doi * 1.15))+
    #  ylab(gsub("\\.", " ", "Depth.of.Invasion"))+xlab("")
    
    # CP score:
    max_cp_score = max(base_table$CP.score, na.rm = TRUE)
    p_cp=base_table%>%
      ggplot(aes(x = Metastasis, y = CP.score)) +geom_boxplot(outliers = FALSE)+
      geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
      stat_summary(fun = median, geom = "crossbar",
                   width = 0.6, fatten = 0, color = "black") +
      theme_bw()+fonts_plot+
      geom_signif(
        comparisons = list(c("Control", "Case")),
        annotations = paste0(
          "Paired p-val = ", format.pval(p_value_CP.score, digits = 2), "\n",
          "Unpaired p-val = ", format.pval(p_value_CP.score_unpaired, digits = 2)
        ),
        textsize = 5,
        vjust = -0.1,
        y_position = log10(max_cp_score) * 1
      ) +scale_y_continuous(trans = "log10") +  
      coord_cartesian(ylim = c(NA,  max_cp_score* 3))+
      ylab("EMC model risk % (log scale)")+xlab("")
    
    # Plot of the differences:
    
    p_cp_diff=base_table%>%
      select("Set.id", "Metastasis", "CP.score")%>% 
      reshape(idvar = "Set.id",
              timevar = "Metastasis",
              direction = "wide")%>%
      mutate(Difference_risk_Case_Control = CP.score.Case - CP.score.Control )%>%
      ggplot(aes(y=Difference_risk_Case_Control))+geom_boxplot()+theme_bw()+xlab("")+
       fonts_plot+ylab("EMC risk difference\n(cases - controls)")+theme(
         axis.text.x = element_text(color = "white"),  # Hide but keep space
         axis.ticks.x = element_line(color = "white")  # Hide ticks but keep space
       )+  annotate("text",
                    x = 0,
                    y = max_cp_score * 1.05,
                    label = paste0("Paired p-val = ", format.pval(p_value_CP.score, digits = 2)),
                    size = 5) 
      

    
    p_diff=base_table%>%
      count(Metastasis, Differentiation) %>%
      group_by(Metastasis) %>%
      mutate(percentage = n / sum(n)) %>%
      ggplot(aes(x = Metastasis, y = percentage, fill = Differentiation)) +
      geom_col(position = "fill") +
      scale_y_continuous(
        labels = scales::percent_format(),
        breaks = c(0, 0.25, 0.5, 0.75, 1)  # Set breaks at 0%, 25%, 50%, 75%, 100%
      ) +
      ylab("Percentage") +
      theme_bw()+fonts_plot+
      geom_signif(
        comparisons = list(c("Control", "Case")),
        annotations = paste0(
          "Paired p-val = ", format.pval(p_value_diff, digits = 2), "\n",
          "Unpaired p-val = ", format.pval(p_value_diff_unpaired, digits = 2)
        ),
        textsize = 5,
        vjust = -0.2,
        y_position = 1
      )+
      coord_cartesian(ylim = c(0, 1.2)) +theme(
        legend.position = "top",  # Move legend to top,
        legend.title = element_text(size = 14),  # Increase legend title size
        legend.text = element_text(size = 12)  
      )+xlab("")
    
    p_tinv=base_table%>%
      count(Metastasis, Tissue.involvement) %>%
      group_by(Metastasis) %>%
      mutate(percentage = n / sum(n)) %>%
      ggplot(aes(x = Metastasis, y = percentage, fill =Tissue.involvement)) +
      geom_col(position = "fill") +
      scale_y_continuous(labels = scales::percent_format()) +
      ylab("Percentage") +
      theme_bw()+fonts_plot+fonts_plot+
      geom_signif(
        comparisons = list(c("Control", "Case")),
        annotations = paste0(
          "Paired p-val = ", format.pval(p_value_tinv, digits = 2), "\n",
          "Unpaired p-val = ", format.pval(p_value_tinv_unpaired, digits = 2)
        ),
        textsize = 5,
        vjust = -0.2,
        y_position = 1
      )+
      coord_cartesian(ylim = c(0, 1.2)) +
      scale_y_continuous(
        labels = scales::percent_format(),
        breaks = c(0, 0.25, 0.5, 0.75, 1)  # Set breaks at 0%, 25%, 50%, 75%, 100%
      )+theme(
        legend.position = "top",  # Move legend to top,
        legend.title = element_text(size = 14),  # Increase legend title size
        legend.text = element_text(size = 12)  
      )+xlab("")
    p_sex=base_table%>%
      count(Metastasis, Sex) %>%
      group_by(Metastasis) %>%
      mutate(percentage = n / sum(n)) %>%
      ggplot(aes(x = Metastasis, y = percentage, fill =Sex)) +
      geom_col(position = "fill") +
      scale_y_continuous(labels = scales::percent_format()) +
      ylab("Percentage") +
      theme_bw()+fonts_plot+fonts_plot+
      geom_signif(
        comparisons = list(c("Control", "Case")),
        annotations = paste0(
          "Paired p-val = ", format.pval(p_value_sex, digits = 2), "\n",
          "Unpaired p-val = ", format.pval(p_value_sex_unpaired, digits = 2)
        ),
        textsize = 5,
        vjust = -0.2,
        y_position = 1
      ) +
      coord_cartesian(ylim = c(0, 1.2))+
      scale_y_continuous(
        labels = scales::percent_format(),
        breaks = c(0, 0.25, 0.5, 0.75, 1)  # Set breaks at 0%, 25%, 50%, 75%, 100%
      )+theme(
        legend.position = "top",  # Move legend to top,
        legend.title = element_text(size = 14),  # Increase legend title size
        legend.text = element_text(size = 12)  
      )+xlab("")
    
    #Age
    max_age=max(base_table$Age, na.rm = TRUE)
    p_age=base_table%>%
      ggplot(aes(x = Metastasis, y = Age)) +geom_boxplot(outliers = FALSE)+
      geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
      stat_summary(fun = median, geom = "crossbar",
                   width = 0.6, fatten = 0, color = "black") +
      theme_bw()+fonts_plot+
      geom_signif(
        comparisons = list(c("Control", "Case")),
        annotations = paste0(
          "Paired p-val = ", format.pval(p_value_age, digits = 2), "\n",
          "Unpaired p-val = ", format.pval(p_value_age_unpaired, digits = 2)
        ),
        textsize = 5,
        vjust = -0.1,
        y_position = max_age * 1
      ) +coord_cartesian(ylim = c(NA, max_age * 1.1))+xlab("")
    
    # Number of cSCC prior:
    max_ncscc = max(base_table$Number.of.cSCC.before.culprit, na.rm = TRUE)
    p_ncscc=base_table%>%
      ggplot(aes(x = Metastasis, y = Number.of.cSCC.before.culprit)) +geom_boxplot(outliers = FALSE)+
      geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
      stat_summary(fun = median, geom = "crossbar",
                   width = 0.6, fatten = 0, color = "black") +
      theme_bw()+fonts_plot+
      geom_signif(
        comparisons = list(c("Control", "Case")),
        annotations = paste0(
          "Paired p-val = ", format.pval(p_value_ncscc, digits = 2), "\n",
          "Unpaired p-val = ", format.pval(p_value_ncscc_unpaired, digits = 2)
        ),
        textsize = 5,
        vjust = -0.1,
        y_position = max_ncscc * 1
      )+coord_cartesian(ylim = c(NA, max_ncscc  * 1.15))+
      ylab(gsub("\\.", " ", "Number.of.cSCC.before.culprit"))+xlab("")
    
    p_pnilvi=base_table%>%
      count(Metastasis, PNI.or.LVI) %>%
      group_by(Metastasis) %>%
      mutate(percentage = n / sum(n)) %>%
      ggplot(aes(x = Metastasis, y = percentage, fill =PNI.or.LVI)) +
      geom_col(position = "fill") +
      scale_y_continuous(labels = scales::percent_format()) +
      ylab("Percentage") +
      theme_bw()+fonts_plot+fonts_plot+
      geom_signif(
        comparisons = list(c("Control", "Case")),
        annotations = paste0(
          "Paired p-val = ", format.pval(p_value_pnilvi, digits = 2), "\n",
          "Unpaired p-val = ", format.pval(p_value_pnilvi_unpaired, digits = 2)
        ),
        textsize = 5,
        vjust = -0.2,
        y_position = 1
      ) +
      coord_cartesian(ylim = c(0, 1.2))+
      scale_y_continuous(
        labels = scales::percent_format(),
        breaks = c(0, 0.25, 0.5, 0.75, 1)  # Set breaks at 0%, 25%, 50%, 75%, 100%
      )+theme(
        legend.position = "top",
        legend.title = element_text(size = 14),  # Increase legend title size
        legend.text = element_text(size = 12)  # Move legend to top
      )+xlab("")
    
    p_tloc=base_table%>%
      count(Metastasis, Tumor.location) %>%
      group_by(Metastasis) %>%
      mutate(percentage = n / sum(n)) %>%
      ggplot(aes(x = Metastasis, y = percentage, fill =Tumor.location)) +
      geom_col(position = "fill") +
      scale_y_continuous(labels = scales::percent_format()) +
      ylab("Percentage") +
      theme_bw()+fonts_plot+fonts_plot+
      geom_signif(
        comparisons = list(c("Control", "Case")),
        annotations = paste0(
          "Paired p-val = ", format.pval(p_value_tloc, digits = 2), "\n",
          "Unpaired p-val = ", format.pval(p_value_tloc_unpaired, digits = 2)
        ),
        textsize = 5,
        vjust = -0.2,
        y_position = 1
      ) +
      coord_cartesian(ylim = c(0, 1.2))+
      scale_y_continuous(
        labels = scales::percent_format(),
        breaks = c(0, 0.25, 0.5, 0.75, 1)  # Set breaks at 0%, 25%, 50%, 75%, 100%
      )+theme(
        legend.position = "top",
        legend.title = element_text(size = 14),  # Increase legend title size
        legend.text = element_text(size = 12)  # Move legend to top
      )+
      ylab(gsub("\\.", " ", "Tumor.location"))+xlab("")
    
    
    # Define the layout matrix
    layout <- rbind(c(1, 2, 3, 4,5),      # First row: 4 plots
                    c(6, 7, 8, 9,NA),     # Second row: 3 plots
                    c(10,11, NA,NA,NA))  
    p_1 = grid.arrange(p_sex, p_age,  p_ncscc, p_tloc,p_tdiam,
                       p_diff, p_tinv, p_pnilvi, 
                        p_bt, p_cp,p_cp_diff,
                       layout_matrix = layout)
    
    ggsave(filename = file.path(dir_results_script, paste0("04_00_01_differences_CP_variables_",imp_info,"_both_tests.pdf")), 
          plot = p_1,
          width = 18, height = 15)
    
    ggsave(filename = file.path(dir_results_script, paste0("04_00_01_differences_CP_variables_",imp_info,"_both_tests_wide.pdf")), 
           plot = p_1,
           width = 40, height = 15)
    
  }
  
#}











