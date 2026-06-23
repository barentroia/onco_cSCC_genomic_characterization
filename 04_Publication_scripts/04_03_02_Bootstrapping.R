#---------------------------------------------------
# Aim: Bootstrapping
# Author: B. Rentroia Pacheco
# Input: Driver mutations
# Output: 
#---------------------------------------------------

# Helper functions:
source(file.path(dir_scripts,"04_Publication_scripts","04_03_02_A1_Bootstrapping_functions.R"))
#-----------------------
# 0. Input features
#-----------------------
# Mutation data:
mut_hits_df = mut_data_drivers %>%
select(Sample,Genes,gene_weight)%>%
  pivot_wider(id_cols = Sample,names_from= Genes,values_from=gene_weight)%>%
  as.data.frame()

# Pathway data:
pathway_hits_df = path_sample %>%
  select(Sample,Pathway,Weighted_burden)%>%
  pivot_wider(id_cols = Sample,names_from= Pathway,values_from=Weighted_burden)%>%
  as.data.frame()

# Copy number data:
cn_split_df = copy_number_split  %>%select(Sample,paste(chr_cols ,"Gain", sep = " - "),paste(chr_cols ,"Loss", sep = " - "))%>%
  as.data.frame()
new_names <- sapply(colnames(cn_split_df), function(x) {
  x <- gsub(" - ", "_", x)
  return(x)
})
colnames(cn_split_df) <- new_names

# Auxilliary variables:
auxilliary_vars = c("Percent_Bait_with_callable_coverage","BWH")

# Add auxilliary variables to dataframes:
model_mut_surv_df = merge(mut_hits_df,sample_summary[,c("Sample.ID", "Metastasis","FU.metastasis.years",auxilliary_vars)],all.x=TRUE,by.x="Sample",by.y="Sample.ID")
model_mut_path_surv_df = merge(pathway_hits_df,sample_summary[,c("Sample.ID", "Metastasis","FU.metastasis.years",auxilliary_vars)],all.x=TRUE,by.x="Sample",by.y="Sample.ID")
model_cnv_surv_df = merge(cn_split_df,sample_summary[,c("Sample.ID", "Metastasis","FU.metastasis.years",auxilliary_vars)],all.x=TRUE,by.x="Sample",by.y="Sample.ID")

#-----------------------
# 1. Bootstrapping
#-----------------------

# Input_features
data_types = c("pathways","cnv")
input_df = list()
if("genes"%in%data_types){
  input_df[["genes"]] = model_mut_surv_df}
if("pathways"%in%data_types){
  input_df[["pathways"]] = model_mut_path_surv_df}
if("cnv"%in%data_types){
  input_df[["cnv"]] = model_cnv_surv_df}

# prefiltering:
pre_filtering = function(df,cutoff_mut,auxilliary_vars,remove_sex){
  
  if(!is.data.frame(df)){
    
    # Find samples that are common accross all samples:
    sample_ids_bs = c()
    for (nr_df in c(1:length(df))){
      sample_ids_bs = c(sample_ids_bs,df[[nr_df]][,"Sample"])
    }
    sample_ids_bs = names(table(sample_ids_bs))[table(sample_ids_bs)==length(df)]
    
    df_y = df[[1]]
    rownames(df_y) = df[[1]][,"Sample"]
    df_y$Metastasis_bin <- ifelse(df_y$Metastasis == "Case", 1, 0)
    df_y = df_y[ sample_ids_bs,]
    # Define response (Surv object)
    y <- Surv(time = df_y$FU.metastasis.years,
              event = df_y$Metastasis_bin)
    
    #Remove the auxilliary variables from the dataframe: 
    gene_cols_all=c()
    for(i in 1:length(df)){
      
      # Only keep the Samples that are common to all dataframes:
      rownames(df[[i]]) = df[[i]]$Sample
      df[[i]] = df[[i]][sample_ids_bs,]
      
      gene_cols = setdiff(colnames(df[[i]]),
                          c("Sample", "Metastasis", "FU.metastasis.years", "Metastasis_bin",auxilliary_vars))
      gene_cols_all <- c(gene_cols_all,gene_cols)
      
      if(length(gene_cols)==1){
        df[[i]] = setNames(data.frame(df[[i]][, gene_cols]), gene_cols)
      }else{
        df[[i]] = as.data.frame(df[[i]][,gene_cols])
      }
      
      # Cutoff on minimum mutations:
      df[[i]] = df[[i]][,colSums(abs(df[[i]])>0)>cutoff_mut[i]]
      
      if(remove_sex){
        #print(colnames(df))
        cols_keep =setdiff(colnames(df[[i]]),colnames(df[[i]])[grepl("chrX|chrY",colnames(df[[i]]))])
        df[[i]] = df[[i]][,cols_keep]
        #print(colnames(df))
      }
    }
  

    df =do.call(cbind, unname(df))
    
    
  }else{
    y <- Surv(time = df$FU.metastasis.years,
              event = df$Metastasis_bin)
    # Select only gene columns (everything between Tumor_Sample_Barcode and Metastasis)
    gene_cols_all <- setdiff(colnames(df),
                             c("Sample", "Metastasis", "FU.metastasis.years", "Metastasis_bin",auxilliary_vars))
    # Select only input variables:
    df = df[, gene_cols_all]
    # Cutoff on minimum mutations:
    df = df[,abs(colSums(df,na.rm=TRUE))>cutoff_mut]
    
    if(remove_sex){
      #print(colnames(df))
      cols_keep =setdiff(colnames(df),colnames(df)[grepl("chrX|chrY",colnames(df))])
      df = df[,cols_keep]
      #print(colnames(df))
    }
  }

  

  
  x <- as.matrix(df)
  x <-x[,unique(colnames(x))]
  return(list(x=x,y=y))
}

# Assessment of different pre-filtring approaches:
if(recompute_BS){
  for (ct_cnv_freq in c(9,25,27)){
    dir_results_script_cutoff = file.path(dir_results_script,paste0("CNV_freq_cutoff",as.character(ct_cnv_freq)))
    pre_processing_dfs = pre_filtering(input_df,c(2,ct_cnv_freq),auxilliary_vars,TRUE)
    x = pre_processing_dfs[["x"]]
    y = pre_processing_dfs[["y"]]
    
    model_list <- list(
      coxnet_a1 = list(method="coxnet","s" = "lambda.1se", "alpha" = 1, "fsel"="none"),
      coxnet_a1_fsel02 = list(method="coxnet","s" = "lambda.1se", "alpha" = 1, "fsel"="cox","p_thresh"=0.2),
      coxnet_a05 = list(method="coxnet","s" = "lambda.1se", "alpha" = 0.5, "fsel"="none"),
      coxnet_a05_fsel02 = list(method="coxnet","s" = "lambda.1se", "alpha" = 0.5, "fsel"="cox","p_thresh"=0.2),
      coxnet_02 = list(method="coxnet","s" = "lambda.1se", "alpha" = 0.2, "fsel"="none"),
      coxnet_a02_fsel02 = list(method="coxnet","s" = "lambda.1se", "alpha" = 0.2, "fsel"="cox","p_thresh"=0.2),
      onco_score = list(method="Onco_score","fsel"="cox","p_thresh"=0.1),
      onco_score_02 = list(method="Onco_score","fsel"="cox","p_thresh"=0.2))
    
    list_x = list("Genes_Pathways_CNV"=x)
    final_results_filtering <- run_bootstrap_all_models(list_x, y, model_list, n_boot = 100, outdir = dir_results_script_cutoff)
    saveRDS(final_results_filtering, file.path(dir_results_script_cutoff,"final_results_differentfiltering.RDS"))
    
  }
  
  # Check genetic algorithms:
  if(ct_cnv_freq==27){
    model_list <- list(
      coxnet_a05_fsel02 = list(method="coxnet","s" = "lambda.1se", "alpha" = 0.5, "fsel"="cox","p_thresh"=0.2),
      onco_score = list(method="Onco_score","fsel"="cox","p_thresh"=0.1),
      onco_score_02 = list(method="Onco_score","fsel"="cox","p_thresh"=0.2),
      onco_score_GA_no_val_noFPen_02 = list(method="Onco_score","fsel"="GA","GA.val"=FALSE,"GA.feat_penalty"=FALSE,"p_thresh"=0.2),
      onco_score_GA_no_val_wFPen_02 = list(method="Onco_score","fsel"="GA","GA.val"=FALSE,"GA.feat_penalty"=TRUE,"p_thresh"=0.2),
      onco_score_GA_with_val_noFPen_02 = list(method="Onco_score","fsel"="GA","GA.val"=TRUE,"GA.feat_penalty"=FALSE,"p_thresh"=0.2),
      onco_score_GA_with_val_wFPen_02 = list(method="Onco_score","fsel"="GA","GA.val"=TRUE,"GA.feat_penalty"=TRUE,"p_thresh"=0.2)
    )
    list_x = list("Genes_Pathways_CNV"=x)
    final_results_2 <- run_bootstrap_all_models(list_x, y, model_list, n_boot = 20, outdir = file.path(dir_results_script_cutoff,"GAs"))
    saveRDS(final_results_2 , file.path(dir_results_script_cutoff,"GAs","final_results_GA.RDS"))
    
  }
  
  # Check RForests:
  model_list <- list(
    coxnet_a05 = list(method="coxnet","s" = "lambda.1se", "alpha" = 0.5, "fsel"="none"),
    coxnet_a05_fsel02 = list(method="coxnet","s" = "lambda.1se", "alpha" = 0.5, "fsel"="cox","p_thresh"=0.2),
    onco_score = list(method="Onco_score","fsel"="cox","p_thresh"=0.1),
    onco_score_02 = list(method="Onco_score","fsel"="cox","p_thresh"=0.2),
    rsf_fsel_02 = list(method="rsf","num.trees" = 1000, "max.depth" = 4, "fsel"="cox","p_thresh"=0.2))
  list_x = list("Genes_Pathways_CNV"=x)
  final_results_3 <- run_bootstrap_all_models(list_x, y, model_list, n_boot = 100, outdir = file.path(dir_results_script_cutoff,"RFs"))
  saveRDS(final_results_3 , file.path(dir_results_script_cutoff,"RFs","final_results_RF.RDS"))
  
}
dir_results_script_cutoff = file.path(dir_results_script,paste0("CNV_freq_cutoff",as.character(27)))
final_results_filtering <- readRDS(file.path(dir_results_script_cutoff,"final_results_differentfiltering.RDS"))
# Save final coefficients:
final_models_coefs = list()
for(i in 1:length(final_results_filtering$results)){
  model_name = names(final_results_filtering$results)[i]
  final_models_coefs[[model_name]]= final_results_filtering$results[[model_name]]$model_training$selected_coefs
}
max_len <- max(lengths(final_models_coefs))

df <- as.data.frame(lapply(final_models_coefs, function(x) {
  length(x) <- max_len
  x
}))



# Define models and parameters:
#model_list <- list(
  #random_coxnet=list(method="coxnet","fsel"="none","s" = "lambda.min", "alpha" = 1, "fsel"="none"),
  #coxnet = list(method="coxnet","s" = "lambda.1se", "alpha" = 1, "fsel"="none","params_filtering" = c(cutoff_mut = c(2,5),auxilliary_vars=auxilliary_vars,remove_sex=TRUE)),
  #coxnet_02 = list(method="coxnet","s" = "lambda.1se", "alpha" = 0.2, "fsel"="none","params_filtering" = c(cutoff_mut = c(2,5),auxilliary_vars=auxilliary_vars,remove_sex=TRUE)),
  #coxnet_05 = list(method="coxnet","s" = "lambda.1se", "alpha" = 0.5, "fsel"="none","params_filtering" = c(cutoff_mut = c(2,5),auxilliary_vars=auxilliary_vars,remove_sex=TRUE)),
  #coxnet_fsel02 = list(method="coxnet","s" = "lambda.1se", "alpha" = 1, "fsel"="cox","p_thresh"=0.4,"params_filtering" = c(cutoff_mut = c(2,5),auxilliary_vars=auxilliary_vars,remove_sex=TRUE)),
  #coxnet_08 = list(method="coxnet","s" = "lambda.1se", "alpha" = 0.25, "fsel"="none","params_filtering" = c(cutoff_mut = c(2,5),auxilliary_vars=auxilliary_vars,remove_sex=TRUE)),
  #coxnet_1se = list(method="coxnet","s" = "lambda.1se", "alpha" = 0.1, "fsel"="none")
  #random_rsf=list(method="rsf","num.trees" = 1000, "max.depth" = 8, "fsel"="none"),
  #rsf = list(method="rsf","num.trees" = 1000, "max.depth" = 8, "fsel"="none"),
  #rsf_fsel_01 = list(method="rsf","num.trees" = 1000, "max.depth" = 5, "fsel"="cox","p_thresh"=0.1),
  #rsf_fsel_02 = list(method="rsf","num.trees" = 1000, "max.depth" = 5, "fsel"="cox","p_thresh"=0.2),
  #onco_score = list(method="Onco_score","fsel"="cox","p_thresh"=0.1),
  #onco_score_02 = list(method="Onco_score","fsel"="cox","p_thresh"=0.2),
  #onco_score_GA_noFPen_02 = list(method="Onco_score","fsel"="GA","GA.val"=TRUE,"GA.feat_penalty"=FALSE,"p_thresh"=0.2),
  #onco_score_GA_noVal_wFPen = list(method="Onco_score","fsel"="GA","GA.val"=FALSE,"GA.feat_penalty"=TRUE),
  #onco_score_GA_noFPen = list(method="Onco_score","fsel"="GA","GA.val"=TRUE,"GA.feat_penalty"=FALSE),
  #onco_score_GA_wFPen = list(method="Onco_score","fsel"="GA","GA.val"=TRUE,"GA.feat_penalty"=TRUE)
)

#list_x = list("Genes_Pathways_CNV"=x)
#final_results_2 <- run_bootstrap_all_models(list_x, y, model_list, n_boot = 20, outdir = dir_results_script)

#-----------------------
# 2. Final model:
#-----------------------
ct_cnv_freq =27 # this corresponds to a 20% cutoff (27/138)
pre_processing_dfs = pre_filtering(input_df,c(2,ct_cnv_freq),auxilliary_vars,TRUE)
x = pre_processing_dfs[["x"]]
y = pre_processing_dfs[["y"]]
list_x = list("Pathways_CNV"=x)

# Coxnet:
set.seed(1)
model_fit=  cv.glmnet(x = as.matrix(list_x[[1]]), y = y, family = "cox", alpha = 0.5,type.measure="C")
coef_full <- coef(model_fit, s = "lambda.1se")
selected_coefs <- rownames(coef_full)[as.numeric(coef_full) != 0]
surv_probs = predict_survival_coxnet(model_fit, "lambda.1se",y,x,x,5)

baseline_surv=summary(survfit(model_fit, s ="lambda.1se",x=as.matrix(list_x[[1]]), y = y),times=5)$surv
lp=predict(model_fit,scale(as.matrix(x),center=TRUE,scale=FALSE),type="link",s="lambda.1se") # Note: the link option does not provide centered lp. 
baseline_surv^exp(lp)-surv_probs$survival_prob # these are very close. I think the difference might come from rounding in the survival/hazard, because the differences are much smaller for t=1 or 3

# Save the model coefficients and baseline survival:
one_off_coxnet_coeffs=data.frame(Coefficient_names = c("Baseline_survival",selected_coefs),Coefficient_values=c(surv_probs$baseline_survival,coef_full[as.numeric(coef_full) != 0]))
write.xlsx(one_off_coxnet_coeffs,file.path(dir_results_script_cutoff,"Final_model","Model_coefficients.xlsx"))

# Chosen model:
chosen_model = "Genes_Pathways_CNV_coxnet_a05"
freq = final_results_filtering$results[[chosen_model]]$gene_selection_freq
n_boot=200
cutoff=0
# Performance:
# AJCC/BWH/CP_model
cp_data = final_imputed[final_imputed$Sample.ID%in%rownames(x),c("Sample.ID","Metastasis","FU.metastasis.years","CP.score","BWH","AJCC.8")]%>%
  mutate(Metastasis_bin = ifelse(Metastasis=="Case",1,0), 
         BWH_number = case_when(
           BWH == "T1"  ~ 1,
           BWH == "T2a" ~ 2,
           BWH == "T2b" ~ 3,
           BWH == "T3"  ~ 4
         ),
         AJCC_number = case_when(
           AJCC.8 == "T1"  ~ 1,
           AJCC.8 == "T2" ~ 2,
           AJCC.8 == "T3" ~ 3,
           AJCC.8 == "T4"  ~ 4
         ))
y_cp = Surv(cp_data$FU.metastasis.years,cp_data$Metastasis_bin)
1-concordance(y_cp ~ cp_data$CP.score)$concordance # 0.59 for CP score
1-concordance(y_cp ~ cp_data$BWH_number)$concordance # 0.65 for CP score
1-concordance(y_cp ~ cp_data$AJCC_number)$concordance # 0.63 for CP score


data_vec <- final_results_filtering$results[[chosen_model]]$C_632plus_all
df_box <- data.frame(value = data_vec)
q1  <- round(quantile(data_vec, 0.25),3)
med <- round(quantile(data_vec, 0.50),3)
q3  <- round(quantile(data_vec, 0.75),3)
p_boxplot = ggplot(df_box, aes(x = "", y = value)) +
  geom_boxplot(outlier.shape = NA, width = 0.4) +
  geom_jitter(width = 0.1, color = "black") +
  annotate("segment", x = 1.22, xend = 1.25, y = q1,  yend = q1,  color = "grey40") +
  annotate("segment", x = 1.22, xend = 1.25, y = med, yend = med, color = "grey40") +
  annotate("segment", x = 1.22, xend = 1.25, y = q3,  yend = q3,  color = "grey40") +
  annotate("text", x = 1.27, y = q1,  label = paste0("Q1 = ", round(q1,  2)), size = 5, hjust = 0, color = "grey20") +
  annotate("text", x = 1.27, y = med, label = paste0("Q2 = ", round(med, 2)), size = 5, hjust = 0, color = "grey20") +
  annotate("text", x = 1.27, y = q3,  label = paste0("Q3 = ", round(q3,  2)), size = 5, hjust = 0, color = "grey20") +
  labs(x = "", y = "C-index") +
  theme_bw(base_size = 14)

ggsave(file.path(dir_results_script_cutoff,"Final_model","C-index.pdf"),p_boxplot,width=5,height=8)

df_box %>%
  mutate(Bootstrap_rep=1:n_boot)%>%
  rename(Bootstrap_632p_estimate=value)%>%
  write.xlsx(file.path(dir_results_script_cutoff,"Final_model","C-index_source_data.xlsx"))

# Feature robustness:
df <- data.frame(
  model = chosen_model,
  feature = names(freq),
  freq = as.numeric(freq),
  prop = as.numeric(freq) / n_boot
) %>% # This is to make it possible to make the labels bold.
  left_join(
    one_off_coxnet_coeffs %>%
      select(Coefficient_names, Coefficient_values),
    by = c("feature" = "Coefficient_names")
  ) %>%
  rename(Coefficient = Coefficient_values)%>%
  mutate(
    direction = case_when(
      is.na(Coefficient) ~ "Absent",
      Coefficient > 0 ~ "Positive",
      Coefficient < 0 ~ "Negative",
      TRUE ~ "Zero"
    )
  )%>%
  mutate(feature=gsub("_"," ",feature))%>%
  mutate(feature_label = ifelse(feature %in% gsub("_"," ",selected_coefs),
                                paste0("**", feature, "**"),
                                feature)) %>%
  arrange(freq) %>%
  mutate(rank = row_number()) 

p <- ggplot(df, aes(x = rank, y = prop*100,fill=direction)) +
  geom_col() +
  coord_flip() +
  scale_x_continuous(breaks = df$rank, labels = df$feature_label) +
  scale_fill_manual(
    values = c(
      "Positive" = "red2",
      "Negative" = "dodgerblue3",
      "Absent" = "grey70",
      "Zero" = "black"
    )
  )+
  labs(
    title = paste0("Features Selected in ≥", cutoff*100, "% of Bootstraps"),
    x = "Feature",
    y = "Bootstrap inclusion frequency (%)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.y = ggtext::element_markdown(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    strip.text.x = element_text(size = 8),
    strip.text.y = element_text(size = 8)
  ) +
  ylim(0, 100)+guides(fill="none")
ggsave(file.path(dir_results_script_cutoff,"Final_model","Top_features.pdf"),p,width=6,height=10)

df%>%mutate(nboot = n_boot,selected_in_training_model = ifelse(feature %in% gsub("_"," ",selected_coefs),"Yes","No"),
            Boostrap_inclusion_freq_pct = prop*100)%>%select(feature,freq,nboot,Boostrap_inclusion_freq_pct ,selected_in_training_model )%>%
  write.xlsx(file.path(dir_results_script_cutoff,"Final_model","Top_features_sourcedata.xlsx"))


# Survival curves:
# Make groups:
make_risk_groups <- function(lp, fraction = 2) {
  # fraction = 2 → median split, 3 → tertiles, etc.
  # cut into quantiles
  q <- quantile(lp, probs = seq(0, 1, length.out = fraction + 1))
  
  # use cut to assign group labels
  groups <- cut(lp, breaks = q, include.lowest = TRUE, labels = FALSE)
  return(groups)
}

# Example: 2 groups (median split)
risk_group <- make_risk_groups(lp, fraction = 2)
# Survival curve for the full cohort:
fit_km <- survfit(y ~ risk_group)
ggsurvplot(
  fit_km,
  data = as.data.frame(list_x[[1]]),
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  palette = c("darkred", "#2E9FDF") # customize for up to 3 groups
)

# Compute most frequent risk group acorss bootstraps, to produce an internally validated survival curve:
# Compute risk cutoff on the bootstrap samples, and apply it to the oob samples:
mat <- final_results_filtering$results[[chosen_model]]$oob_predictions
mat_bs <- final_results_filtering$results[[chosen_model]]$bs_predictions
pct_oob <- matrix(NA,nrow=nrow(final_results_filtering$results[[chosen_model]]$oob_predictions),
                  ncol=ncol(final_results_filtering$results[[chosen_model]]$oob_predictions))
for(j in 1:ncol(pct_oob)){
  risk_cutoff = quantile(mat_bs[, j],0.5,na.rm=TRUE)
  #merged  =dplyr::coalesce(mat_bs[, j], mat[, j])
  if(length(unique(mat_bs[, j]))!=1){
    pct_oob[,j]=ifelse(!is.na(mat[,j])&mat[,j]>risk_cutoff,"low_risk","high_risk")
    pct_oob[is.na(mat[, j]),j]=NA
  }
}

# Extract the most frequent risk group for each sample on the OOB samples:
get_mode <- function(x) {
  x <- x[!is.na(x)]
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Check that we have a similar C-index to what is reported for OOB samples:
pct_oob_estimates = apply(pct_oob, 1, function(x) get_mode(x))
#pct_oob_estimates = apply(pct_oob,1,function(x)median(x,na.rm=TRUE))
names(pct_oob_estimates) = rownames(mat_bs)
concordance(y ~ pct_oob_estimates)$concordance # 0.59, similar to what is reported.

# Compute survival curve:
fit_km <- survfit(y ~ pct_oob_estimates )
surv_curve = ggsurvplot(
  fit_km,
  data = as.data.frame(list_x[[1]]),
  risk.table = TRUE,
  pval = TRUE,
  palette = c("darkred", "#2E9FDF"),
  xlim = c(0, 5),
  conf.int = TRUE,
  break.time.by = 1,
  legend.labs = c("High-risk", "Low-risk"), 
  xlab="Years since CSCC diagnosis"
)
pdf(file.path(dir_results_script_cutoff,"Final_model","Internally_validated_Survival_curve.pdf"), width = 4.5, height = 5.5)
print(surv_curve, newpage = FALSE)
dev.off()


# Corresponding source data:
fit_summary <- summary(fit_km)

source_data <- data.frame(
  time = fit_summary$time,
  survival = fit_summary$surv,
  lower = fit_summary$lower,
  upper = fit_summary$upper,
  n_risk = fit_summary$n.risk,
  n_event = fit_summary$n.event,
  group = fit_summary$strata
)%>%mutate(group=gsub("pct_oob_estimates=","",group))%>%
  mutate(group=ifelse(group=="high_risk","High-risk","Low-risk"))

write.xlsx(source_data, file.path(dir_results_script_cutoff,"Final_model","Internally_validated_Survival_curve_sourcedata.xlsx"), row.names = FALSE)

# Check that it is independent of callable coverage:
df_check = merge(data.frame(Sample_id = names(pct_oob_estimates),lp=pct_oob_estimates),sample_summary[,c("Sample.ID","Percent_Bait_with_callable_coverage","Metastasis","FU.metastasis.years")],by.x="Sample_id",by.y="Sample.ID")
df_check$Metastasis_num = ifelse(df_check$Metastasis=="Case",1,0)
coxph(Surv(FU.metastasis.years,Metastasis_num)~lp+Percent_Bait_with_callable_coverage,data=df_check)
df_check=merge(df_check,cp_data[,c("Sample.ID","CP.score","BWH_number","AJCC_number")],by.x="Sample_id",by.y="Sample.ID")
coxph(Surv(FU.metastasis.years,Metastasis_num)~lp+CP.score,data=df_check)
coxph(Surv(FU.metastasis.years,Metastasis_num)~lp+BWH_number,data=df_check)
coxph(Surv(FU.metastasis.years,Metastasis_num)~lp+AJCC_number,data=df_check)
