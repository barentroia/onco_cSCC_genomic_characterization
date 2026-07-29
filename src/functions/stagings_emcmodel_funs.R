predict_lp <- function(age,
                       sex,
                       n_cscc,
                       t_location,
                       t_diam,
                       t_involvement,
                       differentiation,
                       pni_or_lvi){
  
  #' Implement linear predictor to compute absolute risk
  #'
  #' @param age vector. Vector with age (in decades)
  #' @param sex vector. Vector with gender
  #' @param n_cscc vector. Vector with the number of prior cSCCs
  #' (excluding current cSCC)
  #' @param t_location vector. Vector with the tumor location
  #' @param t_diam vector. Vector with tumor diameter (cm)
  #' @param t_involvement vector. Vector with tissue involvement
  #' @param differentiation vector. Vector with differentiation grade
  #' @param pni_or_lvi vector. Vector with presence of perineural
  #' or lymphovascular invasion
  #' 
  #' @return result of linear predictor
  
  # Adjust variables
  ## Sex
  sex <- ifelse(sex == "Male", 1, 0)
  ## Tumor diameter
  t_diam <- ifelse(as.numeric(t_diam) > 4, 4, t_diam)
  ## Number of cSCC
  n_cscc <- ifelse(n_cscc > 4, 4 , n_cscc)
  ## Tumor location
  t_location_sneck <- ifelse(t_location == "Scalp/neck", 1, 0)
  t_location_face <- ifelse(t_location == "Face", 1, 0)
  ## Tumor involvement
  tinv_sc <- ifelse(t_involvement == "Subcutaneous fat", 1, 0)
  tinv_beysc <- ifelse(t_involvement == "Beyond subcutaneous fat", 1, 0)
  ## Differentiation grade
  differentiation <- ifelse(differentiation == "Poor/undifferentiated", 1, 0)
  ## Perineural or lymphovascular invasion
  pni_or_lvi <- ifelse(pni_or_lvi == "Yes", 1, 0)
  
  # Remove mean values of continuous variables
  ## Age
  age <- age/10 - 7.47
  ## Number of cSCC
  n_cscc <- n_cscc - 0.29 
  ## Tumor diameter
  t_diam <- t_diam - 1.20 
  
  # Coefficients of the linear predictor
  coeffs_model <- c(0.25, 0.51, 0.57, -0.72, 0.52, 0.54, 0.33, 1.4, 1.3, 0.8)
  df <- data.frame(age = age,
                   sex = sex,
                   n_cscc = n_cscc,
                   t_location_sneck = t_location_sneck,
                   t_location_face = t_location_face,
                   t_diam = t_diam,
                   tinv_sc = tinv_sc,
                   tinv_beysc = tinv_beysc,
                   differentiation = differentiation,
                   pni_or_lvi = pni_or_lvi)
  
  # Compute linear predictor
  lp <- apply(as.matrix(df), 1, function(x) sum(as.numeric(t(coeffs_model)%*%x)))
  
  lp
  
}

ajcc_staging <- function(t_diam,
                         pni,
                         t_involvement,
                         depth_inv,
                         inv_bones){
  #' Compute AJCC staging
  #'
  #' @param t_diam vector. Vector with tumor diameter
  #' @param pni vector. Vector with PNI binary
  #' @param t_involvement vector. Vector with tissue involvement
  #' @param depth_inv vector. Vector with depth of invasion
  #' @param inv_bones vector. Vector with invasion of bones
  #'
  #' @return AJCC_staging
  #'
  
  # Compute T-stages
  T1 <- ifelse(t_diam < 2 & ((pni == "No" |(pni == "Yes" & t_involvement == "Dermis")) & t_involvement != "Beyond subcutaneous fat" & depth_inv <= 6), 1, 0)
  T2 <- ifelse(t_diam >=2 & t_diam <4 & ((pni == "No"| (pni == "Yes" & t_involvement == "Dermis")) & t_involvement != "Beyond subcutaneous fat" & depth_inv <=6), 1, 0)
  T3 <- ifelse((t_diam >= 4 | (pni == "Yes" & t_involvement != "Dermis") | t_involvement == "Beyond subcutaneous fat" | depth_inv > 6) & inv_bones == "No", 1, 0)
  T4 <- ifelse(inv_bones == "Yes", 1, 0)
  
  # Assign correct T-stages
  AJCC_staging <- ifelse(T4 == 1, "T4",
                         ifelse(T3 == 1,"T3",
                                ifelse(T2 == 1,"T2",
                                       ifelse(T1 == 1,"T1", NA))))
  
  AJCC_staging
  
}

bwh_staging <- function(t_diam,
                        differentiation,
                        pni,
                        t_involvement,
                        inv_bones){
  #' Compute BWH staging
  #'
  #' @param t_diam vector. Tumor diameter variable
  #' @param differentiation vector. Vector with differentiation grade
  #' @param pni vector. Vector with PNI binary
  #' @param t_involvement vector. Vector with tissue involvement
  #' @param inv_bones vector. Vector with invasion of bones
  #'
  #' @return BWH staging
  #'
  
  # Compute T-stages
  hr_t_diam <- ifelse(t_diam >= 2, 1, 0)
  hr_differentiation <- ifelse(differentiation == "Poor/undifferentiated", 1, 0)
  hr_pni <- ifelse(pni == "Yes", 1, 0)
  hr_t_involvement = ifelse(t_involvement == "Beyond subcutaneous fat", 1, 0)
  sum_hr <- hr_t_diam + hr_differentiation + hr_pni + hr_t_involvement
  sum_hr <-  ifelse(!(is.na(inv_bones)) & inv_bones == "Yes", 4, sum_hr)
  
  # Assign correct T-stages
  BWH_staging = ifelse(sum_hr == 0, "T1" ,
                       ifelse(sum_hr == 1, "T2a", 
                              ifelse(sum_hr %in% c(2,3), "T2b",
                                     ifelse(sum_hr >= 4, "T3", NA))))
  
  BWH_staging
  
}