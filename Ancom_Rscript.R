rm(list = ls())
library(tidyverse)
library(peakRAM)
library(compositions)
source("/fast/AG_Forslund/ChiaYu/Other_tools_for_comparison/ANCOM/ANCOM-Code-Archive-master/programs/ancom.R")
#source("/Users/Jessica/Documents/Lab/Other_tools_for_comparison/ANCOM-Code-Archive-master/programs/ancom.R")
## This is using conda env longdat on cluster
# ANCOM tutorial: https://github.com/FrederickHuangLin/ANCOM-Code-Archive

args<-commandArgs(TRUE)
num_start = args[1]
num_end = args[2]
num_of_covariates = 16

for (i in num_start:num_end) {
  sim <- read.table(paste0("/fast/AG_Forslund/ChiaYu/LongDat_package_simulation_test/Many_random_covariates_cont/16covariate/Many_random_16covariate_sim20_cont/simulated_many_random_covariate_data/Many_random_16covariates_sim20_set", i, ".txt"),
                    sep = "\t", header = T)
   #sim <- read.table("/Users/Jessica/Documents/Lab/Other_tools_for_comparison/ZIBR/Many_random_covariates_100feature/Many_random_16covariates_sim20_set1.txt",
   #                  sep = "\t", header = T)
  # Only the treatment group
  sim <- sim %>% 
    dplyr::filter(group == "Treatment")
  sim$ID <- paste0("p_", sim$ID, "_", sim$time)
  
  sim <- sim %>% 
    mutate(.before = 1,
           Sample_id = paste0(str_split_fixed(ID, pattern = "_", n = 3)[ ,1], "_", str_split_fixed(ID, pattern = "_", n = 3)[ , 2]))
  
  #Extract only 100 features (10 Diff_bug, 90 Non_diff_bug)
  features = sim[ , c((5 + num_of_covariates):(5 + num_of_covariates + 9), (5 + num_of_covariates + 20): (5 + num_of_covariates + 20 + 89))] # Feature

  N = ncol(features)
  
  # Metadata
  if (num_of_covariates != 0) {
    my_meta_data <- sim %>% 
      dplyr::select(c(1:4, 5: (4 + num_of_covariates)))
  } else {
    my_meta_data <- sim %>% 
      dplyr::select(1, 2, 3, 4)
  }
  
  my_meta_data$time <- as.factor(my_meta_data$time)
  
  feature_table = features
  sample_var = "ID" # Unique sample name for each sample
  group_var = NULL # Character. The name of the group indicator.
  out_cut = 0.05 # This is the default value
  zero_cut = (N-6)/N #Numerical fraction between 0 and 1. Taxa with proportion of zeroes greater than zero_cut are not included in the analysis.(N-6/N) means there are at least 6 nonzero values, this is align with longdat nonzero_count_cutoff2
  lib_cut = 0 # Numeric. Samples with library size less than lib_cut are not included in the analysis.
  neg_lb = T
  meta_data = my_meta_data
  
  #Since there's no >= 2 groups in my case (all are Treatment group), no need to run feature_table_pre_process to detect structural zero
  #prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
  #                                   out_cut, zero_cut, lib_cut, neg_lb)
  
  ## Edit the objects to fit ANCOM's requirement
  features_t <- as.data.frame(t(features)) #Data frame representing OTU table with taxa in rows (rownames) and samples in columns (colnames)
  colnames(features_t) <- my_meta_data$ID
  
  
  feature_table = features_t 
  meta_data = my_meta_data 
  struc_zero = NULL # Structural zero info
  
  main_var = "time"
  p_adj_method = "fdr"
  alpha = 0.1
  adj_formula = paste(colnames(my_meta_data)[5:ncol(my_meta_data)], collapse = " + ")
  rand_formula = "~ 1 | Sample_id" #Individual random effect
  
  mem <- peakRAM(res <- ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
                              alpha, adj_formula, rand_formula))
  
  write.table(x = res$out, file =paste0("/fast/AG_Forslund/ChiaYu/Other_tools_for_comparison/ANCOM/Many_random_covariates_100features/16covariate_100feature/Many_random_100feature_16covariate_sim20/ANCOM_result/Many_random_16covariates_sim20_100features_set", i, "_Ancom_out.csv"),
              sep = "\t", row.names = T, col.names = NA, quote = F)
  write.table(x = res$p_data, file =paste0("/fast/AG_Forslund/ChiaYu/Other_tools_for_comparison/ANCOM/Many_random_covariates_100features/16covariate_100feature/Many_random_100feature_16covariate_sim20/ANCOM_result/Many_random_16covariates_sim20_100features_set", i, "_Ancom_p_data.txt"),
              sep = "\t", row.names = T, col.names = NA, quote = F)
  write.table(x = res$q_data, file =paste0("/fast/AG_Forslund/ChiaYu/Other_tools_for_comparison/ANCOM/Many_random_covariates_100features/16covariate_100feature/Many_random_100feature_16covariate_sim20/ANCOM_result/Many_random_16covariates_sim20_100features_set", i, "_Ancom_q_data.txt"),
              sep = "\t", row.names = T, col.names = NA, quote = F)
  write.table(x = mem, file =paste0("/fast/AG_Forslund/ChiaYu/Other_tools_for_comparison/ANCOM/Many_random_covariates_100features/16covariate_100feature/Many_random_100feature_16covariate_sim20/ANCOM_result/Many_random_16covariates_sim20_100features_set", i, "_Ancom_MemTime.txt"), 
              sep = "\t", row.names = F, col.names = T, quote = F)
}
