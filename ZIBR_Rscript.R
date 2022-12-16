rm(list = ls())
library(tidyverse)
library(peakRAM)
library(ZIBR)

args<-commandArgs(TRUE)

## This is using conda env longdat on cluster
# ZIBR tutorial: https://github.com/chvlyl/ZIBR

num_start = args[1]
num_end = args[2]
num_of_covariates = 16

for (i in num_start:num_end) {
  sim <- read.table(paste0("/fast/AG_Forslund/ChiaYu/LongDat_package_simulation_test/Many_random_covariates_cont/16covariate/Many_random_16covariate_sim19_cont/simulated_many_random_covariate_data/Many_random_16covariates_sim19_set", i, ".txt"),
                    sep = "\t", header = T)
   #sim <- read.table("Many_random_covariates/Many_random_16covariates_sim19_set1.txt", sep = "\t", header = T)
  # Only the treatment group
  sim <- sim %>% 
    dplyr::filter(group == "Treatment")
  sim$ID <- paste0("p_", sim$ID, "_", sim$time)
  
  sim <- sim %>% 
    mutate(.before = 1,
           Sample_id = paste0(str_split_fixed(ID, pattern = "_", n = 3)[ ,1], "_", str_split_fixed(ID, pattern = "_", n = 3)[ , 2]))
  
  #Extract only 100 features (10 Diff_bug, 90 Non_diff_bug)
  features = sim[ , c((5 + num_of_covariates):(5 + num_of_covariates + 9), (5 + num_of_covariates + 20): (5 + num_of_covariates + 20 + 89))] # Feature
  features_rel <- features/rowSums(features) #Get relative abundance cuz ZIBR needs it

  N = ncol(features_rel)
  
  # Metadata
  if (num_of_covariates != 0) {
    my_meta_data <- sim %>% 
      dplyr::select(c(3, 5: (4 + num_of_covariates)))
  } else {
    my_meta_data <- sim %>% 
      dplyr::select(3)
  }
  
  ############## ZIBR 
  # https://github.com/chvlyl/ZIBR
  results <- as.data.frame(matrix(nrow = ncol(features_rel), ncol = 3))
  colnames(results) <- c("log_time_p", "beta_time_p", "joint_time_p")
  rownames(results) <- colnames(features_rel)
  rownames(results) <- colnames(features_rel)
  
  mem <- peakRAM(
    for (j in 1:ncol(features_rel)) {
      print(j)
      zibr.fit <- zibr(logistic.cov = my_meta_data, 
                       beta.cov = my_meta_data, 
                       Y = features_rel[ , j], 
                       subject.ind = sim$Sample_id,
                       time.ind = sim$time)
      
      results$log_time_p[j] <- zibr.fit$logistic.est.table[2, 2]
      results$beta_time_p[j] <- zibr.fit$beta.est.table[2, 2]
      results$joint_time_p[j] <- zibr.fit$joint.p[1]
    }
  )
  
  results$joint_time_q <- p.adjust(results$joint_time_p, method = "fdr")
  
  write.table(x = results, file =paste0("/fast/AG_Forslund/ChiaYu/Other_tools_for_comparison/ZIBR/Many_random_covariates_100feature/16covariate_100feature/Many_random_100feature_16covariate_sim19/Zibr_result/Many_random_100feature_16covariates_sim19_set", i, "_ZIBR_result.txt"),
              sep = "\t", row.names = T, col.names = NA, quote = F)
  write.table(x = mem, file =paste0("/fast/AG_Forslund/ChiaYu/Other_tools_for_comparison/ZIBR/Many_random_covariates_100feature/16covariate_100feature/Many_random_100feature_16covariate_sim19/Zibr_result/Many_random_100feature_16covariates_sim19_set", i, "_ZIBR_MemTime.txt"), 
              sep = "\t", row.names = F, col.names = T, quote = F)
}
