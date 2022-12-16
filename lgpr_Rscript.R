rm(list = ls())
library(tidyverse)
library(peakRAM)
library(lgpr)

args<-commandArgs(TRUE)

## This is using conda env microbiomeDASim on cluster
# lgpr tutorial: https://jtimonen.github.io/lgpr-usage/tutorials/basic.html

num_start = args[1]
num_end = args[2]
num_of_covariates = 16

for (k in num_start:num_end) {
  sim <- read.table(paste0("/fast/AG_Forslund/ChiaYu/LongDat_package_simulation_test/Many_random_covariates_cont/16covariate/Many_random_16covariate_sim20_cont/simulated_many_random_covariate_data/Many_random_16covariates_sim20_set", k, ".txt"),
                    sep = "\t", header = T)
   #sim <- read.table("Many_random_covariates/Many_random_16covariates_sim20_set1.txt", sep = "\t", header = T)
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
  
  ############## lgpr 
  result_df_95 <- as.data.frame(matrix(nrow = ncol(features), ncol = 2))
  colnames(result_df_95) <- c("Case_signal", "All_true_ones")
  rownames(result_df_95) <- colnames(features)
  result_df_90 <- as.data.frame(matrix(nrow = ncol(features), ncol = 2))
  colnames(result_df_90) <- c("Case_signal", "All_true_ones")
  rownames(result_df_90) <- colnames(features)
  result_df_85 <- as.data.frame(matrix(nrow = ncol(features), ncol = 2))
  colnames(result_df_85) <- c("Case_signal", "All_true_ones")
  rownames(result_df_85) <- colnames(features)
  result_df_80 <- as.data.frame(matrix(nrow = ncol(features), ncol = 2))
  colnames(result_df_80) <- c("Case_signal", "All_true_ones")
  rownames(result_df_80) <- colnames(features)
  result_df_75 <- as.data.frame(matrix(nrow = ncol(features), ncol = 2))
  colnames(result_df_75) <- c("Case_signal", "All_true_ones")
  rownames(result_df_75) <- colnames(features)
  result_df_70 <- as.data.frame(matrix(nrow = ncol(features), ncol = 2))
  colnames(result_df_70) <- c("Case_signal", "All_true_ones")
  rownames(result_df_70) <- colnames(features)
  result_df_65 <- as.data.frame(matrix(nrow = ncol(features), ncol = 2))
  colnames(result_df_65) <- c("Case_signal", "All_true_ones")
  rownames(result_df_65) <- colnames(features)
  result_df_60 <- as.data.frame(matrix(nrow = ncol(features), ncol = 2))
  colnames(result_df_60) <- c("Case_signal", "All_true_ones")
  rownames(result_df_60) <- colnames(features)
  result_df_55 <- as.data.frame(matrix(nrow = ncol(features), ncol = 2))
  colnames(result_df_55) <- c("Case_signal", "All_true_ones")
  rownames(result_df_55) <- colnames(features)
  result_df_50 <- as.data.frame(matrix(nrow = ncol(features), ncol = 2))
  colnames(result_df_50) <- c("Case_signal", "All_true_ones")
  rownames(result_df_50) <- colnames(features)
  result_df_45 <- as.data.frame(matrix(nrow = ncol(features), ncol = 2))
  colnames(result_df_45) <- c("Case_signal", "All_true_ones")
  rownames(result_df_45) <- colnames(features)
  result_df_40 <- as.data.frame(matrix(nrow = ncol(features), ncol = 2))
  colnames(result_df_40) <- c("Case_signal", "All_true_ones")
  rownames(result_df_40) <- colnames(features)
  result_df_35 <- as.data.frame(matrix(nrow = ncol(features), ncol = 2))
  colnames(result_df_35) <- c("Case_signal", "All_true_ones")
  rownames(result_df_35) <- colnames(features)
  result_df_30 <- as.data.frame(matrix(nrow = ncol(features), ncol = 2))
  colnames(result_df_30) <- c("Case_signal", "All_true_ones")
  rownames(result_df_30) <- colnames(features)
  result_df_25 <- as.data.frame(matrix(nrow = ncol(features), ncol = 2))
  colnames(result_df_25) <- c("Case_signal", "All_true_ones")
  rownames(result_df_25) <- colnames(features)
  result_df_20 <- as.data.frame(matrix(nrow = ncol(features), ncol = 2))
  colnames(result_df_20) <- c("Case_signal", "All_true_ones")
  rownames(result_df_20) <- colnames(features)
  result_df_15 <- as.data.frame(matrix(nrow = ncol(features), ncol = 2))
  colnames(result_df_15) <- c("Case_signal", "All_true_ones")
  rownames(result_df_15) <- colnames(features)
  result_df_10 <- as.data.frame(matrix(nrow = ncol(features), ncol = 2))
  colnames(result_df_10) <- c("Case_signal", "All_true_ones")
  rownames(result_df_10) <- colnames(features)
  result_df_5 <- as.data.frame(matrix(nrow = ncol(features), ncol = 2))
  colnames(result_df_5) <- c("Case_signal", "All_true_ones")
  rownames(result_df_5) <- colnames(features)
  
  
  all_rel <- list()
  all_sel_95 <- list()
  all_sel_90 <- list()
  all_sel_85 <- list()
  all_sel_80 <- list()
  all_sel_75 <- list()
  all_sel_70 <- list()
  all_sel_65 <- list()
  all_sel_60 <- list()
  all_sel_55 <- list()
  all_sel_50 <- list()
  all_sel_45 <- list()
  all_sel_40 <- list()
  all_sel_35 <- list()
  all_sel_30 <- list()
  all_sel_25 <- list()
  all_sel_20 <- list()
  all_sel_15 <- list()
  all_sel_10 <- list()
  all_sel_5 <- list()
  
  mem <- peakRAM(
    for (i in 1:ncol(features)) {
      print(i)
      data_sub <- cbind(sim[ , c(1:(4+num_of_covariates))], features[ , i])
      colnames(data_sub)[ncol(data_sub)] <- "y"
      data_sub$Sample_id <- as.factor(data_sub$Sample_id)
      data_sub$ID <- as.factor(data_sub$ID)
      data_sub$time <- as.numeric(data_sub$time)
      data_sub$group <- as.factor(data_sub$group)
      
      # All continuous variables have Sample_id as random effect
      if (num_of_covariates > 0) {
        my_formula <- as.formula(paste0("y ~ time + time|Sample_id + ", 
                                        paste(colnames(data_sub)[5:(ncol(data_sub) - 1)], collapse = " + ")))
      } else {
        my_formula <- as.formula(paste0("y ~ time + time|Sample_id"))
      }
      model <- lgpr::create_model(formula = my_formula, 
                                  data = data_sub, 
                                  verbose = FALSE)
      #print(model)
      
      fit <- lgp(formula = as.formula(as.character(model@model_formula)),
                 data     = data_sub,
                 iter     = 100,
                 chains   = 4,
                 refresh  = 0, quiet = T,
                 cores = 1)
      #print(fit)
      Relevance <- relevances(fit, reduce = mean)
      rel <- data.frame(Relevance)
      all_rel[[i]] <- rel
      
      result_sub_95 <- lgpr::select(fit, threshold = 0.95)
      result_sub_90 <- lgpr::select(fit, threshold = 0.90)
      result_sub_85 <- lgpr::select(fit, threshold = 0.85)
      result_sub_80 <- lgpr::select(fit, threshold = 0.80)
      result_sub_75 <- lgpr::select(fit, threshold = 0.75) 
      result_sub_70 <- lgpr::select(fit, threshold = 0.70)
      result_sub_65 <- lgpr::select(fit, threshold = 0.65)
      result_sub_60 <- lgpr::select(fit, threshold = 0.60)
      result_sub_55 <- lgpr::select(fit, threshold = 0.55)
      result_sub_50 <- lgpr::select(fit, threshold = 0.50)
      result_sub_45 <- lgpr::select(fit, threshold = 0.45)
      result_sub_40 <- lgpr::select(fit, threshold = 0.40)
      result_sub_35 <- lgpr::select(fit, threshold = 0.35)
      result_sub_30 <- lgpr::select(fit, threshold = 0.30)
      result_sub_25 <- lgpr::select(fit, threshold = 0.25)
      result_sub_20 <- lgpr::select(fit, threshold = 0.20)
      result_sub_15 <- lgpr::select(fit, threshold = 0.15)
      result_sub_10 <- lgpr::select(fit, threshold = 0.10)
      result_sub_5 <- lgpr::select(fit, threshold = 0.05)
      
      all_sel_95[[i]] <- result_sub_95
      all_sel_90[[i]] <- result_sub_90
      all_sel_85[[i]] <- result_sub_85
      all_sel_80[[i]] <- result_sub_80
      all_sel_75[[i]] <- result_sub_75
      all_sel_70[[i]] <- result_sub_70
      all_sel_65[[i]] <- result_sub_65
      all_sel_60[[i]] <- result_sub_60
      all_sel_55[[i]] <- result_sub_55
      all_sel_50[[i]] <- result_sub_50
      all_sel_45[[i]] <- result_sub_45
      all_sel_40[[i]] <- result_sub_40
      all_sel_35[[i]] <- result_sub_35
      all_sel_30[[i]] <- result_sub_30
      all_sel_25[[i]] <- result_sub_25
      all_sel_20[[i]] <- result_sub_20
      all_sel_15[[i]] <- result_sub_15
      all_sel_10[[i]] <- result_sub_10
      all_sel_5[[i]] <- result_sub_5
      
      result_df_95$Case_signal[i] <- result_sub_95$Component[1]
      result_df_95$All_true_ones[i] <- paste(collapse = ",", rownames(result_sub_95)[result_sub_95$Component == TRUE])
      result_df_90$Case_signal[i] <- result_sub_90$Component[1]
      result_df_90$All_true_ones[i] <- paste(collapse = ",", rownames(result_sub_90)[result_sub_90$Component == TRUE])
      result_df_85$Case_signal[i] <- result_sub_85$Component[1]
      result_df_85$All_true_ones[i] <- paste(collapse = ",", rownames(result_sub_85)[result_sub_85$Component == TRUE])
      result_df_80$Case_signal[i] <- result_sub_80$Component[1]
      result_df_80$All_true_ones[i] <- paste(collapse = ",", rownames(result_sub_80)[result_sub_80$Component == TRUE])
      result_df_75$Case_signal[i] <- result_sub_75$Component[1]
      result_df_75$All_true_ones[i] <- paste(collapse = ",", rownames(result_sub_75)[result_sub_75$Component == TRUE])
      result_df_70$Case_signal[i] <- result_sub_70$Component[1]
      result_df_70$All_true_ones[i] <- paste(collapse = ",", rownames(result_sub_70)[result_sub_70$Component == TRUE])
      result_df_65$Case_signal[i] <- result_sub_65$Component[1]
      result_df_65$All_true_ones[i] <- paste(collapse = ",", rownames(result_sub_65)[result_sub_65$Component == TRUE])
      result_df_60$Case_signal[i] <- result_sub_60$Component[1]
      result_df_60$All_true_ones[i] <- paste(collapse = ",", rownames(result_sub_60)[result_sub_60$Component == TRUE])
      result_df_55$Case_signal[i] <- result_sub_55$Component[1]
      result_df_55$All_true_ones[i] <- paste(collapse = ",", rownames(result_sub_55)[result_sub_55$Component == TRUE])
      result_df_50$Case_signal[i] <- result_sub_50$Component[1]
      result_df_50$All_true_ones[i] <- paste(collapse = ",", rownames(result_sub_50)[result_sub_50$Component == TRUE])
      result_df_45$Case_signal[i] <- result_sub_45$Component[1]
      result_df_45$All_true_ones[i] <- paste(collapse = ",", rownames(result_sub_45)[result_sub_45$Component == TRUE])
      result_df_40$Case_signal[i] <- result_sub_40$Component[1]
      result_df_40$All_true_ones[i] <- paste(collapse = ",", rownames(result_sub_40)[result_sub_40$Component == TRUE])
      result_df_35$Case_signal[i] <- result_sub_35$Component[1]
      result_df_35$All_true_ones[i] <- paste(collapse = ",", rownames(result_sub_35)[result_sub_35$Component == TRUE])
      result_df_30$Case_signal[i] <- result_sub_30$Component[1]
      result_df_30$All_true_ones[i] <- paste(collapse = ",", rownames(result_sub_30)[result_sub_30$Component == TRUE])
      result_df_25$Case_signal[i] <- result_sub_25$Component[1]
      result_df_25$All_true_ones[i] <- paste(collapse = ",", rownames(result_sub_25)[result_sub_25$Component == TRUE])
      result_df_20$Case_signal[i] <- result_sub_20$Component[1]
      result_df_20$All_true_ones[i] <- paste(collapse = ",", rownames(result_sub_20)[result_sub_20$Component == TRUE])
      result_df_15$Case_signal[i] <- result_sub_15$Component[1]
      result_df_15$All_true_ones[i] <- paste(collapse = ",", rownames(result_sub_15)[result_sub_15$Component == TRUE])
      result_df_10$Case_signal[i] <- result_sub_10$Component[1]
      result_df_10$All_true_ones[i] <- paste(collapse = ",", rownames(result_sub_10)[result_sub_10$Component == TRUE])
      result_df_5$Case_signal[i] <- result_sub_5$Component[1]
      result_df_5$All_true_ones[i] <- paste(collapse = ",", rownames(result_sub_5)[result_sub_5$Component == TRUE])
    }
  )
  write.table(result_df_95, paste0("/fast/AG_Forslund/ChiaYu/Other_tools_for_comparison/lgpr/Many_random_covariates_100feature/16covariate_100feature/Many_random_100feature_16covariate_sim20/lgpr_result/Many_random_16covariates_sim20_100features_set", k, "_lgpr_result_df_95.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)
  write.table(result_df_90, paste0("/fast/AG_Forslund/ChiaYu/Other_tools_for_comparison/lgpr/Many_random_covariates_100feature/16covariate_100feature/Many_random_100feature_16covariate_sim20/lgpr_result/Many_random_16covariates_sim20_100features_set", k, "_lgpr_result_df_90.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)
  write.table(result_df_85, paste0("/fast/AG_Forslund/ChiaYu/Other_tools_for_comparison/lgpr/Many_random_covariates_100feature/16covariate_100feature/Many_random_100feature_16covariate_sim20/lgpr_result/Many_random_16covariates_sim20_100features_set", k, "_lgpr_result_df_85.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)
  write.table(result_df_80, paste0("/fast/AG_Forslund/ChiaYu/Other_tools_for_comparison/lgpr/Many_random_covariates_100feature/16covariate_100feature/Many_random_100feature_16covariate_sim20/lgpr_result/Many_random_16covariates_sim20_100features_set", k, "_lgpr_result_df_80.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)
  write.table(result_df_75, paste0("/fast/AG_Forslund/ChiaYu/Other_tools_for_comparison/lgpr/Many_random_covariates_100feature/16covariate_100feature/Many_random_100feature_16covariate_sim20/lgpr_result/Many_random_16covariates_sim20_100features_set", k, "_lgpr_result_df_75.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)
  write.table(result_df_70, paste0("/fast/AG_Forslund/ChiaYu/Other_tools_for_comparison/lgpr/Many_random_covariates_100feature/16covariate_100feature/Many_random_100feature_16covariate_sim20/lgpr_result/Many_random_16covariates_sim20_100features_set", k, "_lgpr_result_df_70.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)
  write.table(result_df_65, paste0("/fast/AG_Forslund/ChiaYu/Other_tools_for_comparison/lgpr/Many_random_covariates_100feature/16covariate_100feature/Many_random_100feature_16covariate_sim20/lgpr_result/Many_random_16covariates_sim20_100features_set", k, "_lgpr_result_df_65.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)
  write.table(result_df_60, paste0("/fast/AG_Forslund/ChiaYu/Other_tools_for_comparison/lgpr/Many_random_covariates_100feature/16covariate_100feature/Many_random_100feature_16covariate_sim20/lgpr_result/Many_random_16covariates_sim20_100features_set", k, "_lgpr_result_df_60.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)
  write.table(result_df_55, paste0("/fast/AG_Forslund/ChiaYu/Other_tools_for_comparison/lgpr/Many_random_covariates_100feature/16covariate_100feature/Many_random_100feature_16covariate_sim20/lgpr_result/Many_random_16covariates_sim20_100features_set", k, "_lgpr_result_df_55.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)
  write.table(result_df_50, paste0("/fast/AG_Forslund/ChiaYu/Other_tools_for_comparison/lgpr/Many_random_covariates_100feature/16covariate_100feature/Many_random_100feature_16covariate_sim20/lgpr_result/Many_random_16covariates_sim20_100features_set", k, "_lgpr_result_df_50.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)
  write.table(result_df_45, paste0("/fast/AG_Forslund/ChiaYu/Other_tools_for_comparison/lgpr/Many_random_covariates_100feature/16covariate_100feature/Many_random_100feature_16covariate_sim20/lgpr_result/Many_random_16covariates_sim20_100features_set", k, "_lgpr_result_df_45.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)
  write.table(result_df_40, paste0("/fast/AG_Forslund/ChiaYu/Other_tools_for_comparison/lgpr/Many_random_covariates_100feature/16covariate_100feature/Many_random_100feature_16covariate_sim20/lgpr_result/Many_random_16covariates_sim20_100features_set", k, "_lgpr_result_df_40.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)
  write.table(result_df_35, paste0("/fast/AG_Forslund/ChiaYu/Other_tools_for_comparison/lgpr/Many_random_covariates_100feature/16covariate_100feature/Many_random_100feature_16covariate_sim20/lgpr_result/Many_random_16covariates_sim20_100features_set", k, "_lgpr_result_df_35.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)
  write.table(result_df_30, paste0("/fast/AG_Forslund/ChiaYu/Other_tools_for_comparison/lgpr/Many_random_covariates_100feature/16covariate_100feature/Many_random_100feature_16covariate_sim20/lgpr_result/Many_random_16covariates_sim20_100features_set", k, "_lgpr_result_df_30.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)
  write.table(result_df_25, paste0("/fast/AG_Forslund/ChiaYu/Other_tools_for_comparison/lgpr/Many_random_covariates_100feature/16covariate_100feature/Many_random_100feature_16covariate_sim20/lgpr_result/Many_random_16covariates_sim20_100features_set", k, "_lgpr_result_df_25.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)
  write.table(result_df_20, paste0("/fast/AG_Forslund/ChiaYu/Other_tools_for_comparison/lgpr/Many_random_covariates_100feature/16covariate_100feature/Many_random_100feature_16covariate_sim20/lgpr_result/Many_random_16covariates_sim20_100features_set", k, "_lgpr_result_df_20.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)
  write.table(result_df_15, paste0("/fast/AG_Forslund/ChiaYu/Other_tools_for_comparison/lgpr/Many_random_covariates_100feature/16covariate_100feature/Many_random_100feature_16covariate_sim20/lgpr_result/Many_random_16covariates_sim20_100features_set", k, "_lgpr_result_df_15.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)
  write.table(result_df_10, paste0("/fast/AG_Forslund/ChiaYu/Other_tools_for_comparison/lgpr/Many_random_covariates_100feature/16covariate_100feature/Many_random_100feature_16covariate_sim20/lgpr_result/Many_random_16covariates_sim20_100features_set", k, "_lgpr_result_df_10.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)
  write.table(result_df_5, paste0("/fast/AG_Forslund/ChiaYu/Other_tools_for_comparison/lgpr/Many_random_covariates_100feature/16covariate_100feature/Many_random_100feature_16covariate_sim20/lgpr_result/Many_random_16covariates_sim20_100features_set", k, "_lgpr_result_df_5.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)
  
  write.table(x = mem, file = paste0("/fast/AG_Forslund/ChiaYu/Other_tools_for_comparison/lgpr/Many_random_covariates_100feature/16covariate_100feature/Many_random_100feature_16covariate_sim20/lgpr_result/Many_random_16covariates_sim20_100features_set", k, "_lgpr_MemTime.txt"), sep = "\t",
              col.names = T, quote = F, row.names = F)
}
