library(microbiomeDASim)
library(tidyverse)
library(LongDat)
library(peakRAM)
library(ggpubr)

args<-commandArgs(TRUE)

num_start = args[1]
num_end = args[2]

### First generate 16 covariates, reusing the simulated data
for (i in num_start:num_end) {
  set.seed(i)
  #Load simulated data
  raw <- read.table(paste0("/fast/AG_Forslund/ChiaYu/LongDat_package_simulation_test/microbiomeDASim/simulation_microbiomeDASim_22/simulated_data/Sim_microbiomeDASim_data_22_set", i, ".txt"),
                    sep = "\t", header = T, stringsAsFactors = F)
  
  #Simulate covariates
  #Ref: https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables
  complement <- function(y, rho, x) {
    if (missing(x)) x <- rnorm(length(y)) # Optional: supply a default if `x` is not given
    y.perp <- residuals(lm(x ~ y))
    rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
  }
  
  # Totally 16 covariates
  # The rho are 0.25, 0.5, 0.75, 0.99, randomly pick the 16 rhos from these 4 levels
  rho_levels <- c(0.25, 0.5, 0.75, 0.99)
  
  set.seed(i)
  my_levels <- sample(x = rho_levels, size = 16, replace = T)
 
  a1 = complement(y = raw$time, rho = my_levels[1])
  a2 = complement(y = raw$time, rho = my_levels[2])
  a3 = complement(y = raw$time, rho = my_levels[3])
  a4 = complement(y = raw$time, rho = my_levels[4])
  a5 = complement(y = raw$time, rho = my_levels[5])
  a6 = complement(y = raw$time, rho = my_levels[6])
  a7 = complement(y = raw$time, rho = my_levels[7])
  a8 = complement(y = raw$time, rho = my_levels[8])
  a9 = complement(y = raw$time, rho = my_levels[9])
  a10 = complement(y = raw$time, rho = my_levels[10])
  a11 = complement(y = raw$time, rho = my_levels[11])
  a12 = complement(y = raw$time, rho = my_levels[12])
  a13 = complement(y = raw$time, rho = my_levels[13])
  a14 = complement(y = raw$time, rho = my_levels[14])
  a15 = complement(y = raw$time, rho = my_levels[15])
  a16 = complement(y = raw$time, rho = my_levels[16])
  
  raw <- raw %>%
    mutate(.after = 3,
           Simulated_covariate_a1 = a1,
           Simulated_covariate_a2 = a2,
           Simulated_covariate_a3 = a3,
           Simulated_covariate_a4 = a4,
           Simulated_covariate_a5 = a5,
           Simulated_covariate_a6 = a6,
           Simulated_covariate_a7 = a7,
           Simulated_covariate_a8 = a8,
           Simulated_covariate_a9 = a9,
           Simulated_covariate_a10 = a10,
           Simulated_covariate_a11 = a11,
           Simulated_covariate_a12 = a12,
           Simulated_covariate_a13 = a13,
           Simulated_covariate_a14 = a14,
           Simulated_covariate_a15 = a15,
           Simulated_covariate_a16 = a16,
           
           )
  
write.table(raw, paste0("/fast/AG_Forslund/ChiaYu/LongDat_package_simulation_test/Many_random_covariates_cont/16covariate/Many_random_16covariate_sim22_cont/simulated_many_random_covariate_data/Many_random_16covariates_sim22_set", i, ".txt"), sep = "\t",
              row.names = F, col.names = T, quote = F)
}

### Run longdat_cont
for (i in num_start:num_end) {
  print(i)
  input_table <- read.table(paste0("/fast/AG_Forslund/ChiaYu/LongDat_package_simulation_test/Many_random_covariates_cont/16covariate/Many_random_16covariate_sim22_cont/simulated_many_random_covariate_data/Many_random_16covariates_sim22_set", i, ".txt"),
                            sep = "\t", header = T)
  mem <- peakRAM(cont_result <- longdat_cont(input = input_table,
                                             data_type = "count",
                                             test_var = "time",
                                             variable_col = 20,
                                             fac_var = c(1, 3),
                                             not_used =  c(3),
                                             adjustMethod = "fdr",
                                             model_q = 0.1,
                                             posthoc_q = 0.05,
                                             theta_cutoff = 2^20,
                                             nonzero_count_cutoff1 = 9,
                                             nonzero_count_cutoff2 = 5,
                                             verbose = T))

  write.table(x = cont_result[[1]], file = paste0("/fast/AG_Forslund/ChiaYu/LongDat_package_simulation_test/Many_random_covariates_cont/16covariate/Many_random_16covariate_sim22_cont/longdat_result/Many_random_16covariates_sim22_set", i, "_longdatcont_result_table.txt"), sep = "\t",
             col.names = T, quote = F, row.names = F)
  write.table(x = cont_result[[2]], file = paste0("/fast/AG_Forslund/ChiaYu/LongDat_package_simulation_test/Many_random_covariates_cont/16covariate/Many_random_16covariate_sim22_cont/longdat_result/Many_random_16covariates_sim22_set", i, "_longdatcont_confounder_table.txt"), sep = "\t",
            col.names = T, quote = F, row.names = F)
  write.table(x = mem, file = paste0("/fast/AG_Forslund/ChiaYu/LongDat_package_simulation_test/Many_random_covariates_cont/16covariate/Many_random_16covariate_sim22_cont/longdat_result/Many_random_16covariates_sim22_set", i, "_longdatcont_memory_time.txt"), sep = "\t",
            col.names = T, quote = F, row.names = F)
}
