library(microbiomeDASim)
library(tidyverse)
library(LongDat)
library(peakRAM)
library(permute)

args <- commandArgs(TRUE)

num_start = args[1]
num_end = args[2]

### Run longdat_disc
# The simulated data is the same as the one used in sim1, no need to generate again
for (i in num_start:num_end) {
  print(i)
  
  input_table <- read.table(paste0("/fast/AG_Forslund/ChiaYu/LongDat_package_simulation_test/microbiomeDASim/simulation_microbiomeDASim_2_conf1/simulated_conf_data/Sim2_conf_data_1_set", i, ".txt"),
                            sep = "\t", header = T)
  
  # Shuffle the time variable in stratified manner (within each individual)
  set.seed(i)
  random_order = shuffle(input_table$time, control = how(blocks = input_table$ID))
  input_table <- input_table %>% 
    dplyr::mutate(Shuffled_time = input_table$time[random_order], 
                  .after = 3)
  
  write.table(x = input_table, file = paste0("/fast/AG_Forslund/ChiaYu/LongDat_package_simulation_test/Test_time_shuffled/Time_shuffled_sim2_conf1/time_shuffled_data/Time_shuffled_sim2_conf1_set", i, ".txt"), sep = "\t",
              row.names = T, col.names = NA, quote = F)

  ## Run longdat_cont
  mem <- peakRAM(cont_result <- longdat_cont(input = input_table,
                                             data_type = "count",
                                             test_var = "Shuffled_time",
                                             variable_col = 6,
                                             fac_var = c(1, 5),
                                             not_used =  c(3, 5),
                                             adjustMethod = "fdr",
                                             model_q = 0.1,
                                             posthoc_q = 0.05,
                                             theta_cutoff = 2^20,
                                             nonzero_count_cutoff1 = 9,
                                             nonzero_count_cutoff2 = 5,
                                             verbose = F))

  write.table(x = cont_result[[1]], file = paste0("/fast/AG_Forslund/ChiaYu/LongDat_package_simulation_test/Test_time_shuffled/Time_shuffled_sim2_conf1/longdat_result/Time_shuffled_sim2_conf1_set", i, "_longdatcont_result_table.txt"), sep = "\t",
              row.names = T, col.names = NA, quote = F)
  write.table(x = cont_result[[2]], file = paste0("/fast/AG_Forslund/ChiaYu/LongDat_package_simulation_test/Test_time_shuffled/Time_shuffled_sim2_conf1/longdat_result/Time_shuffled_sim2_conf1_set", i, "_longdatcont_confounder_table.txt"), sep = "\t",
              row.names = T, col.names = NA, quote = F)
  write.table(x = mem, file = paste0("/fast/AG_Forslund/ChiaYu/LongDat_package_simulation_test/Test_time_shuffled/Time_shuffled_sim2_conf1/longdat_result/Time_shuffled_sim2_conf1_set", i, "_longdatcont_memory_time.txt"), sep = "\t",
              row.names = F, col.names = T, quote = F)
}
