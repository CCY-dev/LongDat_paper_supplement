library(tidyverse)
library(LongDat)
library(peakRAM)

for (i in 1:100) {
  print(i)
  input = read.table(paste0("/fast/AG_Forslund/ChiaYu/LongDat_package_simulation_test/microbiomeDASim/simulation_microbiomeDASim_2_conf1/simulated_conf_data/Sim2_conf_data_1_set", i, ".txt"),
                     sep = "\t", header = T)
  mem <- peakRAM(cont_result <- longdat_cont(input = input,
                                             data_type = "count",
                                             test_var = "Simulated_confounder",
                                             variable_col = 5,
                                             fac_var = c(1, 4),
                                             not_used =  c(4),
                                             adjustMethod = "fdr",
                                             model_q = 0.1,
                                             posthoc_q = 0.05,
                                             theta_cutoff = 2^20,
                                             nonzero_count_cutoff1 = 9,
                                             nonzero_count_cutoff2 = 5,
                                             verbose = F))
  #Check memory and time
  mem

  #Check result
  result_table <- cont_result[[1]]
}
