library(tidyverse)
library(LongDat)
library(peakRAM)

### Run longdat_cont
for (i in 1:100) {
  print(i)
  input = read.table(paste0("/fast/AG_Forslund/ChiaYu/LongDat_package_simulation_test/microbiomeDASim/simulation_microbiomeDASim_2/simulated_data/Sim_microbiomeDASim_data_2_set", i, ".txt"),
                     sep = "\t", header = T)
  mem <- peakRAM(cont_result <- longdat_cont(input = input,
                                             data_type = "count",
                                             test_var = "time",
                                             variable_col = 4,
                                             fac_var = c(1, 3),
                                             not_used =  c(3),
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

