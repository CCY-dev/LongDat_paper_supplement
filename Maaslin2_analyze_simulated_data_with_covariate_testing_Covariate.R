library(tidyverse)
library(ggpubr)
library(peakRAM)
library(Maaslin2)

### First generate simulated data using MicrobiomeDASim
for (i in 1:100) {
  metatable <- read.table(paste0("/fast/AG_Forslund/ChiaYu/Maaslin2_for_comparison/NEGBIN_mode/simulation_microbiomeDASim_2_conf1_testVar_confounder/data/Sim2_conf_data_1_set", i, "_meta.txt"),
                          sep = "\t", header = T)
  # N is the number of samples
  N = nrow(metatable)
  
    mem <- peakRAM(fit_data <- Maaslin2(
    input_data = paste0("/fast/AG_Forslund/ChiaYu/Maaslin2_for_comparison/NEGBIN_mode/simulation_microbiomeDASim_2_conf1_testVar_confounder/data/Sim2_conf_data_1_set", i, "_feature.txt"),
    input_metadata = paste0("/fast/AG_Forslund/ChiaYu/Maaslin2_for_comparison/NEGBIN_mode/simulation_microbiomeDASim_2_conf1_testVar_confounder/data/Sim2_conf_data_1_set", i, "_meta.txt"),
    output = paste0("/fast/AG_Forslund/ChiaYu/Maaslin2_for_comparison/NEGBIN_mode/simulation_microbiomeDASim_2_conf1_testVar_confounder_naive_model/Maaslin_result/Sim2_conf_data_1_set", i, "_testVar_confounder_naive_model"),
    min_abundance = 0,
    min_variance = 0,
    min_prevalence = 5/N, #At least 5 samples should have non-zero value
    max_significance = 0.1,
    fixed_effects = "Simulated_confounder",
    random_effects = "ID",
    analysis_method = "NEGBIN",
    correction = "BH",
    standardize = FALSE,
    normalization = 'NONE',
    transform = "NONE",
    plot_heatmap = FALSE,
    plot_scatter = FALSE
  ))
}

