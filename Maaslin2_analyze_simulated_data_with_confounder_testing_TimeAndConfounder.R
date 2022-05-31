library(tidyverse)
library(ggpubr)
library(peakRAM)
library(Maaslin2)

### First generate simulated data using MicrobiomeDASim
for (i in 1:100) {
  sim <- read.table(paste0("/fast/AG_Forslund/ChiaYu/LongDat_package_simulation_test/microbiomeDASim/simulation_microbiomeDASim_2_conf1/simulated_conf_data/Sim2_conf_data_1_set", i, ".txt"),
                    sep = "\t", header = T)
  
  # Select only the treatment group
  sim <- sim %>% 
    dplyr::filter(group == "Treatment")
  
  # Add letter to ID so that maaslin won't take it as number
  sim$ID <- paste0(sim$ID, "a")
  
  meta <- as.data.frame(sim[ , 1:4]) %>%
    dplyr::mutate(Sample = paste0(ID, "_", time)) %>% 
    tibble::column_to_rownames("Sample") %>% 
    dplyr::select(-c(4))

  feature <- as.data.frame(sim[ , c(1, 2, 3, 5:ncol(sim))]) %>%
    dplyr::mutate(Sample = paste0(ID, "_", time)) %>% 
    tibble::column_to_rownames("Sample") %>% 
    dplyr::select(-c(1:3))

  feature_round <- round(feature, 0)

  # N is the number of samples
  N = nrow(meta)
  
  write.table(meta, paste0("/fast/AG_Forslund/ChiaYu/Maaslin2_for_comparison/NEGBIN_mode/simulation_microbiomeDASim_2_conf1_testVar_confounder/data/Sim2_conf_data_1_set", i, "_meta.txt"), sep = "\t",
              row.names = T, col.names = NA, quote = F)
  write.table(feature_round, paste0("/fast/AG_Forslund/ChiaYu/Maaslin2_for_comparison/NEGBIN_mode/simulation_microbiomeDASim_2_conf1_testVar_confounder/data/Sim2_conf_data_1_set", i, "_feature.txt"), sep = "\t",
              row.names = T, col.names = NA, quote = F)
  
  # No fixed effect specified means that all meta columns will be fixed effects
  mem <- peakRAM(fit_data <- Maaslin2(
    input_data = paste0("/fast/AG_Forslund/ChiaYu/Maaslin2_for_comparison/NEGBIN_mode/simulation_microbiomeDASim_2_conf1_testVar_confounder/data/Sim2_conf_data_1_set", i, "_feature.txt"),
    input_metadata = paste0("/fast/AG_Forslund/ChiaYu/Maaslin2_for_comparison/NEGBIN_mode/simulation_microbiomeDASim_2_conf1_testVar_confounder/data/Sim2_conf_data_1_set", i, "_meta.txt"),
    output = paste0("/fast/AG_Forslund/ChiaYu/Maaslin2_for_comparison/NEGBIN_mode/simulation_microbiomeDASim_2_conf1_testVar_confounder/Maaslin_result/Sim2_conf_data_1_set", i, "_testVar_confounder_Maaslin_NEGBIN"),
    min_abundance = 0,
    min_variance = 0,
    min_prevalence = 5/N, #At least 5 samples should have non-zero value
    max_significance = 0.1,
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

