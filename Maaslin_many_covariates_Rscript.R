library(tidyverse)
library(ggpubr)
library(peakRAM)
library(Maaslin2)

### First generate simulated data using MicrobiomeDASim
args <- commandArgs(TRUE)

num_start = args[1]
num_end = args[2]

for (i in num_start:num_end) {
  sim <- read.table(paste0("/fast/AG_Forslund/ChiaYu/LongDat_package_simulation_test/Many_random_covariates_cont/16covariate/Many_random_16covariate_sim19_cont/simulated_many_random_covariate_data/Many_random_16covariates_sim19_set", i, ".txt"),
                    sep = "\t", header = T)
  
  # Only the treatment group
  sim <- sim %>% 
    dplyr::filter(group == "Treatment")
  
  # Add letter to ID so that maaslin won't take it as number
  sim$ID <- paste0(sim$ID, "a")
  
  meta <- as.data.frame(sim[ , 1:19]) %>%
    dplyr::mutate(Sample = paste0(ID, "_", time)) %>% 
    tibble::column_to_rownames("Sample") %>% 
    dplyr::select(-3)

  feature <- as.data.frame(sim[ , c(1, 2, 20:ncol(sim))]) %>%
    dplyr::mutate(Sample = paste0(ID, "_", time)) %>% 
    tibble::column_to_rownames("Sample") %>% 
    dplyr::select(-c(1:2))

  feature_round <- round(feature, 0)
  
  N = nrow(meta)

  write.table(meta, paste0("/fast/AG_Forslund/ChiaYu/Maaslin2_for_comparison/NEGBIN_mode/Many_random_covariates_maaslin_NB/16covariate/Many_random_16covariate_sim19_NB/data/Many_random_16covariates_sim19_set", i, "_meta.txt"), sep = "\t",
              row.names = T, col.names = NA, quote = F)
  write.table(feature_round, paste0("/fast/AG_Forslund/ChiaYu/Maaslin2_for_comparison/NEGBIN_mode/Many_random_covariates_maaslin_NB/16covariate/Many_random_16covariate_sim19_NB/data/Many_random_16covariates_sim19_set", i, "_feature.txt"), sep = "\t",
              row.names = T, col.names = NA, quote = F)
  
  fixed_effects_names <- colnames(meta)[2:length(colnames(meta))]
  
  mem <- peakRAM(fit_data <- Maaslin2(
    input_data = paste0("/fast/AG_Forslund/ChiaYu/Maaslin2_for_comparison/NEGBIN_mode/Many_random_covariates_maaslin_NB/16covariate/Many_random_16covariate_sim19_NB/data/Many_random_16covariates_sim19_set", i, "_feature.txt"),
    input_metadata = paste0("/fast/AG_Forslund/ChiaYu/Maaslin2_for_comparison/NEGBIN_mode/Many_random_covariates_maaslin_NB/16covariate/Many_random_16covariate_sim19_NB/data/Many_random_16covariates_sim19_set", i, "_meta.txt"),
    output = paste0("/fast/AG_Forslund/ChiaYu/Maaslin2_for_comparison/NEGBIN_mode/Many_random_covariates_maaslin_NB/16covariate/Many_random_16covariate_sim19_NB/Maaslin_result/Many_random_16covariates_sim19_set", i, "_Maaslin_NEGBIN"),
    min_prevalence = 5/N,
    min_variance = 0.0,
    min_abundance = 0.0,
    max_significance = 0.1,
    fixed_effects = fixed_effects_names,
    random_effects = "ID",
    analysis_method = "NEGBIN",
    correction = "BH",
    standardize = FALSE,
    normalization = 'NONE',
    transform = "NONE",
    plot_heatmap = FALSE,
    plot_scatter = FALSE
  ))
  write.table(x = mem, file =paste0("/fast/AG_Forslund/ChiaYu/Maaslin2_for_comparison/NEGBIN_mode/Many_random_covariates_maaslin_NB/16covariate/Many_random_16covariate_sim19_NB/Maaslin_result/Many_random_16covariates_sim19_set", i, "_Maaslin_NEGBIN_MemTime.txt"), 
              sep = "\t", row.names = F, col.names = T, quote = F)
}
