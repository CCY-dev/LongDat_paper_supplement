library(microbiomeDASim)
library(tidyverse)
library(LongDat)
library(peakRAM)

### First generate simulated data using MicrobiomeDASim
for (i in 1:100) {
  set.seed(i)
  sim <- gen_norm_microbiome(features=200,
                             diff_abun_features=20,
                             n_control=10,
                             n_treat=300, #The number of treated samples (10, 20, 38, 75, 150, 300)
                             control_mean=0,
                             sigma=10,
                             num_timepoints=2,
                             t_interval=c(1,2),
                             rho=0.8,
                             corr_str="ar1",
                             func_form="linear",
                             beta=c(0, 5), #c(0, 5) for effect size median around 0.2 and c(0, 10) for effect size median around 0.5
                             missing_pct=0.4,
                             missing_per_subject=1,
                             miss_val=0)


  meta <- as.data.frame(sim[[2]]) %>%
    rownames_to_column("Sample")

  feature <- as.data.frame(t(sim[[1]]))

  feature_round <- round(feature, 0) %>%
    rownames_to_column("Sample")

  full_data_round <- inner_join(meta, feature_round, by = "Sample") %>%
    dplyr::select(-Sample_ID)

  # Take only the treatment group
  full_data_round_treatment <- full_data_round %>%
    filter(group == "Treatment") %>%
    dplyr::select(-Sample)

  write.table(full_data_round_treatment, paste0("microbiomeDASim/simulation_microbiomeDASim_3/simulated_data/Sim_microbiomeDASim_data_3_set", i, ".txt"), sep = "\t",
              row.names = F, col.names = T, quote = F)
}

