rm(list = ls())
library(tidyverse)
library(SparseDOSSA2)
setwd("/Users/Jessica/Documents/Lab/SparseDossa2/Simulation_data/sim1")

# The following is for simulating longitudinal data
# with SparseDOSSA 2, using the "Stool" fit. The spiked microbes are associated with a case-control 
# variable, and each is associated with its own Gaussian "random effect" 
# variable that's shared across time points within the same subject.

n_t <- 2 # 2 time points
n_s <- 20 # number of subjects
effect_size = 1.5 # 1.5 corresponds to median being 0.2 rho; 4 corresponds to median being 0.5 rho

# Create a random SparseDossa2 simulated df first
set.seed(101)
random_sim <- SparseDOSSA2(template = "Stool",  new_features = FALSE)
random_sim_df <- random_sim$simulated_data

# Take the microbe names from the random_sim_dif for my real simuation
# SparseDossa2 generates 332 microbes, I want 10% of them are spiked so that's 33
# 17 of them increase and 16 of them decrease
feature_up <- rownames(random_sim_df)[1:17]
feature_down <- rownames(random_sim_df)[18:33]
n_spiked <- length(feature_up) + length(feature_down)

for (i in 1:100) {
  set.seed(i)
  
  # Metadata matrix for spiking in. This is three columns: the first column is
  # case-control variable (here it's time variable). The rests are subject-specific random effects,
  # each corresponding to one spiked microbe.
  mat_metadata <- cbind(
    rep(rep(c(0, 1)), times = n_s), # time variable
    vapply(seq(1, n_spiked), # subject-specific random effects (a metadata column for each feature)
           function(i) {
             rep(rnorm(n = n_s), each = n_t)
           },
           rep(0.0, length = n_s * n_t)) # two random effect variable, each
    # corresponding to one spiked microbe
  )
  colnames(mat_metadata) <- c("time",
                              paste0("subject_feature_", seq(1, n_spiked)))
  
  rownames(mat_metadata) <- paste0("Sample_", rep(1:(n_s), each = n_t), "_", mat_metadata[ , 1])
  
  
  # Create feature-metadata spike-in data frame, per SparseDOSSA 2 standard.
  metadata_spike_1 <- data.frame(feature_spiked = c(feature_up),
                                 effect_size = effect_size)
  metadata_spike_2 <- data.frame(feature_spiked = c(feature_down),
                                 effect_size = -effect_size)
  metadata_spike_3 <-  data.frame(metadata_datum = 1)
  metadata_spike_4 <- tidyr::crossing(metadata_spike_3, metadata_spike_1) #create all combinations
  metadata_spike_5 <- tidyr::crossing(metadata_spike_3, metadata_spike_2) #create all combinations
  metadata_spike_6 <- rbind(metadata_spike_4, metadata_spike_5) # This is associating with "time"
  
  metadata_spike_7 <- data.frame(feature_spiked = c(feature_up, feature_down),
                                 effect_size = 1)
  metadata_spike_8 <-  data.frame(metadata_datum = 2:(n_spiked + 1))
  metadata_spike_9 <- cbind(metadata_spike_8, metadata_spike_7) # This is associating with the subject-specific random effects (a metadata column for each feature)
  
  df_metadata_spike <- rbind(metadata_spike_6, metadata_spike_9) %>% 
    tidyr::crossing(associated_property = c("prevalence", "abundance"))
  
  # Start simulation
  sim_SparseDOSSA2 <- 
    SparseDOSSA2::SparseDOSSA2(
      template = "Stool",
      new_features = FALSE, 
      n_sample = nrow(mat_metadata), 
      spike_metadata = "both",
      metadata_matrix = mat_metadata, 
      feature_metadata_spike_df = df_metadata_spike,
      verbose = FALSE)
  
  simulated_df <- sim_SparseDOSSA2$simulated_data
  spiked_microbe <- sim_SparseDOSSA2[["spike_metadata"]][["feature_metadata_spike_df"]][["feature_spiked"]]
  #n_distinct(spiked_microbe)
  #match(unique(spiked_microbe), rownames(simulated_df))
  
  #sort(unname(rowSums(simulated_df))) # Abundance per feature
  sort(unname(colSums(simulated_df))) #Sequencing depth per sample
  
  simulated_df_t <- as.data.frame(t(simulated_df)) %>% 
    mutate(.before = 1, 
           ID = str_split_fixed(colnames(simulated_df), pattern = "_", n = 3)[ , 2],
           time = str_split_fixed(colnames(simulated_df), pattern = "_", n = 3)[ , 3])
  
  simulated_df_t$ID <- as.factor(simulated_df_t$ID)
  simulated_df_t$time <-  as.numeric(simulated_df_t$time)
  
  # Change the names of the microbes for downstream convenience
  colnames(simulated_df_t)[3:(2+n_spiked)] <- paste0("Diff_Bug", 1:33)
  colnames(simulated_df_t)[(3+n_spiked):ncol(simulated_df_t)] <- paste0("NoDiffBug_", 1:(ncol(simulated_df_t)-n_spiked -2))
  
  # Change the ID names so that they are characters
  simulated_df_t$ID <- paste0("s_", simulated_df_t$ID)
  
  write.table(simulated_df_t, paste0("/Users/Jessica/Documents/Lab/SparseDossa2/Simulation_data/sim1/Raw/Sparsedossa2_sim1_set", i, ".txt"), sep = "\t",
              row.names = F, col.names = T, quote = F)
}


