rm(list = ls())
library(tidyverse)

# Add covariates 1, 4, 16 to all combinations
sim_name <- paste0("sim", 1:12)
type_name <- c("CLR", "CSS", "GMPR", "Rarefied", "Raw", "TMM", "TSS")
type_name_file_tag <-  c("_clr", "_css", "_gmpr", "_rarefied", "", "_tmm", "_tss")

#Simulate covariates with desired rho (between time and the covariate)
#Ref: https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables
complement <- function(y, rho, x) {
  if (missing(x)) x <- rnorm(length(y)) # Optional: supply a default if `x` is not given
  y.perp <- residuals(lm(x ~ y))
  rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
}

for (j in 1:length(sim_name)) {
  this_sim <- sim_name[j]
  print(paste0("j = ", j))
  
  for (n in 1:length(type_name)) {
    this_type <- type_name[n]
    this_type_name_file_tag <- type_name_file_tag[n]
    print(paste("type =", this_type))
    
    for (i in 1:100) {
      print(paste0("i = ", i))
      
      raw <- read.table(paste0("/fast/AG_Forslund/ChiaYu/Benchmark_sparsedossa2/simulation_data/", this_sim, "/", this_type ,"/Sparsedossa2_", this_sim, "_set", i, this_type_name_file_tag, ".txt"),
                            sep = "\t", header = T)
      
      ### Add 1 covariate
      set.seed(i*j*n)
      # The rho are 0.25, 0.5, 0.75, 0.99, randomly pick the 1 rhos from these 4 levels
      rho_levels <- c(0.25, 0.5, 0.75, 0.99)
      my_levels <- sample(x = rho_levels, size = 1, replace = T)
      
      a1 = complement(y = raw$time, rho = my_levels[1])
      
      added_1cov <- raw %>%
        mutate(.after = 2,
               Simulated_covariate_a1 = a1
        )
      
      write.table(added_1cov, paste0("/fast/AG_Forslund/ChiaYu/Benchmark_sparsedossa2/simulation_data/", this_sim, "/", this_type, "/Sparsedossa2_", this_sim, "_set", i, "_", this_type, "_1cov.txt"), sep = "\t",
                  row.names = F, col.names = T, quote = F)
      
      ### Add 4 covariates
      set.seed(i*j*n)
      # The rho are 0.25, 0.5, 0.75, 0.99, randomly pick the 4 rhos from these 4 levels
      rho_levels <- c(0.25, 0.5, 0.75, 0.99)
      my_levels <- sample(x = rho_levels, size = 4, replace = T)
      
      a1 = complement(y = raw$time, rho = my_levels[1])
      a2 = complement(y = raw$time, rho = my_levels[2])
      a3 = complement(y = raw$time, rho = my_levels[3])
      a4 = complement(y = raw$time, rho = my_levels[4])
      
      added_4cov <- raw %>%
        mutate(.after = 2,
               Simulated_covariate_a1 = a1,
               Simulated_covariate_a2 = a2,
               Simulated_covariate_a3 = a3,
               Simulated_covariate_a4 = a4
        )
      
      write.table(added_4cov, paste0("/fast/AG_Forslund/ChiaYu/Benchmark_sparsedossa2/simulation_data/", this_sim, "/", this_type, "/Sparsedossa2_", this_sim, "_set", i, "_", this_type, "_4cov.txt"), sep = "\t",
                  row.names = F, col.names = T, quote = F)
      
      ### Add 16 covariates
      set.seed(i*j*n)
      # The rho are 0.25, 0.5, 0.75, 0.99, randomly pick the 16 rhos from these 4 levels
      rho_levels <- c(0.25, 0.5, 0.75, 0.99)
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
      
      added_16cov <- raw %>%
        mutate(.after = 2,
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
               Simulated_covariate_a16 = a16
        )
      
      write.table(added_16cov, paste0("/fast/AG_Forslund/ChiaYu/Benchmark_sparsedossa2/simulation_data/", this_sim, "/", this_type, "/Sparsedossa2_", this_sim, "_set", i, "_", this_type, "_16cov.txt"), sep = "\t",
                  row.names = F, col.names = T, quote = F)
    }
  }
}

