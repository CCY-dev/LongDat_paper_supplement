library(tidyverse)

for (i in 1:100) {
  print(i)
  #Load simulated data
  raw <- read.table(paste0("microbiomeDASim/simulation_microbiomeDASim_2/simulated_data/Sim_microbiomeDASim_data_2_set", i, ".txt"),
                    sep = "\t", header = T, stringsAsFactors = F)

  #Add a dummy variable correlating with time
  complement <- function(y, rho, x) {
    if (missing(x)) x <- rnorm(length(y)) # Optional: supply a default if `x` is not given
    y.perp <- residuals(lm(x ~ y))
    rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
  }

  y = raw$time
  rho = 0.25 #rho values are 0.25, 0.50, 0.75 or 0.99
  set.seed(111)
  a = complement(y = y, rho = rho)

  raw <- raw %>%
    mutate(Simulated_confounder = a,.before = 2)

  write.table(raw, paste0("microbiomeDASim/simulation_microbiomeDASim_2_conf1/simulated_conf_data/Sim2_conf_data_1_set", i, ".txt"), sep = "\t",
              row.names = F, col.names = T, quote = F)
  }


