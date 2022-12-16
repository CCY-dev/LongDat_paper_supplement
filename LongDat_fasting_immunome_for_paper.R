rm(list = ls())
setwd("~/Documents/Lab/Longdat_paper_writing/Figures/2022End_resubmission/")
library(tidyverse)
library(readxl)
library(patchwork)
library(LongDat)

#################################
raw <- read.table(file = "/Users/Jessica/Documents/Lab/CORONA/masters/immunome.master.r.dose.r", quote = "", sep = "\t", 
                  header = T, check.names = FALSE)
metadata <- raw[ , c(1:23)]
raw_value <- raw[ , -c(1:23)]
  
  
## Divide into percentage and non-percentage (measurement) data
# Only those columns with all values <= 1 are percentage data
percent_or_not <- c()
for (i in 1:ncol(raw_value)) {
  percent_or_not[i] <- ifelse(all(raw_value[ , i] <= 1, na.rm = T),  #all = Are All Values True?
                              yes = T, no = F)
}

raw_value_perc <- raw_value[ , percent_or_not == T]
raw_value_measure <- raw_value[ , percent_or_not == F]

# Because the proportional model needs 0 < y < 1 
# The smallest value in raw_value_perc is 0.0000144
# Replace 0 in raw_value_perc with 0.0000144 * 0.01
# Replace 1 in raw_value_perc with 1 - 0.0000144 * 0.01
raw_value_perc[raw_value_perc == 0] <- 0.0000144 * 0.01
raw_value_perc[raw_value_perc == 1] <- 1 - (0.0000144 * 0.01)


percentage_data <- cbind(metadata, raw_value_perc)
measurement_data <- cbind(metadata, raw_value_measure)

# LongDat percentage data
percent_disc <- longdat_disc(input = percentage_data,
                          data_type = "proportion",
                          test_var = "Case",
                          variable_col = 24,
                          fac_var = c(1:6),
                          not_used =  c(5),
                          adjustMethod = "fdr",
                          model_q = 0.1,
                          posthoc_q = 0.05,
                          theta_cutoff = 2^20,
                          nonzero_count_cutoff1 = 9,
                          nonzero_count_cutoff2 = 5,
                          verbose = T)
#write.table(x = percent_disc[[1]], file = "20221201_fasting_immunome_percentage_Longdat_disc_result_table.txt", sep = "\t",
#            row.names = F, col.names = T, quote = F)
#write.table(x = percent_disc[[2]], file = "20221201_fasting_immunome_percentage_Longdat_disc_covariates.txt", sep = "\t",
#            row.names = F, col.names = T, quote = F)
#saveRDS(percent_disc, "20221201_fasting_immunome_percentage_Longdat_disc.rds")

#cuneiform_plot
#cuneiform_plot_fasting <- cuneiform_plot
#fix(cuneiform_plot_fasting) #Manually add "scale_x_discrete(labels =c("Fasting", "Refeeding", "Study"))"
percent_disc <- readRDS("20221201_fasting_immunome_percentage_Longdat_disc.rds")
(percent_disc_plot <- cuneiform_plot_fasting(result_table = percent_disc[[1]], x_axis_order = 
                                       c("Effect_1_2", "Effect_2_3", "Effect_1_3")))

ggsave("20221201_fasting_immunome_percentage_Longdat_disc_result.pdf", device = "pdf", width = 12, height = 14)


# LongDat measurement data
measure_disc <- longdat_disc(input = measurement_data,
                             data_type = "measurement",
                             test_var = "Case",
                             variable_col = 24,
                             fac_var = c(1:6),
                             not_used =  c(5),
                             adjustMethod = "fdr",
                             model_q = 0.1,
                             posthoc_q = 0.05,
                             theta_cutoff = 2^20,
                             nonzero_count_cutoff1 = 9,
                             nonzero_count_cutoff2 = 5,
                             verbose = T)
write.table(x = measure_disc[[1]], file = "20221201_fasting_immunome_measurement_Longdat_disc_result_table.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)
write.table(x = measure_disc[[2]], file = "20221201_fasting_immunome_measurement_Longdat_disc_covariates.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)
#saveRDS(measure_disc, "20221201_fasting_immunome_measurement_Longdat_disc.rds")

#cuneiform_plot
measure_disc <- readRDS("20221201_fasting_immunome_measurement_Longdat_disc.rds")

(measure_disc_plot <- cuneiform_plot_fasting(result_table = measure_disc[[1]], x_axis_order = 
                                       c("Effect_1_2", "Effect_2_3", "Effect_1_3")))


#ggsave("20221201_fasting_immunome_measurement_Longdat_disc_result.pdf", device = "pdf", width = 10, height = 14)



