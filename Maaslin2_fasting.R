library(Maaslin2)
library(tidyverse)
library(dplyr)
data1 <- read.table("Fasting_species.txt",
                   sep = "\t", header = T)
feature <- data1[ , c(4, 24:ncol(data1))] %>%
  column_to_rownames("SampleID")
meta <- data1[ , c(1:23)] %>%
  column_to_rownames("SampleID")

meta$Patient <- paste0(meta$Patient, "a")
meta$Case <- paste0(meta$Case, "a")
meta$corona.dash <- NULL

write.table(feature, "fasting_feature.txt",
            sep = "\t", quote = F, row.names = T, col.names = NA)
write.table(meta, "fasting_meta.txt",
            sep = "\t", quote = F, row.names = T, col.names = NA)
N = nrow(meta)

### No covariates included in the run
fit_data <- Maaslin2(
  input_data = "fasting_feature.txt",
  input_metadata = "fasting_meta.txt",
  output = "Masslin2_NEGBIN_fasting_species",
  min_prevalence = 5/N,
  min_variance = 0.0,
  min_abundance = 0.0,
  max_significance = 0.1,
  random_effects = "Patient",
  analysis_method = "NEGBIN",
  correction = "BH",
  fixed_effects = "Case",  # "Case" means time point
  reference = "Case,1a",
  standardize = FALSE,
  normalization = 'NONE',
  transform = "NONE",
  plot_heatmap = FALSE,
  plot_scatter = FALSE
  )


### Covariates included in the run
fit_data2 <- Maaslin2(
  input_data = "fasting_feature.txt",
  input_metadata = "fasting_meta.txt",
  output = "Masslin2_NEGBIN_fasting_species_with_covariates",
  min_prevalence = 5/N,
  min_variance = 0.0,
  min_abundance = 0.0,
  max_significance = 0.1,
  random_effects = "Patient",
  analysis_method = "NEGBIN",
  correction = "BH",
  fixed_effects = c("Case", "Status135", "sex", "age", "ACE_INHIBITOR", "ANTIDEPRESSANT",
                    "ANTIDIABETIC", "ANTI_PLATELET", "AT2_BLOCKER", "BETA_BLOCKER",
                    "CALCIUMANTAGONIST", "DIURETIC", "IMIDAZOLE_AGONIST",
                    "METFORMIN", "PPI", "STATIN", "THYROID", "VITAMIN_D", "XANTHAN_OXIDASE_INHIBITOR"), # "Case" means time point
  reference = c("Case,1a", "Status135,non-R"),
  standardize = FALSE,
  normalization = 'NONE',
  transform = "NONE",
  plot_heatmap = FALSE,
  plot_scatter = FALSE
)
