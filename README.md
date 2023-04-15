# LongDat_paper_supplement

This repository holds the supplementary materials for the LongDat paper.

# Description of each supplementary material

| Name                                                     | Description                                                                                      |
|----------------------------------------------------------|--------------------------------------------------------------------------------------------------|
| MicrobiomeDASim_generate_simulated_data.R                | Generate simulated longitudinal data by using microbiomeDASim                                    |
| MicrobiomeDASim_generate_simulated_data_add_covariate.R | Add a dummy variable correlating with the time variable to the data generated by microbiomeDASim |
| LongDat_analyze_simulated_data.R                         | Run LongDat on the data created in "MicrobiomeDASim_generate_simulated_data.R"                   |
| LongDat_analyze_simulated_data_with_covariate.R         | Run LongDat on the data created in "MicrobiomeDASim_generate_simulated_data_add_covariate.R"    |
| LongDat_fasting_example.R                                | Run LongDat on "Fasting_species.txt"                                                             |
| Fasting_species.txt                                      | The fasting gut microbial abundance table at species level. Case means time points in this study |
|immunome.master.r.dose.r                                  |The fasting immunome data. Case means time points in this study|
|LongDat_fasting_immunome_for_paper.R|                      Run LongDat on the immunome data|
| LongDat_fasting_result_table.txt                         | The result table from running LongDat on "Fasting_species.txt"                                   |
| LongDat_fasting_covariate_table.txt                     | The covariate table from running LongDat on "Fasting_species.txt"                               |
| Fasting_immunome_measurement_Longdat_disc_result_table.txt| The result table from running LongDat on the non-percentage data in"immunome.master.r.dose.r"   |
| Fasting_immunome_measurement_Longdat_disc_covariates.txt| The covariate table from running LongDat on the non-percentage data in"immunome.master.r.dose.r"   |
| Fasting_immunome_percentage_Longdat_disc_result_table.txt| The result table from running LongDat on the percentage data in"immunome.master.r.dose.r"   |
| Fasting_immunome_percentage_Longdat_disc_covariates.txt| The covariate table from running LongDat on the percentage data in"immunome.master.r.dose.r"   |
| Maaslin2_analyze_simulated_data.R                        | Run Maaslin2 on the data created in "MicrobiomeDASim_generate_simulated_data.R"                  |
| Maaslin2_analyze_simulated_data_with_covariate_testing_Covariate.R | Run Maaslin2 on the data created in "MicrobiomeDASim_generate_simulated_data_add_covariate.R" and testing the simulated_confounder as fixed effect.|
| Maaslin2_analyze_simulated_data_with_covariate_testing_TimeAndCovariate.R | Run Maaslin2 on the data created in "MicrobiomeDASim_generate_simulated_data_add_covariate.R" and testing both simulated_confounder and time as fixed effect.|
|Maaslin2_fasting.R                                         | Run Maaslin2 on "Fasting_species.txt"                                                           |
|Masslin2_NEGBIN_fasting_species_no_covariate_case_full_result.txt | The result from running Maaslin2_fasting.R. Here no potential covariates were included in the run. Case means time points. |
|Masslin2_NEGBIN_fasting_species_with_covariate_case_full_result.txt |  The result from running Maaslin2_fasting.R. Here potential covariates were included in the run. Case means time points. |
| Maaslin_many_covariates_Rscript.R|          Run Maaslin2 on multiple covariate simulated data|
|Longdat_many_random_covariates_Rscript.R| Run LongDat on multiple covariate simulated data|
|Shuffle_time_run_longdat_Rscript.R| Run LongDat on negative control data with the time shuffled against other variables|
|Maaslin2_TSS_Rscript.R| Run Maaslin2 by using TSS and linear model mode on simulated data|
|lgpr_Rscript.R| Run lgpr on multiple covariate simulated data|
|Ancom_Rscript.R| Run ANCOM on multiple covariate simulated data|
|ZIBR_Rscript.R| Run ZIBR on multiple covariate simulated data|
|Fasting_genus_raw_master.txt| The microbial abundance data at genus level in the fasting study|
|Fasting_genus_rarefied_master.txt| The rarefied microbial abundance data at genus level in the fasting study|
|Fasting_genus_tss_master.txt| The TSS-normalized microbial abundance data at genus level in the fasting study|
|Fasting_genus_clr_master.txt| The CLR-transformed microbial abundance data at genus level in the fasting study|
|Fasting_genus_gmpr_master.txt| The GMPR-normalized microbial abundance data at genus level in the fasting study|
|Fasting_genus_tmm_master.txt| The TMM-normalized microbial abundance data at genus level in the fasting study|
|Fasting_genus_css_master.txt| The CSS-normalized microbial abundance data at genus level in the fasting study|
|Generate_shuffled_fasting_data.R| Shuffles the fasting genus data for suppl. fig. 16|
|SpasreDOSSA2_longitudinal_raw_simulation.R|This script generates the longitudinal simulated data using SparseDOSSA2|
|SpasreDOSSA2_rarefy_and_normalize.R| Rarefy/normalize the simulated data of SparseDOSSA2|
|SpasreDOSSA2_add_covariate.R| Adds covariates to the simulated data of SparseDOSSA2|
|Generate_negative_control_data.R|The R script used to generate negative control data for suppl. fig. 9|




