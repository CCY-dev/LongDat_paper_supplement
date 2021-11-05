library(LongDat)

### Read in the input data
fasting <- read.table("/Users/Jessica/Documents/Lab/CORONA/corona.motu2matrix.rarefied.r.Species.r.master.r.dose_replaced.r",
                      sep = "\t", header = T)

### Run LongDat
set.seed(100)
test_disc <- longdat_disc(input = fasting,
                          data_type = "count",
                          test_var = "Case", #Case is the time variable
                          variable_col = 24,
                          fac_var = c(1:6),
                          not_used =  c(4))


### Plot the results
test_plot_disc <- cuneiform_plot(result_table = test_disc[[1]])

