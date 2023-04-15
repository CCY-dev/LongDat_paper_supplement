# This is executed using the R on cluster conda env longdat
library(tidyverse)
library(permute)

modes <- c("gmpr", "rarefied", "raw", "tmm", "tss", "clr", "css")

for (j in 1:length(modes)) {
  print(modes[j])
  this_mode <- modes[j]

for (i in 1:100) {
  print(i)
  
  input_table <- read.table(paste0("/fast/AG_Forslund/ChiaYu/Fasting_negative_validation/fasting_data/correct/Fasting_genus_", this_mode, "_master.txt"),
                            sep = "\t", header = T)
  
  input_table <- input_table %>% filter(Case != 3)
    
    
  # Shuffle the time variable in stratified manner (within each individual)
  set.seed(i)
  random_order = shuffle(input_table$Case, control = how(blocks = input_table$Patient))
  input_table <- input_table %>% 
    dplyr::mutate(Shuffled_case = input_table$Case[random_order], 
                  .after = 2)
  
  write.table(x = input_table, file = paste0("/fast/AG_Forslund/ChiaYu/Fasting_negative_validation/fasting_data/shuffled/", this_mode, "/Case_shuffled_fasting_genus_",this_mode, "_set", i, ".txt"), sep = "\t",
              row.names = F, col.names = T, quote = F)
}
}