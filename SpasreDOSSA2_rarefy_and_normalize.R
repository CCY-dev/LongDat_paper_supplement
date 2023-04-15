rm(list = ls())
setwd("/Users/Jessica/Documents/Lab/SparseDossa2/Simulation_data/")
library(tidyverse)
library(rtk)
library(compositions)
library(GMPR)
library(edgeR)
library(MetaLonDA)

# Do normalization/rarefaction on sim1 ~ sim12
sim_name <- paste0("sim", 1:12)

for (j in 1:length(sim_name)) {
  this_sim <- sim_name[j]
  print(paste0("j = ", j))
  for (i in 1:100) {
    print(paste0("i = ", i))
    sim_raw <- read.table(paste0("/Users/Jessica/Documents/Lab/SparseDossa2/Simulation_data/", this_sim, "/Raw/Sparsedossa2_", this_sim, "_set", i, ".txt"),
                          sep = "\t", header = T)
    sim_raw_features <- sim_raw[ , -c(1:2)]
    
    #### Rarefaction
    # depth = 0.95*lowest_depth because it's the same as the value on RTK using server
    rarefied_res <- rtk::rtk(sim_raw_features, depth = min(rowSums(sim_raw_features)) * 0.95,
                             ReturnMatrix = 1, seed = 42, margin = 1)
    rarefied_features <- as.data.frame(rarefied_res$raremat)
    #rowSums(rarefied_features)
    rarefied_df <- cbind(sim_raw[ , c(1:2)], rarefied_features)
    
    write.table(rarefied_df, paste0("/Users/Jessica/Documents/Lab/SparseDossa2/Simulation_data/", this_sim, "/Rarefied/Sparsedossa2_", this_sim, "_set", i, "_rarefied.txt"), sep = "\t",
                row.names = F, col.names = T, quote = F)
    
    #### TSS (same as converting to relative abundance)
    tss_features <- sim_raw_features/rowSums(sim_raw_features)
    #rowSums(tss_features)
    tss_df <- cbind(sim_raw[ , c(1:2)], tss_features)
    
    write.table(tss_df, paste0("/Users/Jessica/Documents/Lab/SparseDossa2/Simulation_data/", this_sim, "/TSS/Sparsedossa2_", this_sim, "_set", i, "_tss.txt"), sep = "\t",
                row.names = F, col.names = T, quote = F)
    
    #### CLR 
    # Add a pseudocount 1 before CLR
    clr_features <- as.data.frame(compositions::clr(x = (sim_raw_features + 1)))
    #clr_features2 <- as.data.frame(clr(x = (sim_raw_features + 1)/rowSums(sim_raw_features + 1))) # using relative abundance gives the same result as above
    #rowSums(clr_features)
    clr_df <- cbind(sim_raw[ , c(1:2)], clr_features)
    
    write.table(clr_df, paste0("/Users/Jessica/Documents/Lab/SparseDossa2/Simulation_data/", this_sim, "/CLR/Sparsedossa2_", this_sim, "_set", i, "_clr.txt"), sep = "\t",
                row.names = F, col.names = T, quote = F)
    
    #### GMPR 
    #https://github.com/jchen1981/GMPR/blob/master/GMPR.Example.R
    gmpr_size_factor <- GMPR::GMPR(sim_raw_features)
    gmpr_features <- sim_raw_features/gmpr_size_factor #Counts are normalized by size factors
    gmpr_df <- cbind(sim_raw[ , c(1:2)], gmpr_features)
    
    write.table(gmpr_df, paste0("/Users/Jessica/Documents/Lab/SparseDossa2/Simulation_data/", this_sim, "/GMPR/Sparsedossa2_", this_sim, "_set", i, "_gmpr.txt"), sep = "\t",
                row.names = F, col.names = T, quote = F)
    
    #### TMM 
    #https://www.biostars.org/p/9475236/
    #make the DGEList
    sim_dge <- edgeR::DGEList(as.matrix(t(sim_raw_features)), #row as features, col as samples 
                              samples = as.data.frame(sim_raw[ , c(1:2)]))
    #calculate TMM normalization factors
    sim_dge <- edgeR::calcNormFactors(sim_dge, method = "TMM")
    #get the normalized counts:
    cpms <- edgeR::cpm(sim_dge, log=FALSE, normalized.lib.sizes = T)
    tmm_features <- t(cpms)
    #rowSums(tmm_features)
    tmm_df <- cbind(sim_raw[ , c(1:2)], tmm_features)
    
    write.table(tmm_df, paste0("/Users/Jessica/Documents/Lab/SparseDossa2/Simulation_data/", this_sim, "/TMM/Sparsedossa2_", this_sim, "_set", i, "_tmm.txt"), sep = "\t",
                row.names = F, col.names = T, quote = F)
    
    #### CSS (Cumulative Sum Scaling)
    # https://github.com/aametwally/MetaLonDA
    #https://www.rdocumentation.org/packages/MetaLonDA/versions/1.1.8/topics/normalize
    # Sometimes CSS returns error (don't know why), so uses "try" to make sure if there's error
    # https://stackoverflow.com/questions/11316609/how-can-i-determine-if-try-returned-an-error-or-not
     css_try <- try(t(MetaLonDA::normalize(count = t(sim_raw_features)))) #samples as col, feature as row
     
      if (class(css_try) == "try-error") { #if there's error, and 1 pseudocount
        css_features <- t(MetaLonDA::normalize(count = t(sim_raw_features + 1))) # add a pseudocount
      } else { #if no error, run as usual
        css_features <- t(MetaLonDA::normalize(count = t(sim_raw_features)))
      }
      
    css_df <- cbind(sim_raw[ , c(1:2)], css_features)
    
    write.table(css_df, paste0("/Users/Jessica/Documents/Lab/SparseDossa2/Simulation_data/", this_sim, "/CSS/Sparsedossa2_", this_sim, "_set", i, "_css.txt"), sep = "\t",
                row.names = F, col.names = T, quote = F)
  }
}

