################################### Plot_Real_Data_Analyses.R ##################################### 
## This file obtains the results of the real data analyses and plots the results (see Section    ##
## 6.3-6.4).                                                                                     ##
##                                                                                               ##
## Note that to run the code, the file 'Perform_Real_Data_Analyses.R' should have been           ##
## executed. Alternatively, the content of the file 'Real_Data_Analyses.RDAta' should be loaded. ##
##                                                                                               ##
## In addition, to create the heatmaps, the functions in the files                               ##
## 'Helpfunctions_General.R', 'Helpfunctions_Real_Data.R' and                                    ##
## 'Helpfunctions_Performance_Measures_and_Plots.R' should be loaded.                            ## 
###################################################################################################

library(ggplot2)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(scales)

###################################################################################################
## Compute the entropy R2 of the LC models (Section 6.3.1)                                       ##
###################################################################################################

get_entropy(LC_2016_restricted)
get_entropy(LC_2017_restricted)
get_entropy(LC_2018_restricted)

get_entropy(LC_2016_ok)
get_entropy(LC_2017_ok)
get_entropy(LC_2018_ok)

get_entropy(LC_2016_recoded)
get_entropy(LC_2017_recoded)
get_entropy(LC_2018_recoded)

###################################################################################################
## Compute the standard deviation of tree-MILC's PPEs (Section 6.3.2)                            ##
###################################################################################################

sqrt(get_var_PPEs_treeMILC(treeMILC_2016_restricted))
sqrt(get_var_PPEs_treeMILC(treeMILC_2017_restricted))
sqrt(get_var_PPEs_treeMILC(treeMILC_2018_restricted))

sqrt(get_var_PPEs_treeMILC(treeMILC_2016_ok))
sqrt(get_var_PPEs_treeMILC(treeMILC_2017_ok))
sqrt(get_var_PPEs_treeMILC(treeMILC_2018_ok))

sqrt(get_var_PPEs_treeMILC(treeMILC_2016_recoded))
sqrt(get_var_PPEs_treeMILC(treeMILC_2017_recoded))
sqrt(get_var_PPEs_treeMILC(treeMILC_2018_recoded))

##################################################################################################
## Obtain results for population proportion estimates (Section 6.3.3)                           ##
##################################################################################################

# Fill data frame with model results
results_props <- data.frame(prop1 = NA, prop2 = NA, prop3 = NA, year = NA, model = NA, cov = NA)
ok <- "1. Without missing covariates"
recoded <- "2. Missing covariates using HMM approach"
restricted <- "3. Missing covariates with direct effects"

results_props <- rbind(results_props, c(get_proportions(LC_2016_ok), 2016, "LC", ok))
results_props <- rbind(results_props, c(get_proportions(LC_2017_ok), 2017, "LC", ok))
results_props <- rbind(results_props, c(get_proportions(LC_2018_ok), 2018, "LC", ok))
results_props <- rbind(results_props, c(get_proportions(LC_2016_recoded), 2016, "LC", recoded))
results_props <- rbind(results_props, c(get_proportions(LC_2017_recoded), 2017, "LC", recoded))
results_props <- rbind(results_props, c(get_proportions(LC_2018_recoded), 2018, "LC", recoded))
results_props <- rbind(results_props, c(get_proportions(LC_2016_restricted), 2016, "LC", restricted))
results_props <- rbind(results_props, c(get_proportions(LC_2017_restricted), 2017, "LC", restricted))
results_props <- rbind(results_props, c(get_proportions(LC_2018_restricted), 2018, "LC", restricted))
results_props <- rbind(results_props, c(get_proportions(treeMILC_2016_ok), 2016, "tree-MILC", ok))
results_props <- rbind(results_props, c(get_proportions(treeMILC_2017_ok), 2017, "tree-MILC", ok))
results_props <- rbind(results_props, c(get_proportions(treeMILC_2018_ok), 2018, "tree-MILC", ok))
results_props <- rbind(results_props, c(get_proportions(treeMILC_2016_recoded), 2016, "tree-MILC", recoded))
results_props <- rbind(results_props, c(get_proportions(treeMILC_2017_recoded), 2017, "tree-MILC", recoded))
results_props <- rbind(results_props, c(get_proportions(treeMILC_2018_recoded), 2018, "tree-MILC", recoded))
results_props <- rbind(results_props, c(get_proportions(treeMILC_2016_restricted), 2016, "tree-MILC", restricted))
results_props <- rbind(results_props, c(get_proportions(treeMILC_2017_restricted), 2017, "tree-MILC", restricted))
results_props <- rbind(results_props, c(get_proportions(treeMILC_2018_restricted), 2018, "tree-MILC", restricted))
results_props <- results_props[-1,]

## Add HMM, LFS and ER proportions
results_props$HMM_prop1 <- results_props$HMM_prop2 <- results_props$HMM_prop3 <- NA
results_props$LFS_prop1 <- results_props$LFS_prop2 <- results_props$LFS_prop3 <- NA
results_props$ER_prop1 <- results_props$ER_prop2 <- results_props$ER_prop3 <- NA

year <- c(2016, 2017, 2018)
HMM_prop1 <- c(0.508, 0.515, 0.522) #resp. 2016, 2017 and 2018
HMM_prop2 <- c(0.307, 0.302, 0.292) #,,
HMM_prop3 <- c(0.186, 0.183, 0.186) #,,
LFS_prop1 <- c(0.571, 0.569, 0.579) #,,
LFS_prop2 <- c(0.310, 0.306, 0.295) #,,
LFS_prop3 <- c(0.119, 0.125, 0.126) #,,
ER_prop1 <- c(0.494, 0.495, 0.502)  #,,
ER_prop2 <- c(0.306, 0.302, 0.292)  #,,
ER_prop3 <- c(0.200, 0.203, 0.206)  #,,

for (i in 1:3) {
  results_props[results_props$year == year[i], "HMM_prop1"] <- HMM_prop1[i]
  results_props[results_props$year == year[i], "HMM_prop2"] <- HMM_prop2[i]
  results_props[results_props$year == year[i], "HMM_prop3"] <- HMM_prop3[i]
  results_props[results_props$year == year[i], "LFS_prop1"] <- LFS_prop1[i]
  results_props[results_props$year == year[i], "LFS_prop2"] <- LFS_prop2[i]
  results_props[results_props$year == year[i], "LFS_prop3"] <- LFS_prop3[i]
  results_props[results_props$year == year[i], "ER_prop1"] <- ER_prop1[i]
  results_props[results_props$year == year[i], "ER_prop2"] <- ER_prop2[i]
  results_props[results_props$year == year[i], "ER_prop3"] <- ER_prop3[i]
}

# Convert LC and tree-MILC results to long format
results_LC_tM <- pivot_longer(results_props, cols = c("prop1", "prop2", "prop3"), names_to = "Contract")
results_LC_tM <- results_LC_tM[, -c(4:6)]
results_LC_tM$Contract <- as.factor(results_LC_tM$Contract)
levels(results_LC_tM$Contract) <- c("Permanent", "Other", "Flexible")
results_LC_tM$Contract <- factor(results_LC_tM$Contract, levels=c("Permanent", "Flexible", "Other"))
setnames(results_LC_tM, old = c("value"), new = c("prop"))
results_LC_tM$prop <- as.numeric(results_LC_tM$prop)

# Convert HMM results to long format
results_HMMs <- pivot_longer(results_props, cols = c("HMM_prop1", "HMM_prop2", "HMM_prop3"), names_to = "Contract")
results_HMMs$Contract <- as.factor(results_HMMs$Contract)
levels(results_HMMs$Contract) <- c("Permanent", "Other", "Flexible")
results_HMMs$Contract <- factor(results_HMMs$Contract, levels=c("Permanent", "Flexible", "Other"))
results_HMMs <- results_HMMs[, -c(1:3)]
setnames(results_HMMs, old = c("value"), new = c("HMM"))
results_HMMs$HMM <- as.numeric(results_HMMs$HMM)

# Convert LFS results to long format
results_LFS <- pivot_longer(results_props, cols = c("LFS_prop1", "LFS_prop2", "LFS_prop3"), names_to = "Contract")
results_LFS <- results_LFS[, c("year", "model", "cov", "Contract", "value")]
results_LFS$Contract <- as.factor(results_LFS$Contract)
levels(results_LFS$Contract) <- c("Permanent", "Other", "Flexible")
results_LFS$Contract <- factor(results_LFS$Contract, levels=c("Permanent", "Flexible", "Other"))
setnames(results_LFS, old = c("value"), new = c("prop"))

# Convert ER results to long format
results_ER <- pivot_longer(results_props, cols = c("ER_prop1", "ER_prop2", "ER_prop3"), names_to = "Contract")
results_ER <- results_ER[, c("year", "model", "cov", "Contract", "value")]
results_ER$Contract <- as.factor(results_ER$Contract)
levels(results_ER$Contract) <- c("Permanent", "Other", "Flexible")
results_ER$Contract <- factor(results_ER$Contract, levels=c("Permanent", "Flexible", "Other"))
setnames(results_ER, old = c("value"), new = c("prop"))

# Combine results for LC, tree-MILC, HMM, LFS, and ER
results_props_HMM <- left_join(results_LC_tM, results_HMMs, by = c("year", "model", "cov", "Contract", "LFS_prop1", "LFS_prop2", "LFS_prop3", "ER_prop1", "ER_prop2", "ER_prop3"))
results_props_LFS <- left_join(results_LFS, results_HMMs, by = c("year", "model", "cov", "Contract"))
results_props_ER <- left_join(results_ER, results_HMMs, by = c("year", "model", "cov", "Contract"))
results_HMMs$model <- "HMM"
results_props_LFS$model <- "LFS"
results_props_ER$model <- "ER"
results_HMMs <- distinct(results_HMMs)
results_props_LFS <- distinct(results_props_LFS)
results_props_ER <- distinct(results_props_ER)
results_HMMs$prop <- results_HMMs$HMM
combined_results <- rbind(results_props_HMM, results_HMMs)
combined_results <- rbind(combined_results, results_props_LFS)
combined_results <- rbind(combined_results, results_props_ER)
combined_results$cov <- as.factor(combined_results$cov)
combined_results$prop <- round(combined_results$prop, 3)
combined_results$model <- factor(combined_results$model, levels = c("LFS", "ER", "HMM", "LC", "tree-MILC"))

###################################################################################################
## Plot population proportion estimates (Section 6.3.1 and 6.4)                                  ##
###################################################################################################

group.colors <- c(LFS = "grey85", ER = "grey75", HMM = "grey60", LC = "#F8766D", `tree-MILC` = "#619CFF")

#Plots in Section 6.3: PPEs for (best) approach with direct effects and parameter restrictions
PPE_restricted <- combined_results[combined_results$cov == restricted,] %>%
  ggplot(aes(x = model, y = prop, fill = model)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  facet_grid(Contract ~ year, labeller = label_value, scales = "fixed") +
  scale_fill_manual(name = "", values = group.colors) +
  geom_hline(aes(yintercept = HMM, colour = "black"), colour = "black") +
  labs(y = "Population proportion estimates") + 
  labs(x = "") +
  geom_text(aes(label = prop), size = 3.5, hjust = 0.5, vjust = 2.7, position = "stack") +
  theme(legend.position = "top")

#Plots in Section 6.4: PPEs for additional analyses
PPE_2016 <- combined_results[combined_results$year == 2016,] %>%
  ggplot(aes(x = model, y = prop, fill = model)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  facet_grid(Contract ~ cov, labeller = label_value, scales = "fixed") +
  scale_fill_manual(name = "", values = group.colors) +
  labs(y = "Population proportion estimates") +
  labs(x = "") +
  geom_text(aes(label = prop), size = 3.5, hjust = 0.5, vjust = 2.85) +
  geom_hline(aes(yintercept = HMM, colour = "black"), colour = "black")

PPE_2017 <- combined_results[combined_results$year == 2017,] %>%
  ggplot(aes(x = model, y = prop, fill = model)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  facet_grid(Contract ~ cov, labeller = label_value, scales = "fixed") +
  scale_fill_manual(name = "", values = group.colors) +
  labs(y = "Population proportion estimates") +
  labs(x = "") +
  geom_text(aes(label = prop), size = 3.5, hjust = 0.5, vjust = 2.85) +
  geom_hline(aes(yintercept = HMM, colour = "black"), colour = "black")

PPE_2018 <- combined_results[combined_results$year == 2018,] %>%
  ggplot(aes(x = model, y = prop, fill = model)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  facet_grid(Contract ~ cov, labeller = label_value, scales = "fixed") +
  scale_fill_manual(name = "", values = group.colors) +
  labs(y = "Population proportion estimates") +
  labs(x = "") +
  geom_text(aes(label = prop), size = 3.5, hjust = 0.5, vjust = 2.85) +
  geom_hline(aes(yintercept = HMM, colour = "black"), colour = "black")


####################################################################################################
## Define measurement error probability matrices as estimated by the HMMs (Section 6.3.4 and 6.4) ##
####################################################################################################

ER.2016 <- matrix(c(0.932, 0.001, 0.067, 0.004, 0.996, 0, 0.107, 0.001, 0.892), ncol = 3, byrow = TRUE)
LFS.2016 <- matrix(c(0.985, 0.01, 0.005, 0.04, 0.934, 0.026, 0.332, 0.071, 0.597), ncol = 3, byrow = TRUE)
HMM_ME.2016 <- list(ER.2016, LFS.2016)

ER.2017 <- matrix(c(0.929, 0.001, 0.07, 0.001, 0.998, 0.001, 0.089, 0, 0.911), ncol = 3, byrow = TRUE)
LFS.2017 <- matrix(c(0.982, 0.014, 0.005, 0.042, 0.928, 0.03, 0.296, 0.071, 0.633), ncol = 3, byrow = TRUE)
HMM_ME.2017 <- list(ER.2017, LFS.2017)

ER.2018 <- matrix(c(0.931, 0.002, 0.068, 0.004, 0.996, 0, 0.083, 0.001, 0.916), ncol = 3, byrow = TRUE)
LFS.2018 <- matrix(c(0.982, 0.011, 0.007, 0.04, 0.933, 0.027, 0.302, 0.074, 0.623), ncol = 3, byrow = TRUE)
HMM_ME.2018 <- list(ER.2018, LFS.2018)

###################################################################################################
## Obtain ME probability estimates for the LC and the tree-MILC models (Section 6.3.4 and 6.4)   ##                                                          ##
###################################################################################################

## Note that 1 = Permanent, 2 = Other, and 3 = Flexible.

# ME probability matrices for (best) approach with direct effects and parameter restrictions (Section 6.1-6.3)
get_ME(LC_2016_restricted)
get_ME(LC_2017_restricted)
get_ME(LC_2018_restricted)
get_ME(treeMILC_2016_restricted)
get_ME(treeMILC_2017_restricted)
get_ME(treeMILC_2018_restricted)

# ME probability matrices for additional approach without missing covariates (Section 6.4)
get_ME(LC_2016_ok)
get_ME(LC_2017_ok)
get_ME(LC_2018_ok)
get_ME(treeMILC_2016_ok)
get_ME(treeMILC_2017_ok)
get_ME(treeMILC_2018_ok)

# ME probability matrices for additional approach without HMM recodings (Section 6.4)
get_ME(LC_2016_recoded)
get_ME(LC_2017_recoded)
get_ME(LC_2018_recoded)
get_ME(treeMILC_2016_recoded)
get_ME(treeMILC_2017_recoded)
get_ME(treeMILC_2018_recoded)

###################################################################################################
## Plot differences in ME probability estimates (Section 6.3.4 and 6.4)                          ##
###################################################################################################

# Differences in ME probability estimates for (best) approach with direct effects and parameter restrictions (Section 6.1-6.3)
HM_LC_2016_restricted <- get_ME_heatmap(get_ME(LC_2016_restricted), HMM_ME.2016) 
HM_LC_2017_restricted <- get_ME_heatmap(get_ME(LC_2017_restricted), HMM_ME.2017) 
HM_LC_2018_restricted <- get_ME_heatmap(get_ME(LC_2018_restricted), HMM_ME.2018)
HM_treeMILC_2016_restricted <- get_ME_heatmap(get_ME(treeMILC_2016_restricted), HMM_ME.2016)  
HM_treeMILC_2017_restricted <- get_ME_heatmap(get_ME(treeMILC_2017_restricted), HMM_ME.2017) 
HM_treeMILC_2018_restricted <- get_ME_heatmap(get_ME(treeMILC_2018_restricted), HMM_ME.2018) 

# Differences in ME probability estimates for additional approach without missing covariates (Section 6.4)
HM_LC_2016_ok <- get_ME_heatmap(get_ME(LC_2016_ok), HMM_ME.2016)
HM_LC_2017_ok <- get_ME_heatmap(get_ME(LC_2017_ok), HMM_ME.2017) 
HM_LC_2018_ok <- get_ME_heatmap(get_ME(LC_2018_ok), HMM_ME.2018) 
HM_treeMILC_2016_ok <- get_ME_heatmap(get_ME(treeMILC_2016_ok), HMM_ME.2016) 
HM_treeMILC_2017_ok <- get_ME_heatmap(get_ME(treeMILC_2017_ok), HMM_ME.2017) 
HM_treeMILC_2018_ok <- get_ME_heatmap(get_ME(treeMILC_2018_ok), HMM_ME.2018)

# Differences in ME probability estimates for additional approach without HMM recodings (Section 6.4)
HM_LC_2016_recoded <- get_ME_heatmap(get_ME(LC_2016_recoded), HMM_ME.2016) 
HM_LC_2017_recoded <- get_ME_heatmap(get_ME(LC_2017_recoded), HMM_ME.2017)
HM_LC_2018_recoded <- get_ME_heatmap(get_ME(LC_2018_recoded), HMM_ME.2018) 
HM_treeMILC_2016_recoded <- get_ME_heatmap(get_ME(treeMILC_2016_recoded), HMM_ME.2016) 
HM_treeMILC_2017_recoded <- get_ME_heatmap(get_ME(treeMILC_2017_recoded), HMM_ME.2017)
HM_treeMILC_2018_recoded <- get_ME_heatmap(get_ME(treeMILC_2018_recoded), HMM_ME.2018)
