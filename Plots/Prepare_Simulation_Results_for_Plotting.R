######################## Prepare_Simulation_Results_for_Plotting.R ############################### 
## This file contains the code that is required to prepare the results of both simulation       ##
## studies (and the initial analyses) in Chapters 4-5 for plotting.                             ##
##                                                                                              ##
## Note that the code in this file requires that the code in the file                           ##
## 'Perform_Simulation_Study_1.R' or 'Perform_Simulation_Study_2.R' has been executed.          ##
##                                                                                              ##
## In other words, the code assumes that the following objects are in the global environment:   ##
##    - all_results (data.frame): Data frame containing the results of the simulation study     ##
##    - true_proportions (numeric): Vector containing the true proportions                      ##
##    - ME_matrix1 (matrix): Matrix with the true ME probabilities for 10% ME                   ##
##    - ME_matrix2 (matrix): Matrix with the true ME probabilities for 20% ME                   ##
##    - ME_matrix3 (matrix): Matrix with the true ME probabilities for 30% ME                   ##
##    - ME_matrix4 (matrix): Matrix with the true ME probabilities for a realistic 7% ME        ##
##    - LC_models (list): List containing the (detailed) model results for LC analysis          ##
##    - LCT_models (list): List containing the (detailed) model results for LCT analysis        ##
##    - treeMILC_models (list): List containing the (detailed) model results for tree-MILC      ## 
##    - LC_results (data.frame): Data frame containing the model results for LC analysis        ##
##    - LCT_results (data.frame): Data frame containing the model results for LC analysis       ##
##    - treeMILC_results (data.frame): Data frame containing the model results for treeMILC     ##
##    - LC_summary (data.frame): Data frame containing the model results for LC analysis        ##
##    - LCT_summary (data.frame): Data frame containing the model results for LC analysis       ##
##    - treeMILC_summary (data.frame): Data frame containing the model results for treeMILC     ##
##    - get_rmse_ME (function): Function that computes the RMSE of the ME probabilities         ##
##################################################################################################

# Load packages
library(ggplot2)
library(RColorBrewer)
library(scales)
library(dplyr)

##################################################################################################
## Prepare data for plotting in general                                                         ##
##################################################################################################

# Rename the values of cov_problem to Z1, Z2 etc.
all_results$cov_problem <- factor(all_results$cov_problem) 
all_results$cov_problem <- relevel(all_results$cov_problem, "null")
levels(all_results$cov_problem) <- c("None","Z1","Z1 and Z2")

# Change column names to label the variables correctly in the plots
names(all_results)[which(names(all_results) == "cov_problem")] <- "C" 

##################################################################################################
## Prepare data to plot population proportion estimates (PPEs) (i.e. expected value and RMSE)   ##                                                                    ##
##################################################################################################

# Create data frames for each contract type
perm <- all_results[,c(names(all_results)[1:6],"prop1","rmse_prop1","sd_prop1","type","ME1.1","ME2.1","sd_ME1.1","sd_ME2.1","diag","entropy")]
other <- all_results[,c(names(all_results)[1:6],"prop2","rmse_prop2","sd_prop2","type","ME1.2","ME2.2","sd_ME1.2","sd_ME2.2","diag","entropy")]
flex <- all_results[,c(names(all_results)[1:6],"prop3","rmse_prop3","sd_prop3","type","ME1.3","ME2.3","sd_ME1.3","sd_ME2.3","diag","entropy")]
names(perm) <- names(flex) <- names(other) <- c(names(all_results)[1:6],"prop","rmse","sd","type","ME1","ME2","sd_ME1","sd_ME2","diag","entropy")
perm$Contract <- "Permanent"
other$Contract <- "Other"
flex$Contract <- "Flexible"

# Add true proportions
perm$hline <- true_proportions[1]
other$hline <- true_proportions[2]
flex$hline <- true_proportions[3]

# Combine data frames 
plot_df <- rbind(perm, flex, other)
rm(perm, flex, other)
plot_df$Contract <- as.factor(plot_df$Contract)
plot_df$Contract <- factor(plot_df$Contract, levels = c("Permanent", "Flexible", "Other"))

# Get variables into the right format
plot_df$ME <- as.factor(plot_df$ME)
plot_df$ME <- relevel(plot_df$ME, 4)
levels(plot_df$ME) <- c("Realistic 7%", "10%", "20%", "30%")

# Add true ME probabilities (10%, 20%, and 30%)
plot_df$mline <- NA
plot_df[plot_df$Contract == "Permanent" & plot_df$ME == "10%", ]$mline <- ME_matrix1[1, 1] 
plot_df[plot_df$Contract == "Other" & plot_df$ME == "10%", ]$mline <- ME_matrix1[2, 2] 
plot_df[plot_df$Contract == "Flexible" & plot_df$ME == "10%", ]$mline <- ME_matrix1[3, 3] 
plot_df[plot_df$Contract == "Permanent" & plot_df$ME == "20%", ]$mline <- ME_matrix2[1, 1] 
plot_df[plot_df$Contract == "Other" & plot_df$ME == "20%", ]$mline <- ME_matrix2[2, 2] 
plot_df[plot_df$Contract == "Flexible" & plot_df$ME == "20%", ]$mline <- ME_matrix2[3, 3] 
plot_df[plot_df$Contract == "Permanent" & plot_df$ME == "30%", ]$mline <- ME_matrix3[1, 1] 
plot_df[plot_df$Contract == "Other" & plot_df$ME == "30%", ]$mline <- ME_matrix3[2, 2] 
plot_df[plot_df$Contract == "Flexible" & plot_df$ME == "30%", ]$mline <- ME_matrix3[3, 3]

##################################################################################################
## Prepare data to plot expected value of measurement error probability estimates (MEPEs)       ##
##################################################################################################

# Create data frame that contains the expected value of MEPEs for Y1 and Y3 in conditions with a realistic 7% ME
plot_rmse_ME2 <- plot_df[!is.na(plot_df$ME2), c("indicator","cov_ok","C","N","ME","id","type","ME2","sd_ME2","Contract","mline")]
levels(plot_rmse_ME2$ME)[1] <- "Real. 7% (ind. 2)"
plot_rmse_ME2[plot_rmse_ME2$Contract == "Permanent" & plot_rmse_ME2$ME == "Real. 7% (ind. 2)", ]$mline <- ME_matrix4b[1, 1] 
plot_rmse_ME2[plot_rmse_ME2$Contract == "Other" & plot_rmse_ME2$ME == "Real. 7% (ind. 2)", ]$mline <- ME_matrix4b[2, 2] 
plot_rmse_ME2[plot_rmse_ME2$Contract == "Flexible" & plot_rmse_ME2$ME == "Real. 7% (ind. 2)", ]$mline <- ME_matrix4b[3, 3] 
names(plot_rmse_ME2)[8:9] <- c("ME1","sd_ME1")

# Create data frame that contains the expected value of MEPEs for Y2 and Y4 in conditions with a realistic 7% ME
plot_rmse_ME1 <- plot_df[, c("indicator","cov_ok","C","N","ME","id","type","ME1","sd_ME1","Contract","mline")]
levels(plot_rmse_ME1$ME)[1] <- "Real. 7% (ind. 1)"
plot_rmse_ME1[plot_rmse_ME1$Contract == "Permanent" & plot_rmse_ME1$ME == "Real. 7% (ind. 1)", ]$mline <- ME_matrix4a[1, 1] 
plot_rmse_ME1[plot_rmse_ME1$Contract == "Other" & plot_rmse_ME1$ME == "Real. 7% (ind. 1)", ]$mline <- ME_matrix4a[2, 2] 
plot_rmse_ME1[plot_rmse_ME1$Contract == "Flexible" & plot_rmse_ME1$ME == "Real. 7% (ind. 1)", ]$mline <- ME_matrix4a[3, 3] 

# Combine data frames
plot_ME <- rbind(plot_rmse_ME2, plot_rmse_ME1)
rm(plot_rmse_ME1, plot_rmse_ME2)
plot_ME$ME <- as.factor(plot_ME$ME)
plot_ME$ME <- relevel(plot_ME$ME, "Real. 7% (ind. 2)")
plot_ME$ME <- relevel(plot_ME$ME, "Real. 7% (ind. 1)")

##################################################################################################
## Prepare data to plot RMSE of measurement error probability estimators (MEPEs)                ##
##################################################################################################

# Compute RMSE of ME probability estimators
LC_ME_rmse <- get_rmse_ME(LC_models, LC_results) # Note that this takes approx. 5 minutes to run
LCT_ME_rmse <- get_rmse_ME(LCT_models, LCT_results) # Note that this takes approx. 5 minutes to run
treeMILC_ME_rmse <- get_rmse_ME(treeMILC_models, treeMILC_results) # Note that this takes approx. 5 minutes to run
LC_ME_rmse$rmse <- as.numeric(LC_ME_rmse$rmse)
LCT_ME_rmse$rmse <- as.numeric(LCT_ME_rmse$rmse)
treeMILC_ME_rmse$rmse <- as.numeric(treeMILC_ME_rmse$rmse)
LC_ME_rmse$type <- "LC"
LCT_ME_rmse$type <- "LCT"
treeMILC_ME_rmse$type <- "tree-MILC"

# Combine results
plot_rmse_ME <- rbind(LC_ME_rmse, LCT_ME_rmse, treeMILC_ME_rmse)
plot_rmse_ME$ME <- as.numeric(plot_rmse_ME$ME)
plot_rmse_ME$ind <- as.numeric(plot_rmse_ME$ind)

# Rename the values of cov_problem to Z1, Z2 etc.
plot_rmse_ME$cov_problem <- factor(plot_rmse_ME$cov_problem) 
plot_rmse_ME$cov_problem <- relevel(plot_rmse_ME$cov_problem, "null")
levels(plot_rmse_ME$cov_problem) <- c("None","Z1","Z1 and Z2")

# Change column names to label the variables correctly in the plots
names(plot_rmse_ME)[which(names(plot_rmse_ME) == "cov_problem")] <- "C" 

# Combine indicators Y1/Y3 and Y2/Y4
plot_rmse_ME[plot_rmse_ME$ME == 4 & (plot_rmse_ME$ind == 1 | plot_rmse_ME$ind == 3), ]$ind <- 1
plot_rmse_ME[plot_rmse_ME$ME == 4 & (plot_rmse_ME$ind == 2 | plot_rmse_ME$ind == 4), ]$ind <- 2

# Compute mean RMSE over all indicators
plot_rmse_ME_not4 <- plot_rmse_ME[plot_rmse_ME$ME != 4, ] 
plot_rmse_ME_not4 <- as.data.frame(plot_rmse_ME_not4 %>% group_by(indicator, cov_ok, C, N, ME, id, Contract, type) %>% summarise_at(c("rmse"), mean))
plot_rmse_ME_not4$ind <- NA
plot_rmse_ME_4 <- plot_rmse_ME[plot_rmse_ME$ME == 4, ] 
plot_rmse_ME_4 <- as.data.frame(plot_rmse_ME_4 %>% group_by(indicator, cov_ok, C, N, ME, id, ind, Contract, type) %>% summarise_at(c("rmse"), mean))
plot_rmse_ME <- rbind(plot_rmse_ME_4, plot_rmse_ME_not4)
rm(plot_rmse_ME_4, plot_rmse_ME_not4)

# Get variables into the right format
plot_rmse_ME$Contract <- as.factor(plot_rmse_ME$Contract)
plot_rmse_ME$ME <- as.character(plot_rmse_ME$ME)
plot_rmse_ME[plot_rmse_ME$ME == 4 & plot_rmse_ME$ind == 1, ]$ME <- "Real. 7% (ind. 1)"
plot_rmse_ME[plot_rmse_ME$ME == 4 & plot_rmse_ME$ind == 2, ]$ME <- "Real. 7% (ind. 2)"
plot_rmse_ME$ME <- as.factor(plot_rmse_ME$ME)
plot_rmse_ME$ME <- relevel(plot_rmse_ME$ME, "Real. 7% (ind. 2)")
plot_rmse_ME$ME <- relevel(plot_rmse_ME$ME, "Real. 7% (ind. 1)")
levels(plot_rmse_ME$ME)[3:5] <- c("10%", "20%", "30%")
levels(plot_rmse_ME$Contract) <- c("Permanent", "Other", "Flexible")
plot_rmse_ME$Contract <- factor(plot_rmse_ME$Contract, levels = c("Permanent", "Flexible", "Other"))
