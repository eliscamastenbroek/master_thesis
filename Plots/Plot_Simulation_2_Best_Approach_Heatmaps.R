############################# Plot_Simulation_2_BestApproach.R ################################### 
## This file contains the code that is required to create the plots for the simulation study    ##
## in Chapter 5 using the best approach (i.e. with direct effects and parameter restrictions).  ##
##                                                                                              ##
## To create most plots, the following objects are required:                                    ##
##    - plot_df (data.frame): Data frame with results to plot the PPEs, the entropy R2, and     ##
##      the mean summed bias                                                                    ##
##    - plot_ME (data.frame): Data frame with results to plot the expected value of the MEPEs   ##
##    - plot_rmse_ME (data.frame): Data frame with results to plot the RMSE of the MEPEs        ##
##                                                                                              ##
## These objects can be obtained by either:                                                     ##
##    1) Performing the simulation study in the file 'Perform_Simulation2.R' and preparing the  ##
##       data for plotting using the file 'Prepare_Simulation_Results_for_Plotting.R', or:      ##
##    2) Loading the RData in the file 'Simulation2_Reduced.RData'.                             ##
##                                                                                              ##
## Note that to create the heatmaps, the following objects are required in addition:            ##
##    - LC_models (list): List containing the (detailed) model results for LC analysis          ##
##    - LCT_models (list): List containing the (detailed) model results for LCT analysis        ##
##    - treeMILC_models (list): List containing the (detailed) model results for tree-MILC      ## 
## --- toevoegen van RData ---                                                                  ##
##################################################################################################

# Load packages
library(ggplot2)
library(RColorBrewer)
library(scales)
library(dplyr)

# Define which colours to use in the plots
colours <- c(LC = "#F8766D", LCT = "#00BA38", "tree-MILC" = "#619CFF") 

# Prepare data for plotting
all_results$cov_problem <- factor(all_results$cov_problem) # Rename the values of cov_problem to Z1, Z2 etc.
all_results$cov_problem <- relevel(all_results$cov_problem, "null")
levels(all_results$cov_problem) <- c("None","Z1","Z1 and Z2")
names(plot_df)[which(names(plot_df) == "cov_problem")] <- "C" # Change column names to label the variables correctly in the plots
names(plot_ME)[which(names(plot_ME) == "cov_problem")] <- "C"
names(plot_rmse_ME)[which(names(plot_rmse_ME) == "cov_problem")] <- "C"
plot_rmse_ME$C <- factor(plot_rmse_ME$C) # Rename the values of cov_problem to Z1, Z2 etc.
plot_rmse_ME$C <- relevel(plot_rmse_ME$C, "null")
levels(plot_rmse_ME$C) <- c("None","Z1","Z1 and Z2")

##################################################################################################
## Create plots for PPEs (i.e. expected value and RMSE)                                         ##
##################################################################################################

# Adjust n to create plots for either n=1,000 or n=10,000
n = 1000

# Expected value (n=n)
plot_df[plot_df$N==n,] %>% ggplot(aes(x=indicator, y=prop, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="Model", values=colours) +
  labs(y= "Expected value", x = "Number of indicators") +
  facet_grid(Contract + C ~ ME, labeller=label_both) +
  geom_hline(aes(yintercept=hline), color="black", linewidth=0.3) +
  geom_errorbar(aes(ymin=prop-sd, ymax=prop+sd), width=.25, linewidth=0.3, position=position_dodge(.9)) +
  theme(legend.position="top")

# RMSE (n=n)
plot_df[plot_df$N==n,] %>% ggplot(aes(x=indicator, y=rmse, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="Model", values=colours) +
  labs(y= "RMSE", x = "Number of indicators") +
  facet_grid(Contract + C ~ ME, labeller=label_both) +
  theme(legend.position="top")

# Expected value (in-text; permanent and n=1,000)
plot_df[plot_df$N==1000 & plot_df$Contract=="Permanent",] %>% ggplot(aes(x=indicator, y=prop, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="", values=colours) +
  labs(y= "Expected\nvalue", x = "Number of indicators") +
  facet_grid(C ~ ME,labeller=label_both) +
  geom_hline(aes(yintercept=hline), color="black", size=0.3) +
  geom_errorbar(aes(ymin=prop-sd, ymax=prop+sd), width=.25, size=0.3, position=position_dodge(.9)) +
  theme(legend.position="top")

# RMSE (in-text; permanent and n=1,000)
plot_df[plot_df$N==1000 & plot_df$Contract=="Permanent",] %>% ggplot(aes(x=indicator, y=rmse, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="", values=colours) +
  labs(y= "RMSE", x = "Number of indicators") +
  facet_grid(C ~ ME,labeller=label_both) +
  theme(legend.position="top") +
  guides(fill = guide_legend(direction = "horizontal"))  +
  theme(legend.position="top")

##################################################################################################
## Create plots for MEPEs (i.e. expected value and RMSE)                                        ##
##################################################################################################

# Expected value (n=n)
plot_ME[plot_ME$N==n,] %>% ggplot(aes(x=indicator, y=ME1, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="Model", values=colours) +
  labs(y= "(Mean) expected value", x = "Number of indicators") +
  geom_errorbar(aes(ymin=ME1-sd_ME1, ymax=ME1+sd_ME1), width=.25, size=0.3, position=position_dodge(.9)) +
  facet_grid(Contract + C ~ ME,labeller = label_both) +
  geom_hline(aes(yintercept=mline), color="black", size=0.3) +
  theme(legend.position="top")

# RMSE (n=n)
plot_rmse_ME[plot_rmse_ME$N==n,] %>% ggplot(aes(x=indicator, y=rmse, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="Model", values=colours) +
  labs(y= "RMSE", x = "Number of indicators") +
  facet_grid(Contract + C ~ ME, labeller = label_both) +
  theme(legend.position="top")

# Expected value (in-text; permanent and n=1,000)
plot_ME[plot_ME$N==1000 & plot_ME$Contract=="Permanent",] %>% ggplot(aes(x=indicator, y=ME1, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="", values=colours) +
  labs(y= "(Mean) expected value", x = "Number of indicators") +
  geom_errorbar(aes(ymin=ME1-sd_ME1, ymax=ME1+sd_ME1), width=.25, size=0.3, position=position_dodge(.9)) +
  guides(fill = guide_legend(direction = "horizontal")) +
  facet_grid(C ~ ME,labeller=label_both) +
  geom_hline(aes(yintercept=mline), color="black", size=0.3) +
  theme(legend.position="top")

# RMSE (in-text; permanent and n=1,000)
plot_rmse_ME[plot_rmse_ME$N==1000 & plot_rmse_ME$Contract=="Permanent",] %>% ggplot(aes(x=indicator, y=rmse, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="Model", values=colours) +
  labs(y= "(Mean) RMSE", x = "Number of indicators") +
  guides(fill = guide_legend(direction = "horizontal")) +
  theme(legend.position="top") +
  facet_grid(C ~ ME,labeller=label_both)

##################################################################################################
## Plot entropy R2                                                                              ##
##################################################################################################

#Entropy R2 (all)
plot_df[plot_df$type!="tree-MILC",] %>% ggplot(aes(x=indicator, y=entropy, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="Model", values=colours) +
  labs(y= "Mean entropy R-Squared", x = "Number of indicators") +
  facet_grid(ME ~ N + C, labeller=label_both) + 
  theme(legend.position="top")

#Entropy R2 (in-text; n=1,000 and 30% ME)
plot_df[plot_df$type!="tree-MILC" & plot_df$N==1000 & plot_df$ME=="30%",] %>% ggplot(aes(x=indicator, y=entropy, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="Model", values=colours) +
  labs(y= "Mean entropy R-Squared", x = "Number of indicators") +
  facet_wrap(~ C, ncol=4, labeller=label_both) +
  theme(legend.position="top")

##################################################################################################
## Plot mean summed bias                                                                        ##
##################################################################################################

plot_df %>% ggplot(aes(x=indicator, y=diag, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="", values=colours) +
  labs(y= "Mean summed bias", x = "Number of indicators") +
  facet_grid(N + C ~ ME,labeller=label_both) +
  theme(legend.position="top")

##################################################################################################
## Create plots for heatmaps                                                                    ##
##################################################################################################

ind <- 2:3
cov_problem <- c(NULL,"baanduur", "baanduur-SBIgroep")

for (j in c(1000, 10000)) {
  for (m in cov_problem) {
    for (k in ind) {
      if (k == 2) {
        cov_ok <- "q"
      } else {
        cov_ok <- NULL
      }
      
      # Plot bias and variance for 10%, 20%, and 30% ME
      bias1 <- get_ME_heatmap(LC_models, LCT_models, treeMILC_models, LC_results, 1, k, cov_ok, m, j)
      variance1 <- get_ME_heatmap(LC_models, LCT_models, treeMILC_models, LC_results, 2, k, cov_ok, m, j)
      
      # Plot bias and variance for realistic 7% ME
      bias2 <- get_ME_heatmap_realistic(LC_models, LCT_models, treeMILC_models, LC_results, 1, k, cov_ok, m, j)
      variance2 <- get_ME_heatmap_realistic(LC_models, LCT_models, treeMILC_models, LC_results, 1, k, cov_ok, m, j)
      
    }  
  }
}
