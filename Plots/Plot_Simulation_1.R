#################################### Plot_Simulation_1.R ######################################### 
## This file contains the code that is required to create the plots for the simulation study    ##
## in Chapter 4. To create most plots, the following objects are required to be in the global   ##
## environment:                                                                                 ##
##    - plot_df (data.frame): Data frame with results to plot the PPEs, the entropy R2, and     ##
##      the mean summed bias                                                                    ##
##    - plot_ME (data.frame): Data frame with results to plot the expected value of the MEPEs   ##
##    - plot_rmse_ME (data.frame): Data frame with results to plot the RMSE of the MEPEs        ##
## These objects can be created by:                                                             ##
##    1) Performing the simulation study in the file 'Perform_Simulation1.R' and preparing the  ##
##       data for plotting using the file 'Prepare_Simulation_Results_for_Plotting.R', or:      ##
##    2) Loading the RData in the file 'Simulation1_Reduced.RData'.                             ##
## Note that to create the heatmaps, the following objects are also required:                   ##
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
colours <- c(LC = "#F8766D", LCT = "#00BA38", "tree-MILC" ="#619CFF") 

##################################################################################################
## Create plots for PPEs (i.e. expected value and RMSE)                                         ##
##################################################################################################

# Expected value (all)
plot_df %>% ggplot(aes(x=indicator, y=prop, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="Model", values=colours) +
  labs(y= "Expected value", x = "Number of indicators") +
  facet_grid(Contract + N ~ ME, labeller=label_both) +
  geom_hline(aes(yintercept=hline), color="black", linewidth=0.3) +
  geom_errorbar(aes(ymin=prop-sd, ymax=prop+sd), width=.25, size=0.3, position=position_dodge(.9)) +
  theme(legend.position="top")
  
# RMSE (all)
plot_df %>% ggplot(aes(x=indicator, y=rmse, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="Model", values=colours) +
  labs(y= "RMSE", x = "Number of indicators") +
  facet_grid(Contract + N ~ ME, labeller=labeller(.rows=label_value, .cols=label_both)) +
  theme(legend.position="top")
  
# Expected value (permanent and n=1,000)
plot_df[plot_df$N==1000 & plot_df$Contract=="Permanent",] %>% ggplot(aes(x=indicator, y=prop, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="", values=colours) +
  labs(y= "Expected\nvalue", x = "Number of indicators") +
  facet_wrap(~ME, ncol=4, labeller=label_value) +
  geom_hline(aes(yintercept=hline), color="black", size=0.3) +
  geom_errorbar(aes(ymin=prop-sd, ymax=prop+sd), width=.25, size=0.3, position=position_dodge(.9)) +
  theme(legend.position="top")
  
# RMSE (permanent and n=1,000)
plot_df[plot_df$N==1000 & plot_df$Contract=="Permanent",] %>% ggplot(aes(x=indicator, y=rmse, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="", values=colours) +
  labs(y= "RMSE", x = "Number of indicators") +
  facet_wrap(~ME, ncol=4, labeller=label_value) +
  theme(legend.position="top") +
  guides(fill = guide_legend(direction = "horizontal"))  +
  theme(legend.position="top")

##################################################################################################
## Create plots for MEPEs (i.e. expected value and RMSE)                                        ##
##################################################################################################

# Expected value (all)
plot_ME %>% ggplot(aes(x=indicator, y=ME1, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="Model", values=colours) +
  labs(y= "(Mean) expected value", x = "Number of indicators") +
  geom_errorbar(aes(ymin=ME1-sd_ME1, ymax=ME1+sd_ME1), width=.25, size=0.3, position=position_dodge(.9)) +
  facet_grid(Contract + N ~ ME, labeller=label_both) +
  geom_hline(aes(yintercept=mline), color="black", size=0.3) +
  theme(legend.position="top")

# RMSE (all)
plot_rmse_ME %>% ggplot(aes(x=indicator, y=rmse, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="Model", values=colours) +
  labs(y= "RMSE", x = "Number of indicators") +
  facet_grid(Contract + N ~ ME, labeller=labeller(.rows=label_both, .cols=label_both)) +
  theme(legend.position="top")

# Expected value (permanent and n=1,000)
plot_ME[plot_ME$N==1000 & plot_ME$Contract=="Permanent",] %>% ggplot(aes(x=indicator, y=ME1, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="", values=colours) +
  labs(y= "(Mean) expected value", x = "Number of indicators") +
  geom_errorbar(aes(ymin=ME1-sd_ME1, ymax=ME1+sd_ME1), width=.25, size=0.3, position=position_dodge(.9)) +
  guides(fill = guide_legend(direction = "horizontal")) +
  facet_wrap(~ME, ncol=5, labeller=label_value) +
  geom_hline(aes(yintercept=mline), color="black", size=0.3) +
  theme(legend.position="top")

# RMSE (permanent and n=1,000)
plot_rmse_ME[plot_rmse_ME$N==1000 & plot_rmse_ME$Contract=="Permanent",] %>% ggplot(aes(x=indicator, y=rmse, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="Model", values=colours) +
  labs(y= "(Mean) RMSE", x = "Number of indicators") +
  guides(fill = guide_legend(direction = "horizontal")) +
  theme(legend.position="top") +
  facet_wrap(~ ME, ncol=5, labeller=label_value)

##################################################################################################
## Plot entropy R2                                                                              ##
##################################################################################################

#Entropy R2 (all)
plot_df[plot_df$type!="tree-MILC",] %>% ggplot(aes(x=indicator, y=entropy, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="Model", values=colours) +
  labs(y= "Mean entropy R-Squared", x = "Number of indicators") +
  facet_grid(rows = vars(N), cols=vars(ME), labeller=label_both) +
  theme(legend.position="top")

#Entropy R2 (n=1,000)
plot_df[plot_df$type!="tree-MILC" & plot_df$N==1000,] %>% ggplot(aes(x=indicator, y=entropy, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="Model", values=colours) +
  labs(y= "Mean entropy R-Squared", x = "Number of indicators") +
  facet_wrap(~ME, ncol=4, labeller=label_both) +
  theme(legend.position="top")

##################################################################################################
## Plot mean summed bias                                                                        ##
##################################################################################################

plot_df %>% ggplot(aes(x=indicator, y=diag, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="", values=colours) +
  labs(y= "Mean summed bias", x = "Number of indicators") +
  facet_wrap(~ ME, ncol=5, labeller=label_value) +
  theme(legend.position="top")

##################################################################################################
## Create plots for heatmaps                                                                    ##
##################################################################################################

# Simulation study 1
ind <- 2:4
cov_problem <- "null"

# Simulation study 2
ind <- 2:3
cov_problem <- c("baanduur", "SBIgroep", "")

for (j in c(1000, 10000)) {
  for (m in cov_problem) {
    for (k in ind) {
      if (k == 2) {
        cov_ok <- "q"
      } else {
        cov_ok <- NULL
      }
      # Plot bias and variance for 10%, 20%, and 30% ME
      bias1 <- get_ME_heatmap(LC_models, LCT_models, treeMILC_models, LC_results, 1, k, cov_ok, cov_problem, j)
      variance1 <- get_ME_heatmap(LC_models, LCT_models, treeMILC_models, LC_results, 2, k, cov_ok, cov_problem, j)
      
      # Plot bias and variance for realistic 7% ME
      bias2 <- get_ME_heatmap_realistic(LC_models, LCT_models, treeMILC_models, LC_results, 1, k, cov_ok, cov_problem, j)
      variance2 <- get_ME_heatmap_realistic(LC_models, LCT_models, treeMILC_models, LC_results, 1, k, cov_ok, cov_problem, j)

    }  
  }
}