############################# Plot_Simulation_2_Best_Approach.R ################################## 
## This file contains the code that is required to create the plots (except for the heatmaps)   ##
## for the simulation study in Chapter 5 using the best approach (i.e. with direct effects and  ##
## parameter restrictions).                                                                     ##
##                                                                                              ##
## To run the code, the following objects are required:                                         ##
##    - plot_df (data.frame): Data frame with results to plot the PPEs, the entropy R2, and     ##
##      the mean summed bias                                                                    ##
##    - plot_ME (data.frame): Data frame with results to plot the expected value of the MEPEs   ##
##    - plot_rmse_ME (data.frame): Data frame with results to plot the RMSE of the MEPEs        ##
##                                                                                              ##
## These objects can be obtained by either:                                                     ##
##    1) Performing the simulation study in the file 'Perform_Simulation2.R' and preparing the  ##
##       data for plotting using the file 'Prepare_Simulation_Results_for_Plotting.R', or:      ##
##    2) Loading the RData in the file 'Simulation2_Reduced.RData'.                             ##
##################################################################################################

# Load packages
library(ggplot2)
library(dplyr)

# Define which colours to use in the plots
colours <- c(LC = "#F8766D", LCT = "#00BA38", "tree-MILC" = "#619CFF") 

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
