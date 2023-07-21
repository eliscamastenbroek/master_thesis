######################## Plot_Simulation_2_Less_Optimal_Approach.R ############################### 
## This file contains the code that is required to create the plot (except for the heatmaps)    ##
## for the simulation study in Chapter 5 using the less optimal approach (see Section 5.1.1).   ##
##                                                                                              ##
## To run the code, the following object is required:                                           ##
##    - plot_df (data.frame): Data frame with results to plot the PPEs                          ##
##                                                                                              ##
## This object can be obtained by either:                                                       ##
##    1) Performing the simulation study in the file 'Perform_Simulation_2.R' and preparing the ##
##       data for plotting using the file 'Prepare_Simulation_Results_for_Plotting.R', or:      ##
##    2) Loading the RData in the file 'Simulation_2_Less_Optimal_Approach_Reduced.RData'.      ##
##################################################################################################

# Load packages
library(ggplot2)
library(dplyr)

# Define which colours to use in the plot
colours <- c(LC = "#F8766D", LCT = "#00BA38", "tree-MILC" = "#619CFF") 

# Expected value (in-text; permanent and n=1,000)
plot_df[plot_df$N==1000 & plot_df$Contract=="Permanent",] %>% ggplot(aes(x=indicator, y=prop, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="", values=colours) +
  labs(y= "Expected\nvalue", x = "Number of indicators") +
  facet_grid(C ~ ME,labeller=label_both) +
  geom_hline(aes(yintercept=hline), color="black", size=0.3) +
  geom_errorbar(aes(ymin=prop-sd, ymax=prop+sd), width=.25, size=0.3, position=position_dodge(.9)) +
  theme(legend.position="top")
