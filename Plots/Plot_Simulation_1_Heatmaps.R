################################### Plot_Simulation_1_Heatmaps.R ################################# 
## This file contains the code that is required to create the heatmaps for the simulation study ##
## in Chapter 4.                                                                                ##
##                                                                                              ##
## Note that to create the heatmaps, the functions in the files 'Helpfunctions_General.R',      ##
## 'Helpfunctions_Simulations.R' and 'Helpfunctions_Performance_Measures_and_Plots.R' are       ## 
## required. In addition, the following objects are required:                                   ##
##    - plot_df (data.frame): Data frame with results to plot the PPEs, the entropy R2, and     ##
##      the mean summed bias                                                                    ##
##    - plot_ME (data.frame): Data frame with results to plot the expected value of the MEPEs   ##
##    - plot_rmse_ME (data.frame): Data frame with results to plot the RMSE of the MEPEs        ##
##    - LC_models (list): List containing the (detailed) model results for LC analysis          ##
##    - LCT_models (list): List containing the (detailed) model results for LCT analysis        ##
##    - treeMILC_models (list): List containing the (detailed) model results for tree-MILC      ## 
##                                                                                              ##
## These objects can be obtained by either:                                                     ##
##    1) Performing the simulation study in the file 'Perform_Simulation_1.R' and preparing the ##
##       data for plotting using the file 'Prepare_Simulation_Results_for_Plotting.R', or:      ##
##    2) Loading the RData in the file 'Simulation_1_Complete.RData'. Note that due to          ##
##       storage limits in GitHub, this file is only available in the folder 'Stage_Elisca'     ##
##       at CBS.                                                                                ##
##                                                                                              ##
## Finally, to create the heatmaps, a working directory should be specified in line 33. The     ##
## plots will be saved in this directory.                                                       ##
##################################################################################################

# Load packages
library(ggplot2)
library(RColorBrewer)
library(scales)
library(dplyr)

setwd("your_working_directory_here")

ind <- 2:4

for (j in c(1000, 10000)) {
    for (k in ind) {
      if (k == 2) {
        cov_ok <- "q"
      } else {
        cov_ok <- NULL
      }

      # Model name
      name = paste(k, cov_ok, j, sep="-") 
      
      # Plot bias for models with 10%, 20%, and 30% ME
      heatmap <- get_ME_heatmap(LC_models, LCT_models, treeMILC_models, LC_results, 1, k, cov_ok, NULL, j)
      filename <- paste0(name, "_Bias_10-20-30.pdf")
      ggsave(plot = heatmap, width = 3.9, height = 3, dpi = 300, filename = filename)

      # Plot variance for models with 10%, 20%, and 30% ME
      heatmap <- get_ME_heatmap(LC_models, LCT_models, treeMILC_models, LC_results, 2, k, cov_ok, NULL, j)
      filename <- paste0(name, "_Variance_10-20-30.pdf")
      ggsave(plot = heatmap, width = 3.9, height = 3, dpi = 300, filename = filename)
    
      # Plot bias for models with a realistic 7% ME
      heatmap <- get_ME_heatmap_realistic(LC_models, LCT_models, treeMILC_models, LC_results, 1, k, cov_ok, NULL, j)
      filename <- paste0(name, "_Bias_7.pdf")
      ggsave(plot = heatmap, width = 3.9, height = 3, dpi = 300, filename = filename)

      # Plot variance for models with a realistic 7% ME
      heatmap <- get_ME_heatmap_realistic(LC_models, LCT_models, treeMILC_models, LC_results, 1, k, cov_ok, NULL, j)
      filename <- paste0(name, "_Variance_7.pdf")
      ggsave(plot = heatmap, width = 3.9, height = 3, dpi = 300, filename = filename)

    }  
}
