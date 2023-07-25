############################### Perform_Initial_Analysis_Approach_1.R ############################
## This file contains the code that is required to perform the initial analysis using the       ##
## first approach to include missing covariates as described in Chapter 5.                      ##
##                                                                                              ##
## To run the code:                                                                             ##
##    - Specify a working directory in line 23.                                                 ##
##    - Make sure that the files 'exampleDat_1000.txt' and 'exampleDat_10000.txt' are in        ##
##      the working directory.                                                                   ##
##    - Specify the argument 'folder' in lines 26 and 27.                                       ##
##    - Make sure that the functions in the files 'Helpfunctions_General.R',                    ##
##      'Helpfunctions_Simulations.R', 'Helpfunctions_Performance_Measures_and_Plots.R',        ##
##      'Methods_Initial_Analysis_Approach_1', and 'Simulate_Data_2.R' are loaded.              ##
##################################################################################################

# Load packages
library(dplyr)
library(data.table)

# Ignore redundant warnings from dplyr 
options(dplyr.summarise.inform = FALSE) 

# Set working directory
setwd("your_working_directory_here") 

# Perform LCT and tree-MILC analysis on a data set with 10\% ME, three indicators and two missing covariates
perform_lct(9, 3, NULL, c("baanduur","SBIgroep"), 1000, 1, folder = "your_folder_here")
perform_treeMILC(9, 3, NULL, c("baanduur","SBIgroep"), 1000, 1, folder = "your_folder_here")
