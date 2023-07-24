############################# Perform_Real_Data_Analyses.R ####################################### 
## This file reads in data from the ER and the LFS, creates the right subsets, and performs     ##  
## LC and tree-MILC analysis on the data (see Chapter 6).                                       ## 
##                                                                                              ##
## To run the code:                                                                             ##
##    - Set a working directory in line 19.                                                     ##
##    - Specify the 'folder ' arguments in lines 82-103 (i.e, a path to where the               ##
##      model results should be stored).                                                        ##
##    - Make sure that the functions in the files 'Helpfunctions_General.R',                    ##
##      'Helpfunctions_Real_Data.R', 'Helpfunctions_Performance_Measures_and_Plots.R', and      ##
##      'Methods_Real_Data.R' should be loaded.                                                 ##
##################################################################################################

## Load packages and set working directory
library(haven)     
library(data.table) 

## Set working directory and ignore dplyr warnings
setwd("add_working_directory_here")
options(dplyr.summarise.inform=FALSE)

##################################################################################################
## Data initialisation                                                                           ## 
##################################################################################################

## Read in data sets for 2016, 2017 and 2018
data2016 <- read_spss('HMMdata2016_2017_voor_Frederick.sav')
data2017 <- read_spss('HMMdata2017_2018_voor_Frederick.sav')
data2018 <- read_spss('HMMdata2018_2019_voor_Frederick.sav')
setDT(data2016)
setDT(data2017)
setDT(data2018)

## Select records for which both the ER and LFS are recorded in the first quarter of each year and in the right age group
data2016 <- data2016[ !is.na(data2016$contract) & !is.na(data2016$contractEBB) & (data2016$LftKlasse==3 | data2016$LftKlasse==4) & data2016$cohort==1, ]
data2017 <- data2017[ !is.na(data2017$contract) & !is.na(data2017$contractEBB) & (data2017$LftKlasse==3 | data2017$LftKlasse==4) & data2017$cohort==1, ]
data2018 <- data2018[ !is.na(data2018$contract) & !is.na(data2018$contractEBB) & (data2018$LftKlasse==3 | data2018$LftKlasse==4) & data2018$cohort==1, ]
data2016$year = 2016
data2017$year = 2017
data2018$year = 2018

## Select first observations per respondent
data2016 <- data2016[, .SD[1], by = persnr]
data2017 <- data2017[, .SD[1], by = persnr]
data2018 <- data2018[, .SD[1], by = persnr]

## Combine subsets
combined <- rbind(data2016, data2017, data2018)

## Recode to 1 = permanent, 2 = other, 3 = flexible to match the cluster names in the simulation studies
combined$contract[combined$contract == 3] <- 4
combined$contractEBB[combined$contractEBB == 3] <- 4
combined$contractEBB[combined$contractEBB == 2] <- 3
combined$contract[combined$contract == 2] <- 3
combined$contract[combined$contract == 4] <- 2
combined$contractEBB[combined$contractEBB == 4] <- 2

## Create subsets for the models in Section 6.1-6.3
data2016_original <- as.data.frame(combined[combined$year == 2016, ])
data2017_original <- as.data.frame(combined[combined$year == 2017, ])
data2018_original <- as.data.frame(combined[combined$year == 2018, ])

## Create subsets for the models in Section 6.4 with HMM recodings
combined_recoded <- combined
combined_recoded$SBIgroep[combined_recoded$SBIgroep == 10] <- 1 
combined_recoded$BaanduurKlasse[combined_recoded$BaanduurKlasse == 7] <- 1
combined_recoded$grootteklasse[combined_recoded$grootteklasse == 4] <- 3
combined_recoded$softwarecluster[combined_recoded$softwarecluster == 6] <- 5
data2016_recoded <- combined_recoded[combined_recoded$year == 2016, ]
data2017_recoded <- combined_recoded[combined_recoded$year == 2017, ]
data2018_recoded <- combined_recoded[combined_recoded$year == 2018, ]

##################################################################################################
## After initialising the data, LC and tree-MILC analysis can be performed.                     ## 
##################################################################################################

## Define covariates
cov_ok = c("geslacht", "land", "proxy", "opleiding")
cov_problem = c("Contracturen", "SBIgroep", "BaanduurKlasse", "grootteklasse", "softwarecluster")

## Perform LC and tree-MILC analyiss with missing covariates using direct effects and parameter restrictions (Section 6.1-6.3)
LC_2016_restricted <- perform_lc(cov_ok = c(cov_ok, "Contracturen"), cov_problem=c("SBIgroep", "BaanduurKlasse", "grootteklasse", "softwarecluster"), dat = data2016_original, folder = "add_folder_here")
LC_2017_restricted <- perform_lc(cov_ok = c(cov_ok, "Contracturen"), cov_problem=c("SBIgroep", "BaanduurKlasse", "grootteklasse", "softwarecluster"), dat = data2017_original, folder = "add_folder_here")
LC_2018_restricted <- perform_lc(cov_ok = c(cov_ok, "Contracturen"), cov_problem=c("SBIgroep", "BaanduurKlasse", "grootteklasse", "softwarecluster"), dat = data2018_original, folder = "add_folder_here")
treeMILC_2016_restricted <- perform_treeMILC(cov_ok = c(cov_ok, "Contracturen"), cov_problem=c("SBIgroep", "BaanduurKlasse", "grootteklasse", "softwarecluster"), dat = data2016_original, folder = "add_folder_here")
treeMILC_2017_restricted <- perform_treeMILC(cov_ok = c(cov_ok, "Contracturen"), cov_problem=c("SBIgroep", "BaanduurKlasse", "grootteklasse", "softwarecluster"), dat = data2017_original, folder = "add_folder_here")
treeMILC_2018_restricted <- perform_treeMILC(cov_ok = c(cov_ok, "Contracturen"), cov_problem=c("SBIgroep", "BaanduurKlasse", "grootteklasse", "softwarecluster"), dat = data2018_original, folder = "add_folder_here")

## Perform LC and tree-MILC analysis without missing covariates (Section 6.4)
LC_2016_ok <- perform_lc(cov_ok = c(cov_ok), cov_problem=NULL, dat = data2016_original, folder = "add_folder_here")
LC_2017_ok <- perform_lc(cov_ok = c(cov_ok), cov_problem=NULL, dat = data2017_original, folder = "add_folder_here")
LC_2018_ok <- perform_lc(cov_ok = c(cov_ok), cov_problem=NULL, dat = data2018_original, folder = "add_folder_here")
treeMILC_2016_ok <- perform_treeMILC(cov_ok = c(cov_ok), cov_problem=NULL, dat = data2016_original, folder = "add_folder_here")
treeMILC_2017_ok <- perform_treeMILC(cov_ok = c(cov_ok), cov_problem=NULL, dat = data2017_original, folder = "add_folder_here")
treeMILC_2018_ok <- perform_treeMILC(cov_ok = c(cov_ok), cov_problem=NULL, dat = data2018_original, folder = "add_folder_here")

## Perform LC and tree-MILC analysis with missing covariates using HMM recodings (Section 6.4)               
LC_2016_recoded <- perform_lc(cov_ok = c(cov_ok, cov_problem), cov_problem=NULL, dat = data2016_recoded, folder = "add_folder_here")
LC_2017_recoded <- perform_lc(cov_ok = c(cov_ok, cov_problem), cov_problem=NULL, dat = data2017_recoded, folder = "add_folder_here")
LC_2018_recoded <- perform_lc(cov_ok = c(cov_ok, cov_problem), cov_problem=NULL, dat = data2018_recoded, folder = "add_folder_here")
treeMILC_2016_recoded <- perform_treeMILC(cov_ok = c(cov_ok, cov_problem), cov_problem=NULL, dat = data2016_recoded, folder = "add_folder_here")
treeMILC_2017_recoded <- perform_treeMILC(cov_ok = c(cov_ok, cov_problem), cov_problem=NULL, dat = data2017_recoded, folder = "add_folder_here")
treeMILC_2018_recoded <- perform_treeMILC(cov_ok = c(cov_ok, cov_problem), cov_problem=NULL, dat = data2018_recoded, folder = "add_folder_here")
