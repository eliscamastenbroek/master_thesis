############################### Initialisation_RealData.R ######################################## 
## This file reads in data from the ER and the LFS and creates the right subsets for performing ## 
## the analyses in Chapter 6.                                                                   ##
##################################################################################################

## Load packages and set working directory
library(haven)     
library(data.table) 
library(dplyr)
setwd('add_working_directory_here')

##################################################################################################
## 1. Select the correct records                                                                ##
##################################################################################################

## Read data sets for 2016, 2017 and 2018
data2016 <- read_spss('HMMdata2016_2017_voor_Frederick.sav')
data2017 <- read_spss('HMMdata2017_2018_voor_Frederick.sav')
data2018 <- read_spss('HMMdata2018_2019_voor_Frederick.sav')
setDT(data2016)
setDT(data2017)
setDT(data2018)

## Select records 1) for which both the ER and LFS are recorded, 2) in the first quarter of each year, 3) in the right age group
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
combined = rbind(data2016, data2017, data2018)

## Recode such that 1=permanent, 2=other, 3=flexible (i.e. to match the cluster names in the simulation studies)
combined[combined$contract==3,]$contract=4
combined[combined$contractEBB==3,]$contractEBB=4
combined[combined$contractEBB==2,]$contractEBB=3
combined[combined$contract==2,]$contract=3
combined[combined$contract==4,]$contract=2
combined[combined$contractEBB==4,]$contractEBB=2

##################################################################################################
## 2. Create subsets in original form (see Section 6.1)                                         ##
##################################################################################################

data2016_original = as.data.frame(combined[combined$year == 2016,])
data2017_original = as.data.frame(combined[combined$year == 2017,])
data2018_original = as.data.frame(combined[combined$year == 2018,])

##################################################################################################
## 3. Create subsets with HMM recodings (see Section 6.4)                                       ##
##################################################################################################

## Replace missing covariate categories to existing covariates categories
combined_recoded = combined

combined_recoded[combined_recoded$SBIgroep==10,]$SBIgroep = 1
combined_recoded[combined_recoded$BaanduurKlasse==7,]$BaanduurKlasse = 1
combined_recoded[combined_recoded$grootteklasse==4,]$grootteklasse = 3
combined_recoded[combined_recoded$softwarecluster==6,]$softwarecluster = 5

data2016_recoded = combined_recoded[combined_recoded$year==2016,]
data2017_recoded = combined_recoded[combined_recoded$year==2017,]
data2018_recoded = combined_recoded[combined_recoded$year==2018,]