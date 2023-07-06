################################### Analyses_RealData.R ########################################## 
## After initialising the data using 'Initialisation_RealData.R', this file can be used to      ##
## perform LC and tree-MILC analysis on data from the ER and the LFS (see Chapter 6).           ##
##################################################################################################

## Set working directory and ignore dplyr warnings
setwd("add_working_directory_here")
options(dplyr.summarise.inform=FALSE)

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
