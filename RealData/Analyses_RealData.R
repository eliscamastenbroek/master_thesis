################################### Analyses_RealData.R ########################################## 
## In this file, LC and tree-MILC analysis is performed on data from the ER and the LFS (see    ##
## Chapter 6).                                                                                  ##
##################################################################################################

# Set working directory
setwd("Y:/Elisca/EchteData")
options(dplyr.summarise.inform=FALSE)

# Define covariates
cov_ok = c("geslacht", "land", "proxy", "opleiding")
cov_problem = c("Contracturen", "SBIgroep", "BaanduurKlasse", "grootteklasse", "softwarecluster")

##################################################################################################
## 1. (All) covariates with direct effects and parameter restrictions (Section 6.1-6.3)         ##
##################################################################################################

LC_2016_restricted <- perform_lc(cov_ok = c(cov_ok, "Contracturen"), cov_problem=c("SBIgroep", "BaanduurKlasse", "grootteklasse", "softwarecluster"), dat = data2016_original, folder = "Y:/Elisca/EchteData/2016/LC/CovariatesRestricted/")
LC_2017_restricted <- perform_lc(cov_ok = c(cov_ok, "Contracturen"), cov_problem=c("SBIgroep", "BaanduurKlasse", "grootteklasse", "softwarecluster"), dat = data2017_original, folder = "Y:/Elisca/EchteData/2017/LC/CovariatesRestricted/")
LC_2018_restricted <- perform_lc(cov_ok = c(cov_ok, "Contracturen"), cov_problem=c("SBIgroep", "BaanduurKlasse", "grootteklasse", "softwarecluster"), dat = data2018_original, folder = "Y:/Elisca/EchteData/2018/LC/CovariatesRestricted/")
treeMILC_2016_restricted <- perform_treeMILC(cov_ok = c(cov_ok, "Contracturen"), cov_problem=c("SBIgroep", "BaanduurKlasse", "grootteklasse", "softwarecluster"), dat = data2016_original, folder = "Y:/Elisca/EchteData/2016/treeMILC/CovariatesRestricted/")
treeMILC_2017_restricted <- perform_treeMILC(cov_ok = c(cov_ok, "Contracturen"), cov_problem=c("SBIgroep", "BaanduurKlasse", "grootteklasse", "softwarecluster"), dat = data2017_original, folder = "Y:/Elisca/EchteData/2017/treeMILC/CovariatesRestricted/")
treeMILC_2018_restricted <- perform_treeMILC(cov_ok = c(cov_ok, "Contracturen"), cov_problem=c("SBIgroep", "BaanduurKlasse", "grootteklasse", "softwarecluster"), dat = data2018_original, folder = "Y:/Elisca/EchteData/2018/treeMILC/CovariatesRestricted/")

##################################################################################################
## 2. No missing covariates (Section 6.4)                                                       ##
##################################################################################################

LC_2016_ok <- perform_lc(cov_ok = c(cov_ok), cov_problem=NULL, dat = data2016_original, folder = "Y:/Elisca/EchteData/2016/LC/CovariatesOK/")
LC_2017_ok <- perform_lc(cov_ok = c(cov_ok), cov_problem=NULL, dat = data2017_original, folder = "Y:/Elisca/EchteData/2017/LC/CovariatesOK/")
LC_2018_ok <- perform_lc(cov_ok = c(cov_ok), cov_problem=NULL, dat = data2018_original, folder = "Y:/Elisca/EchteData/2018/LC/CovariatesOK/")
treeMILC_2016_ok <- perform_treeMILC(cov_ok = c(cov_ok), cov_problem=NULL, dat = data2016_original, folder = "Y:/Elisca/EchteData/2016/treeMILC/CovariatesOK/")
treeMILC_2017_ok <- perform_treeMILC(cov_ok = c(cov_ok), cov_problem=NULL, dat = data2017_original, folder = "Y:/Elisca/EchteData/2017/treeMILC/CovariatesOK/")
treeMILC_2018_ok <- perform_treeMILC(cov_ok = c(cov_ok), cov_problem=NULL, dat = data2018_original, folder = "Y:/Elisca/EchteData/2018/treeMILC/CovariatesOK/")

##################################################################################################
## 3. Covariates with HMM recodings (Section 6.4)                                               ##
##################################################################################################

LC_2016_recoded <- perform_lc(cov_ok = c(cov_ok, cov_problem), cov_problem=NULL, dat = data2016_recoded, folder = "Y:/Elisca/EchteData/2016/LC/CovariatesRecoded/")
LC_2017_recoded <- perform_lc(cov_ok = c(cov_ok, cov_problem), cov_problem=NULL, dat = data2017_recoded, folder = "Y:/Elisca/EchteData/2017/LC/CovariatesRecoded/")
LC_2018_recoded <- perform_lc(cov_ok = c(cov_ok, cov_problem), cov_problem=NULL, dat = data2018_recoded, folder = "Y:/Elisca/EchteData/2018/LC/CovariatesRecoded/")
treeMILC_2016_recoded <- perform_treeMILC(cov_ok = c(cov_ok, cov_problem), cov_problem=NULL, dat = data2016_recoded, folder = "Y:/Elisca/EchteData/2016/treeMILC/CovariatesRecoded/")
treeMILC_2017_recoded <- perform_treeMILC(cov_ok = c(cov_ok, cov_problem), cov_problem=NULL, dat = data2017_recoded, folder = "Y:/Elisca/EchteData/2017/treeMILC/CovariatesRecoded/")
treeMILC_2018_recoded <- perform_treeMILC(cov_ok = c(cov_ok, cov_problem), cov_problem=NULL, dat = data2018_recoded, folder = "Y:/Elisca/EchteData/2018/treeMILC/CovariatesRecoded/")
