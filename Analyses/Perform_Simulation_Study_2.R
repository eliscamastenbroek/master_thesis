#################################### Perform_Simulation_Study_2.R ################################
## This file contains the code that is required to perform the simulation study with missing    ##
## covariates as described in Chapter 5. In the first part, the simulation study is performed.  ##
## In the second part, the results of the simulation study are obtained.                        ##
##                                                                                              ##
## Note that depending on whether the functions from the file                                   ##
## 'Methods_Initial_Analysis_Approach_2.R' or the file 'Methods_Simulation_Studies.R' are       ##
## loaded, missing covariates are included using either the second approach (see Section 5.1.2) ##
## or the third (and best) approach (see Section 5.1.3).                                        ##
##                                                                                              ##
## To run the code:                                                                             ##
##    - In line 32, a working directory should be set.                                          ##
##    - In lines 236, 240, 244, the argument 'folder' should be specified.                      ##
##    - The functions in the files 'Helpfunctions_General.R', 'Helpfunctions_Simulations.R',    ##
##      'Helpfunctions_Performance_Measures_and_Plots.R', and 'Simulate_Data_2.R' should be     ##
##      loaded.                                                                                 ##
##    - The functions in the file 'Methods_LessOptimalApproach.R' OR 'Methods_BestApproach.R'   ##
##      should be loaded.                                                                       ##
##    - The files 'exampleDat_1000.dat' and 'exampleDat_10000.dat' are required to be in the    ##
##      working directory.                                                                      ##
##################################################################################################

# Initialisations
library(dplyr)
library(data.table)
library(stringr)

# Ignore redundant warnings from dplyr
options(dplyr.summarise.inform = FALSE)  

# Set working directory
setwd("your_working_directory_here") 

#####################################################################################################
## 2. Perform simulation study 2                                                                   ##
#####################################################################################################

## Specification of simulation conditions
ind <- c(2, 3)
N <- c(1000,10000)
ME <- c(1:4)
iteration <- c(1:50)
cov_problem <- c("NULL", "baanduur", "baanduur-SBIgroep")

## Create lists to store results in
LC_models <- list()
LCT_models <- list()
treeMILC_models <- list()

## Keep track of number of models
m <- 0 

## Keep track of potential errors
errors <- 0
error_vec <- vector()

## Execute simulations
for(l in iteration){
  for(k in ME){
    for(i in ind){
      if(i == 2){
        cov_ok <- "q"
      } else {
        cov_ok <- NULL
      }
      
      for(j in N){
        for(n in cov_problem){
          m <- m + 1
          
          ## Create model name
          name <- paste0(i, "-", ifelse(is.null(cov_ok), "NULL", cov_ok), "-", n, "-", j, "-", k)
          print(paste0("---------------", m, "/", length(iteration) * length(cov_problem) * length(ind) * length(N) * length(ME) , "------------------"))
          
          ## Get parameters in right format
          if(n == "NULL"){
            n <- NULL
          } else if(n == "baanduur-SBIgroep"){
            n <- c("baanduur", "SBIgroep")
          }
          
          ## Perform LC, LCT and tree-MILC
          LC <- perform_lc(l, i, cov_ok, n, j, k, folder="your_folder_here")
          LC_models <- append(LC_models, list(LC))
          print(paste0("LC model ", name, " complete."))
          
          LCT <- perform_lct(l, i, cov_ok, n, j, k, folder="your_folder_here")
          LCT_models <- append(LCT_models, list(LCT))
          print(paste0("LCT model ", name, " complete."))
          
          treeMILC <- perform_treeMILC(l, i, cov_ok, n, j, k, folder="your_folder_here")
          treeMILC_models  <- append(treeMILC_models , list(treeMILC))
          
          #Keep track of errors
          if((LC[[3]] == "Error")|(LCT[[3]] == "Error")|(treeMILC[[3]] == "Error")){
            errors <- errors + 1
            error_vec <- c(error_vec, LC[[1]]$id)
          }
          
          print(paste0("treeMILC model ", name, " complete."))
          print(paste0("Number of errors caught: ", errors,  " in ", error_vec))
        }
      }
    }
    
    #If a certain data set does not exist, create it
    if(exists(paste0("simDat", k, "_iteration", l, "_1000"))){
      rm(list=paste0("simDat", k, "_iteration", l, "_1000"))
    }
    if(exists(paste0("simDat", k, "_iteration", l, "_10000"))){
      rm(list=paste0("simDat", k, "_iteration", l, "_10000"))
    }
  }
}

##################################################################################################
## Get results of the simulation study                                                          ##
##################################################################################################

## True proportions in simulated data
true_proportions <- c(0.6178, 0.2473, 0.1349)
names(true_proportions) <- c("Permanent", "Other", "Flexible")

## Create (true) ME matrices
ME_matrix1 <- create_ME_matrix(-log(18), -log(18), log(324), log(324), log(18), log(18)) #10% ME
ME_matrix2 <- create_ME_matrix(-3 * log(2), -3 * log(2), 6 * log(2), 6 * log(2), 3 * log(2), 3 * log(2)) #20% ME
ME_matrix3 <- create_ME_matrix(-1.54045, -1.54045, 3.0809, 3.0809, 1.54045, 1.54045) #30% ME
ME_matrix4a <- create_ME_matrix(-4.4917, -4.94368, 7.64123, 5.6678, 3.09876, 4.55064) #Realistic 7% ME
ME_matrix4b <- create_ME_matrix(-6.14311, -2.63157, 11.9482, 5.03275, 1.72427, 1.53296) #Realistic 7% ME
ME_matrix4 <- list(ME_matrix4a, ME_matrix4b) #Realistic 7% ME

#Get results
LC_results <- get_results(LC_models)
LCT_results <- get_results(LCT_models)
treeMILC_results <- get_results(treeMILC_models)

#Get summary
LC_summary <- get_summary("LC", LC_results)
LCT_summary <- get_summary("LCT", LCT_results)
treeMILC_summary <- get_summary("tree-MILC", treeMILC_results)

#Merge summaries
all_results <- rbind(LC_summary, LCT_summary, treeMILC_summary)
