################################ Perform_Simulation_1.R ##########################################
## This file contains the code that is required to perform the simulation study without missing ##
## covariates as described in Chapter 4. In the first part, the simulation study is performed.  ##
## In the second part, the results of the simulation study are obtained.                        ##
##                                                                                              ##
## To run the code:                                                                             ##
##    - In line 22, a working directory should be specified.                                    ##
##    - In lines 66, 71, 76, the argument 'folder' should be specified.                         ##
##    - The functions in the files 'Helpfunctions_General.R', 'Helpfunctions_Simulations.R',    ##
##      'Helpfunctions_Performance_Measures_and_Plots.R', and 'Methods_BestApproach.R'          ##
##      should be loaded.                                                                       ##
##################################################################################################

# Load packages
library(dplyr)
library(data.table)

# Ignore redundant warnings from dplyr 
options(dplyr.summarise.inform = FALSE) 

# Set working directory
setwd("your_working_directory_here") 

##################################################################################################
## Perform simulation study 1                                                                   ##
##################################################################################################

# Specification of simulation conditions
ind <- c(2:4)
N <- c(1000,10000)
ME <- c(1:4)
iteration <- c(1:50)

# Create lists to store results in
LC_models <- list()
LCT_models <- list()
treeMILC_models <- list()

# Keep track of number of models
m <- 0 

# Keep track of potential errors
errors <- 0
error_vec <- vector()

# Execute simulations
for(l in iteration){
  for(k in ME){
    for(i in ind){
      
      # Add covariate if the number of indicators is two
      if(i == 2){
        cov_ok <- "q"
      } else {
        cov_ok <- NULL
      }
      
      for(j in N){
        # Create model name
        m <- m + 1
        name <- paste0(i, "-", cov_ok, "-", j, "-", k)
        row <- c(l, i, cov, j, k, name)
        print(paste0("---------------", m, "/", length(iteration) * length(ind) * length(N) * length(ME), "------------------"))
        
        # LC
        LC <- perform_lc(l, i, cov_ok, cov_problem = NULL, j, k, folder="add_your_folder_here")
        LC_models <- append(LC_models, list(LC))
        print(paste0("LC model ", name, " complete."))
        
        # LCT
        LCT <- perform_lct(l, i, cov_ok, cov_problem = NULL, j, k, folder="add_your_folder_here")
        LCT_models <- append(LCT_models, list(LCT))
        print(paste0("LCT model ", name, " complete."))
        
        # tree-MILC
        treeMILC <- perform_treeMILC(l, i, cov_ok, cov_problem = NULL, j, k, folder="add_your_folder_here")
        treeMILC_models <- append(treeMILC_models, list(treeMILC))
        print(paste0("tree-MILC model ", name, " complete."))
        
        # Keep track of errors
        if((LC[[3]] == "Error") | (LCT[[3]] == "Error") | (treeMILC[[3]] == "Error")){
          errors <- errors + 1
          error_vec <- c(error_vec, LC[[1]]$id)
        }
        print(paste0("Number of errors caught: ", errors, " in ", error_vec))
        
      }
    }
    
    # Remove data set from global environment if not needed anymore   
    if(exists(paste0("simDat", k, "_iteration", l))){
      rm(list=paste0("simDat", k, "_iteration", l))
    }
  }
}

##################################################################################################
## Get results of the simulation study                                                          ##
##################################################################################################

# True proportions in simulated data
true_proportions <- c(0.6061153, 0.2577143, 0.1361704)
names(true_proportions) <- c("Permanent", "Other", "Flexible")

# Create (true) ME matrices
ME_matrix1 <- create_ME_matrix(-log(18), -log(18), log(324), log(324), log(18), log(18)) #10% ME
ME_matrix2 <- create_ME_matrix(-3*log(2), -3*log(2), 6*log(2), 6*log(2), 3*log(2), 3*log(2)) #20% ME
ME_matrix3 <- create_ME_matrix(-1.54045, -1.54045, 3.0809, 3.0809, 1.54045, 1.54045) #30% ME
ME_matrix4a <- create_ME_matrix(-4.4917, -4.94368, 7.64123, 5.6678, 3.09876, 4.55064) #Realistic 7% ME
ME_matrix4b <- create_ME_matrix(-6.14311, -2.63157, 11.9482, 5.03275, 1.72427, 1.53296) #Realistic 7% ME
ME_matrix4 <- list(ME_matrix4a, ME_matrix4b) #Realistic 7% ME

# Get results
LC_results <- get_results(LC_models)
LCT_results <- get_results(LCT_models)
treeMILC_results <- get_results(treeMILC_models)

# Get summary
LC_summary <- get_summary("LC", LC_results)
LCT_summary <- get_summary("LCT", LCT_results)
treeMILC_summary <- get_summary("tree-MILC", treeMILC_results)

# Merge summaries
all_results <- rbind(LC_summary, LCT_summary, treeMILC_summary)
