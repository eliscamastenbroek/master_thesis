##################################################################################################
## This file contains the code that is required to perform the simulation study without missing ##
## covariates as described in Chapter 5. The file is divided into three parts:                  ##
##   1. Functions 'simulate_data' and 'create_subset' to simulate data.                         ##
##   2. Perform the simulation study.                                                           ##
##   3. Get the results of the simulation study.                                                ## 
##                                                                                              ##
## Note that depending on whether the functions from the file 'Methods_LessOptimalApproach.R'   ##
## or the file 'Methods_BestApproach.R' are loaded, missing covariates are included using       ##
## either the less optimal approach (see Section 5.1.1) or the best approach                    ##
## (see Section 5.1.2) with direct effects and parameter restrictions.                          ##
##                                                                                              ##
## To run the code:                                                                             ##
##    - In line 32, a working directory should be set.                                          ##
##    - In lines 236, 240, 244, the argument 'folder' should be specified.                      ##
##    - The functions in the files 'Helpfunctions_General.R', 'Helpfunctions_Simulations.R',    ##
##      'Helpfunctions_Performance_Measures_and_Plots.R' should be loaded.                      ##
##    - The functions in the file 'Methods_LessOptimalApproach.R' OR 'Methods_BestApproach.R'   ##
##      should be loaded.                                                                       ##
##    - The files 'exampleDat_1000.dat' and 'exampleDat_10000.dat' should be in the             ##
##      working directory.                                                                      ##
##################################################################################################

# Initialisations
library(dplyr)
library(data.table)

# Ignore redundant warnings from dplyr
options(dplyr.summarise.inform = FALSE)  

# Set working directory
setwd("your_working_directory_here") 

##################################################################################################
## Simulate a data set                                                                          ##
## @param seed (int): Seed                                                                      ##
## @param N (int): Size of the data set                                                         ##
## @param ME (int): Amount of measurement error (1=10%, 2=20%, 3=30%, 4=realistic)              ##
## @param folder (string): Folder to save files in                                              ##
## @returns (data.frame): A simulated data set (n=10,000)                                       ##
##################################################################################################

simulate_data <- function(seed, N, ME) {
  
  # Select the correct exampleDat.dat file depending on N
  if(N == 1000){
    filepath_input <- paste0("exampleDat_1000.dat")
  } else {
    filepath_input <- paste0("exampleDat_10000.dat")
  }
  
  # Create Latent Gold script for LC 
  script_part1 <- paste0("version = 6.0\ninfile '", filepath_input, "' \n\nmodel
    title 'simulation", ME, "';
    options
    algorithm
        tolerance=1e-08 emtolerance=0.01 emiteration=500 nriterations=500;
    startvalues
        seed=1 sets=100 tolerance=1e-05 iterations=100;
    montecarlo
        seed=1 replicates=500 tolerance=1e-008;
    bayes
        categorical=1 variances=1 latent=1 poisson=1;
    missing includeall;
    output       
    	parameters=first standarderrors profile reorderclasses iterationdetails;\n")
  
  # Specify parameter estimates
  contract <- c(0.8252, -0.3192)
  q <- c(-0.1506, -0.1584) 
  SBIgroep <- c(-2.5392, 2.2426)
  baanduur <- c(-4.6105, -3.4073)
  data_param <- gsub(",", "", toString(c(contract, q, SBIgroep, baanduur)))
  
  # Parameters for different levels of measurement error
  if (ME == 1) {
    outfile_name <- paste0("simDat1_iteration", seed, ".dat")
    a2 <- a3 <- -log(18)  # Coefficients for measurement error matrix
    b22 <- b33 <- log(324)
    b32 <- b23 <- log(18)
    ME_coefs <- c(a2, a3, b22, b32, b23, b33, "\n")
    ME_coefs <- gsub(",", "", toString(rep(ME_coefs, 3)))
    parameters <- paste("\n", data_param, "\n", ME_coefs, "}\nend model")
  } else if (ME == 2) {
    outfile_name <- paste0("simDat2_iteration", seed, ".dat")
    a2 <- a3 <- -3 * log(2)  # Coefficients for measurement error matrix
    b22 <- b33 <- 6 * log(2)
    b32 <- b23 <- 3 * log(2)
    ME_coefs <- c(a2, a3, b22, b32, b23, b33, "\n")
    ME_coefs <- gsub(",", "", toString(rep(ME_coefs, 3)))
    parameters <- paste("\n", data_param, "\n", ME_coefs, "}\nend model")
  } else if (ME == 3) {
    outfile_name <- paste0("simDat3_iteration", seed, ".dat")
    a2 <- a3 <- -1.54045  # Coefficients for measurement error matrix
    b22 <- b33 <- 3.0809
    b32 <- b23 <- 1.54045
    ME_coefs <- c(a2, a3, b22, b32, b23, b33, "\n")
    ME_coefs <- gsub(",", "", toString(rep(ME_coefs, 3)))
    parameters <- paste("\n", data_param, "\n", ME_coefs, "}\nend model")
  } else if (ME == 4) {
    outfile_name <- paste0("simDat4_iteration", seed, ".dat")
    Y1_a2 <- -4.4917
    Y1_a3 <- -4.94368
    Y1_b22 <- 7.64123
    Y1_b32 <- 4.55064
    Y1_b23 <- 3.09876
    Y1_b33 <- 5.6678
    Y2_a2 <- -6.14311
    Y2_a3 <- -2.63157
    Y2_b22 <- 11.9482
    Y2_b32 <- 1.53296
    Y2_b23 <- 1.72427
    Y2_b33 <- 5.03275
    ME_coefs_Y1 <- c(Y1_a2, Y1_a3, Y1_b22, Y1_b32, Y1_b23, Y1_b33, "\n")
    ME_coefs_Y2 <- c(Y2_a2, Y2_a3, Y2_b22, Y2_b32, Y2_b23, Y2_b33, "\n")
    ME_coefs <- gsub(",", "", toString(rep(c(ME_coefs_Y1, ME_coefs_Y2, ME_coefs_Y1), 1))) 
    parameters <- paste(data_param, "\n", ME_coefs, "}\nend model")
  }
  
  script_part2 <- paste0("\toutfile '", outfile_name, "' simulation=1 seed=", seed, ";
    variables
         caseid id;
         caseweight w;
         dependent Y1 nominal 3, Y2 nominal 3, Y3 nominal 3;
         independent q nominal, SBIgroep nominal, baanduur nominal;
         latent cluster nominal 3;
     equations
         cluster <- 1 + q + SBIgroep + baanduur;			
         Y1      <- 1 + cluster;	
         Y2      <- 1 + cluster;
         Y3      <- 1 + cluster;
{ ")
  
  # Combine parts of the script
  script <- paste0(script_part1, script_part2, parameters)
  writeLines(script, paste0("simDat", ME, "_iteration", seed, "_script.lgs"))
  
  # Execute Latent Gold script
  shell(paste0('"C:/Program Files/LatentGOLDnet6.0/lg60.exe" ', paste0("simDat", ME, "_iteration", seed, "_script.lgs"), ' /b'))
  
  # Import simulated data set
  simDat <- as.data.frame(fread(outfile_name, dec = ","))
  
  # Add extra 'problem covariate' category
  id_other_Y1 <- simDat[simDat$Y1 == 2, ]$id   
  id_add_extra_cat <- sample(id_other_Y1, (0.9 * length(id_other_Y1))) 
  ncat1 <- length(levels(factor(simDat[, which(names(simDat) == "baanduur")]))) 
  ncat2 <- length(levels(factor(simDat[, which(names(simDat) == "SBIgroep")]))) 
  simDat[id_add_extra_cat, which(names(simDat) == "baanduur")] <- ncat1 + 1  
  simDat[id_add_extra_cat, which(names(simDat) == "SBIgroep")] <- ncat2 + 1   
  
  return(simDat)
}

##################################################################################################
## Function to create a subset from a simulated data set                                        ##
## @param iteration (int): Iteration number                                                     ## 
## @param ind (int): Number of indicators                                                       ##
## @param cov_ok (vector): Vector with the names of non-missing covariates                      ##
## @param cov_problem (vector): Vector with the names of missing covariates                     ##
## @param N (int): Size of the data set                                                         ##
## @param ME (int): Amount of measurement error (1=10%, 2=20%, 3=30%, 4=realistic)              ##
## @returns (data.frame): A subset                                                              ##
##################################################################################################

create_subset <- function(iteration, ind, cov_ok, cov_problem, N, ME) {
  
  # If a certain data set does not exist, create it  
  if(!exists(paste0("simDat", ME, "_iteration", iteration,"_", N))){
    assign(paste0("simDat", ME, "_iteration", iteration, "_", N), simulate_data(iteration, N, ME), envir=globalenv())
  }
  
  data <- get(paste0("simDat", ME, "_iteration", iteration, "_", N))
  
  # Set seed to get the same data set for every model within each iteration
  set.seed(iteration) 
  select_cases <- sample(1:nrow(data), N, replace = FALSE)
  
  # Remove redundant columns
  all_ind <- c("Y1", "Y2", "Y3", "Y4")
  ind <- all_ind[1:ind]
  subset <- data[select_cases, c("id", ind, cov_ok, cov_problem)]
  
  return(subset)
}

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
## 3. Get results of the simulation study                                                       ##
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
