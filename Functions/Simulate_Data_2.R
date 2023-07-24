####################################### Simulate_Data_2.R ########################################
## This file contains two functions that are required to perform the second simulation study    ##
## (including the initial analyses) (see Chapter 5):                                            ##
##    - simulate_data: Function to simulate a data set                                          ##
##    - create_subset: Function to create a subset from a simulated data set                    ##
##################################################################################################

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
