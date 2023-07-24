##################################### Simulate_Data_1.R ##########################################
## This file contains two functions that are required to perform the first simulation study     ##
## (see Chapter 4):                                                                             ##
##    - simulate_data: Function to simulate a data set                                          ##
##    - create_subset: Function to create a subset from a simulated data set                    ##
##################################################################################################

##################################################################################################
## Function to simulate a data set                                                              ##
## @param seed (int): Seed                                                                      ##
## @param ME (int): Amount of measurement error (1=10%, 2=20%, 3=30%, 4=realistic)              ##
## @returns (data.frame): A simulated data set (n=10,000)                                       ##
##################################################################################################

simulate_data <- function(seed, ME) {
  
  # Create data structure file required by Latent GOLD
  filepath_input <- "exampleData.dat"
  exampleData <- "id q Y1 Y2 Y3 Y4 n
1 1 1 1 1 1 2500
2 2 2 2 2 2 2500
3 1 3 3 3 3 2500
4 2 3 3 3 3 2500"
  writeLines(exampleData, filepath_input)
  
  # Create Latent Gold script
  script_part1 <- paste0("version = 6.0\ninfile '", filepath_input, "' \n\nmodel
    title 'simulation", ME, "';
    options
    algorithm
        tolerance=1e-08 emtolerance=0.01 emiteration=1000 nriterations=1000;
    startvalues
        seed=1 sets=100 tolerance=1e-05 iterations=100;
    montecarlo
        seed=1 replicates=500 tolerance=1e-008;
    bayes
        categorical=1 variances=1 latent=1 poisson=1;
    missing includeall;
    output       
    	parameters=first standarderrors profile reorderclasses iterationdetails;\n")
  outfile_path <- paste0("simDat", ME, "_iteration", seed, ".dat")
  
  # Logit parameters for P(X=k) based on real data
  data_a2 <- -0.9662
  data_a3 <- -1.6260
  data_b22 <- 0.22
  data_b32 <- 0.2607
  data_param <- paste(data_a2, data_a3, data_b22, data_b32)
  
  # Logit parameters for P(Yj=yj|X=k) for different levels of measurement error (ME)
  if (ME == 1) {
    a2 <- a3 <- -log(18)  
    b22 <- b33 <- log(324)
    b32 <- b23 <- log(18)
    ME_coefs <- c(a2, a3, b22, b32, b23, b33, "\n")
    ME_coefs <- gsub(",", "", toString(rep(ME_coefs, 4)))
    parameters <- paste("\n", data_param, "\n", ME_coefs, "}\nend model")
  } else if (ME == 2) {
    a2 <- a3 <- -3 * log(2) 
    b22 <- b33 <- 6 * log(2)
    b32 <- b23 <- 3 * log(2)
    ME_coefs <- c(a2, a3, b22, b32, b23, b33, "\n")
    ME_coefs <- gsub(",", "", toString(rep(ME_coefs, 4)))
    parameters <- paste("\n", data_param, "\n", ME_coefs, "}\nend model")
  } else if (ME == 3) {
    a2 <- a3 <- -1.54045 
    b22 <- b33 <- 3.0809
    b32 <- b23 <- 1.54045
    ME_coefs <- c(a2, a3, b22, b32, b23, b33, "\n")
    ME_coefs <- gsub(",", "", toString(rep(ME_coefs, 4)))
    parameters <- paste("\n", data_param, ME_coefs, "}\nend model")
  } else if (ME == 4) {
    # Coefficients for indicator 1 (and 3)
    Y1_a2 <- -4.4917 
    Y1_a3 <- -4.94368
    Y1_b22 <- 7.64123
    Y1_b32 <- 4.55064
    Y1_b23 <- 3.09876
    Y1_b33 <- 5.6678
    # Coefficients for indicator 2 (and 4)
    Y2_a2 <- -6.14311 
    Y2_a3 <- -2.63157
    Y2_b22 <- 11.9482
    Y2_b32 <- 1.53296
    Y2_b23 <- 1.72427
    Y2_b33 <- 5.03275
    ME_coefs_Y1 <- c(Y1_a2, Y1_a3, Y1_b22, Y1_b32, Y1_b23, Y1_b33, "\n")
    ME_coefs_Y2 <- c(Y2_a2, Y2_a3, Y2_b22, Y2_b32, Y2_b23, Y2_b33, "\n")
    ME_coefs <- gsub(",", "", toString(rep(c(ME_coefs_Y1, ME_coefs_Y2), 2))) 
    parameters <- paste(data_param, "\n", ME_coefs, "}\nend model")
  }
  
  script_part2 <- paste0("\toutfile '", outfile_path, "' simulation=1 seed=", seed, ";
    variables
         caseid id;
         caseweight n;
         dependent Y1 nominal 3, Y2 nominal 3, Y3 nominal 3, Y4 nominal 3;
         independent q nominal;
         latent cluster nominal 3;
     equations
         cluster <- 1 + q;			
         Y1      <- 1 + cluster;	
         Y2      <- 1 + cluster;
         Y3      <- 1 + cluster;
         Y4      <- 1 + cluster;
{ ")
  
  # Combine parts of script
  script <- paste0(script_part1, script_part2, parameters)
  script_path <- paste0("simDat", ME, "_iteration", seed, "_script.lgs")
  writeLines(script, script_path)
  
  # Execute Latent Gold script
  shell(paste0('"C:/Program Files/LatentGOLDnet6.0/lg60.exe" ', script_path, ' /b'))
  
  # Import simulated data set
  simDat <- read.delim(outfile_path, sep = "\t", dec = ",")
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
  if (!exists(paste0("simDat", ME, "_iteration", iteration))) {
    assign(paste0("simDat", ME, "_iteration", iteration), simulate_data(iteration, ME), envir = globalenv())
  }
  
  data <- get(paste0("simDat", ME, "_iteration", iteration))
  
  # Set seed to get the same data set for every model within each iteration
  set.seed(iteration) 
  select_cases <- sample(1:nrow(data), N, replace = FALSE)
  
  # Remove redundant columns
  all_ind <- c("Y1", "Y2", "Y3", "Y4")
  ind <- all_ind[1:ind]
  subset <- data[select_cases, c("id", ind, cov_ok, cov_problem)]
  
  return(subset)
}
