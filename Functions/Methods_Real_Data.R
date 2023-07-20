################################### Methods_RealData.R ###########################################
## This file contains the functions that are required to perform LC and tree-MILC analysis on   ##
## real data from the ER and the LFS (see Chapter 6):                                           ##
##    - perform_lc                                                                              ##
##    - perform_treeMILC                                                                        ##
##################################################################################################

##################################################################################################
## Perform LC analysis                                                                          ##
## @param cov_ok (vector): Vector with the names of non-missing covariates                      ##
## @param cov_problem (vector): Vector with the names of missing covariates                     ##
## @param dat (data.frame) (optional): Data frame to perform LC analysis on                     ## 
## @param dat_path (character) (optional): Path to data set to perform LC analyis on            ##
## @param folder (string): Folder to save files in                                              ##
## @returns (list): A list that consists of:                                                    ## 
##  [[1]] Data frame with an overview of model parameters (iteration, ind, cov, N, ME)          ##
##  [[2]] Data frame with model results (posterior probabilities for each observation)          ##
##  [[3]] Character indicating whether the model is "Good" or whether there was an "Error"      ## 
##################################################################################################

perform_lc <- function(cov_ok, cov_problem, dat = NULL, dat_path = NULL, folder) {
  
  # Store model information
  to_return <- store_model_info(cov_ok, cov_problem)
  model_name <- to_return[[1]]$id
  
  # Check if a data set or a filepath is provided
  if (is.null(dat) && is.null(dat_path)) {
    stop("Input required!")
  } else if (is.null(dat)) {
    dat_path <- dat_path
    dat <- as.data.frame(fread(dat_path, sep = "\t", dec = ","))
    dat <- dat[, which(colnames(dat) %in% c("contract", "contractEBB", "persnr", cov_ok, cov_problem))]
  } else {
    dat <- dat[, which(colnames(dat) %in% c("contract", "contractEBB", "persnr", cov_ok, cov_problem))]
    dat_path <- paste0(folder, "LC_", model_name, "_data.dat")
    fwrite(dat, file = dat_path, sep = "\t")
  }
  
  # Create Latent Gold script for LC
  filepath_output <- paste0(folder, paste0("LC_", model_name, "_output.dat"))
  script <- generate_script("LC", cov_ok, cov_problem, model_name, dat = dat, dat_path, filepath_output)
  script_path <- paste0(folder, "LC_", model_name, "_script.lgs")
  writeLines(script, script_path)
  
  # Execute Latent Gold script
  shell(paste0('"C:/Program Files/LatentGOLDnet6.0/lg60.exe" ', script_path, ' /b'))
  
  # Read model output
  model_output <- as.data.frame(fread(filepath_output, dec = ","))
  
  # Rename some columns and add data frame to return list
  setnames(model_output, old = c("Cluster#1", "Cluster#2", "Cluster#3", "Cluster#"), new = c("p1", "p2", "p3", "cluster"))
  to_return <- append(to_return, list(model_output))
  
  # Make sure clusters are assigned the right names
  to_return <- fix_cluster_assignment(type = "LC", results = to_return)
  
  # Check if no warning was given
  model_lst <- paste(readLines(paste0(folder, "LC_", model_name, "_script.lst")), collapse = "\n")
  
  if (grepl("WARNING", model_lst, fixed = TRUE, useBytes = TRUE)) {
    to_return <- append(to_return, list("Error"))
  } else {
    to_return <- append(to_return, list("Good"))
  }
  
  return(to_return)
}

##################################################################################################
## Perform tree-MILC analysis                                                                   ##
## @param cov_ok (vector): Vector with the names of non-missing covariates                      ##
## @param cov_problem (vector): Vector with the names of missing covariates                     ##
## @param dat (data.frame) (optional): Data frame to perform LC analysis on                     ## 
## @param dat_path (character) (optional): Path to data set to perform LC analyis on            ##
## @param folder (string): Folder to save files in                                              ##
## @param M (int): Number of bootstrap samples                                                  ##
## @returns (list): A list that consists of:                                                    ## 
##  [[1]] Data frame with an overview of model parameters (iteration, ind, cov, N, ME)          ##
##  [[2]] Data frame with model results (posterior probabilities for each observation)          ##
##  [[3]] Character indicating whether the model is "Good" or whether there was an "Error"      ## 
##################################################################################################

perform_treeMILC <- function(cov_ok, cov_problem, dat = NULL, dat_path = NULL, folder, M = 5) {
  
  # Store model information
  to_return <- store_model_info(cov_ok, cov_problem)
  model_name <- to_return[[1]]$id
  
  # Create help vector to combine results later
  all_indicators <- c("contract", "contractEBB")
  by_vector <- c(all_indicators, cov_ok, cov_problem)
  
  # Check if a data set or a filepath is provided
  if (is.null(dat) && is.null(dat_path)) {
    stop("Input required!")
  }
  
  # If a filepath is provided, obtain data set and keep only relevant columns
  if (is.null(dat)) {
    dat_org <- as.data.frame(fread(dat_path, sep = "\t", dec = ","))
    dat_org <- dat_org[, which(colnames(dat_org) %in% c(by_vector, "persnr"))]
  } else {
    dat_org <- dat[, which(colnames(dat) %in% c(by_vector, "persnr"))]
    dat_path <- paste0(folder, "LC_", model_name, "_data.dat")
    fwrite(dat, file = dat_path, sep = "\t")
  }
  
  # Count combinations of indicators and covariates in the original data set
  count_dat <- as.data.frame(dat_org[, which(colnames(dat_org) %in% by_vector)] %>% group_by_all() %>% summarise(COUNT = n()))
  count_dat <- count_dat[, -ncol(count_dat)]
  
  # Create list to store the results for each bootstrap sample
  bootstrap_results <- list()
  
  # For each bootstrap sample
  for (i in 1:M) {
    dat <- dat_org
    set.seed(10 + i) # Set seed to get different bootstrap samples
    sample_ids <- sample(1:nrow(dat), nrow(dat), replace = TRUE)
    sample <- dat[sample_ids, ]
    sample$persnr <- row.names(sample) # Change id to capture duplicates
    boot_name <- paste0("tree_MILC", model_name, "_boot", i)
    boot_output <- perform_lc(cov_ok, cov_problem, dat = sample, folder = folder) # Perform LC with 3 classes and combine posterior probabilities
    
    # If model contains an error, add error to model output
    if (boot_output[[3]] == "Error") {
      to_return <- append(to_return, list("Error"))
      to_return <- append(to_return, list("Error"))
      return(to_return)
    }
    
    boot_output <- boot_output[[2]]
    
    # Count combinations of indicators and covariates in the bootstrap sample
    count_boot <- as.data.frame(boot_output[, -which(colnames(boot_output) %in% c("persnr", "cluster"))] %>% group_by_all() %>% summarise(COUNT = n()))
    count_boot <- count_boot[, -ncol(count_boot)]
    
    # If not all combinations are present in the bootstrap sample
    if (nrow(count_boot) != nrow(count_dat)) {
      # Get data frame with parameters
      filename_lc <- paste0(folder, "LC_", model_name, "_script.lst")
      lc_lst <- readChar(filename_lc, file.info(filename_lc)$size)
      lc_lst <- strsplit(lc_lst, split = "Regression Parameters")
      lc_lst <- strsplit(lc_lst[[1]][2], split = "Paired Comparisons")
      lc_lst <- lc_lst[[1]][1]
      parameters_path <- paste0(folder, "tree_MILC_", model_name, "_boot", i, "_parameters.dat")
      writeLines(lc_lst, parameters_path)
      parameters <- suppressWarnings(fread(parameters_path, sep = "\t", dec = ","))[, 1:6]
      names(parameters) <- c("term1", "term2", "term3", "term4", "term5", "coef")
      parameters$coef <- as.numeric(parameters$coef)
      remove <- which(grepl("Cluster(1)", parameters$term1, fixed = TRUE, useBytes = TRUE) | grepl("(1)", parameters$term3, fixed = TRUE, useBytes = TRUE) | grepl("contractEBB(1)", parameters$term1, fixed = TRUE, useBytes = TRUE) | (grepl("Cluster(1)", parameters$term5, fixed = TRUE, useBytes = TRUE) & (grepl("contract(1)", parameters$term1, fixed = TRUE, useBytes = TRUE))) | (grepl("Cluster(2)", parameters$term5, fixed = TRUE, useBytes = TRUE) & (grepl("contract(1)", parameters$term1, fixed = TRUE, useBytes = TRUE))) | (grepl("Cluster(3)", parameters$term5, fixed = TRUE, useBytes = TRUE) & (grepl("contract(1)", parameters$term1, fixed = TRUE, useBytes = TRUE))))
      parameters <- parameters[-remove, ]
      parameters$term1 <- do.call(paste0, parameters[, c("term1", "term2", "term3", "term4", "term5")])
      parameters <- parameters[, c(1, 6)]
      parameters <- as.data.frame(parameters)
      
      # Estimate extra LC model with obtained parameters as starting values
      output_path <- paste0(folder, "tree_MILC_", model_name, "_boot", i, "_dat_org_posteriors.dat")
      script <- generate_script_treeMILC_extra(cov_ok, cov_problem, model_name, dat = dat, filepath_input = dat_path, filepath_output = output_path, parameters)
      script_path <- paste0(folder, "tree_MILC_", model_name, "_boot", i, "_dat_org_posteriors.lgs")
      writeLines(script, script_path)
      shell(paste0('"C:/Program Files/LatentGOLDnet6.0/lg60.exe" ', script_path, ' /b'))
      
      # Read output of extra LC model that contains posterior probabilities for every observation in the original data set
      dat <- as.data.frame(fread(output_path, sep = "\t", dec = ","))
      setnames(dat, old = c("Cluster#1", "Cluster#2", "Cluster#3", "Cluster#"), new = c("p1", "p2", "p3", "cluster"))
      dat <- fix_cluster_assignment(type = "LC", list(1, dat))
      dat <- dat[[2]]
      dat <- dat[, -which(colnames(dat) == "cluster")]
    } else {
      dat <- left_join(x = dat, y = count_boot, by = by_vector)
    }
    
    # Sample from obtained posterior membership probabilities
    dat$imp1 <- apply(dat, 1, impute_value_step1)
    dat$persnr <- as.character(dat$persnr)
    
    # Create subsets of cases with imputed values of 1 (with original indicator values```R
    # Create subsets of cases with imputed values of 1 (with original indicator values) and write to file
    subset_ids <- dat[dat$imp1 == 1, ]$persnr
    subset <- dat[dat$persnr %in% subset_ids, ]
    filepath_subset <- paste0(folder, "tree_MILC", model_name, "_boot", i, "_subset.dat")
    fwrite(subset, file = filepath_subset, sep = "\t")
    
    # Create and execute Latent Gold script for the second model
    filepath_output <- paste0(folder, boot_name, "_step2_output.dat")
    script <- generate_script("treeMILC", cov_ok, cov_problem, model_name, dat = subset, filepath_subset, filepath_output)
    script_path <- paste0(folder, boot_name, "_step2_script.lgs")
    writeLines(script, script_path)
    shell(paste0('"C:/Program Files/LatentGOLDnet6.0/lg60.exe" ', script_path, ' /b'))
    
    # Import model output
    model2_output <- as.data.frame(fread(paste0(folder, boot_name, "_step2_output.dat"), dec = ","))
    
    # Check if a warning was given
    model_lst <- paste(readLines(paste0(folder, boot_name, "_step2_script.lst")), collapse = "\n")
    if (grepl("WARNING", model_lst, fixed = TRUE, useBytes = TRUE)) {
      to_return <- append(to_return, list("Error"))
      to_return <- append(to_return, list("Error"))
      return(to_return)
    }
    
    # Add the results from the second model to the results from the first model
    count_step2 <- as.data.frame(model2_output[, -which(colnames(boot_output) %in% c("persnr", "Cluster#"))] %>% group_by_all() %>% summarise(COUNT = n()))
    count_step2 <- count_step2[, -ncol(count_step2)]
    dat <- left_join(x = dat, y = count_step2, by = c(all_indicators, cov_ok, cov_problem))
    
    # Sample again from obtained posterior membership probabilities
    dat$imp2 <- apply(dat, 1, impute_value_step2)
    
    #Combine imputations from both models 
    dat$cluster <- dat$imp1
    dat[dat$imp1 == 1 & ((!is.na(dat$imp2) & dat$imp2 == 1)),]$cluster <- 1
    dat[dat$imp1 == 1 & ((!is.na(dat$imp2) & dat$imp2 == 2)),]$cluster <- 3
    
    #Select and rename columns
    new_names <- c("step2_p1","step2_p3")
    setnames(dat,old = c("Cluster#1","Cluster#2"),new = new_names)
    dat <- dat[, -which(colnames(dat) == "Cluster#")]
    
    #Fix cluster assignment (if necessary)
    dat <- fix_cluster_bootstrap("real", dat)
    bootstrap_results <- append(bootstrap_results, list(dat))
  }
  
  #Add results to final output list
  to_return <- append(to_return, list(bootstrap_results))
  return(to_return)
}
