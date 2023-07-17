################################## Methods_BestApproach.R ########################################
## This file contains the functions required to perform LC, LCT and tree-MILC analysis using    ##
## the best approach for including missing covariates from Chapter 5 (see Section               ##
## 5.1.2). The current file contains the following functions:                                   ##
##    - generate_script                                                                         ##
##    - perform_lc                                                                              ##
##    - perform_lct                                                                             ##
##    - perform_treeMILC                                                                        ##
##################################################################################################

##################################################################################################
## Perform LC analysis                                                                          ##
## @param iteration (int): Iteration number                                                     ## 
## @param ind (int): Number of indicators                                                       ## 
## @param cov_ok (vector): Vector with the names of non-missing covariates                      ##
## @param cov_problem (vector): Vector with the names of missing covariates                     ##
## @param N (int): Size of data set                                                             ## 
## @param ME (int): Amount of measurement error (1=10%, 2=20%, 3=30%, 4=realistic)              ##
## @param folder (string): Folder to save files in                                              ##
## @returns (list): A list that consists of:                                                    ## 
##  [[1]] Data frame with an overview of model parameters (iteration, ind, cov, N, ME)          ##
##  [[2]] Data frame with model results (posterior probabilities for each observation)          ##
##  [[3]] Character indicating whether the model is "Good" or whether there was an "Error"     ## 
##################################################################################################

perform_lc <- function(iteration, ind, cov_ok, cov_problem, N, ME, dat = NULL, folder) {
  
  # Store model information
  to_return <- store_model_info(iteration, ind, cov_ok, cov_problem, N, ME)
  model_name <- to_return[[1]]$id
  
  # Write data set to file
  if (is.null(dat)) {
    dat <- create_subset(iteration, ind, cov_ok, cov_problem, N, ME)
  } else {
    dat <- dat
  }
  
  filepath_input <- paste0(folder, "LC_", model_name, "_data.dat")
  fwrite(dat, file = filepath_input, sep = "\t")
  
  # Create Latent Gold script for LC
  filepath_output <- paste0(folder, paste0("LC_", model_name, "_output.dat"))
  script <- generate_script("LC", ind, cov_ok, cov_problem, N, dat = dat, model_name, filepath_input, filepath_output)
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
## Perform LCT analysis                                                                         ##
## @param iteration (int): Iteration number                                                     ## 
## @param ind (int): Number of indicators                                                       ## 
## @param cov_ok (vector): Vector with the names of non-missing covariates                      ##
## @param cov_problem (vector): Vector with the names of missing covariates                     ##
## @param N (int): Size of data set                                                             ## 
## @param ME (int): Amount of measurement error (1=10%, 2=20%, 3=30%, 4=realistic)              ##
## @param folder (string): Folder to save files in                                              ##
## @returns (list): A list that consists of:                                                    ## 
##  [[1]] Data frame with an overview of model parameters (iteration, ind, cov, N, ME)          ##
##  [[2]] Data frame with model results (posterior probabilities for each observation)          ##
##  [[3]] Character indicating whether the model is "Good" or whether there was an "Error"     ## 
##################################################################################################

perform_lct <- function(iteration, ind, cov_ok, cov_problem, N, ME, folder) {
  
  # Store model information
  to_return <- to_return1 <- store_model_info(iteration, ind, cov_ok, cov_problem, N, ME)
  model_name <- to_return[[1]]$id
  
  # Write data set to use to file
  dat <- create_subset(iteration, ind, cov_ok, cov_problem, N, ME)
  dat_name <- paste0("LCT_", model_name, "_step1_data.dat")
  fwrite(dat, file = paste0(folder, dat_name), sep = "\t")
  
  # Perform LC with 3 classes
  model1 <- perform_lc(iteration, ind, cov_ok, cov_problem, N, ME, folder = folder) # Since the same seed is used, the data set used in perform_lc should be identical
  model1_output <- model1[[2]] 
  
  # Check if no warning was given
  if (model1[[3]] == "Error") {
    to_return <- append(to_return, list("Error"))
    to_return <- append(to_return, list("Error"))
    return(to_return)
  }
  
  # Combine posterior probabilities for the classes 'permanent' and 'flexible' in the model output
  model1_output$p1 <- 1 - model1_output$p2 # Note that this is actually p1 + p3, but we call it p1 for the 'generate_script' function to work
  
  # Write data set to use in second step to file
  dat_name <- paste0("LCT_", model_name, "_step1_output.dat")
  fwrite(model1_output, file = paste0(folder, dat_name), sep = "\t")
  
  # Create Latent Gold script for second LC model
  filepath_input <- paste0(folder, dat_name)
  filepath_output <- paste0(folder, "LCT_", model_name, "_step2_output.dat")
  script <- generate_script("LCT", ind, cov_ok, cov_problem, N, model_name, dat, filepath_input, filepath_output)
  script_path <- paste0(folder, "LCT_", model_name, "_step2_script.lgs")
  writeLines(script, script_path)
  
  # Execute script in Latent Gold and read model output
  shell(paste0('"C:/Program Files/LatentGOLDnet6.0/lg60.exe" ', script_path, ' /b'))
  model2_output <- as.data.frame(fread(filepath_output, dec = ","))
  
  # Combine results from both models to compute final posterior probabilities and fix column names
  all_indicators <- c("Y1", "Y2", "Y3", "Y4")
  by_vector <- c("id", all_indicators[1:ind], c(cov_ok, cov_problem))
  combined_output <- left_join(x = model1_output, y = model2_output, by = by_vector)
  new_names <- c("p1.1", "p2.1", "p2.2", "cluster.2", "cluster.1")
  setnames(combined_output, old = c("p1.y", "Cluster#1", "Cluster#2", "Cluster#", "cluster"), new = new_names)
  combined_output <- combined_output[, c(by_vector, new_names)]
  combined_output$p2.2 <- fix_number_notation(combined_output$p2.2)
  combined_output$p2.1 <- fix_number_notation(combined_output$p2.1)
  
  if (!is.numeric(combined_output$cluster.2)) {
    combined_output[combined_output$cluster.2 == ".", ]$cluster.2 <- 2
    combined_output$cluster.2 <- as.numeric(combined_output$cluster.2)
  }
  
  # Compute posterior probabilities and combine cluster assignments from step 1 and step 2
  combined_output$p1 <- combined_output$p2 <- combined_output$p3 <- combined_output$cluster <- NA
  
  for (i in 1:nrow(combined_output)) {
    combined_output[i, ]$p2 <- 1 - as.numeric(combined_output[i, ]$p1.1)
    combined_output[i, ]$p1 <- as.numeric(combined_output[i, ]$p1.1) * as.numeric(combined_output[i, ]$p2.1)
    combined_output[i, ]$p3 <- as.numeric(combined_output[i, ]$p1.1) * as.numeric(combined_output[i, ]$p2.2)
    
    if (combined_output[i, ]$cluster.1 == 2) {
      combined_output[i, ]$cluster <- 2
    } else if ((combined_output[i, ]$cluster.1 == 1 & combined_output[i, ]$cluster.2 == 1) | (combined_output[i, ]$cluster.1 == 3 & combined_output[i, ]$cluster.2 == 1)) {
      combined_output[i, ]$cluster <- 1
    } else if ((combined_output[i, ]$cluster.1 == 1 & combined_output[i, ]$cluster.2 == 2) | (combined_output[i, ]$cluster.1 == 3 & combined_output[i, ]$cluster.2 == 2)) {
      combined_output[i, ]$cluster <- 3
    }
  }
  
  # Remove redundant columns
  remove_cols <- c("cluster.1", "cluster.2", "p2.1", "p2.2", "p1.1")
  combined_output <- combined_output[, -which(names(combined_output) %in% remove_cols)]
  to_return <- append(to_return, list(combined_output))
  
  # Make sure clusters are assigned the right names
  to_return <- fix_cluster_assignment(type = "LCT", results = to_return)
  
  # Check if no warning was given
  model_lst <- paste(readLines(paste0(folder, "LCT_", model_name, "_step2_script.lst")), collapse = "\n")
  if (grepl("WARNING", model_lst, fixed = TRUE, useBytes = TRUE)) {
    to_return <- append(to_return, list("Error"))
  } else {
    to_return <- append(to_return, list("Good"))
  }
  
  return(to_return)
}

##################################################################################################
## Perform tree-MILC analysis                                                                   ##
## @param iteration (int): Iteration number                                                     ## 
## @param ind (int): Number of indicators                                                       ## 
## @param cov_ok (vector): Vector with the names of non-missing covariates                      ##
## @param cov_problem (vector): Vector with the names of missing covariates                     ##
## @param N (int): Size of data set                                                             ## 
## @param ME (int): Amount of measurement error (1=10%, 2=20%, 3=30%, 4=realistic)              ##
## @param M (int): Number of bootstrap samples                                                  ##
## @param folder (string): Folder to save files in                                              ##
## @returns (list): A list that consists of:                                                    ## 
##  [[1]] Data frame with an overview of model parameters (iteration, ind, cov, N, ME)          ##
##  [[2]] Data frame with model results (posterior probabilities for each observation)          ##
##  [[3]] Character indicating whether the model is "Good" or whether there was an "Error"     ## 
##################################################################################################

perform_treeMILC <- function(iteration, ind, cov_ok, cov_problem, N, ME, M = 5, folder) {
  
  # Store model information
  to_return <- store_model_info(iteration, ind, cov_ok, cov_problem, N, ME)
  model_name <- to_return[[1]]$id
  
  # Create subset and write to file
  dat_org <- create_subset(iteration, ind, cov_ok, cov_problem, N, ME)
  dat_org_path <- paste0(folder, "tree_MILC_", model_name, "_dat_org.dat")
  fwrite(dat_org, file = dat_org_path, sep = "\t")
  
  # Count combinations of indicators and covariates in the original data set
  count_dat <- as.data.frame(dat_org[, -which(colnames(dat_org) %in% c("id"))] %>%
                               group_by_all() %>%
                               summarise(COUNT = n()))
  count_dat <- count_dat[, -ncol(count_dat)]
  
  # Create help vector to combine results later
  all_indicators <- c("Y1", "Y2", "Y3", "Y4")
  by_vector <- c(all_indicators[1:ind], cov_ok, cov_problem)
  
  # Create list to store the results for each bootstrap sample
  bootstrap_results <- list()
  
  # For each bootstrap sample
  for (i in 1:M) {
    dat <- dat_org
    set.seed(i)  # Set seed to get different bootstrap samples
    sample <- dat[sample(1:N, N, replace = TRUE), ]
    boot_name <- paste0("tree_MILC", model_name, "_boot", i)
    
    # Perform LC with 3 classes and combine posterior probabilities
    boot_output <- perform_lc(iteration, ind, cov_ok, cov_problem, N, ME, dat = sample, folder = folder)
    
    # If model contains an error, add error to model output
    if (boot_output[[3]] == "Error") {
      to_return <- append(to_return, list("Error"))
      to_return <- append(to_return, list("Error"))
      return(to_return)
    }
    
    boot_output <- boot_output[[2]]
    
    # Count combinations of indicators and covariates in the bootstrap sample
    count_boot <- as.data.frame(boot_output[, -which(colnames(boot_output) %in% c("id", "cluster"))] %>%
                                  group_by_all() %>%
                                  summarise(COUNT = n()))
    count_boot <- count_boot[, -ncol(count_boot)]
    
    # If not all combinations are present in the bootstrap sample
    if (nrow(count_boot) != nrow(count_dat)) {
      # Get data frame with parameters
      filename_lc <- paste0(folder, "LC_", model_name, "_script.lst")
      lc_lst <- readChar(filename_lc, file.info(filename_lc)$size)
      lc_lst <- strsplit(lc_lst, split = "Regression Parameters")
      lc_lst <- strsplit(lc_lst[[1]][2], split = "Paired Comparisons")
      parameters_path <- paste0(folder, "tree_MILC_", model_name, "_boot", i, "_parameters.dat")
      writeLines(lc_lst[[1]][1], parameters_path)
      parameters <- suppressWarnings(fread(parameters_path, sep = "\t", dec = ","))[, 1:6]
      names(parameters) <- c("term1", "term2", "term3", "term4", "term5", "coef")
      parameters$coef <- as.numeric(parameters$coef)
      remove <- which(grepl("Cluster(1)", parameters$term1, fixed = TRUE, useBytes = TRUE) |
                        grepl("Y2(1)", parameters$term1, fixed = TRUE, useBytes = TRUE) |
                        grepl("Y3(1)", parameters$term1, fixed = TRUE, useBytes = TRUE) |
                        grepl("Y4(1)", parameters$term1, fixed = TRUE, useBytes = TRUE) |
                        grepl("(1)", parameters$term3, fixed = TRUE, useBytes = TRUE) |
                        (grepl("Cluster(1)", parameters$term5, fixed = TRUE, useBytes = TRUE) &
                           (grepl("Y1(1)", parameters$term1, fixed = TRUE, useBytes = TRUE))) |
                        (grepl("Cluster(2)", parameters$term5, fixed = TRUE, useBytes = TRUE) &
                           (grepl("Y1(1)", parameters$term1, fixed = TRUE, useBytes = TRUE))) |
                        (grepl("Cluster(3)", parameters$term5, fixed = TRUE, useBytes = TRUE) &
                           (grepl("Y1(1)", parameters$term1, fixed = TRUE, useBytes = TRUE))))
      parameters <- as.data.frame(parameters[-remove, ])
      
      # Estimate extra LC model with obtained parameters as starting values
      output_path <- paste0(folder, "tree_MILC_", model_name, "_boot", i, "_dat_org_posteriors.dat")
      script <- generate_script_treeMILC_extra(ind, cov_ok, cov_problem, model_name, dat = dat,
                                               filepath_input = dat_org_path, filepath_output = output_path,
                                               parameters)
      script_path <- paste0(folder, "tree_MILC_", model_name, "_boot", i, "_dat_org_posteriors.lgs")
      writeLines(script, script_path)
      shell(paste0('"C:/Program Files/LatentGOLDnet6.0/lg60.exe" ', script_path, ' /b'))
      
      # Read model output with posterior probabilities for every observation in the original data set
      dat <- as.data.frame(fread(output_path, sep = "\t", dec = ","))
      setnames(dat, old = c("Cluster#1", "Cluster#2", "Cluster#3", "Cluster#"), new = c("p1", "p2", "p3", "cluster"))
      dat <- fix_cluster_assignment(type = "LC", list(1, dat))
      dat <- dat[[2]]
      dat <- dat[, -which(colnames(dat) == "cluster")] 
    } else {
      # If all combinations are present, add posterior probabilities to the observations in the original data set
      dat <- left_join(x = dat, y = count_boot, by = by_vector)
    }
    
    # Sample from obtained posterior membership probabilities
    dat$imp1 <- apply(dat, 1, impute_value_step1)
    dat$id <- as.character(dat$id)
    
    # Create subsets of cases with imputed values of 1 (with original indicator values) and write to file
    subset_ids <- dat[dat$imp1 == 1, ]$id
    subset <- dat[dat$id %in% subset_ids, ]
    filepath_subset <- paste0(folder, "tree_MILC", model_name, "_boot", i, "_subset.dat")
    fwrite(subset, file = filepath_subset, sep = "\t")
    
    # Create and execute Latent Gold script for second model
    filepath_output <- paste0(folder, boot_name, "_step2_output.dat")
    script <- generate_script("treeMILC", ind, cov_ok, cov_problem, N, model_name, subset, filepath_subset, filepath_output)
    script_path <- paste0(folder, boot_name, "_step2_script.lgs", sep = "")
    writeLines(script, script_path)
    shell(paste0('"C:/Program Files/LatentGOLDnet6.0/lg60.exe" ', script_path, ' /b'))
    
    # Import model output
    model2_output <- as.data.frame(fread(paste0(folder, boot_name, "_step2_output.dat", sep = ""), dec = ","))
    
    # Check if no warning was given
    model_lst <- paste(readLines(paste0(folder, boot_name, "_step2_script.lst")), collapse = "\n")
    if (grepl("WARNING", model_lst, fixed = TRUE, useBytes = TRUE)) {
      to_return <- append(to_return, list("Error"))
      to_return <- append(to_return, list("Error"))
      return(to_return)
    }
    
    # Add the results from the second model to the results from the first model
    count_step2 <- as.data.frame(model2_output[, -which(colnames(boot_output) %in% c("id", "Cluster#"))] %>%
                                   group_by_all() %>% summarise(COUNT = n()))
    count_step2 <- count_step2[, -ncol(count_step2)]
    dat <- left_join(x = dat, y = count_step2, by = c(all_indicators[1:ind], cov_ok, cov_problem))
    
    # Sample again from obtained posterior membership probabilities
    dat$imp2 <- apply(dat, 1, impute_value_step2)
    
    # Combine imputations from both models 
    dat$cluster <- dat$imp1
    dat[dat$imp1 == 1 & ((!is.na(dat$imp2) & dat$imp2 == 1)), ]$cluster <- 1
    dat[dat$imp1 == 1 & ((!is.na(dat$imp2) & dat$imp2 == 2)), ]$cluster <- 3
    
    # Select and rename columns
    new_names <- c("step2_p1", "step2_p3")
    setnames(dat, old = c("Cluster#1", "Cluster#2"), new = new_names)
    dat <- dat[, -which(colnames(dat) == "Cluster#")]
    
    # Fix cluster assignment (if necessary)
    dat <- fix_cluster_bootstrap("sim", dat)
    bootstrap_results <- append(bootstrap_results, list(dat))
  }
  
  # Add results to final output list
  to_return <- append(to_return, list(bootstrap_results))
  to_return <- append(to_return, list("Good"))
  return(to_return)
}
