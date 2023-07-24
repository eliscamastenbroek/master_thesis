############################ Methods_Initial_Analysis_Approach_1.R ###############################
## This file contains the functions required to perform the initial analysis in Section 5.1.1.  ##
## This analysis involves performing LCT and tree-MILC analysis using the first approach to     ## 
## include missing covariates in Chapter 5. This file contains the following functions:         ##
##    - generate_script                                                                         ##
##    - perform_lc                                                                              ##
##    - perform_lct                                                                             ##
##    - perform_treeMILC                                                                        ##
##################################################################################################

##################################################################################################
## Generate a Latent Gold script for to estimate an LC model. This function can be used to      ##  
## estimate regular LC models, as well as LC models in LCT and tree-MILC analysis.              ##                         
## @param type (string): What type of model (e.g. "LC" for regular LC and "LCT" for LCT step 2) ## 
## @param ind (int): Number of indicators                                                       ##
## @param cov_ok (vector): Vector with the names of non-missing covariates                      ##
## @param cov_problem (vector): Vector with the names of missing covariates                     ##
## @param N (int): Size of data set                                                             ##
## @param model_name (string): Name of the model                                                ##
## @param filepath_input (string): Path to input data set                                       ##
## @param filepath_output (string): Path where model output should be stored                    ##
## @param start_val (string): Starting values to obtain posterior probabilities in tree-MILC    ##
## @returns (string): A string containing a Latent Gold script                                  ##
##################################################################################################

generate_script <- function(type, ind, cov_ok, cov_problem, N, model_name, filepath_input, filepath_output, start_val = NULL) {
  
  script_part1 <- paste0("version = 6.0\ninfile '", filepath_input, "' \n\nmodel title '")
  
  # Let the number of sets of starting values depend on the size of the data set
  if (N < 10000) {
    sets <- 3200
  } else {
    sets <- 100
  }
  
  script_part2 <- paste0("';
    options
    algorithm
        tolerance=1e-08 emtolerance=0.01 emiterations=1000 nriterations=1000;
    startvalues
        seed=1 sets=", sets, " tolerance=1e-05 iterations=100;
    bayes
        categorical=1 variances=1 latent=1 poisson=1;
    missing includeall;
    output
    	parameters=")
  
  script_part3 <- paste0("first standarderrors profile reorderclasses iterationdetails;
    	outfile '", filepath_output, "' classification keep=id;
    variables\n") 
  
  # Adjust number of clusters depending on what type of analysis is performed
  if (type == "LC") {
    latent_var <- paste0("\n\tlatent Cluster nominal 3;
    equations\n")
  } else {
    latent_var <- paste0("\n\tlatent Cluster nominal 2;
    equations\n")
  }
  
  # For LCT step 2, use posterior probabilities from step 1 as weights 
  if (type == "LCT2") {
    caseweight <- "caseweight p1;\n"
  } else {
    caseweight <- ""
  }
  
  # Adjust equations depending on the number of indicators
  if (ind == 2) {
    dep_ind <- "\tdependent Y1 nominal, Y2 nominal;"
    dep_ind_eq <- "\tY1 <- 1 + Cluster;\n\tY2 <- 1 + Cluster;"
  } else if (ind == 3) {
    dep_ind <- "\tdependent Y1 nominal, Y2 nominal, Y3 nominal;"
    dep_ind_eq <- "\tY1 <- 1 + Cluster;\n\tY2 <- 1 + Cluster;\n\tY3 <- 1 + Cluster;"
  } else if (ind == 4) {
    dep_ind <- "\tdependent Y1 nominal, Y2 nominal, Y3 nominal, Y4 nominal;"
    dep_ind_eq <- "\tY1 <- 1 + Cluster;\n\tY2 <- 1 + Cluster;\n\tY3 <- 1 + Cluster;\n\tY4 <- 1 + Cluster;"
  }
  
  # Adjust equations depending on which covariates to include (if any)
  cov <- c(cov_ok,cov_problem)
  
  if (is.null(cov)) {
    dep_cov <- ""
    latent_var_eq <- "\tCluster <- 1;\n"
  } else {
    dep_cov <- "\n\tindependent"
    latent_var_eq <- "\tCluster <- 1"
    for (i in 1:length(cov)) {
      if(i == length(cov)){
        dep_cov <- paste(dep_cov, cov[i], "nominal;")
        latent_var_eq <- paste0(latent_var_eq, " + ", cov[i], ";\n")
      } else {
        dep_cov <- paste(dep_cov, cov[i], "nominal,")
        latent_var_eq <- paste0(latent_var_eq, " + ", cov[i])
      }
    }
  }
  
  if(!is.null(start_val)){
    start_val = paste("\n{", start_val, "}")
  }
  
  # Combine all parts of the script
  script <- paste0(script_part1, model_name, script_part2, script_part3, caseweight, dep_ind, dep_cov,
                   latent_var, latent_var_eq, dep_ind_eq, start_val, "\nend model")
  return(script)
}

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

perform_lct <- function(iteration, ind, cov_ok = NULL, cov_problem = NULL, N, ME, folder) {
  
  # Store information about the model in the output list 
  to_return <- store_model_info(iteration, ind, cov_ok, cov_problem, N, ME)
  model_name <- to_return[[1]]$id

  # Write data set to use to file
  dat <- create_subset(iteration, ind, cov_ok, cov_problem, N, ME)
  dat_name <- paste0("LCT_", model_name, "_step1_data.dat")
  write.table(x = dat, file = paste0(folder, dat_name), row.names = FALSE, quote = FALSE)
  
  # Perform LC with 3 classes without problem covariate(s). 
  model1_output <- perform_lc(iteration, ind, cov_ok, cov_problem = NULL, N, ME, folder = folder)[[2]]
  
  # Combine posterior probabilities for the classes 'permanent' and 'flexible' in the model output
  model1_output$p1 <- 1 - model1_output$p2
  
  # The output of model 1 does not contain the problem covariates, so add them here:
  cov_problem_indexes <- which(names(dat) == cov_problem)
  for (i in cov_problem_indexes) {
    model1_output[, names(dat)[i]] <- dat[, i]
  }
  
  # Write data set to use in the second step to file
  dat_name <- paste0("LCT_", model_name, "_step1_output.dat")
  write.table(x = model1_output, file = paste0(folder, dat_name), row.names = FALSE, quote = FALSE)
  
  # Create Latent Gold script for the second LC model with problem covariate(s)
  filepath_input <- paste0(folder, dat_name)
  filepath_output <- paste0(folder, "LCT_", model_name, "_step2_output.dat")
  script <- generate_script("LCT2", ind, cov_ok = cov_ok,cov_problem, N, model_name, filepath_input, filepath_output)
  script_path <- paste0(folder, "LCT_", model_name, "_step2_script.lgs")
  writeLines(script, script_path)
  
  # Execute script in Latent Gold and read model output
  shell(paste0('"C:/Program Files/LatentGOLDnet6.0/lg60.exe" ', script_path, ' /b'))
  model2_output <- read.delim(filepath_output, sep = "\t", dec = ",")

  # Combine results from both models to compute final posterior probabilities and fix column names
  all_indicators <- c("Y1", "Y2", "Y3", "Y4")
  by_vector <- c("id", all_indicators[1:ind], c(cov_ok, cov_problem))
  combined_output <- left_join(x = model1_output, y = model2_output, by = by_vector)
  new_names <- c("p1.1", "p2.1", "p2.2", "cluster.2", "cluster.1")
  setnames(combined_output, old = c("p1.y", "Cluster.1", "Cluster.2", "Cluster.", "cluster"), new = new_names)
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
    boot_output <- perform_lc(iteration, ind, cov_ok, cov_problem = NULL, N, ME, dat = sample, folder = folder)
    model_name2 <- boot_output[[1]]$id
    
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
      
      # Get starting values from previous model and put them in a string
      filename_lc <- paste0(folder, "LC_", model_name2, "_script.lst")
      lc_lst <- readChar(filename_lc, file.info(filename_lc)$size)
      lc_lst <- strsplit(lc_lst, split = "Regression Parameters")
      lc_lst <- strsplit(lc_lst[[1]][2], split = "Paired Comparisons")
      parameters_path <- paste0(folder, "tree_MILC_", model_name2, "_boot", i, "_parameters.dat")
      writeLines(lc_lst[[1]][1], parameters_path)
      parameters <- as.data.frame(suppressWarnings(fread(parameters_path, sep = "\t", dec = ",")))[1:5]
      names(parameters) <- c("term1", "term2", "term3", "coef", "SE")
      parameters$coef <- as.numeric(parameters$coef)
      parameters <- as.data.frame(parameters[parameters$SE!=".",])
      start_val <- paste(parameters$coef[!is.na(parameters$coef)], collapse = " ")
      
      # Estimate extra LC model with obtained parameters as starting values
      output_path <- paste0(folder, "tree_MILC_", model_name2, "_boot", i, "_dat_org_posteriors.dat")
      script <- generate_script("LC",ind, cov_ok, cov_problem = NULL, N, model_name2,
                                filepath_input = dat_org_path, filepath_output = output_path,
                                start_val = start_val)
      script_path <- paste0(folder, "tree_MILC_", model_name2, "_boot", i, "_dat_org_posteriors.lgs")
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
    
    # The output of model 1 does not contain the problem covariates, so add them here:
    cov_problem_indexes <- which(names(dat_org) == cov_problem)
    for (i in cov_problem_indexes) {
      dat[, names(dat_org)[i]] <- dat_org[, i]
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
    script <- generate_script("treeMILC", ind = ind, cov_ok = cov_ok, cov_problem = cov_problem, N = N, model_name = model_name, filepath_input = filepath_subset, filepath_output = filepath_output)
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
    count_step2 <- as.data.frame(model2_output[, -which(colnames(model2_output) %in% c("id", "Cluster#"))] %>%
                                   group_by_all() %>% summarise(COUNT = n()))
    count_step2 <- count_step2[, -ncol(count_step2)]

    print(head(dat))
    print(head(count_step2))
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

    # Fix cluster assignment (if necessary)
    dat <- fix_cluster_bootstrap("sim", dat)
    bootstrap_results <- append(bootstrap_results, list(dat))
  }
  
  # Add results to final output list
  to_return <- append(to_return, list(bootstrap_results))
  to_return <- append(to_return, list("Good"))
  return(to_return)
}
