################################# Methods_Helpfunctions.R ######################################## 
## This file contains the functions required to perform LC, LCT and tree-MILC analysis without  ##
## missing covariates (see Chapter 4) and with missing covariates with direct effects and       ##
## parameter restrictions (see Chapter 5). Note that in this file, the term 'clusters' is used  ## 
## instead of the term 'latent class'.                                                          ##
##    - create_subset                                                                           ## 
##    - generate_script                                                                         ##
##    - generate_script_treeMILC_extra                                                          ##
##    - fix_cluster_assignment                                                                  ##
##    - fix_cluster_bootstrap                                                                   ##
##    - fix_number_notation                                                                     ##         
##    - fix_extra_covcat                                                                        ##
##    - store_model_info                                                                        ##
##    - impute_value_step1                                                                      ##
##    - impute_value_step2                                                                      ##
##    - get_ME                                                                                  ##
##    - get_ME_help                                                                             ##
##################################################################################################

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
## @returns (string): A string containing a Latent Gold script                                  ##
##################################################################################################

generate_script <- function(type, ind, cov_ok, cov_problem, N, model_name, dat, filepath_input, filepath_output) {
  
  # Create vectors with characters to use for restrictions later on
  letters <- c("aa", "bb", "cc", "dd", "ee", "ff", "gg", "hh", "ii", "jj")
  letters2 <- c("kk", "ll", "mm", "nn", "oo", "pp", "qq", "rr", "ss", "tt")
  letters_count <- letters2_count <- 0 # Keep track of which letters have been used
  
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
        tolerance=1e-08 emtolerance=0.01 emiterations=10000000 nriterations=10000000;
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
  
  # Adjust the number of clusters depending on what type of analysis is performed
  if (type == "LC") {
    latent_var <- paste0("\n\tlatent Cluster nominal 3;
    equations\n")
  } else {
    latent_var <- paste0("\n\tlatent Cluster nominal 2;
    equations\n")
  }
  
  # For LCT step 2, use posterior probabilities from step 1 as weights 
  if (type == "LCT") {
    caseweight <- "caseweight p1;\n"
  } else {
    caseweight <- ""
  }
  
  # Adjust equations depending on the number of indicators
  if (ind == 2) {
    dep_ind <- "\tdependent Y1 nominal, Y2 nominal;"
    dep_ind_eq2 <- "\n\tY2 <- 1 | Cluster;"
  } else if (ind == 3) {
    dep_ind <- "\tdependent Y1 nominal, Y2 nominal, Y3 nominal;"
    dep_ind_eq2 <- "\n\tY2 <- 1 | Cluster;\n\tY3 <- 1 | Cluster;"
  } else if (ind == 4) {
    dep_ind <- "\tdependent Y1 nominal, Y2 nominal, Y3 nominal, Y4 nominal;"
    dep_ind_eq2 <- "\n\tY2 <- 1 | Cluster;\n\tY3 <- 1 | Cluster;\n\tY4 <- 1 | Cluster;"
  }
  
  # Adjust equations depending on which covariates to include (if any)
  cov <- c(cov_ok, cov_problem)
  restrictions1 <- restrictions2 <- ""
  dep_cov <- ifelse(is.null(cov), "", "\n\tindependent ")
  dep_ind_eq1 <- "\n\tY1 <- 1 | Cluster"
  
  if (!is.null(cov)) {
    for (i in 1:length(cov)) {
      if (i == 1) {
        dep_cov <- paste0(dep_cov, " ", cov[i], " nominal")
      } else {
        dep_cov <- paste0(dep_cov, ", ", cov[i], " nominal")
      }
    }
    dep_cov <- paste0(dep_cov, ";")
  }
  latent_var_eq <- "\tCluster <- 1"
  
  if (!is.null(cov_ok)) {
    for (i in 1:length(cov_ok)) {
      latent_var_eq <- paste(latent_var_eq, "+", cov_ok[i])
    }
  }
  
  # Specify restrictions
  if (!is.null(cov_problem)) {
    for (i in 1:length(cov_problem)) {
      letters_count <- letters_count + 1
      which_col <- which(colnames(dat) == cov_problem[i])
      num_cats <- length(levels(factor(dat[, which_col])))
      dep_ind_eq1 <- paste0(dep_ind_eq1, " + (", letters[i], "~ful) 1 | ", cov_problem[i])
      latent_var_eq <- paste0(latent_var_eq, " + (", letters2[i], ") ", cov_problem[i])
      
      # Specify the first restriction (see Section 5.1)
      for (j in 1:num_cats) {
        if (j != num_cats) {
          restrictions2 <- paste0(restrictions2, "\n\t", letters[i], "[", j, ",] = 0;")
        } else {
          restrictions2 <- paste0(restrictions2, "\n\t", letters[i], "[", j, ",1] = -100;")
          restrictions2 <- paste0(restrictions2, "\n\t", letters[i], "[", j, ",2] = 0;")
          restrictions2 <- paste0(restrictions2, "\n\t", letters[i], "[", j, ",3] = -100;")
        }
      }   
      
      # Specify the second restriction (see Section 5.1)
      if (length(cov_problem) > 1 & (i > 1)) {
        par2 <- (num_cats - 1) * 2
        par1 <- par2 - 1
        if (type == "LC") {
          restrictions1 <- paste0(restrictions1, paste0("\n\t", letters2[i], "[1,", par1, "] = 0; ", letters2[i], "[1,", par2, "] = 0;"))
        } else {
          restrictions1 <- paste0(restrictions1, paste0("\n\t", letters2[i], "[1,", num_cats - 1, "] = 0;"))
        }
      }
    }
  }
  
  latent_var_eq <- paste0(latent_var_eq, ";")
  dep_ind_eq1 <- paste0(dep_ind_eq1, ";")
  
  # Combine all parts of the script
  script <- paste0(script_part1, model_name, script_part2, script_part3, caseweight, dep_ind, dep_cov,
                   latent_var, latent_var_eq, dep_ind_eq1, dep_ind_eq2, restrictions1, restrictions2, "\nend model")
  
  return(script)
}

##################################################################################################
## Generate a Latent Gold script to estimate an extra model in tree-MILC analysis to obtain     ##
## posterior probabilities for response patterns that were not present in the                   ##
## bootstrap sample. This is done by estimating an LC model on the original data set while      ##
## using the logit parameters from the previous model as starting values and setting the number ## 
## of EM- and NR-iterations to 0.                                                               ## 
## @param ind (int): Number of indicators                                                       ##
## @param cov_ok (vector): Vector with the names of non-missing covariates                      ##
## @param cov_problem (vector): Vector with the names of missing covariates                     ##
## @param N (int): Size of data set                                                             ## 
## @param model_name (string): Name of the model                                                ##
## @param filepath_input (string): Path to input data set                                       ##
## @param filepath_output (string): Path where model output should be stored                    ##
## @param par (data.frame): Data frame with logit parameters as estimated by the previous model ## 
## @returns (string): A string containing a Latent Gold script                                  ##
##################################################################################################

generate_script_treeMILC_extra <- function(ind, cov_ok, cov_problem, model_name, dat, filepath_input, filepath_output, par) {
  
  # Create vectors with characters to use for restrictions later on
  pars <- c("aa", "bb", "cc", "dd", "ee", "ff", "gg", "hh", "ii", "jj")
  pars_count <- 0
  
  script_part1 <- paste0("version = 6.0\ninfile '", filepath_input, "' \n\nmodel title '")
  script_part2 <- paste0("';
    options
    algorithm
        tolerance=1e-08 emtolerance=0.01 emiterations=0 nriterations=0;
    startvalues
        seed=1 sets=100 tolerance=1e-05 iterations=100;
    bayes
        categorical=1 variances=1 latent=1 poisson=1;
    missing includeall;
    output
    	parameters=")
  script_part3 <- paste0("first standarderrors profile reorderclasses iterationdetails;
    	outfile '", filepath_output, "' classification keep=id;
    variables\n")
  
  latent_var <- paste0("\n\tlatent Cluster nominal 3;
    equations\n")

  # Adjust equations depending on which covariates to include (if any)
  cov <- c(cov_ok, cov_problem)
  restrictions <- restrictions2 <- ""
  latent_var_eq <- "\tCluster <- (zz) 1"
  dep_cov <- ifelse(is.null(cov), "", "\n\tindependent ")
  dep_ind_eq1 <- "\n\tY1 <- (jj) 1 | Cluster"
  
  #Specify restrictions
  if (!is.null(cov_problem)) {
    for (i in 1:length(cov_problem)) {
      pars_count <- pars_count + 1
      which_col <- which(colnames(dat) == cov_problem[i])
      num_cats <- length(levels(factor(dat[, which_col])))
      dep_ind_eq1 <- paste0(dep_ind_eq1, " + (", pars[i], "~ful) 1 | ", cov_problem[i]) #Specify second restriction (see Section 5.1)
      
      for (j in 1:num_cats) {
        if (j != num_cats) {
          restrictions2 <- paste0(restrictions2, "\n\t", pars[i], "[", j, ",] = 0;") #Specify first restriction (see Section 5.1)
        } else {
          restrictions2 <- paste0(restrictions2, "\n\t", pars[i], "[", j, ",1] = -100;")
          restrictions2 <- paste0(restrictions2, "\n\t", pars[i], "[", j, ",2] = 0;")
          restrictions2 <- paste0(restrictions2, "\n\t", pars[i], "[", j, ",3] = -100;")
        }
      }   
    }
  }
  
  if (!is.null(cov)) {
    for (i in 1:length(cov)) {
      if (i == 1) {
        dep_cov <- paste0(dep_cov, " ", cov[i], " nominal")
      } else {
        dep_cov <- paste0(dep_cov, ", ", cov[i], " nominal")
      }
    }
    dep_cov <- paste0(dep_cov, ";")
  }  
  
  #Specify start values
  restrictions <- paste0(restrictions, "\n\tzz[1,1] ~= ", par[1, ]$coef, ";\n")
  restrictions <- paste0(restrictions, "\tzz[1,2] ~= ", par[2, ]$coef, ";\n")
  df_counter <- 2
  if (!is.null(cov)) {
    for (i in 1:length(cov)) {
      which_col <- which(colnames(dat) == cov[i])
      num_cats <- length(levels(factor(dat[, which_col])))
      pars_count <- pars_count + 1
      latent_var_eq <- paste0(latent_var_eq, " + (", pars[pars_count], ") ", cov[i])
      par2 <- (num_cats - 1) * 2
      for (j in 1:par2) {
        df_counter <- df_counter + 1
        if (par[df_counter, ]$coef == 0) {
          restrictions <- paste0(restrictions, "\t", pars[pars_count], "[1,", j, "] = ", par[df_counter, ]$coef, ";\n")
        } else {
          restrictions <- paste0(restrictions, "\t", pars[pars_count], "[1,", j, "] ~= ", par[df_counter, ]$coef, ";\n")
        }
      }
    }
  }
  
  # Adjust equations depending on the number of indicators
  if (ind == 2) {
    dep_ind <- "\tdependent Y1 nominal, Y2 nominal;"
    dep_ind_eq2 <- "\n\tY2 <- (kk) 1 | Cluster;"
  } else if (ind == 3) {
    dep_ind <- "\tdependent Y1 nominal, Y2 nominal, Y3 nominal;"
    dep_ind_eq2 <- "\n\tY2 <- (kk) 1 | Cluster;\n\tY3 <- (ll) 1 | Cluster;"
  } else if (ind == 4) {
    dep_ind <- "\tdependent Y1 nominal, Y2 nominal, Y3 nominal, Y4 nominal;"
    dep_ind_eq2 <- "\n\tY2 <- (kk) 1 | Cluster;\n\tY3 <- (ll) 1 | Cluster;\n\tY4 <- (mm) 1 | Cluster;"
  }
  
  #Continue specifying start values
  all_indicators <- c("Y1", "Y2", "Y3", "Y4")
  indicators_parname <- c("jj", "kk", "ll", "mm")
  for (i in 1:ind) {
    df_counter <- grep(all_indicators[i], par$term1)[1]
    restrictions <- paste0(restrictions, "\t", indicators_parname[i], "[1,1] ~= ", par[df_counter, ]$coef, ";\n")
    restrictions <- paste0(restrictions, "\t", indicators_parname[i], "[1,2] ~= ", par[df_counter + 1, ]$coef, ";\n")
    restrictions <- paste0(restrictions, "\t", indicators_parname[i], "[2,1] ~= ", par[df_counter + 2, ]$coef, ";\n")
    restrictions <- paste0(restrictions, "\t", indicators_parname[i], "[2,2] ~= ", par[df_counter + 3, ]$coef, ";\n")
    restrictions <- paste0(restrictions, "\t", indicators_parname[i], "[3,1] ~= ", par[df_counter + 4, ]$coef, ";\n")
    restrictions <- paste0(restrictions, "\t", indicators_parname[i], "[3,2] ~= ", par[df_counter + 5, ]$coef, ";\n")
  }
  
  latent_var_eq <- paste0(latent_var_eq, ";")
  dep_ind_eq1 <- paste0(dep_ind_eq1, ";")
  
  # Combine all parts of the script
  script <- paste0(script_part1, model_name, script_part2, script_part3, dep_ind, dep_cov,
                   latent_var, latent_var_eq, dep_ind_eq1, dep_ind_eq2, restrictions, restrictions2, "\nend model")
  return(script)
}

##################################################################################################
## Fix cluster assignments for LC and LCT (i.e. assign the right names to the right clusters)   ##                                                                       ##
## @param type (string): Type of model (only relevant is "LC" for LC)                           ##
## @param results (list): Results of one particular LC, LCT or tree-MILC model                  ##
## @returns (list): Same as input, but with corrected cluster assignments (if necessary)        ##
##################################################################################################

fix_cluster_assignment <- function(type = NULL, results) {
  
  # Get list of measurement error matrices (i.e. one per indicator)
  ME_list <- get_ME(results)
  
  # Create one matrix that contains the average values of all matrices
  summed_matrices <- ME_list[[1]]
  for (j in 2:length(ME_list)) {
    summed_matrices <- summed_matrices + ME_list[[j]]
  }
  mean_matrix <- summed_matrices / length(ME_list)
  
  # Find all possible diagonal combinations and store each sum in a data frame
  all_diagonals <- data.frame()
  
  if (type == "LC") { 
    for (i in 1:3) {
      for (j in 1:3) {
        for (k in 1:3) {
          if (length(unique(c(i, j, k))) == 3) {
            all_diagonals <- rbind(all_diagonals, c(i, j, k, mean_matrix[i, 1], mean_matrix[j, 2], mean_matrix[k, 3],
                                                    sum(mean_matrix[i, 1], mean_matrix[j, 2], mean_matrix[k, 3])))
          }
        }
      }
    }
  } else {
    # For non-LC models: only find combinations for 1 and 3, because 2 is already correct
    all_diagonals <- rbind(all_diagonals, c(1, 2, 3, mean_matrix[1, 1], mean_matrix[2, 2], mean_matrix[3, 3], sum(mean_matrix[1, 1], mean_matrix[2, 2], mean_matrix[3, 3])))
    all_diagonals <- rbind(all_diagonals, c(3, 2, 1, mean_matrix[3, 1], mean_matrix[2, 2], mean_matrix[1, 3], sum(mean_matrix[3, 1], mean_matrix[2, 2], mean_matrix[1, 3])))
  }
  colnames(all_diagonals) <- c("1", "2", "3", "d1", "d2", "d3", "sum")
  
  # Find out which combination of diagonal values yields the highest sum 
  which_max <- which.max(as.vector(all_diagonals[, 7]))
  
  # Find out how to reassign clusters
  max_1 <- all_diagonals[which_max, 1]
  max_2 <- all_diagonals[which_max, 2]
  max_3 <- all_diagonals[which_max, 3]
  reassignment <- c(as.numeric(max_1), as.numeric(max_2), as.numeric(max_3))
  
  # If clusters need not to be reassigned, return original results
  if (identical(reassignment, c(1, 2, 3))) {
    return(results)
  } else {
    # Create a list to store the corrected results in (i.e. in the same format as the input)
    to_return <- list(results[[1]])
    results <- results[[2]]
    
    # Assign posterior probabilities to the right clusters
    posteriors <- data.frame(p1 = results$p1, p2 = results$p2, p3 = results$p3)
    results$p1 <- posteriors[, max_1]
    results$p2 <- posteriors[, max_2]
    results$p3 <- posteriors[, max_3]
    to_return <- append(to_return, list(results))
    return(to_return)
  }
}

##################################################################################################
## Fix cluster assignments per bootstrap sample in tree-MILC                                     ## 
## @param boot_results (data.frame): Results of one bootstrap sample                            ##
## @returns (data.frame): Same as input, but with corrected cluster assignments if necessary    ##
##################################################################################################

fix_cluster_bootstrap <- function(boot_results) {
  
  # Compute ME probability matrix (i.e. compute a contingency table for the values in indicator Y2 and the imputed values and normalise by rows)
  table <- table(boot_results$Y2, boot_results$cluster)
  prop_table <- prop.table(table, margin = 1)
  
  # Find out which combination of diagonals yields the highest sum of diagonal values (assume that 2 is already correct)
  all_diagonals <- data.frame()
  all_diagonals <- rbind(all_diagonals, c(1, 2, 3, prop_table[1, 1], prop_table[2, 2], prop_table[3, 3], sum(prop_table[1, 1], prop_table[2, 2], prop_table[3, 3])))
  all_diagonals <- rbind(all_diagonals, c(3, 2, 1, prop_table[3, 1], prop_table[2, 2], prop_table[1, 3], sum(prop_table[3, 1], prop_table[2, 2], prop_table[1, 3])))
  colnames(all_diagonals) <- c("1", "2", "3", "d1", "d2", "d3", "sum")
  which_max <- which.max(as.vector(all_diagonals[, 7]))
  
  # Find out how to reassign cluster names
  max_1 <- all_diagonals[which_max, 1]
  max_3 <- all_diagonals[which_max, 3]
  reassignment <- c(max_1, max_3)
  
  # If clusters need not to be reassigned, return original results
  if (identical(reassignment, c(1, 3))) {
    return(boot_results)
  } else { # Reassign cluster names
    boot_results$new_cluster <- NA
    boot_results[boot_results$cluster == max_1, ]$new_cluster <- 1
    boot_results[boot_results$cluster == 2, ]$new_cluster <- 2
    boot_results[boot_results$cluster == max_3, ]$new_cluster <- 3
    boot_results$cluster <- boot_results$new_cluster
    boot_results <- boot_results[, -which(colnames(boot_results) == "new_cluster")]
  }
  
  return(boot_results)
}

##################################################################################################
## Fix number notation in Latent GOLD output (i.e. convert '1e-02' to '.001')                   ##
## @param vector (vector): Input vector (potentially) containing numbers in wrong format        ## 
## @returns (vector): Output vector containing corrected numbers                                ##
##################################################################################################

fix_number_notation = function(vector) {
  
  # If the vector is numeric, all values are already in the correct format
  if (is.numeric(vector)) {
    return(vector)
  } else {
    # Create vector to store corrected results in
    return_vec <- rep(NA, length(vector))
    
    for (i in 1:length(vector)) {
      # Remove spaces from the string and split string by 'e-'
      removed_spaces <- gsub(" ", "", vector[i])  
      split_vec <- str_split(removed_spaces, "e-") 
      
      # If the string did contain 'e-'
      if (length(split_vec[[1]]) > 1) {
        # Convert number to the right format and add to corrected results vector
        nominator <- as.numeric(gsub(",", ".", split_vec[[1]][1]))  # Replace comma with dot
        denominator <- as.numeric(split_vec[[1]][2])  
        final_number <- nominator / (10^denominator)  
        return_vec[i] <- final_number
      } else {
        temp_string <- gsub(",", ".", split_vec[[1]][1])  # Replace comma with dot
        
        # In some cases, the output was '.', which meant that the probability was 0
        if (temp_string == ".") {
          return_vec[i] <- 0  
        } else { 
          return_vec[i] <- as.numeric(gsub(",", ".", split_vec[[1]][1]))  # Replace comma with dot
        }
      }
    }
    
    return(return_vec)
  }
}

##################################################################################################
## Assign observations with a missing covariate to an existing covariate category               ##
## @param data (data.frame): Data frame that contains observations with a missing covariate     ## 
## @param cov_index (int): Integer denoting the column number of a missing covariate            ##
## @returns (data.frame): Data frame with observations assigned to existing covariate category  ## 
##################################################################################################

fix_extra_covcat <- function(data, cov_index) {
  # Find extra category
  ncat <- length(levels(factor(data[, cov_index])))
  
  # Identify the largest group
  maxcat <- which.max(summary(factor(data[, cov_index])))
  
  # Assign the largest group as the new category for the observations with the extra category
  data[data[, cov_index] == ncat, cov_index] <- maxcat
  
  return(data)
}

##################################################################################################
## Get a list that consists of a data frame with model information                              ##
## @param iteration (int): Iteration number                                                     ## 
## @param ind (int): Number of indicators                                                       ## 
## @param cov_ok (vector): Vector with the names of non-missing covariates                      ##
## @param cov_problem (vector): Vector with the names of missing covariates                     ##
## @param N (int): Size of data set                                                             ## 
## @param ME (int): Amount of measurement error (1=10%, 2=20%, 3=30%, 4=realistic)              ##
## @returns (list): A list that consists of:                                                    ##
##  [[1]]: A data frame with model information                                                  ##
##################################################################################################

store_model_info <- function(iteration, ind, cov_ok, cov_problem, N, ME) {
  if (is.null(cov_ok)) {
    cov_ok1 <- "null"
  } else {
    cov_ok1 <- paste(cov_ok, collapse = "-")
  }
  
  if (is.null(cov_problem)) {
    cov_problem1 <- "null"
  } else {
    cov_problem1 <- paste(cov_problem, collapse = "-")
  }
  
  cov <- c(cov_ok1, cov_problem1)
  model_name <- paste(iteration, ind, cov_ok1, cov_problem1, N, ME, sep = "-")
  model_info1 <- data.frame(iteration = iteration, ind = ind, cov_ok = cov_ok1, cov_problem = cov_problem1)
  model_info2 <- data.frame(N = N, ME = ME, id = model_name)
  model_info <- cbind(model_info1, model_info2)
  model_info <- list(model_info)
  
  return(model_info)
}

##################################################################################################
## Obtain first round of imputations in tree-MILC analysis                                      ##
## @param x (vector): Vector that contains the posterior probabilities p1, p2, and p3           ##
## @returns (int): Imputed value                                                                ##
##################################################################################################

impute_value_step1 = function(x) {
  return(which(rmultinom(1, 1, c(as.numeric(x["p1"]) + as.numeric(x["p3"]), as.numeric(x["p2"]))) == 1))
}

##################################################################################################
## Obtain second round of imputations in tree-MILC analysis                                     ##
## @param x (vector): Vector that contains the posterior probabilities p1 and p2                ##
## @returns (int): Imputed value                                                                ##
##################################################################################################

impute_value_step2 = function(x) {
  if (is.na(x["Cluster#1"])) {
    return(NA)
  } else {
    return(which(rmultinom(1, 1, c(as.numeric(x["Cluster#1"]), as.numeric(x["Cluster#2"]))) == 1))
  }
}

##################################################################################################
## Compute measurement error matrices (one per indicator) for an LC, LCT or tree-MILC model     ##
## @param results (list): Results of an LC, LCT or tree-MILC model                              ##
## @returns (list): List of measurement error matrices (one per indicator)                      ##
##################################################################################################

get_ME = function(results){
  
  #For LC and LCT: Compute a ME probability matrix for each indicator using the get_ME_help function
  if (class(results[[2]]) != "list") {
    results = results[[2]] #Ignore first data frame with model information
    
    which_ind = results[, colnames(results) %in% c("Y1", "Y2", "Y3", "Y4")]
    to_return = list()
    for (i in 1:ncol(which_ind)) {
      to_return = append(to_return, list(get_ME_help(results, which_ind[, i])))
    }
    return(to_return)
  }
  
  #For tree-MILC: Compute ME matrix for each indicator based on the imputed values
  else {
    
    all_indicators = c("Y1", "Y2", "Y3", "Y4")
    num_of_ind = results[[1]]$ind
    results_2 = results[[2]]
    to_return = list() 
    
    for (i in 1:num_of_ind) {
      cluster_index = which(names(results_2[[1]]) == "cluster")
      ind_index = which(colnames(results_2[[1]]) == all_indicators[i])
      
      #Compute averages of the ME probability matrices for each bootstrap sample (per indicator) 
      summed_matrix = prop.table(table(results_2[[1]][, cluster_index], results_2[[1]][, ind_index]), 1)
      for (j in 2:length(results_2)) {
        cluster_index = which(names(results_2[[j]]) == "cluster")
        ind_index = which(colnames(results_2[[j]]) == all_indicators[i])
        summed_matrix = summed_matrix + prop.table(table(results_2[[j]][, cluster_index], results_2[[j]][, ind_index]), 1)
      }
      mean_matrix = summed_matrix / length(results_2)
      to_return = append(to_return, list(mean_matrix))
    }
    
    return(to_return)
  }
}

##################################################################################################
## Helpfunction to compute a measurement error matrix for one indicator                         ##
## @param results (data.frame): Data frame with posterior probabilities from one LC or LCT model##  
## @param ind_vec (vector): Vector containing posterior probabilities for one indicator         ##
## @returns (matrix): Measurement error matrix for one indicator                                ##
##################################################################################################

get_ME_help = function(results, ind_vec) {
  
  results$yx = ind_vec
  indicator_matrix = data.frame()
  
  for (i in 1:3) {
    #Compute the sum of p1, p2, and p3 for each indicator value
    indicator_matrix = rbind(indicator_matrix, sum(results[results$yx == i, ]$p1))
    indicator_matrix = rbind(indicator_matrix, sum(results[results$yx == i, ]$p2))
    indicator_matrix = rbind(indicator_matrix, sum(results[results$yx == i, ]$p3))
  }
  
  #Create a matrix and normalise the rows
  indicator_matrix = matrix(as.vector(indicator_matrix[, 1]), nrow = 3, ncol = 3, byrow = FALSE)
  indicator_matrix = prop.table(indicator_matrix, margin = 1)
  
  return(indicator_matrix)
}
