########################## Performance_Measures_and_Plots.R ###################################### 
## This file contains the functions needed to compute the performance measures and creating the ## 
## plots in Chapters 4 and 5:                                                                   ##
##    - create_ME_matrix                                                                        ##                                                                          
##    - get_entropy                                                                             ## 
##    - get_proportions                                                                         ##
##    - get_mean_summed_bias                                                                    ##
##    - get_diff_matrix                                                                         ##
##    - get_mean_matrix                                                                         ##
##    - get_ME_per_contract                                                                     ## 
##    - get_var_ME                                                                              ##
##################################################################################################

##################################################################################################
## Create a matrix with the measurement error probabilities used to simulate data               ##
## @param a2,...b32 (int): Logit parameters used to simulate the data                           ##
## @returns (matrix): Matrix with measurement error probabilities                               ##
##################################################################################################

create_ME_matrix <- function(a2, a3, b22, b33, b23, b32) {
  
  # Compute rows separately
  row1 <- c(1 / (1 + exp(a2) + exp(a3)), exp(a2) / (1 + exp(a2) + exp(a3)), exp(a3) / (1 + exp(a2) + exp(a3)))
  row2 <- c(1 / (1 + exp(a2 + b22) + exp(a3 + b32)), exp(a2 + b22) / (1 + exp(a2 + b22) + exp(a3 + b32)), exp(a3 + b32) / (1 + exp(a2 + b22) + exp(a3 + b32)))
  row3 <- c(1 / (1 + exp(a2 + b23) + exp(a3 + b33)), exp(a2 + b23) / (1 + exp(a2 + b23) + exp(a3 + b33)), exp(a3 + b33) / (1 + exp(a2 + b23) + exp(a3 + b33)))
  
  # Combine rows into a matrix
  matrix <- matrix(c(row1, row2, row3), nrow = 3, ncol = 3, byrow = TRUE)
  
  return(matrix)
}

##################################################################################################
## Calculate entropy and entropy R2                                                             ##
## @param results (list): Results of one particular LC or LCT model                             ##
## @returns (vector): A vector containing the entropy and the entropy R2                        ##
##################################################################################################

get_entropy <- function(results) {
  
  # Ignore first data frame with model information
  results <- results[[2]] 
  
  # For LC and LCT
  if (class(results) != "list") { 
    lc1 <- sum((results$p1) * log(results$p1))
    lc1 <- replace(lc1, lc1 == "NaN", 0) # In case some posterior probabilities are 0
    lc2 <- sum((results$p2) * log(results$p2))
    lc2 <- replace(lc2, lc2 == "NaN", 0)
    lc3 <- sum((results$p3) * log(results$p3))
    lc3 <- replace(lc3, lc3 == "NaN", 0)
    entropy <- -(lc1 + lc2 + lc3)
    entropy_squared <- 1 - (entropy / (nrow(results) * log(3)))
  }
  
  # For tree-MILC
  else { 
    entropy <- entropy_squared <- NA
  }
  
  entropy_vector <- c(entropy_squared, entropy)
  names(entropy_vector) <- c("Entropy R-squared", "Entropy")
  
  return(entropy_vector)
}

##################################################################################################
## Compute population proportion estimates                                                      ##
## @param results (list): Results of one particular LC, LCT or tree-MILC model                  ##
## @returns (vector): A vector containing the estimated proportions per cluster                 ##
##################################################################################################

get_proportions <- function(results) {
  results <- results[[2]]   # Ignore first data frame with model information
  
  # For LC or LCT: Compute proportions based on posterior probabilities
  if (class(results) != "list") {        
    prop1 <- sum(results$p1)
    prop2 <- sum(results$p2)
    prop3 <- sum(results$p3)
    proportions_vector <- c(prop1, prop2, prop3) / nrow(results)
    names(proportions_vector) <- 1:3
    return(proportions_vector)
  } 
  
  # For tree-MILC: Compute proportions based on imputations
  else {                    
    # Compute proportions per bootstrap
    prop_per_bootstrap <- list()
    for (i in 1:length(results)) {
      prop_per_bootstrap <- append(prop_per_bootstrap, list(summary(factor(results[[i]]$cluster)) / nrow(results[[i]])))
    }
    
    # Pool results
    pooled_proportions <- colMeans(bind_rows(prop_per_bootstrap))   
    return(pooled_proportions)
  }
}

##################################################################################################
## Get mean summed bias (see Section 4.4.2)                                                     ##
## @param results (list): Results of one particular LC, LCT or tree-MILC model                  ##
## @returns (int): Mean summed bias                                                             ##
##################################################################################################

get_mean_summed_bias <- function(results) {
  
  # Get estimated measurement error probability matrices
  ME_results <- get_ME(results)

  # Get true measurement error probability matrices
  ME <- results[[1]]$ME
  ME_matrix <- get(paste0("ME_matrix", ME))

  # If the amount of measurement error is realistic
  if (length(ME_matrix) == 2) {
    diff_sum <- 0 # Create integer to store the mean summed bias in
    for (j in 1:length(ME_results)) { # Iterate over the number of indicators
      # Define k to represent whether we need to use the ME probability matrix for Y1 and Y3 or the ME probability matrix for Y2 and Y4
      k <- j 
      if (j > 2) {
        k <- k - 2
      }
      
      # Compute the bias of each element on the diagonal
      for (i in 1:3) { 
        diff_sum <- diff_sum + (ME_results[[j]][i, i] - ME_matrix[[k]][i, i])
      }
    }
    
    return(diff_sum / length(ME_results)) # Return mean summed bias

  }
  
  # If the amount of measurement error is 10%, 20%, or 30%
  else {
    diff_sum <- 0 # Create integer to store the mean summed bias in
    for (j in 1:length(ME_results)) { # Iterate over the number of indicators 
      # Compute the bias of each element on the diagonal
      for (i in 1:3) { 
        diff_sum <- diff_sum + (ME_results[[j]][i, i] - ME_matrix[i, i])
      }
    }
    
    # Return mean summed bias
    return(diff_sum / length(ME_results))
  }
}

##################################################################################################
## Get for every indicator a matrix with the differences between the input ME matrix            ##
## and the estimated ME matrix                                                                  ##
## @param results (list): Results of one particular LC, LCT or tree-MILC model                  ##
## @returns (list): List of matrices with for every indicator the differences between the       ## 
##                  input and estimated ME matrix                                               ##
##################################################################################################

get_diff_matrix <- function(results) {
  
  # Get true measurement error probability matrix
  ME <- results[[1]]$ME
  ind <- results[[1]]$ind
  ME_results <- get_ME(results)
  ME_matrix <- get(paste0("ME_matrix", ME))
  
  # Create list to store the results
  to_return <- list()
  
  # If the amount of measurement error is 10%, 20%, or 30%
  if (ME != 4) {
    # Get for every indicator the differences between the input and the estimated ME matrix and add to return list
    for (i in 1:ind) {
      to_return <- append(to_return, list((ME_results[[i]] - ME_matrix))) 
    }
  }
  
  # If the amount of measurement error is realistic
  if (ME == 4) {
    # Get differences between input and estimated ME matrix for Y1 and Y2 and add to return list
    for (i in 1:2) {
      to_return <- append(to_return, list((ME_results[[i]] - ME_matrix[[i]])))  
    }
    # Get differences between input and estimated ME matrix for Y3 and add to return list 
    if (ind > 2) {
      to_return <- append(to_return, list((ME_results[[3]] - ME_matrix[[1]])))
    }
    # Get differences between input and estimated ME matrix for Y4 and add to return list
    if (ind > 3) {
      to_return <- append(to_return, list((ME_results[[4]] - ME_matrix[[2]])))
    }
  }
  
  return(to_return)
}

##################################################################################################
## Takes in a list of (estimated) ME probability matrices for models with the same parameters   ##  
## but different iterations. This function computes for each indicator the mean (estimated) ME  ##  
## probability matrix and returns these matrices in a list.                                     ##
## @param results (list): A list of a list of matrices                                          ##
## @returns (list): List of the mean ME probability matrix per indicator                        ##
##################################################################################################

get_mean_matrix <- function(list) {
  num_matrices <- length(list) # Number of lists
  num_ind <- length(list[[1]]) # Number of matrices per list
  to_return <- list() # List to store results in
  
  # Compute mean of the matrices per indicator
  for (i in 1:num_ind) {
    temp_matrix <- list[[1]][[i]]
    for (j in 2:num_matrices) {
      temp_matrix <- temp_matrix + list[[j]][[i]]
    }
    temp_matrix <- temp_matrix / num_matrices
    to_return <- append(to_return, list(temp_matrix))
  }
  return(to_return)
}

##################################################################################################
## Compute (mean) ME probability estimate for one particular contract type.                     ##
## @param results (model): Results of one particular LC, LCT or tree-MILC model                 ##
## @param contract (int): Contract type (1=permanent, 2=other, 3=flexible)                      ##
## @returns (int): Mean ME probability estimate for one particular contract type.               ##
##################################################################################################

get_ME_per_contract <- function(model, contract){

  # If the amount of measurement error is realistic
  ME <- model[[1]]$ME
  
  if(ME == 4){

    # Get (estimated) ME probability matrix
    ME_estimated <- get_ME(model)
    ind <- length(ME_estimated)
    
    first <- ME_estimated[seq(1, ind, 2)] # ME probability matrices for Y1 and Y3
    second <- ME_estimated[seq(2, ind, 2)] # ME probability matrices for Y2 and Y4
    
    # If Y1 and Y3, compute the mean ME probability estimate for both indicators
    if(length(first) > 1){
      x <- get_mean_matrix(first)
      first <- matrix(unlist(x), ncol = 3, byrow = TRUE)
      first <- first[contract, contract]
    # If only Y1, get the ME probability estimate
    } else {
      first <- first[[1]][contract, contract]
    }
    
    # If Y2 and Y4, compute the mean ME probability estimate for both indicators
    if(length(second) > 1){
      y <- get_mean_matrix(second)
      second <- matrix(unlist(y), ncol = 3, byrow = TRUE)
      second <- second[contract, contract]
    # If only Y2, get the ME probability estimate	
    } else {
      second = second[[1]][contract,contract]
    }

    return(c(first, second))
  # If the amount of measurement error is 10%, 20%, or 30%   
  } else {
    # Compute the mean of ME probability matrices for all indicators
    ME_estimated <- get_mean_matrix(get_ME(model))
    ME_estimated <- matrix(unlist(ME_estimated), ncol = 3, byrow = TRUE)
    return(ME_estimated[contract, contract])
  }
}
##################################################################################################
## Compute the variance of the ME probability estimates                                         ##
## -------------------------------------------------------------------------------------------- ##
## @param results (models): List of LC, LCT or tree-MILC models                                 ##
## @param ind (int): Number of indicators                                                       ##
## @param cov_ok (vector): Vector with the names of non-missing covariates                      ##
## @param cov_problem (vector): Vector with the names of missing covariates                     ##
## @param N (int): Size of data set                                                             ##
## @param ME (int): Amount of measurement error (1=10%, 2=20%, 3=30%, 4=realistic)              ##
## @param ind_type (int): In case of realistic amount of ME, which indicator type (1=Y1 and Y3, ##
## 2=Y2 and Y4)                                                                                 ##
## @returns (list): List of matrices with for every indicator the variance of the ME            ## 
## probability estimates                                                                        ##
##################################################################################################

get_var_ME <- function(models, ind, cov_ok, cov_problem, N, ME, ind_type = NULL) {
  
  # Define which indicators to include
  if (ME == 4) {
    select_inds <- seq(ind_type, ind, 2)  # Include only a subset of indicators in case of a realistic amount of ME
  } else {
    select_inds <- 1:ind
  }
  
  # Get list of ME probability matrices for all models with the same parameters
  ME_list <- list()
  for (i in 1:length(models)) {
    if ((models[[i]][[1]]$ind == ind) && 
        (models[[i]][[1]]$cov_ok == cov_ok) && 
        (models[[i]][[1]]$cov_problem == cov_problem) && 
        (models[[i]][[1]]$N == N) && 
        (models[[i]][[1]]$ME == ME)) {
      ME_list <- append(ME_list, list(get_ME(models[[i]])))
    }
  }
  
  # Create a list to store results in
  to_return <- list()
  
  # For every indicator
  for (i in select_inds) {
    
    # Create a vector to store variances in
    var_vec <- vector()
    
    # Iterate over every position in a ME probability matrix
    for (m in 1:3) {
      for (n in 1:3) {
        
        # For every position, store ME probability estimates of all models in a temporary vector
        vec <- vector()
        
        for (o in 1:length(ME_list)) {
          vec <- c(vec, ME_list[[o]][[i]][m, n])
        }
        
        # Compute the variance of all items in the temporary vector and add the result to the variance vector
        var_vec <- c(var_vec, var(vec))
      }
    }
    
    # Convert variance vector to a matrix and add to the return list
    temp_matrix <- matrix(var_vec, ncol = 3, byrow = TRUE)
    to_return <- append(to_return, list(temp_matrix))
  }
  
  return(to_return)
}

##################################################################################################
## Get data frame with the results of all LC, LCT or tree-MILC models                           ##
## -------------------------------------------------------------------------------------------- ##
## @param models (list): List of models                                                         ##
## @returns (data.frame): Data frame with the results for each model                            ##
##################################################################################################

get_results <- function(models) {
  
  # Create data frame to store results
  results_df <- data.frame(iteration = NA,
                           indicator = NA,
                           cov_ok = NA,
                           cov_problem = NA,
                           N = NA,
                           ME = NA,
                           id = NA,
                           full_id = NA,
                           prop1 = NA,
                           prop2 = NA,
                           prop3 = NA,
                           entropy = NA,
                           ME1.1 = NA,
                           ME1.2 = NA,
                           ME1.3 = NA,
                           ME2.1 = NA,
                           ME2.2 = NA,
                           ME2.3 = NA,
                           diag = NA)
  
  for (i in 1:length(models)) {
    print(i)
    model <- models[[i]]
    model1 <- model[[1]] # This object contains the model information
    model_name <- paste(as.numeric(model1$ind), model1$cov_ok, model1$cov_problem, as.numeric(model1$N), as.numeric(model1$ME), sep = "-")
    
    # For every contract type, get the estimated ME probability of observing correctly (i.e., element on the diagonal in the matrix)
    if (model1$ME == 4) { # In case of a realistic amount of ME, split for indicators (Y1 and Y3) and (Y2 and Y4)
      ME1.1 <- get_ME_per_contract(model, 1)[1]
      ME2.1 <- get_ME_per_contract(model, 1)[2]
      ME1.2 <- get_ME_per_contract(model, 2)[1]
      ME2.2 <- get_ME_per_contract(model, 2)[2]
      ME1.3 <- get_ME_per_contract(model, 3)[1]
      ME2.3 <- get_ME_per_contract(model, 3)[2]
    } else { 
      ME1.1 <- get_ME_per_contract(model, 1)
      ME1.2 <- get_ME_per_contract(model, 2)
      ME1.3 <- get_ME_per_contract(model, 3)
      ME2.1 <- ME2.2 <- ME2.3 <- NA
    }
    
    # Compute performance measures and add to data frame
    row <- c(as.numeric(model1$iteration),
             as.numeric(model1$ind),
             model1$cov_ok,
             model1$cov_problem,
             as.numeric(model1$N),
             as.numeric(model1$ME),
             model_name,
             full_id = model1$id,
             get_proportions(model)[1],
             get_proportions(model)[2],
             get_proportions(model)[3],
             ifelse(length(model[[2]]) != 5, get_entropy(model)[1], NA),
             ME1.1,
             ME1.2,
             ME1.3,
             ME2.1,
             ME2.2,
             ME2.3,
             get_mean_summed_bias(model))
    results_df[i, ] <- row
  }
  
  return(results_df)
}

##################################################################################################
## Get data frame with a summary of the results of the LC, LCT or tree-MILC models              ##
## -------------------------------------------------------------------------------------------- ##
## @param results (data.frame): Data frame with results                                         ##
## @returns (data.frame): Data frame with summary of results (e.g. averaged over iterations)    ##
##################################################################################################

get_summary <- function(type, model_results) {
  
  # Convert columns to numeric
  convert <- c("prop1", "prop2", "prop3", "entropy", "ME1.1", "ME1.2", "ME1.3", "ME2.1", "ME2.2", "ME2.3", "diag")
  model_results[, convert] <- apply(model_results[, convert], 2, function(x) as.numeric(as.character(x)))
  
  # Compute the mean of estimates and performance measures for models with the same parameters 
  summary1 <- as.data.frame(model_results %>% 
                              group_by(indicator, cov_ok, cov_problem, N, ME, id) %>%   
                              summarise_at(c("prop1", "prop2", "prop3", "entropy", "ME1.1", "ME1.2", "ME1.3", "ME2.1", "ME2.2", "ME2.3", "diag"), mean))
  
  # Compute RMSE
  summary1$rmse_prop1 <- summary1$rmse_prop2 <- summary1$rmse_prop3 <- NA
  split_df <- split(model_results, model_results$id) # Split the data frame into smaller lists that only contain models with the same parameters
  for (i in 1:length(split_df)) {
    ME <- split_df[[i]]$ME[1]
    nsim <- nrow(split_df[[i]])
    ME_matrix <- get(paste0("ME_matrix", ME))
    summary1[summary1$id == split_df[[i]]$id[1], ]$rmse_prop1 <- sqrt(sum((split_df[[i]]$prop1 - true_proportions[1])^2) / nsim)
    summary1[summary1$id == split_df[[i]]$id[1], ]$rmse_prop2 <- sqrt(sum((split_df[[i]]$prop2 - true_proportions[2])^2) / nsim)
    summary1[summary1$id == split_df[[i]]$id[1], ]$rmse_prop3 <- sqrt(sum((split_df[[i]]$prop3 - true_proportions[3])^2) / nsim)
  }
  
  # Compute the standard deviation of estimates and performance measures for models with the same parameters 
  summary2 <- as.data.frame(model_results %>% 
                              group_by(indicator, cov_ok, cov_problem, N, ME, id) %>% 
                              summarise_at(c("prop1", "prop2", "prop3", "entropy", "ME1.1", "ME1.2", "ME1.3", "ME2.1", "ME2.2", "ME2.3"), sd))[,-c(1:6)]
  names(summary2) <- c("sd_prop1", "sd_prop2", "sd_prop3", "sd_entropy", "sd_ME1.1", "sd_ME1.2", "sd_ME1.3", "sd_ME2.1", "sd_ME2.2", "sd_ME2.3")
  
  # Combine results and add model type
  summary <- cbind(summary1, summary2)
  summary$type <- type
  
  return(summary)
}

##################################################################################################
## Get RMSE of ME probability estimates (on the diagonals)                                      ##
## -------------------------------------------------------------------------------------------- ##
## @param models (list): List of models                                                         ##
## @param model_results (data.frame): Data frame with model results                             ##
## @returns (data.frame): Data frame with RMSE of ME probability estimates (on the diagonals)   ##
##################################################################################################

get_rmse_ME <- function(models, model_results) {
  
  # Create data frame to store results
  id <- unique(model_results$id)
  df <- data.frame(indicator = NA, cov_ok = NA, cov_problem = NA, N = NA, ME = NA, id = NA, rmse = NA, ind = NA, Contract = NA)
  
  # For every contract type, compute estimated ME probability using the get_rmse_ME_help function
  for (i in id) {
    model_sub <- models[as.numeric(rownames(model_results[model_results$id == i,]))]
    for (j in 1:num_ind) {
      for (k in 1:3) {
        row <- c(model_sub[[1]][[1]]$ind, model_sub[[1]][[1]]$cov_ok, model_sub[[1]][[1]]$cov_problem, model_sub[[1]][[1]]$N, model_sub[[1]][[1]]$ME, i, get_rmse_ME_help(model_sub, j, k), j, k)
        df <- rbind(df, row)
      }
    }
  }
  return(df[-1, ])
} 

##################################################################################################
## Get RMSE of estimated ME probability per indicator and per contract type                     ##
## -------------------------------------------------------------------------------------------- ##
## @param model_list (list): List of models                                                     ##
## @param indicator (int): Integer indicating which indicator to look at                        ##
## @param contract (int): Integer indicating which contract type to look at                     ##
## @returns (int): RMSE                                                                         ##
##################################################################################################

get_rmse_ME_help <- function(model_list, indicator, contract) {
  
  # Store estimated measurement error probabilities in a vector
  ME_vec <- vector()
  for (i in model_list) {
    ME_vec <- c(ME_vec, get_ME(i)[[indicator]][contract, contract])
  }
  
  # Get true ME probability matrix  
  ME_matrix <- get(paste0("ME_matrix", as.numeric(model_list[[1]][[1]]$ME)))
  
  # For a realistic amount of ME, get the correct true value (for Y1/Y3 or Y2/Y4)
  if (model_list[[1]][[1]]$ME == 4) { 
    if (indicator == 1 | indicator == 3) {
      true_value <- ME_matrix[[1]][contract, contract]
    } else {
      true_value <- ME_matrix[[2]][contract, contract]
    }
  # For 10%, 20%, and 30%, get the true value
  } else {
    true_value <- ME_matrix[contract, contract]
  }
  
  # Compute RMSE
  rmse <- sqrt(sum((ME_vec - true_value)^2) / length(model_list))
  return(rmse)
}

##################################################################################################
## Create a grid of ME heatmaps for LC, LCT and tree-MILC models with 10%, 20% and 30% ME       ##
## -------------------------------------------------------------------------------------------- ##
## @param LC_models (list): List of LC models                                                   ##
## @param LCT_models (list): List of LCT models                                                 ##
## @param treeMILC_models (list): List of treeMILC models                                       ##
## @param model_results (data.frame): Results of the LC, LCT or tree-MILC models                ##
## @param type (int): Integer denoting whether to show bias or variance (1=bias, 2=variance)    ##
## @param ind (int): Number of indicators                                                       ##
## @param cov_ok (vector): Vector with the names of non-missing covariates                      ##
## @param cov_problem (vector): Vector with the names of missing covariates                     ##
## @param N (int): Size of data set                                                             ##
## @returns (ggplot object): ME heatmap                                                         ##
##################################################################################################

get_ME_heatmap <- function(LC_models, LCT_models, treeMILC_models, model_results, type, ind, cov_ok, cov_problem, N) {
  
  # Factor ME and rename 1, 2, 3 to 10%, 20%, and 30% for aesthetic purposes (i.e., for in the plots)
  model_results$ME <- as.factor(model_results$ME)
  levels(model_results$ME)[1] <- "10%"
  levels(model_results$ME)[2] <- "20%"
  levels(model_results$ME)[3] <- "30%"
  
  # Create data frame to store results in
  df_to_plot <- data.frame()

  # Create temporary variables with "null" instead of NULL for other functions to work
  cov_ok1 <- cov_ok
  cov_problem1 <- cov_problem
  
  if (is.null(cov_ok)) {
    cov_ok1 <- "null"
  }
  
  if (is.null(cov_problem)) {
    cov_problem1 <- "null"
  }

  # Repeat for every amount of ME:
  for (k in levels(model_results$ME)[1:3]) {
    
    # Select indices of models with the correct number of indicators, covariates, and data set size
    indexes <- as.numeric(row.names(model_results[(model_results$ind == ind) & 
                                                   (model_results$cov_ok == cov_ok1) & 
                                                   (model_results$cov_problem == cov_problem1) &  
                                                   (model_results$N == N) & 
                                                   (model_results$ME == k),]))
    LC_diff_list <- LCT_diff_list <- treeMILC_diff_list <- list()
    N <- as.numeric(N)
    
    # To plot the bias:
    if (type == 1) {
      
      # Create list (for every model) of lists (for every indicator) with matrices that contain the bias
      for (i in indexes) {
        LC_diff_list <- append(LC_diff_list, list(get_diff_matrix(LC_models[[i]])))
        LCT_diff_list <- append(LCT_diff_list, list(get_diff_matrix(LCT_models[[i]])))
        treeMILC_diff_list <- append(treeMILC_diff_list, list(get_diff_matrix(treeMILC_models[[i]]))) 
      }     
      
      # Compute for every indicator the mean of the biases
      LC_mean_diff <- get_mean_matrix(LC_diff_list)
      LCT_mean_diff <- get_mean_matrix(LCT_diff_list)
      treeMILC_mean_diff <- get_mean_matrix(treeMILC_diff_list)
      
      # To plot the variance:
    } else if (type == 2) {
      
      # Create list (for every model) of lists (for every indicator) with matrices that contain the variance
      LC_mean_diff <- get_var_ME(LC_models, ind, cov_ok1, cov_problem1, N, which(levels(model_results$ME) == k))
      LCT_mean_diff <- get_var_ME(LCT_models, ind, cov_ok1, cov_problem1, N, which(levels(model_results$ME) == k))
      treeMILC_mean_diff <- get_var_ME(treeMILC_models, ind, cov_ok1, cov_problem1, N, which(levels(model_results$ME) == k)) 
    }
    
    # Take sum of matrices for every indicator
    LC_df_to_plot <- as.matrix(LC_mean_diff[[1]])
    LCT_df_to_plot <- as.matrix(LCT_mean_diff[[1]])
    treeMILC_df_to_plot <- as.matrix(treeMILC_mean_diff[[1]])
    
    for (i in 2:ind) {
      LC_df_to_plot <- LC_df_to_plot + as.matrix(LC_mean_diff[[i]])
      LCT_df_to_plot <- LCT_df_to_plot + as.matrix(LCT_mean_diff[[i]])
      treeMILC_df_to_plot <- treeMILC_df_to_plot + as.matrix(treeMILC_mean_diff[[i]])
    }
    
    # Compute mean of matrices (s.t. matrices are not separate anymore for every indicator) and convert to data frame
    LC_df_to_plot <- as.data.frame(as.table(LC_df_to_plot / ind))
    LC_df_to_plot$type <- "LC"
    LCT_df_to_plot <- as.data.frame(as.table(LCT_df_to_plot / ind))
    LCT_df_to_plot$type <- "LCT"
    treeMILC_df_to_plot <- as.data.frame(as.table(treeMILC_df_to_plot / ind))
    treeMILC_df_to_plot$type <- "tree-MILC"
    
    # Merge results for LC, LCT, and tree-MILC into one data frame
    temp <- rbind(LC_df_to_plot, LCT_df_to_plot, treeMILC_df_to_plot)
    temp$ME <- k
    colnames(temp) <- c("Model", "Indicator", "Bias", "type", "ME")
    df_to_plot <- rbind(df_to_plot, temp)
    colnames(df_to_plot) <- c("Model", "Indicator", "Bias", "type", "ME")
    
    # Get values in the right format (P = Permanent, F = Flexible, and O = Other)
    df_to_plot$Model <- as.character(df_to_plot$Model)
    df_to_plot$Indicator <- as.character(df_to_plot$Indicator)
    df_to_plot[df_to_plot$Model == "A" | df_to_plot$Model == "1", ]$Model <- "P"
    df_to_plot[df_to_plot$Model == "B" | df_to_plot$Model == "2", ]$Model <- "O"
    df_to_plot[df_to_plot$Model == "C" | df_to_plot$Model == "3", ]$Model <- "F"
    df_to_plot[df_to_plot$Indicator == "A" | df_to_plot$Indicator == "1", ]$Indicator <- "P"
    df_to_plot[df_to_plot$Indicator == "B" | df_to_plot$Indicator == "2", ]$Indicator <- "O"
    df_to_plot[df_to_plot$Indicator == "C" | df_to_plot$Indicator == "3", ]$Indicator <- "F"
    df_to_plot$Model <- factor(df_to_plot$Model, levels = c("O", "F", "P"))
    df_to_plot$Indicator <- factor(df_to_plot$Indicator, levels = c("P", "F", "O"))
  }
  
  # To plot the bias:
  if (type == 1) {
    # Multiply values on the diagonal by -1 to ensure consistent interpretation
    df_to_plot[df_to_plot$Model == "P" & df_to_plot$Indicator == "P", ]$Bias <- -df_to_plot[df_to_plot$Model == "P" & df_to_plot$Indicator == "P", ]$Bias
    df_to_plot[df_to_plot$Model == "O" & df_to_plot$Indicator == "O", ]$Bias <- -df_to_plot[df_to_plot$Model == "O" & df_to_plot$Indicator == "O", ]$Bias
    df_to_plot[df_to_plot$Model == "F" & df_to_plot$Indicator == "F", ]$Bias <- -df_to_plot[df_to_plot$Model == "F" & df_to_plot$Indicator == "F", ]$Bias
    
    # Compute ranges for the y-axis
    max_abs <- max(abs(df_to_plot$Bias))
    rng <- c(-max_abs, max_abs)
    
    # Plot heatmap
    g <- ggplot(df_to_plot, aes(Indicator, Model)) + 
      geom_tile(aes(fill = Bias)) + 
      scale_fill_gradientn(colors = rev(brewer.pal(9, "RdBu")), values = rescale(c(rng[1], 0, rng[2])), limits = c(rng[1], rng[2])) + 
      facet_grid(rows = vars(type), cols = vars(ME), labeller = labeller(.rows = label_value, .cols = label_both)) + 
      labs(fill = "(Mean) Bias")
    
    # To plot the variance:
  } else { 
    # Compute ranges for the y-axis
    rng <- c(min(df_to_plot$Bias), max(df_to_plot$Bias))
    df_to_plot$Variance <- df_to_plot$Bias
    
    # Plot heatmap
    g <- ggplot(df_to_plot, aes(Indicator, Model)) + 
      geom_tile(aes(fill = Variance)) +  
      scale_fill_gradientn(colors = brewer.pal(9, "Reds"), values = rescale(c(0, rng[2])), limits = c(0, rng[2])) + 
      facet_grid(rows = vars(type), cols = vars(ME), labeller = labeller(.rows = label_value, .cols = label_both)) + 
      labs(fill = "(Mean) Variance")
  }
  return(g)
}

##################################################################################################
## Create a grid of ME heatmaps for LC, LCT and tree-MILC models with a realistic 7% ME         ##
## -------------------------------------------------------------------------------------------- ##
## @param LC_models (list): List of LC models                                                   ##
## @param LCT_models (list): List of LCT models                                                 ##
## @param treeMILC_models (list): List of treeMILC models                                       ##
## @param model_results (data.frame): Results of the LC, LCT or tree-MILC models                ##
## @param type (int): Integer denoting whether to show bias or variance (1=bias, 2=variance)    ##
## @param ind (int): Number of indicators                                                       ##
## @param cov_ok (vector): Vector with the names of non-missing covariates                      ##
## @param cov_problem (vector): Vector with the names of missing covariates                     ##
## @param N (int): Size of data set                                                             ##
## @returns (ggplot object): ME heatmap                                                         ##
##################################################################################################

get_ME_heatmap_realistic <- function(LC_models, LCT_models, treeMILC_models, model_results, type, ind, cov_ok, cov_problem, N) {
  
  # Create data frame to store results in
  df_to_plot <- data.frame()

  # Create temporary variables with "null" instead of NULL for other functions to work
  cov_ok1 <- cov_ok
  cov_problem1 <- cov_problem
  
  if (is.null(cov_ok)) {
    cov_ok1 <- "null"
  }
  
  if (is.null(cov_problem)) {
    cov_problem1 <- "null"
  }
  
  # Select indices of models with the correct number of indicators, covariates, and data set size
  indexes <- as.numeric(row.names(model_results[(model_results$ind == ind) & (model_results$cov_ok == cov_ok1) & (model_results$cov_problem == cov_problem1) & (model_results$N == N) & (model_results$ME == 4), ]))
  N <- as.numeric(N)
  
  # Repeat for indicators (Y1 and Y3) and (Y2 and Y3)
  for (k in 1:2) {
    
    LC_diff_list <- LCT_diff_list <- treeMILC_diff_list <- list()
    
    # To plot the bias:
    if (type == 1) {
      
      # Create list (for every model) of lists (for every indicator) with matrices that contain the bias
      for (i in indexes) {
        LC_diff_list <- append(LC_diff_list, list(get_diff_matrix(LC_models[[i]])[seq(k, ind, 2)]))
        LCT_diff_list <- append(LCT_diff_list, list(get_diff_matrix(LCT_models[[i]])[seq(k, ind, 2)]))
        treeMILC_diff_list <- append(treeMILC_diff_list, list(get_diff_matrix(treeMILC_models[[i]])[seq(k, ind, 2)]))
      }
      
      # Get mean (per indicator) of these matrices
      LC_mean_diff <- get_mean_matrix(LC_diff_list)
      LCT_mean_diff <- get_mean_matrix(LCT_diff_list)
      treeMILC_mean_diff <- get_mean_matrix(treeMILC_diff_list)
      
      # To plot the variance:
    } else if (type == 2) {
      
      # Create list (for every model) of lists (for every indicator) with matrices that contain the variance
      LC_mean_diff <- get_var_ME(LC_models, ind, cov_ok1, cov_problem1, ME = 4, N, k)
      LCT_mean_diff <- get_var_ME(LCT_models, ind, cov_ok1, cov_problem1, N, ME = 4, k)
      treeMILC_mean_diff <- get_var_ME(treeMILC_models, ind, cov_ok1, cov_problem1, N, ME = 4, k)
    }
    
    # If there is only one indicator in the model (e.g. only Y1 (and not Y3)), convert to data frame
    if (length(seq(k, ind, 2)) == 1) {
      LC_df_to_plot <- as.data.frame(matrix(unlist(LC_mean_diff), ncol = 3, byrow = FALSE))
      LCT_df_to_plot <- as.data.frame(matrix(unlist(LCT_mean_diff), ncol = 3, byrow = FALSE))
      treeMILC_df_to_plot <- as.data.frame(matrix(unlist(treeMILC_mean_diff), ncol = 3, byrow = FALSE))
      
      # If there are two indicators, get the mean of the indicators and convert to data frame
    } else {
      LC_df_to_plot <- get_mean_matrix(LC_mean_diff)
      LC_df_to_plot <- as.data.frame(matrix(unlist(LC_df_to_plot), ncol = 3, byrow = FALSE))
      LCT_df_to_plot <- get_mean_matrix(LCT_mean_diff)
      LCT_df_to_plot <- as.data.frame(matrix(unlist(LCT_df_to_plot), ncol = 3, byrow = FALSE))
      treeMILC_df_to_plot <- get_mean_matrix(treeMILC_mean_diff)
      treeMILC_df_to_plot <- as.data.frame(matrix(unlist(treeMILC_df_to_plot), ncol = 3, byrow = FALSE))
    }
    
    # Merge results for LC, LCT, and tree-MILC into one data frame
    LC_df_to_plot$Model <- rownames(LC_df_to_plot)
    LC_df_to_plot <- melt(setDT(LC_df_to_plot), id.vars = c("Model"), variable.name = "Indicator")
    LCT_df_to_plot$Model <- rownames(LCT_df_to_plot)
    LCT_df_to_plot <- melt(setDT(LCT_df_to_plot), id.vars = c("Model"), variable.name = "Indicator")
    treeMILC_df_to_plot$Model <- rownames(treeMILC_df_to_plot)
    treeMILC_df_to_plot <- melt(setDT(treeMILC_df_to_plot), id.vars = c("Model"), variable.name = "Indicator")
    LC_df_to_plot$type <- "LC"
    LCT_df_to_plot$type <- "LCT"
    treeMILC_df_to_plot$type <- "tree-MILC"
    temp <- rbind(LC_df_to_plot, LCT_df_to_plot, treeMILC_df_to_plot)
    temp$indicator <- ifelse(k == 1, "Ind. 1", "Ind. 2")
    colnames(temp) <- c("Model", "Indicator", "Bias", "type", "ME")
    df_to_plot <- rbind(df_to_plot, temp)
    
  }
  
  # Get values in the right format (P = Permanent, F = Flexible, and O = Other)
  colnames(df_to_plot) <- c("Model", "Indicator", "Bias", "type", "ME")
  df_to_plot$Model <- as.character(df_to_plot$Model)
  df_to_plot$Indicator <- as.character(df_to_plot$Indicator)
  df_to_plot[df_to_plot$Model == "V1" | df_to_plot$Model == "1", ]$Model <- "P"
  df_to_plot[df_to_plot$Model == "V2" | df_to_plot$Model == "2", ]$Model <- "O"
  df_to_plot[df_to_plot$Model == "V3" | df_to_plot$Model == "3", ]$Model <- "F"
  df_to_plot[df_to_plot$Indicator == "V1" | df_to_plot$Indicator == "1", ]$Indicator <- "P"
  df_to_plot[df_to_plot$Indicator == "V2" | df_to_plot$Indicator == "2", ]$Indicator <- "O"
  df_to_plot[df_to_plot$Indicator == "V3" | df_to_plot$Indicator == "3", ]$Indicator <- "F"
  df_to_plot$Model <- factor(df_to_plot$Model, levels = c("O", "F", "P"))
  df_to_plot$Indicator <- factor(df_to_plot$Indicator, levels = c("P", "F", "O"))
  df_to_plot$ME <- as.factor(df_to_plot$ME)
  df_to_plot$Variance <- df_to_plot$Bias
  
  # To plot the bias
  if (type == 1) {
    
    # Multiply values on the diagonal by -1 to ensure consistent interpretation
    df_to_plot[df_to_plot$Model == "P" & df_to_plot$Indicator == "P", ]$Bias <- -df_to_plot[df_to_plot$Model == "P" & df_to_plot$Indicator == "P", ]$Bias
    df_to_plot[df_to_plot$Model == "O" & df_to_plot$Indicator == "O", ]$Bias <- -df_to_plot[df_to_plot$Model == "O" & df_to_plot$Indicator == "O", ]$Bias
    df_to_plot[df_to_plot$Model == "F" & df_to_plot$Indicator == "F", ]$Bias <- -df_to_plot[df_to_plot$Model == "F" & df_to_plot$Indicator == "F", ]$Bias
    
    # Compute ranges for y-axis
    max_abs <- max(abs(df_to_plot$Bias))
    rng <- c(-max_abs, max_abs)
    
    # Plot heatmap
    g <- ggplot(df_to_plot, aes(Indicator, Model)) + geom_tile(aes(fill = Bias)) + scale_fill_gradientn(colors = rev(brewer.pal(9, "RdBu")), values = rescale(c(rng[1], 0, rng[2])), limits = c(rng[1], rng[2])) + facet_grid(rows = vars(type), cols = vars(ME), labeller = labeller(.rows = label_value, .cols = label_value)) +
      labs(fill = "Bias")
    
    # To plot the variance
  } else if (type == 2) {
    
    # Compute ranges for y-axis
    rng <- c(min(df_to_plot$Variance), max(df_to_plot$Variance))
    
    # Plot heatmap
    g <- ggplot(df_to_plot, aes(Indicator, Model)) + geom_tile(aes(fill = Variance)) + scale_fill_gradientn(colors = brewer.pal(9, "Reds"), values = rescale(c(0, rng[2])), limits = c(0, rng[2])) + facet_grid(rows = vars(type), cols = vars(ME), labeller = labeller(.rows = label_value, .cols = label_value)) +
      labs(fill = "Variance")
  }
  
  return(g)
}
