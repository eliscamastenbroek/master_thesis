################################## Helpfunctions_Simulations.R ###################################
## This file contains a number of (help) functions that are required to perform the             ##
## simulations without missing covariates (see Chapter 4) and with missing covariates (see      ##
## Chapter 5) (all three approaches). This file contains the following functions:               ##                                                                              ##
##    - fix_extra_covcat                                                                        ##
##    - store_model_info                                                                        ##
##    - get_ME                                                                                  ##
################################################################################################## 

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
## Compute measurement error matrices (one per indicator) for an LC, LCT or tree-MILC model     ##
## @param results (list): Results of an LC, LCT or tree-MILC model                              ##
## @returns (list): List of measurement error matrices (one per indicator)                      ##
##################################################################################################

get_ME <- function(results){
  
  # For LC and LCT: Compute a ME probability matrix for each indicator using the get_ME_help function
  if (class(results[[2]]) != "list") {
    results <- results[[2]] # Ignore first data frame with model information
    
    which_ind <- results[, colnames(results) %in% c("Y1", "Y2", "Y3", "Y4")]
    to_return <- list()
    for (i in 1:ncol(which_ind)) {
      to_return <- append(to_return, list(get_ME_help(results, which_ind[, i])))
    }
    return(to_return)
  }
  
  # For tree-MILC: Compute ME matrix for each indicator based on the imputed values
  else {
    
    all_indicators <- c("Y1", "Y2", "Y3", "Y4")
    num_of_ind <- results[[1]]$ind
    results_2 <- results[[2]]
    to_return <- list() 
    
    for (i in 1:num_of_ind) {
      cluster_index <- which(names(results_2[[1]]) == "cluster")
      ind_index <- which(colnames(results_2[[1]]) == all_indicators[i])
      
      # Compute averages of the ME probability matrices for each bootstrap sample (per indicator) 
      summed_matrix <- prop.table(table(results_2[[1]][, cluster_index], results_2[[1]][, ind_index]), 1)
      for (j in 2:length(results_2)) {
        cluster_index <- which(names(results_2[[j]]) == "cluster")
        ind_index <- which(colnames(results_2[[j]]) == all_indicators[i])
        summed_matrix <- summed_matrix + prop.table(table(results_2[[j]][, cluster_index], results_2[[j]][, ind_index]), 1)
      }
      mean_matrix <- summed_matrix / length(results_2)
      to_return <- append(to_return, list(mean_matrix))
    }
    
    return(to_return)
  }
}
