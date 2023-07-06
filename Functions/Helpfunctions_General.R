################################# Methods_Helpfunctions.R ######################################## 
## This file contains a number of (help) functions that are required to perform LC, LCT and     ##
## tree-MILC analysis on simulated data (Chapters 4-5) and real data from the ER and the LFS    ##
## (Chapter 6). This file contains the following functions:                                     ## 
##    - fix_cluster_assignment                                                                  ##
##    - fix_cluster_bootstrap                                                                   ##
##    - fix_number_notation                                                                     ##         
##    - impute_value_step1                                                                      ##
##    - impute_value_step2                                                                      ##
##    - get_ME_help                                                                             ##
##################################################################################################

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