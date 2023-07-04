##################################################################################################
## Function to create a matrix with the measurement error probabilities used to simulate data   ##
## @param a2,...b32 (int): Logit parameters used to simulate the data                           ##
## @returns (matrix): Matrix with measurement error probabilities                               ##
##################################################################################################

create_ME_matrix = function(a2, a3, b22, b33, b23, b32) {
  
  #Compute rows separately
  row1 = c(1 / (1 + exp(a2) + exp(a3)), exp(a2) / (1 + exp(a2) + exp(a3)), exp(a3) / (1 + exp(a2) + exp(a3)))
  row2 = c(1 / (1 + exp(a2 + b22) + exp(a3 + b32)), exp(a2 + b22) / (1 + exp(a2 + b22) + exp(a3 + b32)), exp(a3 + b32) / (1 + exp(a2 + b22) + exp(a3 + b32)))
  row3 = c(1 / (1 + exp(a2 + b23) + exp(a3 + b33)), exp(a2 + b23) / (1 + exp(a2 + b23) + exp(a3 + b33)), exp(a3 + b33) / (1 + exp(a2 + b23) + exp(a3 + b33)))
  
  #Combine rows into a matrix
  matrix = matrix(c(row1, row2, row3), nrow = 3, ncol = 3, byrow = TRUE)
  
  return(matrix)
}

##################################################################################################
## Function to calculate entropy R2                                                             ##
## @param results (list): Results of one particular LC or LCT model                             ##
## @returns (vector): A vector containing the entropy and the entropy R2                        ##
##################################################################################################

get_entropy = function(results){
  
  # Ignore first data frame with model information
  results = results[[2]] 
  
  # For LC and LCT
  if(class(results) != "list"){ 
    lc1 = sum((results$p1) * log(results$p1))
    lc1 = replace(lc1, lc1 == "NaN", 0) # In case some posterior probabilities are 0
    lc2 = sum((results$p2) * log(results$p2))
    lc2 = replace(lc2, lc2 == "NaN", 0)
    lc3 = sum((results$p3) * log(results$p3))
    lc3 = replace(lc3, lc3 == "NaN", 0)
    entropy = -(lc1 + lc2 + lc3)
    entropy_squared = 1 - (entropy / (nrow(results) * log(3)))
  }
  
  # For tree-MILC
  else { 
    entropy = entropy_squared = NA
  }
  
  entropy_vector = c(entropy_squared, entropy)
  names(entropy_vector) = c("Entropy R-squared", "Entropy")
  
  return(entropy_vector)
}

##################################################################################################
## Function to compute the population proportion estimates                                      ##
## @param results (list): Results of one particular LC, LCT or tree-MILC model                  ##
## @returns (vector): A vector containing the estimated proportions per cluster                 ##
##################################################################################################

get_proportions = function(results){
  results = results[[2]]   # Ignore first data frame with model information
  
  # For LC or LCT: Compute proportions based on posterior probabilities
  if(class(results) != "list"){        
    prop1 = sum(results$p1)
    prop2 = sum(results$p2)
    prop3 = sum(results$p3)
    proportions_vector = c(prop1, prop2, prop3) / nrow(results)
    names(proportions_vector) = 1:3
    return(proportions_vector)
  } 
  
  # For tree-MILC: Compute proportions based on imputations
  else {                    
    # Compute proportions per bootstrap
    prop_per_bootstrap = list()
    for(i in 1:length(results)){
      prop_per_bootstrap = append(prop_per_bootstrap, list(summary(factor(results[[i]]$cluster)) / nrow(results[[i]])))
    }
    
    # Pool results
    pooled_proportions = colMeans(bind_rows(prop_per_bootstrap))   
    return(pooled_proportions)
  }
}
