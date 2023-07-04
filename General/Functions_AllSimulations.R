
##################### FUNCTIONS_PERFORMANCE_MEASURES_AND_PLOTS.R ################################# 
## This file contains the functions needed to compute the performance measures and creating the ## 
## plots in Chapters 4 and 5:                                                                   ##
##    - create_ME_matrix                                                                        ##                                                                          
##    - get_entropy                                                                             ## 
##    - get_proportions                                                                         ##
##    - get_mean_summed_bias                                                                    ##
##    - get_diff_matrix                                                                         ##
##    - get_mean_matrix                                                                         ##
##    - get_ME_per_contract                                                                     ##                                                                        
##################################################################################################

##################################################################################################
## Create a matrix with the measurement error probabilities used to simulate data               ##
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
## Calculate entropy and entropy R2                                                             ##
## @param results (list): Results of one particular LC or LCT model                             ##
## @returns (vector): A vector containing the entropy and the entropy R2                        ##
##################################################################################################

get_entropy = function(results){
  
  #Ignore first data frame with model information
  results = results[[2]] 
  
  #For LC and LCT
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
  
  #For tree-MILC
  else { 
    entropy = entropy_squared = NA
  }
  
  entropy_vector = c(entropy_squared, entropy)
  names(entropy_vector) = c("Entropy R-squared", "Entropy")
  
  return(entropy_vector)
}

##################################################################################################
## Compute population proportion estimates                                                      ##
## @param results (list): Results of one particular LC, LCT or tree-MILC model                  ##
## @returns (vector): A vector containing the estimated proportions per cluster                 ##
##################################################################################################

get_proportions = function(results){
  results = results[[2]]   #Ignore first data frame with model information
  
  #For LC or LCT: Compute proportions based on posterior probabilities
  if(class(results) != "list"){        
    prop1 = sum(results$p1)
    prop2 = sum(results$p2)
    prop3 = sum(results$p3)
    proportions_vector = c(prop1, prop2, prop3) / nrow(results)
    names(proportions_vector) = 1:3
    return(proportions_vector)
  } 
  
  #For tree-MILC: Compute proportions based on imputations
  else {                    
    #Compute proportions per bootstrap
    prop_per_bootstrap = list()
    for(i in 1:length(results)){
      prop_per_bootstrap = append(prop_per_bootstrap, list(summary(factor(results[[i]]$cluster)) / nrow(results[[i]])))
    }
    
    #Pool results
    pooled_proportions = colMeans(bind_rows(prop_per_bootstrap))   
    return(pooled_proportions)
  }
}

##################################################################################################
## Get mean summed bias (see Section 4.4.2)                                                     ##
## @param results (list): Results of one particular LC, LCT or tree-MILC model                  ##
## @returns (int): Mean summed bias                                                             ##
##################################################################################################

get_mean_summed_bias = function(results){
  
  #Get estimated measurement error probability matrices
  ME_results = get_ME(results)

  #Get true measurement error probability matrices
  ME = results[[1]]$ME
  ME_matrix = get(paste0("ME_matrix",ME))
  
  #If the amount of measurement error is realistic
  if(length(ME_matrix) == 2){
    diff_sum = 0 #Create integer to store the mean summed bias in
    for(j in 1:length(ME_results)){ #Iterate over the number of indicators
      #Define k to represent whether we need to use the ME probability matrix for Y1 and Y3 or the ME probability matrix for Y2 and Y4
      k = j 
      if(j > 2){
        k = k - 2
      }
      
      #Compute the bias of each element on the diagonal
      for(i in 1:3){ 
        diff_sum = diff_sum + (ME_results[[j]][i,i] - ME_matrix[[k]][i,i])
      }
    }
    
    return(diff_sum / length(ME_results)) #Return mean summed bias

  }
  
  #If the amount of measurement error is 10%, 20% or 30%
  else {
    diff_sum = 0 #Create integer to store the mean summed bias in
    for(j in 1:length(ME_results)){ #Iterate over the number of indicators 
      #Compute the bias of each element on the diagonal
      for(i in 1:3){ 
        diff_sum = diff_sum + (ME_results[[j]][i,i] - ME_matrix[[k]][i,i])
      }
    }
    
    #Return mean summed bias
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

get_diff_matrix = function(results){
  
  #Get true measurement error probability matrix
  ME = results[[1]]$ME
  ind = results[[1]]$ind
  ME_results = get_ME(results)
  ME_matrix = get(paste0("ME_matrix", ME))
  
  #Create list to store the results
  to_return = list()
  
  #If the amount of measurement error is 10%, 20%, or 30%
  if (ME != 4){
    #Get for every indicator the differences between the input and the estimated ME matrix and add to return list
    for (i in 1:ind){
      to_return = append(to_return, list((ME_results[[i]] - ME_matrix))) 
    }
  }
  
  #If the amount of measurement error is realistic
  if (ME == 4){
    #Get differences between input and estimated ME matrix for Y1 and Y2 and add to return list
    for (i in 1:2){
      to_return = append(to_return, list((ME_results[[i]] - ME_matrix[[i]])))  
    }
    #Get differences between input and estimated ME matrix for Y3 and add to return list 
    if (ind > 2){
      to_return = append(to_return, list((ME_results[[3]] - ME_matrix[[1]])))
    }
    #Get differences between input and estimated ME matrix for Y4 and add to return list
    if (ind > 3){
      to_return = append(to_return, list((ME_results[[4]] - ME_matrix[[2]])))
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

get_mean_matrix = function(list) {
  num_matrices = length(list) #Number of lists
  num_ind = length(list[[1]]) #Number of matrices per list
  to_return = list() #List to store results in
  
  #Compute mean of the matrices per indicator
  for (i in 1:num_ind) {
    temp_matrix = list[[1]][[i]]
    for (j in 2:num_matrices) {
      temp_matrix = temp_matrix + list[[j]][[i]]
    }
    temp_matrix = temp_matrix / num_matrices
    to_return = append(to_return, list(temp_matrix))
  }
  return(to_return)
}

##################################################################################################
## Compute (mean) ME probability estimate for one particular contract type.                     ##
## @param results (model): Results of one particular LC, LCT or tree-MILC model                 ##
## @param contract (int): Contract type (1=permanent, 2=other, 3=flexible)                      ##
## @returns (int): Mean ME probability estimate for one particular contract type.               ##
##################################################################################################

get_ME_per_contract = function(model, contract){
  
  #If the amount of measurement error is realistic
  ME = model[[1]]$ME
  if(ME==4){
    #Get (estimated) ME probability matrix
    ME_estimated = get_ME(model)
    ind = length(ME_estimated)
    
    first = ME_estimated[seq(1,ind,2)]  #ME probability matrices for Y1/Y3
    second = ME_estimated[seq(2,ind,2)] #ME probability matrices for Y2/Y4
    
    #If Y1 and Y3, compute mean ME probability estimate for both indicators
    if(length(first)>1){
      first = mean(first[[1]][contract,contract]+first[[2]][contract,contract])
    #If only Y1, get ME probability estimate 
    } else {
      first = first[[1]][contract,contract]
    }
    
    #If Y2 and Y4, compute mean ME probability estimate for both indicators
    if(length(second)>1){
      second = mean(second[[1]][contract,contract]+second[[2]][contract,contract])
    #If only Y2, get ME probability estimate
    } else {
      second = second[[1]][contract,contract]
    }
    return(c(first, second))
  
  #If the amount of measurement error is 10%, 20% or 30%
  } else {
    #Compute mean of ME probability matrices for all indicators
    ME_estimated = get_mean_matrix(get_ME(model))
    ME_estimated = matrix(unlist(ME_estimated), ncol = 3, byrow = TRUE)
    return(ME_estimated[contract,contract])
  }
}
