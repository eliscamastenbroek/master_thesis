## Load required packages
# library(stringr)
# library(dplyr)
# library(tidyr)
# library(data.table)
# library(ExPosition)
# library(gridExtra)
# library(scales)
# library(ggplot2)
# library(RColorBrewer)

##################################################################################################
## Function to generate a data set                                                              ##
## @param seed (int): Seed for generating the data set                                          ##
## @param ME (int): Amount of measurement error (1=0.1, 2=0.2, 3=0.3, 4=realistic)              ##
## @param folder (string): Folder to save files in                                              ##
## @returns (data.frame): A simulated data set of size n=20,000                                 ##
##################################################################################################

simulate_data = function(seed, ME, folder){

  #Create data structure file required by Latent GOLD
  filepath_input = paste0(folder,"exampleData.dat")
  exampleData = paste0("id q Y1 Y2 Y3 Y4 w20000
1 1 1 1 1 1 2500
2 2 2 2 2 2 2500
3 1 3 3 3 3 2500
4 2 3 3 3 3 2500")
  writeLines(exampleData, filepath_input)
  
  #Create Latent Gold script
  script_part1 = paste0("version = 6.0\ninfile '", filepath_input, "' \n\nmodel
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
  outfile_path = paste0(folder, "simDat", ME, "_iteration", seed, ".dat")
  
  #Logit parameters for P(X=k) based on real data
  data_a2 = -0.9662
  data_a3 = -1.6260
  data_b22 = 0.22
  data_b32 = 0.2607
  data_param = paste(data_a2, data_a3, data_b22, data_b32)
  
  #Logit parameters for P(Yj=yj|X=k) for 10% measurement error
  if(ME == 1){
    a2 = a3 = -log(18)  
    b22 = b33 = log(324)
    b32 = b23 = log(18)
    ME_coefs = c(a2, a3, b22, b32, b23, b33, "\n")
    ME_coefs = gsub(",", "", toString(rep(ME_coefs, 4)))
    parameters = paste("\n", data_param, "\n", ME_coefs, "}\nend model")
  }
  
  #Logit parameters for P(Yj=yj|X=k) for 20% measurement error
  else if(ME == 2){
    a2 = a3 = -3*log(2) 
    b22 = b33 = 6*log(2)
    b32 = b23 = 3*log(2)
    ME_coefs = c(a2, a3, b22, b32, b23, b33, "\n")
    ME_coefs = gsub(",", "", toString(rep(ME_coefs, 4)))
    parameters = paste("\n", data_param, "\n", ME_coefs, "}\nend model")
  }
  
  #Logit parameters for P(Yj=yj|X=k) for 30% measurement error
  else if(ME == 3){
    a2 = a3 = -1.54045 
    b22 = b33 = 3.0809
    b32 = b23 = 1.54045
    ME_coefs = c(a2, a3, b22, b32, b23, b33, "\n")
    ME_coefs = gsub(",", "", toString(rep(ME_coefs, 4)))
    parameters = paste("\n", data_param, "\n", ME_coefs, "}\nend model")
  }
  
  #Logit parameters for P(Yj=yj|X=k) for realistic amount of measurement error
  else if(ME == 4){
    
    #Coefficients for indicator 1 (and 3)
    Y1_a2 = -4.4917 
    Y1_a3 = -4.94368
    Y1_b22 = 7.64123
    Y1_b32 = 4.55064
    Y1_b23 = 3.09876
    Y1_b33 = 5.6678
    
    #Coefficients for indicator 2 (and 4)
    Y2_a2 = -6.14311 
    Y2_a3 = -2.63157
    Y2_b22 = 11.9482
    Y2_b32 = 1.53296
    Y2_b23 = 1.72427
    Y2_b33 = 5.03275
    
    ME_coefs_Y1 = c(Y1_a2, Y1_a3, Y1_b22, Y1_b32, Y1_b23, Y1_b33, "\n")
    ME_coefs_Y2 = c(Y2_a2, Y2_a3, Y2_b22, Y2_b32, Y2_b23, Y2_b33, "\n")
    ME_coefs = gsub(",", "", toString(rep(c(ME_coefs_Y1, ME_coefs_Y2), 2))) 
    parameters = paste(data_param, "\n", ME_coefs, "}\nend model")
  }
  
  script_part2 = paste0("\toutfile '", outfile_path, "' simulation=1 seed=", seed, ";
    variables
         caseid id;
         caseweight w20000;
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
  
  #Combine parts of script
  script = paste0(script_part1, script_part2, parameters)
  script_path = paste0(folder, "simDat", ME, "_iteration", seed, "_script.lgs")
  writeLines(script, script_path)
  
  #Execute Latent Gold script
  shell(paste0('"C:/Program Files/LatentGOLDnet6.0/lg60.exe" ', script_path, ' /b'))
  
  #Import simulated data set
  simDat = read.delim(outfile_path, sep="\t", dec=",")
  return(simDat)
}

##################################################################################################
## Function to create a subset of the right size for each analysis                              ##
## @param iteration (int): Iteration number                                                     ## 
## @param ind (int): Number of indicators                                                       ##
## @param cov (int): Number of covariates                                                       ##
## @param N (int): Size of the data set                                                         ##
## @param ME (int): Amount of measurement error (1=0.1, 2=0.2, 3=0.3, 4=realistic)              ##
## @returns (data.frame): A subset                                                              ##
##################################################################################################

create_subset = function(iteration, ind, cov, N, ME) {
  
  #If a certain data set does not exist in the global environment, create it 
  if (!exists(paste0("simDat", ME, "_iteration", iteration))) {
    assign(paste0("simDat", ME, "_iteration", iteration), simulate_data(iteration, ME), envir = globalenv())
  }
  
  #Get simulated data set
  data = get(paste0("simDat", ME, "_iteration", iteration))
  
  #Set seed to get the same data set for every model within each iteration
  set.seed(iteration) 
  select_cases = sample(1:nrow(data), N, replace = FALSE)
  
  #Remove some redundant columns
  if (cov == 1) {
    subset = data[select_cases, c(2:(3 + ind))]
  } else {
    subset = data[select_cases, c(2, 4:(3 + ind))]
  }
  
  return(subset)
}

##################################################################################################
## Function to generate a Latent Gold script for estimating an LC model (i.e. for LC,           ## 
## LCT step 1, LCT step 2, tree-MILC step 1, tree-MILC step 2)                                  ##
## @param type (string): What type of script ("LC" for regular LC and "LCT" for LCT step 2)     ## 
## @param ind (int): Number of indicators                                                       ##
## @param cov (int): Number of covariates                                                       ##
## @param N (int): Size of data set                                                             ##
## @param model_name (string): Name of the model                                                ##
## @param filepath_input (string): Path to input data set                                       ##
## @param filepath_output (string): Path where model output should be stored                    ##
## @returns (string): A string containing a Latent Gold script                                  ##
##################################################################################################

generate_script = function(type, ind, cov, N, model_name, filepath_input, filepath_output){
  
  script_part1 = paste0("version = 6.0\ninfile '", filepath_input, "' \n\nmodel title '")
  
  #Let the number of sets of starting values depend on the size of the data set
  if(N < 10000){
    sets = 3200
  } else {
    sets = 100
  }
  
  script_part2 = paste0("';
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
  
  script_part3 = paste0("first profile standarderrors reorderclasses;
    	outfile '", filepath_output, "' classification keep=id;
    variables\n") 
  
  #Adjust some parameters depending on what type of analysis is performed
  if(type == "LC"){
    latent_var = paste0("\n\tlatent Cluster nominal 3;
    equations\n")
  } else {
    latent_var = paste0("\n\tlatent Cluster nominal 2;
    equations\n")
  }
  
  #For LCT step 2, use posterior probabilities from step 1 as weights 
  if(type == "LCT"){
    caseweight = "caseweight p1;\n"
  } else {
    caseweight = ""
  }
  
  #For LCT step 2 and four indicators, add restriction that 
  #Y2(2) < Cluster(2) = 0 and Y4(2) < Cluster(2) = 0 to avoid boundary solutions 
  if((type == "LCT") & ind == 4){
    restrictions = "\n\tx[1,1] = 0; y[1,1] = 0;"
  } else { 
    restrictions = ""
  }
  
  #Adjust equations depending on the number of indicators
  if(ind == 2){
    dep_ind = "\tdependent Y1 nominal, Y2 nominal;"
    dep_ind_eq = "\tY1 <- 1 + Cluster;\n\tY2 <- 1 + (x) Cluster;"
  } else if(ind ==3){
    dep_ind = "\tdependent Y1 nominal, Y2 nominal, Y3 nominal;"
    dep_ind_eq = "\tY1 <- 1 + Cluster;\n\tY2 <- 1 + (x) Cluster;\n\tY3 <- 1 + Cluster;"
  } else if(ind ==4){
    dep_ind = "\tdependent Y1 nominal, Y2 nominal, Y3 nominal, Y4 nominal;"
    dep_ind_eq = "\tY1 <- 1 + Cluster;\n\tY2 <- 1 + (x) Cluster;\n\tY3 <- 1 + Cluster;\n\tY4 <- 1 + (y) Cluster;"
  }
  
  #Adjust equations depending on whether to include a covariate or not
  if(cov == 1){
    dep_cov = "\n\tindependent q nominal;"
    latent_var_eq = "\tCluster <- 1 + q;\n"
  }
  else if(cov == 0){
    dep_cov = ""
    latent_var_eq = "\tCluster <- 1;\n"
  }
  
  #Combine all parts of the script
  script = paste0(script_part1, model_name, script_part2, script_part3, caseweight, dep_ind, dep_cov, 
                latent_var, latent_var_eq, dep_ind_eq, restrictions, "\nend model")
  return(script)
}

##################################################################################################
## Function to generate a Latent Gold script to obtain posterior probabilities for observations ## 
## in tree-MILC with response patterns that were not present in the bootstrap sample. This is   ##
## done by estimating an LC model on the original data set while using the logit parameters     ##
## from the previous model as starting values and setting the number of EM- and NR-iterations   ## 
## to 0.                                                                                        ##
## @param type (string): What type of script ("LC" for regular LC and "LCT" for LCT step 2)     ## 
## @param ind (int): Number of indicators                                                       ##
## @param cov (int): Number of covariates                                                       ##
## @param N (int): Size of data set                                                             ## 
## @param model_name (string): Name of the model                                                ##
## @param filepath_input (string): Path to input data set                                       ##
## @param filepath_output (string): Path where model output should be stored                    ##
## @param par (data.frame): Data frame with logit parameters as estimated by the previous model ## 
## @returns (string): A string containing a Latent Gold script                                  ##
##################################################################################################

generate_script_treeMILC_step1 = function(type, ind, cov, N, model_name, filepath_input, filepath_output, par){
  
  script_part1 = paste0("version = 6.0\ninfile '",filepath_input,"' \n\nmodel title '")
  
  #Let the number of sets of starting values depend on the size of the data set
  if(N < 10000){
    sets = 3200
  } else {
    sets = 100
  }
  
  script_part2 = paste0("';
    options
    algorithm
        tolerance=1e-08 emtolerance=0.01 emiterations=0 nriterations=0;
    startvalues
        seed=1 sets=",sets," tolerance=1e-05 iterations=100;
    bayes
        categorical=1 variances=1 latent=1 poisson=1;
    missing includeall;
    output
    	parameters=")
  
  script_part3 = paste0("first profile;
    	outfile '",filepath_output,"' classification keep=id;
    variables\n") 
  
  #Adjust some parameters depending on what type of analysis is performed
  latent_var = paste0("\n\tlatent Cluster nominal 3;\nequations\n")
  
  #Adjust equations depending on the number of indicators
  if(ind==2){
    dep_ind = "\tdependent Y1 nominal, Y2 nominal;"
    dep_ind_eq = "\tY1 <- (c) 1 + (d) Cluster;\n\tY2 <- (e) 1 + (f) Cluster;"
  } else if(ind==3){
    dep_ind = "\tdependent Y1 nominal, Y2 nominal, Y3 nominal;"
    dep_ind_eq = "\tY1 <- (c) 1 + (d) Cluster;\n\tY2 <- (e) 1 + (f) Cluster;\n\tY3 <- (g) 1 + (h) Cluster;"
  } else if(ind==4){
    dep_ind = "\tdependent Y1 nominal, Y2 nominal, Y3 nominal, Y4 nominal;"
    dep_ind_eq = "\tY1 <- (c) 1 + (d) Cluster;\n\tY2 <- (e) 1 + (f) Cluster;\n\tY3 <- (g) 1 + (h) Cluster;\n\tY4 <- (i) 1 + (j) Cluster;"
  }
  
  #Create a string where the starting values are stored
  start_val="\n"
  
  #Adjust equations depending on whether to include a covariate or not
  if(cov==1){
    dep_cov = "\n\tindependent q nominal;"
    latent_var_eq = "\tCluster <- (a) 1 + (b) q;\n"
    start_val = paste0(start_val,"\tb[1,1] ~= ",par[3,]$coef,";\n")
    start_val = paste0(start_val,"\tb[1,2] ~= ",par[4,]$coef,";\n")
  }
  else if(cov==0){
    dep_cov = ""
    latent_var_eq = "\tCluster <- (a) 1;\n"
  }
  
  #Specify starting values here
  start_val = paste0(start_val,  "\ta[1,1] ~= ", par[1, ]$coef, ";\n")
  start_val = paste0(start_val, "\ta[1,2] ~= ", par[2, ]$coef, ";\n")
  start_val = paste0(start_val, "\tc[1,1] ~= ", par[4 + cov + cov, ]$coef, ";\n")
  start_val = paste0(start_val, "\tc[1,2] ~= ", par[5 + cov + cov, ]$coef, ";\n")
  start_val = paste0(start_val, "\td[1,1] ~= ", par[6 + cov + cov, ]$coef, ";\n")
  start_val = paste0(start_val, "\td[1,2] ~= ", par[7 + cov + cov, ]$coef, ";\n")
  start_val = paste0(start_val, "\td[1,3] ~= ", par[8 + cov + cov, ]$coef, ";\n")
  start_val = paste0(start_val, "\td[1,4] ~= ", par[9 + cov + cov, ]$coef, ";\n")
  start_val = paste0(start_val, "\te[1,1] ~= ", par[11 + cov + cov, ]$coef, ";\n")
  start_val = paste0(start_val, "\te[1,2] ~= ", par[12 + cov + cov, ]$coef, ";\n")
  start_val = paste0(start_val, "\tf[1,1] ~= ", par[13 + cov + cov, ]$coef, ";\n")
  start_val = paste0(start_val, "\tf[1,2] ~= ", par[14 + cov + cov, ]$coef, ";\n")
  start_val = paste0(start_val, "\tf[1,3] ~= ", par[15 + cov + cov, ]$coef, ";\n")
  start_val = paste0(start_val, "\tf[1,4] ~= ", par[16 + cov + cov, ]$coef, ";\n")
  
  if(ind >= 3){
    start_val = paste0(start_val, "\tg[1,1] ~= ", par[18 + cov + cov, ]$coef, ";\n")
    start_val = paste0(start_val, "\tg[1,2] ~= ", par[19 + cov + cov, ]$coef, ";\n")
    start_val = paste0(start_val, "\th[1,1] ~= ", par[20 + cov + cov, ]$coef, ";\n")
    start_val = paste0(start_val, "\th[1,2] ~= ", par[21 + cov + cov, ]$coef, ";\n")
    start_val = paste0(start_val, "\th[1,3] ~= ", par[22 + cov + cov, ]$coef, ";\n")
    start_val = paste0(start_val, "\th[1,4] ~= ", par[23 + cov + cov, ]$coef, ";\n")
  } 
  if(ind == 4){
    start_val = paste0(start_val, "\ti[1,1] ~= ", par[25 + cov + cov, ]$coef, ";\n")
    start_val = paste0(start_val, "\ti[1,2] ~= ", par[26 + cov + cov, ]$coef, ";\n")
    start_val = paste0(start_val, "\tj[1,1] ~= ", par[27 + cov + cov, ]$coef, ";\n")
    start_val = paste0(start_val, "\tj[1,2] ~= ", par[28 + cov + cov, ]$coef, ";\n")
    start_val = paste0(start_val, "\tj[1,3] ~= ", par[29 + cov + cov, ]$coef, ";\n")
    start_val = paste0(start_val, "\tj[1,4] ~= ", par[30 + cov + cov, ]$coef, ";\n")
  }
  
  #Combine all parts of the script
  script = paste0(script_part1, model_name, script_part2, script_part3, dep_ind, dep_cov, 
                latent_var, latent_var_eq, dep_ind_eq, start_val, "\nend model")
  return(script)
}

##################################################################################################
## Function to fix cluster assignments for LC and LCT models (i.e. assign the right names       ## 
## to the right clusters)                                                                       ##
## @param type (string): What type of model; only relevant is "LC" for LC)                      ##
## @param results (list): Results of one particular LC, LCT or tree-MILC model                  ##
## @returns (list): Same object as input, but with corrected cluster assignments if necessary   ##
##################################################################################################

fix_cluster_assignment = function(type = NULL, results) {
  
  #Get a list of matrices with measurement error per indicator
  ME_list = get_ME(results)
  
  #Create one matrix that contains the average values of all matrices
  summed_matrices = ME_list[[1]]
  for (j in 2:length(ME_list)) {
    summed_matrices = summed_matrices + ME_list[[j]]
  }
  mean_matrix = summed_matrices / length(ME_list)
  
  #Find diagonal combinations of cluster names and store them in a data frame
  all_diagonals = data.frame()
  
  if (type == "LC") { #Find all possible diagonal combinations
    for (i in 1:3) {
      for (j in 1:3) {
        for (k in 1:3) {
          if (length(unique(c(i, j, k))) == 3) {
            all_diagonals = rbind(all_diagonals, c(i, j, k, mean_matrix[i, 1], mean_matrix[j, 2], mean_matrix[k, 3],
                                                   sum(mean_matrix[i, 1], mean_matrix[j, 2], mean_matrix[k, 3])))
          }
        }
      }
    }
  } else {
    #For non-LC models: only find combinations for switched 1 and 3's, because 2 is already correct
    all_diagonals = rbind(all_diagonals, c(1, 2, 3, mean_matrix[1, 1], mean_matrix[2, 2], mean_matrix[3, 3], sum(mean_matrix[1, 1], mean_matrix[2, 2], mean_matrix[3, 3])))
    all_diagonals = rbind(all_diagonals, c(3, 2, 1, mean_matrix[3, 1], mean_matrix[2, 2], mean_matrix[1, 3], sum(mean_matrix[3, 1], mean_matrix[2, 2], mean_matrix[1, 3])))
  }
  
  colnames(all_diagonals) = c("1", "2", "3", "d1", "d2", "d3", "sum")
  
  #Find out which combination of diagonals yields the highest sum of diagonal values
  which_max = which.max(as.vector(all_diagonals[, 7]))
  
  #Find out how to reassign clusters
  max_1 = all_diagonals[which_max, 1]
  max_2 = all_diagonals[which_max, 2]
  max_3 = all_diagonals[which_max, 3]
  reassignment = c(as.numeric(max_1), as.numeric(max_2), as.numeric(max_3))
  
  #If clusters need not to be reassigned, return original results
  if (identical(reassignment, c(1, 2, 3))) {
    return(results)
  } else {
    #Create a list to store the corrected results (to ensure the output has the same format as the input)
    to_return = list(results[[1]])
    results = results[[2]]
    
    #Assign the posterior probabilities to the right cluster names
    posteriors = data.frame(p1 = results$p1, p2 = results$p2, p3 = results$p3)
    results$p1 = posteriors[, max_1]
    results$p2 = posteriors[, max_2]
    results$p3 = posteriors[, max_3]
    to_return = append(to_return, list(results))
    return(to_return)
  }
}


##################################################################################################
## Function to fix cluster assignments per tree-MILC bootstrap sample                           ## 
## @param boot_results (data.frame): Results of one bootstrap sample                            ##
## @returns (data.frame): Same as input, but with corrected cluster assignments if necessary    ##
##################################################################################################

fix_cluster_bootstrap = function(boot_results) {
  
  #Compute a contingency table for the value as observed by indicator Y2 and the imputed value
  table = table(boot_results$Y2, boot_results$cluster)
  prop_table = prop.table(table, margin = 1)
  
  #Find out which combination of diagonals yields the highest sum of diagonal values (assume that 2 is already correct)
  all_diagonals = data.frame()
  all_diagonals = rbind(all_diagonals, c(1, 2, 3, prop_table[1, 1], prop_table[2, 2], prop_table[3, 3], sum(prop_table[1, 1], prop_table[2, 2], prop_table[3, 3])))
  all_diagonals = rbind(all_diagonals, c(3, 2, 1, prop_table[3, 1], prop_table[2, 2], prop_table[1, 3], sum(prop_table[3, 1], prop_table[2, 2], prop_table[1, 3])))
  colnames(all_diagonals) = c("1", "2", "3", "d1", "d2", "d3", "sum")
  which_max = which.max(as.vector(all_diagonals[, 7]))
  
  #Find out how to reassign clusters
  max_1 = all_diagonals[which_max, 1]
  max_3 = all_diagonals[which_max, 3]
  reassignment = c(max_1, max_3)
  
  if (identical(reassignment, c(1, 3))) {
    return(boot_results)
  } else {
    boot_results$new_cluster = NA
    boot_results[boot_results$cluster == max_1, ]$new_cluster = 1
    boot_results[boot_results$cluster == 2, ]$new_cluster = 2
    boot_results[boot_results$cluster == max_3, ]$new_cluster = 3
    boot_results$cluster = boot_results$new_cluster
    boot_results = boot_results[, -which(colnames(boot_results) == "new_cluster")]
  }
  
  return(boot_results)
}


##################################################################################################
## Helpfunction to fix number notation in Latent GOLD output (i.e. convert '1e-02' to '.001')   ##
## @param vector (vector): Input vector (potentially) containing numbers in wrong format        ## 
## @returns (vector): Output vector containing corrected numbers                                ##
##################################################################################################

fix_number_notation = function(vector) {
  
  #If the vector is numeric, all values are already in the correct format
  if (is.numeric(vector)) {
    return(vector)
    
  #If the vector is not numeric, some values need to be corrected
  } else {
    #Create vector to store corrected results in
    return_vec = rep(NA, length(vector))
    
    for (i in 1:length(vector)) {
      #Remove spaces from the string and split string by 'e-'
      removed_spaces = gsub(" ", "", vector[i])  
      split_vec = str_split(removed_spaces, "e-") 
      
      #If the string did contain 'e-'
      if (length(split_vec[[1]]) > 1) {
        #Convert number to the right format and add to corrected results vector
        nominator = as.numeric(gsub(",", ".", split_vec[[1]][1]))  #Replace comma with dot
        denominator = as.numeric(split_vec[[1]][2])  
        final_number = nominator / (10^denominator)  
        return_vec[i] = final_number
        
      #If the string did not contain 'e'
      } else {
        temp_string = gsub(",", ".", split_vec[[1]][1])  #Replace comma with dot
        
        #In some cases, the output was '.', meaning that the probability was 0
        if (temp_string == ".") {
          return_vec[i] = 0  
        #In other cases, simply convert to numberic
        } else { 
          return_vec[i] = as.numeric(gsub(",", ".", split_vec[[1]][1]))  #Replace comma with dot
        }
      }
    }
    return(return_vec)
  }
}


##################################################################################################
## Function to perform LC                                                                       ##
## @param iteration (int): Iteration number                                                     ## 
## @param ind (int): Number of indicators                                                       ## 
## @param cov (int): Number of covariates                                                       ## 
## @param N (int): Size of data set                                                             ## 
## @param ME (int): Amount of measurement error (1=0.1, 2=0.2, 3=0.3, 4=realistic)              ##
## @param folder (string): Folder to save files in                                              ##
## @returns (list): A list that consists of:                                                    ## 
##  [[1]] Data frame with an overview of model parameters (iteration, ind, cov, N, ME)          ##
##  [[2]] Data frame with model results (posterior probabilities for each observation)          ##
##################################################################################################

perform_lc = function(iteration, ind, cov, N, ME, dat = NULL, folder) {
  
  #Store information about the model 
  model_name = paste(iteration, ind, cov, N, ME, sep = "-")
  model_info = data.frame(iteration = iteration, ind = ind, cov = cov, N = N, ME = ME, id = model_name)
  to_return = list(model_info) 
  
  #Write data set to use to file
  if (is.null(dat)) {
    dat = create_subset(iteration, ind, cov, N, ME)
  } else {
    dat = dat
  }
  
  filepath_input = paste0(folder, "LC_", model_name, "_data.dat")
  write.table(x = dat, file = filepath_input, row.names = FALSE, quote = FALSE)
  
  #Create Latent Gold script for LC
  filepath_output = paste0(folder, paste0("LC_", model_name, "_output.dat"))
  script = generate_script("LC", ind, cov, N, ME, model_name, filepath_input, filepath_output)
  script_path = paste0(folder, "LC_", model_name, "_script.lgs")
  writeLines(script, script_path)
  
  #Execute Latent Gold script
  shell(paste0('"C:/Program Files/LatentGOLDnet6.0/lg60.exe" ', script_path, ' /b'))
  
  #Read model output
  model_output = read.delim(filepath_output, sep = "\t", dec = ",")
  
  #Check if no warning was given
  model_lst = paste(readLines(paste0(folder, "LC_", model_name, "_script.lst")), collapse = "\n")
  if (grepl("WARNING", model_lst, fixed = TRUE)) {
    stop("Error")
  } 
  
  #Rename some columns and add data frame to return list
  setnames(model_output, old = c("Cluster.1", "Cluster.2", "Cluster.3", "Cluster."), new = c("p1", "p2", "p3", "cluster"))
  to_return = append(to_return, list(model_output))
  
  #Make sure clusters are assigned the right names
  to_return = fix_cluster_assignment(type = "LC", results = to_return)
  return(to_return)
}

##################################################################################################
## Function to perform LCT                                                                      ##
## @param iteration (int): Iteration number                                                     ## 
## @param ind (int): Number of indicators                                                       ## 
## @param cov (int): Number of covariates                                                       ## 
## @param N (int): Size of data set                                                             ## 
## @param ME (int): Amount of measurement error (1=0.1, 2=0.2, 3=0.3, 4=realistic)              ##
## @param folder (string): Folder to save files in                                              ##
## @returns (list): A list that consists of:                                                    ## 
##  [[1]] Data frame with an overview of model parameters (iteration, ind, cov, N, ME)          ##
##  [[2]] Data frame with model results (posterior probabilities for each observation)          ##
##################################################################################################

perform_lct <- function(iteration, ind, cov, N, ME, folder) {
  
  # Store information about the model 
  model_name <- paste(iteration, ind, cov, N, ME, sep = "-")
  model_info <- data.frame(iteration = iteration, ind = ind, cov = cov, N = N, ME = ME, id = model_name)
  to_return <- list(model_info) 
  
  #Write data set to use to file
  dat <- create_subset(iteration, ind, cov, N, ME)
  dat_name <- paste0("LCT_", model_name, "_step1_data.dat")
  write.table(x = dat, file = paste0(folder, dat_name), row.names = FALSE, quote = FALSE)
  
  #Perform LC with 3 classes and combine posterior probabilities for permanent & flexible
  model1_output <- perform_lc(iteration, ind, cov, N, ME, folder = folder)[[2]]
  model1_output$p1 <- 1 - model1_output$p2 #This is actually p1 + p3, but call it p1 for generate_script function 
  dat_name <- paste0("LCT_", model_name, "_step1_output.dat")
  write.table(x = model1_output, file = paste0(folder, dat_name), row.names = FALSE, quote = FALSE)
  
  #Create Latent Gold script for second LC model
  filepath_input <- paste0(folder, dat_name)
  filepath_output <- paste0(folder, "LCT_", model_name, "_step2_output.dat")
  script <- generate_script("LCT", ind, cov, N, ME, model_name, filepath_input, filepath_output)
  script_path <- paste0(folder, "LCT_", model_name, "_step2_script.lgs")
  writeLines(script, script_path)
  
  #Execute script in Latent Gold and read model output
  shell(paste0('"C:/Program Files/LatentGOLDnet6.0/lg60.exe" ', script_path, ' /b'))
  model2_output <- read.delim(filepath_output, sep = "\t", dec = ",")
  
  #Check if model is valid (i.e. does not contain a warning)
  model_lst <- paste(readLines(paste0(folder, "LCT_", model_name, "_step2_script.lst")), collapse = "\n")
  if (grepl("WARNING", model_lst, fixed = TRUE)) {
    stop("Error in Model 2")
  }
  
  #Combine results from both models to compute final posterior probabilities and fix column names
  combined_output <- left_join(x = model1_output, y = model2_output, by = c("id", "Y1", "Y2"))
  possible_colnames <- c("Y1", "Y2", "id", "cluster", "q.y", "Y3.y", "Y4.y", "p1.y", "Cluster.1", "Cluster.2", "Cluster.")
  neat_colnames <- c("Y1", "Y2", "id", "cluster.1", "q", "Y3", "Y4", "p1.1", "p2.1", "p2.2", "cluster.2")
  real_colnames <- names(combined_output)
  select_cols <- which(real_colnames %in% possible_colnames)
  combined_output <- combined_output[, select_cols]
  
  to_return <- c(to_return, combined_output)
  return(to_return)
}

##################################################################################################
## Function to obtain first round of imputations in tree-MILC                                   ##
## @param x (vector): Vector that contains the posterior probabilities p1, p2, and p3           ##
## @returns (int): Imputed value                                                                ##
##################################################################################################

impute_value_step1 <- function(x) {
  return(which(rmultinom(1, 1, c(as.numeric(x["p1"]) + as.numeric(x["p3"]), as.numeric(x["p2"]))) == 1))
}

##################################################################################################
## Function to obtain second round of imputations in tree-MILC                                  ##
## @param x (vector): Vector that contains the posterior probabilities p1 and p2                ##
## @returns (int): Imputed value                                                                ##
##################################################################################################

impute_value_step2 <- function(x) {
  if (is.na(x["Cluster#1"])) {
    return(NA)
  } else {
    return(which(rmultinom(1, 1, c(as.numeric(x["Cluster#1"]), as.numeric(x["Cluster#2"]))) == 1))
  }
}

##################################################################################################
## Function to perform tree-MILC                                                                ##
## @param iteration (int): Iteration number                                                     ## 
## @param ind (int): Number of indicators                                                       ## 
## @param cov (int): Number of covariates                                                       ## 
## @param N (int): Size of data set                                                             ## 
## @param ME (int): Amount of measurement error (1=0.1, 2=0.2, 3=0.3, 4=realistic)              ##
## @param M (int): Number of bootstrap samples                                                  ##
## @param folder (string): Folder to save files in                                              ##
## @returns (list): A list that consists of:                                                    ## 
##  [[1]] Data frame with an overview of model parameters (iteration, ind, cov, N, ME)          ##
##  [[2]] Data frame with model results (posterior probabilities for each observation)          ##
##################################################################################################

perform_treeMILC <- function(iteration, ind, cov, N, ME, M=5, folder) {
  
  # Store information about the model
  model_name = paste(iteration, ind, cov, N, ME, sep = "-")
  model_info = data.frame(iteration = iteration, ind = ind, cov = cov, N = N, ME = ME, id = model_name)
  to_return = list(model_info)
  
  # Create original data set and write to file
  dat_org = create_subset(iteration, ind, cov, N, ME)
  dat_org_path = paste0(folder, "tree_MILC_", model_name, "_dat_org.dat")
  fwrite(dat_org, file = dat_org_path, sep = "\t")
  
  # Count combinations of indicators (+covariate) in the original data set
  count_dat = as.data.frame(dat_org[, -which(colnames(dat_org) %in% c("id"))] %>% group_by_all() %>% summarise(COUNT = n()))
  count_dat = count_dat[, -ncol(count_dat)]
  
  # Create help vector to combine results later
  all_indicators = c("Y1", "Y2", "Y3", "Y4")
  by_vector = c(all_indicators[1:ind])
  if (cov == 1) {
    by_vector = c(by_vector, "q")
  }
  
  # Create list to store the results for each bootstrap sample in
  bootstrap_results = list()
  
  # For each bootstrap sample
  for (i in 1:M) {
    dat = dat_org
    set.seed(i) # Set seed to get different bootstrap samples
    sample_ids = sample(1:N, N, replace = TRUE)
    sample = dat[sample_ids, ]
    sample$id = row.names(sample) # Change id to make sure it captures duplicates
    boot_name = paste0("tree_MILC_", model_name, "_boot", i)
    
    # Perform LC with 3 classes
    boot_output_step1 = perform_lc(iteration, ind, cov, N, ME, dat = sample, folder = folder)[[2]]
    
    # Count combinations of indicators (+covariate) in the bootstrap sample
    count_boot = as.data.frame(boot_output_step1[, -which(colnames(boot_output_step1) %in% c("id", "cluster"))] %>% group_by_all() %>% summarise(COUNT = n()))
    count_boot = count_boot[, -ncol(count_boot)]
    
    # If not all combinations are present in the bootstrap sample
    if (nrow(count_boot) != nrow(count_dat)) {
      
      # Get data frame with parameters
      filename_lc = paste0(folder, "LC_", iteration, "-", ind, "-", cov, "-", N, "-", ME, "_script.lst")
      lc_lst = readChar(filename_lc, file.info(filename_lc)$size)
      lc_lst = strsplit(lc_lst, split = "Regression Parameters")
      lc_lst = strsplit(lc_lst[[1]][2], split = "Paired Comparisons")
      parameters_path = paste0(folder, "tree_MILC_", model_name, "_boot", i, "_parameters.dat")
      writeLines(lc_lst[[1]][1], parameters_path)
      parameters = suppressWarnings(fread(parameters_path, sep = "\t", dec = ",")[, 1:4])
      names(parameters) = c("term1", "term2", "term3", "coef")
      parameters$coef = as.numeric(parameters$coef)
      remove = which(grepl("(1)", parameters$term1, fixed = TRUE) | grepl("(1)", parameters$term3, fixed = TRUE)) # Remove dummies that are 0
      parameters = parameters[-remove, ]
      
      # Estimate extra LC model with obtained parameters as starting values
      output_path = paste0(folder, "tree_MILC_", model_name, "_boot", i, "_dat_org_posteriors.dat")
      script = generate_script_treeMILC_step1("LC", ind, cov, N, ME, model_name, filepath_input = dat_org_path, filepath_output = output_path, parameters)
      script_path = paste0(folder, "tree_MILC_", model_name, "_boot", i, "_dat_org_posteriors.lgs")
      writeLines(script, script_path)
      shell(paste0('"C:/Program Files/LatentGOLDnet6.0/lg60.exe" ', script_path, ' /b'))
      
      model_lst = paste(readLines(script_path), collapse = "\n")
      if (grepl("WARNING", model_lst, fixed = TRUE)) {
        stop("Error in Posterior Model")
      }
      
      # Read output of extra LC model that contains posterior probabilities for every observation in the original data set
      dat = as.data.frame(fread(output_path, sep = "\t", dec = ","))
      setnames(dat, old = c("Cluster#1", "Cluster#2", "Cluster#3", "Cluster#"), new = c("p1", "p2", "p3", "cluster"))
      dat = dat[, -which(colnames(dat) == "cluster")]
      
    } else {
      # If all combinations are present in the bootstrap sample, add the obtained posterior probabilities to the observations in the original data set
      dat = left_join(x = dat, y = count_boot, by = by_vector)
    }
    
    # Sample from obtained posterior membership probabilities
    dat$imp1 = apply(dat, 1, impute_value_step1)
    
    # Create subsets of cases with imputed values of 1 and write to file
    subset_ids = dat[dat$imp1 == 1, ]$id
    subset = dat[dat$id %in% subset_ids, ]
    filepath_subset = paste0(folder, "tree_MILC_", model_name, "_boot", i, "_subset.dat")
    fwrite(subset, file = filepath_subset, sep = "\t")
    
    # Create and execute Latent Gold script for the second model
    filepath_output = paste0(folder, boot_name, "_step2_output.dat")
    script = generate_script("treeMILC", ind, cov, N, ME, model_name, filepath_subset, filepath_output)
    script_path = paste0(folder, boot_name, "_step2_script.lgs", sep = "")
    writeLines(script, script_path)
    shell(paste0('"C:/Program Files/LatentGOLDnet6.0/lg60.exe" ', script_path, ' /b'))
    
    # Import model output
    model2_output = as.data.frame(fread(paste0(folder, boot_name, "_step2_output.dat", sep = ""), dec = ","))
    
    # Check if model is valid (i.e. does not contain a warning)
    model_lst = paste(readLines(paste0(folder, boot_name, "_step2_script.lst")), collapse = "\n")
    if (grepl("WARNING", model_lst, fixed = TRUE)) {
      stop("Error in Model 2")
    }
    
    # Add the results from the second model to the results from the first model
    count_step2 = as.data.frame(model2_output[, -which(colnames(model2_output) %in% c("id", "Cluster#"))] 
                                %>% group_by_all() %>% summarise(COUNT = n()))
    count_step2 = count_step2[, -ncol(count_step2)]
    dat = left_join(x = dat, y = count_step2, by = c(by_vector))
    
    # Sample again from obtained posterior membership probabilities
    dat$imp2 = apply(dat, 1, impute_value_step2)
    
    # Combine imputations from both models 
    dat$cluster = dat$imp1
    dat[dat$imp1 == 1 & ((!is.na(dat$imp2) & dat$imp2 == 1)), ]$cluster = 1
    dat[dat$imp1 == 1 & ((!is.na(dat$imp2) & dat$imp2 == 2)), ]$cluster = 3
    
    # Select and rename columns
    new_names = c("step2_p1", "step2_p3")
    setnames(dat, old = c("Cluster#1", "Cluster#2"), new = new_names)
    dat = dat[, -which(colnames(dat) == "Cluster#")]
    
    # Fix cluster assignment (if necessary)
    dat = fix_cluster_bootstrap(dat)
    bootstrap_results = append(bootstrap_results, list(dat))
  }
  
  # Add results to final output list
  to_return = append(to_return, list(bootstrap_results))
  return(to_return)
}


##################################################################################################
## Compute measurement error matrices (one per indicator) for an LC, LCT or tree-MILC model     ##
## @param results (list): Results of an LC, LCT or tree-MILC model                              ##
## @returns (list): List of measurement error matrices (one per indicator)                      ##
##################################################################################################

get_ME = function(results){
  
  # For LC and LCT: Compute a ME probability matrix for each indicator using the get_ME_help function
  if (class(results[[2]]) != "list") {
    results = results[[2]] # Ignore first data frame with model information
    
    which_ind = results[, colnames(results) %in% c("Y1", "Y2", "Y3", "Y4")]
    to_return = list()
    for (i in 1:ncol(which_ind)) {
      to_return = append(to_return, list(get_ME_help(results, which_ind[, i])))
    }
    return(to_return)
  }
  
  # For tree-MILC: Compute ME matrix for each indicator based on the imputed values
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
## Helpfunction to compute s measurement error matrix for one indicator                         ##
## @param results (data.frame): Data frame with posterior probabilities from one LC or LCT model##  
## @param ind_vec (vector): Vector containing posterior probabilities for one indicator         ##
## @returns (matrix): Measurement error matrix for one indicator                                ##
##################################################################################################

get_ME_help = function(results, ind_vec) {
  
  results$yx = ind_vec
  indicator_matrix = data.frame()
  
  for (i in 1:3) {
    # Compute the sum of p1, p2, and p3 for each indicator value
    indicator_matrix = rbind(indicator_matrix, sum(results[results$yx == i, ]$p1))
    indicator_matrix = rbind(indicator_matrix, sum(results[results$yx == i, ]$p2))
    indicator_matrix = rbind(indicator_matrix, sum(results[results$yx == i, ]$p3))
  }
  
  # Create a matrix and normalise the rows
  indicator_matrix = matrix(as.vector(indicator_matrix[, 1]), nrow = 3, ncol = 3, byrow = FALSE)
  indicator_matrix = prop.table(indicator_matrix, margin = 1)
  
  return(indicator_matrix)
}

