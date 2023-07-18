################################# Methods_Helpfunctions.R ######################################## 
## This file contains a number of (help) functions that are required to perform LC and          ##
## tree-MILC analysis on real data from the ER and the LFS (see Chapter 6) and to plot the      ##
## results. This file contains the following functions:                                         ##
##    - generate_script                                                                         ##
##    - generate_script_treeMILC_extra                                                          ##
##    - store_model_info                                                                        ##
##    - get_ME                                                                                  ##
##    - get_ME_heatmap                                                                          ##
##    - get_within_variance_treeMILC                                                            ##
##    - get_variance_treeMILC_PPEs                                                              ##
##################################################################################################

##################################################################################################
## Generate a Latent Gold script for to estimate an LC model. This function can be used to      ##  
## estimate regular LC models, as well as LC models in tree-MILC analysis.                      ##                         
## @param type (string): What type of model (e.g. "LC" for regular LC and "LCT" for LCT step 2) ## 
## @param cov_ok (vector): Vector with the names of non-missing covariates                      ##
## @param cov_problem (vector): Vector with the names of missing covariates                     ##
## @param model_name (string): Name of the model                                                ##
## @param dat (data.frame): Data set to perform the analysis on                                 ##
## @param filepath_input (string): Path to input data set                                       ##
## @param filepath_output (string): Path where model output should be stored                    ##
## @returns (string): A string containing a Latent Gold script                                  ##
##################################################################################################

generate_script <- function(type, cov_ok, cov_problem, model_name, dat, filepath_input, filepath_output) {
  
  # Create vectors with characters to use for restrictions later on
  pars <- c("aa", "bb", "cc", "dd", "ee", "ff", "gg", "hh", "ii", "jj")
  pars2 <- c("kk", "ll", "mm", "nn", "oo", "pp", "qq", "rr", "ss", "tt")
  pars_count <- 0
  pars2_count <- 0
  
  script_part1 <- paste0("version = 6.0\ninfile '", filepath_input, "' \n\nmodel title '")
  
  script_part2 <- paste0("';
    options
    algorithm
        tolerance=1e-08 emtolerance=0.01 emiterations=10000000 nriterations=10000000;
    startvalues
        seed=1 sets=100 tolerance=1e-05 iterations=100;
    bayes
        categorical=1 variances=1 latent=1 poisson=1;
    missing includeall;
    output
    	parameters=")
  
  script_part3 <- paste0("first standarderrors profile reorderclasses iterationdetails;
    	outfile '", filepath_output, "' classification keep=persnr;
    variables\n") 
  
  # Adjust the number of clusters depending on the type of analysis performed
  if (type == "LC") {
    latent_var <- paste0("\n\tlatent Cluster nominal 3;
    equations\n")
  } else {
    latent_var <- paste0("\n\tlatent Cluster nominal 2;
    equations\n")
  }
  
  dep_ind <- "\tdependent contract nominal, contractEBB nominal;"
  dep_ind_eq1 <- "\n\tcontract <- 1 | Cluster"
  dep_ind_eq2 <- "\n\tcontractEBB <- 1 | Cluster;"
  
  # Adjust equations and restrictions depending on the included covariates
  cov <- c(cov_ok, cov_problem)
  restrictions1 <- restrictions2 <- ""
  dep_cov <- "\n\tindependent "
  
  for (i in 1:length(cov)) {
    if (i == 1) {
      dep_cov <- paste0(dep_cov, " ", cov[i], " nominal")
    } else {
      dep_cov <- paste0(dep_cov, ", ", cov[i], " nominal")
    }
  }
  
  dep_cov <- paste0(dep_cov, ";")
  latent_var_eq <- "\tCluster <- 1"
  
  if (!is.null(cov_ok)) {
    for (i in 1:length(cov_ok)) {
      latent_var_eq <- paste(latent_var_eq, "+", cov_ok[i])
    }
  }
  
  # Specify restrictions
  if (!is.null(cov_problem)) {
    for (i in 1:length(cov_problem)) {
      pars_count <- pars_count + 1
      which_col <- which(colnames(dat) == cov_problem[i])
      num_cats <- length(levels(factor(dat[, which_col])))
      dep_ind_eq1 <- paste0(dep_ind_eq1, "+ (", pars[i], "~ful) 1 | ", cov_problem[i])
      latent_var_eq <- paste0(latent_var_eq, "+ (", pars2[i], ")", cov_problem[i])
      
      # Specify the first restriction (see Section 5.1)
      for (j in 1:num_cats) {
        if (j != num_cats) {
          restrictions2 <- paste0(restrictions2, "\n\t", pars[i], "[", j, ",] = 0;")
        } else{
          restrictions2 <- paste0(restrictions2, "\n\t", pars[i], "[", j, ",1] = -100;")
          restrictions2 <- paste0(restrictions2, "\n\t", pars[i], "[", j, ",2] = 0;")
          restrictions2 <- paste0(restrictions2, "\n\t", pars[i], "[", j, ",3] = -100;")
        }
      }   
      
      # Specify the second restriction (see Section 5.1)
      if (length(cov_problem) > 1) {
        par2 <- (num_cats - 1) * 2
        par1 <- par2 - 1
        if (type == "LC") {
          restrictions1 <- paste0(restrictions1, paste0("\n\t", pars2[i], "[1,", par1, "] = 0; ", pars2[i], "[1,", par2, "] = 0;"))
        } else {
          restrictions1 <- paste0(restrictions1, paste0("\n\t", pars2[i], "[1,", num_cats - 1, "] = 0;"))
        }
      }
    }
  }
  
  latent_var_eq <- paste0(latent_var_eq, ";")
  dep_ind_eq1 <- paste0(dep_ind_eq1, ";")
  
  # Combine all parts of the script
  script <- paste0(script_part1, model_name, script_part2, script_part3, dep_ind, dep_cov,
                   latent_var, latent_var_eq, dep_ind_eq1, dep_ind_eq2, restrictions1, restrictions2, "\nend model")
  return(script)
}

##################################################################################################
## Generate a Latent Gold script to estimate an extra model in tree-MILC analysis to obtain     ##
## posterior probabilities for response patterns that were not present in the                   ##
## bootstrap sample. This is done by estimating an LC model on the original data set while      ##
## using the logit parameters from the previous model as starting values and setting the number ## 
## of EM- and NR-iterations to 0.                                                               ## 
## @param cov_ok (vector): Vector with the names of non-missing covariates                      ##
## @param cov_problem (vector): Vector with the names of missing covariates                     ##
## @param model_name (string): Name of the model                                                ##
## @param dat (data.frame): Data set to perform the analysis on                                 ##
## @param filepath_input (string): Path to input data set                                       ##
## @param filepath_output (string): Path where model output should be stored                    ##
## @param par (data.frame): Data frame with logit parameters as estimated by the previous model ## 
## @returns (string): A string containing a Latent Gold script                                  ##
##################################################################################################

generate_script_treeMILC_extra <- function(cov_ok, cov_problem, model_name, dat, filepath_input, filepath_output, par){
    
    # Create vectors with characters to use for restrictions later on
    pars <- c("aa", "bb", "cc", "dd", "ee", "ff", "gg", "hh", "ii", "ll", "mm", "nn", "oo", "pp", "qq", "rr", "ss", "tt", "uu", "vv", "ww")
    pars_count  <- 0
    
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
    	outfile '", filepath_output, "' classification keep=persnr;
    variables\n") 
    
    latent_var <- paste0("\n\tlatent Cluster nominal 3;
    equations\n")
    
    # Adjust equations (and restrictions) depending on which covariates to include (if any)
    cov <- c(cov_ok, cov_problem)
    restrictions <- restrictions2 <- ""
    latent_var_eq <- "\tCluster <- (zz) 1"
    dep_cov <- ifelse(is.null(cov), "", "\n\tindependent ")
    dep_ind_eq1 <- "\n\tcontract <- (jj) 1 | Cluster"
    
    # Specify restrictions
    if(!is.null(cov_problem)){
      for(i in 1:length(cov_problem)){
        pars_count <- pars_count + 1
        which_col <- which(colnames(dat) == cov_problem[i])
        num_cats <- length(levels(factor(dat[, which_col])))
        dep_ind_eq1 <- paste0(dep_ind_eq1, " + (", pars[i], "~ful) 1 | ", cov_problem[i]) 
        
        for(j in 1:num_cats){
          if(j != num_cats){
            restrictions2 <- paste0(restrictions2, "\n\t", pars[i], "[", j, ",] = 0;")  #Specify second restriction (see Section 5.1)
          } else{
            restrictions2 <- paste0(restrictions2, "\n\t", pars[i], "[", j, ",1] = -100;")  #Specify first restriction (see Section 5.1)
            restrictions2 <- paste0(restrictions2, "\n\t", pars[i], "[", j, ",2] = 0;")
            restrictions2 <- paste0(restrictions2, "\n\t", pars[i], "[", j, ",3] = -100;")
          }
        }   
      }
    }
    
    if(!is.null(cov)){
      for(i in 1:length(cov)){
        if(i == 1){
          dep_cov <- paste0(dep_cov, " ", cov[i], " nominal")
        } else {
          dep_cov <- paste0(dep_cov, ",  ", cov[i], " nominal")
        }
      }
      dep_cov <- paste0(dep_cov, ";")
    }  
    
    # Specify start values
    restrictions <- paste0(restrictions, "\n\tzz[1,1] ~= ", par[1, ]$coef, ";\n")
    restrictions <- paste0(restrictions, "\tzz[1,2] ~= ", par[2, ]$coef, ";\n")
    
    df_counter <- 2
    
    if(!is.null(cov)){
      for(i in 1:length(cov)){
        which_col <- which(colnames(dat) == cov[i])
        num_cats <- length(levels(factor(dat[, which_col])))
        pars_count <- pars_count + 1
        latent_var_eq <- paste0(latent_var_eq, " + (", pars[pars_count], ") ", cov[i])
        par2 <- (num_cats - 1) * 2
        
        for(j in 1:par2){
          df_counter <- df_counter + 1
          if(par[df_counter, ]$coef == 0){
            restrictions <- paste0(restrictions, "\t", pars[pars_count], "[1,", j, "] = ", par[df_counter, ]$coef, ";\n")
          } else {
            restrictions <- paste0(restrictions, "\t", pars[pars_count], "[1,", j, "] ~= ", par[df_counter, ]$coef, ";\n")
          }
        }
      }
    }
    
    dep_ind <- "\tdependent contract nominal, contractEBB nominal;"
    dep_ind_eq2 <- "\n\tcontractEBB <- (kk) 1 | Cluster;"
    
    df_counter <- grep("contractEBB", par$term1)[1]
    restrictions <- paste0(restrictions, "\tkk[1,1] ~= ", par[df_counter, ]$coef, ";\n")
    restrictions <- paste0(restrictions, "\tkk[1,2] ~= ", par[df_counter+1, ]$coef, ";\n")
    restrictions <- paste0(restrictions, "\tkk[2,1] ~= ", par[df_counter+2, ]$coef, ";\n")
    restrictions <- paste0(restrictions, "\tkk[2,2] ~= ", par[df_counter+3, ]$coef, ";\n")
    restrictions <- paste0(restrictions, "\tkk[3,1] ~= ", par[df_counter+4, ]$coef, ";\n")
    restrictions <- paste0(restrictions, "\tkk[3,2] ~= ", par[df_counter+5, ]$coef, ";\n")
    
    df_counter <- grep("contract", par$term1)[1]
    restrictions <- paste0(restrictions, "\tjj[1,1] ~= ", par[df_counter, ]$coef, ";\n")
    restrictions <- paste0(restrictions, "\tjj[1,2] ~= ", par[df_counter+1, ]$coef, ";\n")
    restrictions <- paste0(restrictions, "\tjj[2,1] ~= ", par[df_counter+2, ]$coef, ";\n")
    restrictions <- paste0(restrictions, "\tjj[2,2] ~= ", par[df_counter+3, ]$coef, ";\n")
    restrictions <- paste0(restrictions, "\tjj[3,1] ~= ", par[df_counter+4, ]$coef, ";\n")
    restrictions <- paste0(restrictions, "\tjj[3,2] ~= ", par[df_counter+5, ]$coef, ";\n")
    
    latent_var_eq <- paste0(latent_var_eq, ";")
    dep_ind_eq1 <- paste0(dep_ind_eq1, ";")
    
    
    #Combine all parts of the script
    script=paste0(script_part1, model_name, script_part2, script_part3, dep_ind, dep_cov, 
                  latent_var, latent_var_eq, dep_ind_eq1, dep_ind_eq2, restrictions, restrictions2, "\nend model")
    return(script)
  }
  

##################################################################################################
## Get a list that consists of a data frame with model information                              ##
## @param cov_ok (vector): Vector with the names of non-missing covariates                      ##
## @param cov_problem (vector): Vector with the names of missing covariates                     ##
## @returns (list): A list that consists of:                                                    ##
##  [[1]]: A data frame with model information                                                  ##
##################################################################################################

store_model_info <- function(cov_ok, cov_problem) {
  
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
  model_name <- paste(cov_ok1, cov_problem1, sep = "-")
  model_info1 <- data.frame(cov_ok = cov_ok1, cov_problem = cov_problem1)
  model_info2 <- data.frame(id = model_name)
  model_info <- cbind(model_info1, model_info2)
  
  to_return <- list(model_info) 
  
  return(to_return)
}

##################################################################################################
## Compute measurement error matrices (one per indicator) for an LC, LCT or tree-MILC model     ##
## @param results (list): Results of an LC, LCT or tree-MILC model                              ##
## @returns (list): List of measurement error matrices (one per indicator)                      ##
##################################################################################################

get_ME <- function(results) {
  
  # For LC and LCT: Compute ME for each indicator using the get_ME_help function 
  if (class(results[[2]]) != "list") {
    results <- results[[2]] # Ignore first data frame with model information
    
    which_ind <- results[, colnames(results) %in% c("contract", "contractEBB")]
    to_return <- list()
    for (i in 1:ncol(which_ind)) {
      to_return <- append(to_return, list(get_ME_help(results, which_ind[, i])))
    }
    return(to_return)
  }
  
  # For tree-MILC: Compute ME matrix for each indicator based on the imputed values
  else {
    all_indicators <- c("contract", "contractEBB")    
    num_of_ind <- 2
    results <- results[[2]]
    to_return <- list()  # Create list to store the averages of all bootstrap sample matrices for each indicator 
    
    for (i in 1:num_of_ind) {
      cluster_index <- which(names(results[[1]]) == "cluster")
      ind_index <- which(names(results[[1]]) == all_indicators[i])
      summed_matrix <- prop.table(table(results[[1]][, cluster_index], results[[1]][, ind_index]), 1)
      for (j in 2:length(results)) {
        summed_matrix <- summed_matrix + prop.table(table(results[[j]][, cluster_index], results[[j]][, ind_index]), 1)
      }
      
      mean_matrix <- summed_matrix / length(results)
      to_return <- append(to_return, list(mean_matrix))
    }
    
    return(to_return)
  }
}

###################################################################################################
## This function creates a heatmap that shows the difference between the ME probability matrices ## 
## as estimated by an HMM and as estimated by an LC or a tree-MILC model.                        ##
## @param ME_estimated (list): List of ME probability matrices (one for each indicator) as       ## 
## estimated by an LC or a tree-MILC model.                                                      ##
## @param ME_HMM (list): List of ME probability matrices (one for each indicator) as             ## 
## estimated by an HMM.                                                                          ##    
###################################################################################################

get_ME_heatmap <- function(ME_estimated, ME_HMM) {
  
  # Compute the difference between the two matrices for each indicator
  diff_ER <- ME_estimated[[1]] - ME_HMM[[1]]
  diff_LFS <- ME_estimated[[2]] - ME_HMM[[2]]
  
  # Convert matrices to data frame (in long format)
  ER_to_plot <- as.data.frame(as.table(diff_ER))
  ER_to_plot$type <- "ER"
  LFS_to_plot <- as.data.frame(as.table(diff_LFS))
  LFS_to_plot$type <- "LFS"
  df_to_plot <- rbind(ER_to_plot, LFS_to_plot)
  colnames(df_to_plot) <- c("Model", "Indicator", "Difference", "type")
  df_to_plot$Model <- as.character(df_to_plot$Model)
  df_to_plot$Indicator <- as.character(df_to_plot$Indicator)
  
  # Rename (automatically assigned) values to P, F, and O (permanent, flexible, and other)
  df_to_plot[df_to_plot$Model %in% c("A", "1"), "Model"] <- "P"
  df_to_plot[df_to_plot$Model %in% c("B", "2"), "Model"] <- "O"
  df_to_plot[df_to_plot$Model %in% c("C", "3"), "Model"] <- "F"
  df_to_plot[df_to_plot$Indicator %in% c("A", "1"), "Indicator"] <- "P"
  df_to_plot[df_to_plot$Indicator %in% c("B", "2"), "Indicator"] <- "O"
  df_to_plot[df_to_plot$Indicator %in% c("C", "3"), "Indicator"] <- "F"
  
  # Get P, F, and O into the correct order
  df_to_plot$Model <- factor(df_to_plot$Model, levels = c("O", "F", "P"))
  df_to_plot$Indicator <- factor(df_to_plot$Indicator, levels = c("P", "F", "O"))
  df_to_plot$Difference <- round(df_to_plot$Difference, 3)
  
  # Define ranges for plots
  max_abs <- max(abs(df_to_plot$Difference))
  rng <- c(-0.15, 0.15)
  
  plot <- ggplot(df_to_plot, aes(Indicator, Model)) +
    geom_tile(aes(fill = Difference)) +
    geom_text(aes(label = Difference)) +
    scale_fill_gradientn(
      colors = rev(brewer.pal(9, "RdBu")),
      values = rescale(c(rng[1], 0, rng[2])),
      limits = c(rng[1], rng[2]),
      n.breaks = 8
    ) +
    facet_wrap(~type, labeller = label_value) +
    labs(fill = "Difference")
  
  return(plot)
}

##################################################################################################
## This is a help function to compute the within variance of tree-MILC estimates                ##
## @param prop_boot (numeric):  Estimate of interest based on a particular bootstrap sample     ##
## @nsize (int): Number of observations in the boostrap sample                                  ##
## @returns (nueric): The fraction in Equation (2.22)                                           ##
##################################################################################################

get_within_variance_treeMILC <- function(prop_boot, nsize) (prop_boot * (1 - prop_boot)) / nsize

##################################################################################################
## This is function to compute the variance of tree-MILC's PPEs                                 ##
## @param results (list): Results of one tree-MILC model                                        ##
## @returns (vector): A vector that contains the variance of tree-MILC's PPEs (for permanent,   ##
## flexible and other)                                                                          ##
##################################################################################################

get_variance_treeMILC_PPEs <- function(results){ 
  
  # Ignore first data frame with model information
  results <- results[[2]]   
  
  n <- nrow(results[[2]])
  prop_per_bootstrap <- list()
  
  for(i in 1 : length(results)){
    prop_per_bootstrap <- append(prop_per_bootstrap, list(summary(factor(results[[i]]$cluster)) / nrow(results[[i]])))
  }
  
  #Pool results
  pooled_proportions <- bind_rows(prop_per_bootstrap) 
  pooled_proportions <- as.data.frame(pooled_proportions)
  pooled_proportions <- split(pooled_proportions, seq(5))
  
  res_var <- lapply(pooled_proportions, fun_var, nsize = n) # within variance
  tvarmat <- colMeans(bind_rows(res_var)) + apply(as.matrix(bind_rows(pooled_proportions)), 2, var) * (1 + 1/5)
  
  return(tvarmat)
}
