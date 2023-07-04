#Load required packages
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)
library(ExPosition)
library(gridExtra)
library(scales)
library(ggplot2)
library(RColorBrewer)

##################################################################################################
## Function to create a matrix with the measurement error probabilities used to simulate data   ##
## @param a2,...b32 (int): Logit parameters used to simulate the data                           ##
## @returns matrix (matrix): Matrix with measurement error probabilities                        ##
##################################################################################################
create_ME_matrix = function(a2, a3, b22, b33, b23, b32){
  row1 = c(1/(1+exp(a2)+exp(a3)), exp(a2)/(1+exp(a2)+exp(a3)), exp(a3)/(1+exp(a2)+exp(a3)))
  row2 = c(1/(1+exp(a2+b22)+exp(a3+b32)), exp(a2+b22)/(1+exp(a2+b22)+exp(a3+b32)), exp(a3+b32)/(1+exp(a2+b22)+exp(a3+b32)))
  row3 = c(1/(1+exp(a2+b23)+exp(a3+b33)), exp(a2+b23)/(1+exp(a2+b23)+exp(a3+b33)), exp(a3+b33)/(1+exp(a2+b23)+exp(a3+b33)))
  matrix = matrix(c(row1, row2, row3), nrow=3, ncol=3, byrow=TRUE)
  return(matrix)
}

##################################################################################################
## Function to generate a data set                                                              ##
## @param seed (int): Seed for generating the data set                                          ##
## @param ME (int): The amount of measurement error (1=0.2, 2=0.3, 3=0.4, 4=realistic)          ##
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
## @param ME (int): Amount of measurement error (1=0.2, 2=0.3, 3=0.5, 4=realistic)              ##
## @returns (data.frame): A subset                                                              ##
##################################################################################################
create_subset = function(iteration, ind, cov, N, ME){
  
  #If a certain data set does not exist in the global environment, create it 
  if(!exists(paste0("simDat",ME,"_iteration",iteration))){
    assign(paste0("simDat",ME,"_iteration",iteration),simulate_data(iteration,ME),envir=globalenv())
  }
  
  #Get simulated data set
  data = get(paste0("simDat",ME,"_iteration",iteration))
  
  #Set seed to get the same data set for every model within each iteration
  set.seed(iteration) 
  select_cases = sample(1:nrow(data),N,replace=FALSE)
  
  #Remove some redundant columns
  if(cov ==1){
    subset = data[select_cases,c(2:(3+ind))]
  } else {
    subset = data[select_cases,c(2,4:(3+ind))]
  }
  
  return(subset)
}

##################################################################################################
## Function to generate a Latent Gold script                                                    ## 
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
  
  #Create a string where the restrictions are stored
  restrictions="\n"
  
  #Adjust equations depending on whether to include a covariate or not
  if(cov==1){
    dep_cov = "\n\tindependent q nominal;"
    latent_var_eq = "\tCluster <- (a) 1 + (b) q;\n"
    restrictions = paste0(restrictions,"\tb[1,1] ~= ",par[3,]$coef,";\n")
    restrictions = paste0(restrictions,"\tb[1,2] ~= ",par[4,]$coef,";\n")
  }
  else if(cov==0){
    dep_cov = ""
    latent_var_eq = "\tCluster <- (a) 1;\n"
  }
  
  restrictions = paste0(restrictions,"\ta[1,1] ~= ",par[1,]$coef,";\n")
  restrictions = paste0(restrictions,"\ta[1,2] ~= ",par[2,]$coef,";\n")
  
  restrictions = paste0(restrictions,"\tc[1,1] ~= ",par[4+cov+cov,]$coef,";\n")
  restrictions = paste0(restrictions,"\tc[1,2] ~= ",par[5+cov+cov,]$coef,";\n")
  restrictions = paste0(restrictions,"\td[1,1] ~= ",par[6+cov+cov,]$coef,";\n")
  restrictions = paste0(restrictions,"\td[1,2] ~= ",par[7+cov+cov,]$coef,";\n")
  restrictions = paste0(restrictions,"\td[1,3] ~= ",par[8+cov+cov,]$coef,";\n")
  restrictions = paste0(restrictions,"\td[1,4] ~= ",par[9+cov+cov,]$coef,";\n")
  
  restrictions = paste0(restrictions,"\te[1,1] ~= ",par[11+cov+cov,]$coef,";\n")
  restrictions = paste0(restrictions,"\te[1,2] ~= ",par[12+cov+cov,]$coef,";\n")
  restrictions = paste0(restrictions,"\tf[1,1] ~= ",par[13+cov+cov,]$coef,";\n")
  restrictions = paste0(restrictions,"\tf[1,2] ~= ",par[14+cov+cov,]$coef,";\n")
  restrictions = paste0(restrictions,"\tf[1,3] ~= ",par[15+cov+cov,]$coef,";\n")
  restrictions = paste0(restrictions,"\tf[1,4] ~= ",par[16+cov+cov,]$coef,";\n")
  
  
  if(ind>=3){
    restrictions = paste0(restrictions,"\tg[1,1] ~= ",par[18+cov+cov,]$coef,";\n")
    restrictions = paste0(restrictions,"\tg[1,2] ~= ",par[19+cov+cov,]$coef,";\n")
    restrictions = paste0(restrictions,"\th[1,1] ~= ",par[20+cov+cov,]$coef,";\n")
    restrictions = paste0(restrictions,"\th[1,2] ~= ",par[21+cov+cov,]$coef,";\n")
    restrictions = paste0(restrictions,"\th[1,3] ~= ",par[22+cov+cov,]$coef,";\n")
    restrictions = paste0(restrictions,"\th[1,4] ~= ",par[23+cov+cov,]$coef,";\n")
  } 
  if(ind==4){
    restrictions = paste0(restrictions,"\ti[1,1] ~= ",par[25+cov+cov,]$coef,";\n")
    restrictions = paste0(restrictions,"\ti[1,2] ~= ",par[26+cov+cov,]$coef,";\n")
    restrictions = paste0(restrictions,"\tj[1,1] ~= ",par[27+cov+cov,]$coef,";\n")
    restrictions = paste0(restrictions,"\tj[1,2] ~= ",par[28+cov+cov,]$coef,";\n")
    restrictions = paste0(restrictions,"\tj[1,3] ~= ",par[29+cov+cov,]$coef,";\n")
    restrictions = paste0(restrictions,"\tj[1,4] ~= ",par[30+cov+cov,]$coef,";\n")
  }
  
  #Combine all parts of the script
  script=paste0(script_part1,model_name,script_part2,script_part3,dep_ind,dep_cov,
                latent_var,latent_var_eq,dep_ind_eq,restrictions,"\nend model")
  return(script)
}

# t = generate_script_treeMILC("LC", 3,0,1000,1, "model_name", "F:/Documents/Thesis/Simulatie/Voorbeeld/simDat_onlyRemoved.dat", "filepath_output.dat",par=test.dat_withoutdummy)
# writeLines(t,"F:/Documents/Thesis/Simulatie/Voorbeeld/script.lgs")
# 







