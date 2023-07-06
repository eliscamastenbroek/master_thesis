##################################### Simulation2.R ##############################################
## This file contains the code that is required to perform the simulation study with missing    ##
## covariates (i.e. direct effects and parameter restrictions) as described in Chapter 5.       ##
## The file is divided into three parts:                                                        ##
##   1. Function 'simulate_data' to simulate data for this specific simulation study.           ##
##   2. Perform the simulation study.                                                           ##
##   3. Get the results of the simulation study.                                                ## 
## To run the code, the functions in the files 'Methods_BestApproach.R' and                     ## 
## 'Methods_Helpfunction.R' are also required.                                                  ##
##################################################################################################

## Load required packages and initialisations
## Initialisations
library(dplyr)
library(data.table)
options(dplyr.summarise.inform = FALSE) # Ignore redundant warnings from dplyr 
setwd("F:/Documents/Thesis/Simulatie/Simulatie_8_met_twee_cov") # Set working directory

##################################################################################################
## Simulate a data set                                                                          ##
## @param seed (int): Seed                                                                      ##
## @param ME (int): Amount of measurement error (1=10%, 2=20%, 3=30%, 4=realistic)              ##
## @param folder (string): Folder to save files in                                              ##
## @returns (data.frame): A simulated data set (n=10,000)                                       ##
##################################################################################################

simulate_data = function(seed, ME){
  
  #Create Latent Gold script for LC
  filepath_input = paste0("exampleDat.dat")
  
  script_part1 = paste0("version = 6.0\ninfile '",filepath_input,"' \n\nmodel
    title 'simulation",ME,"';
    options
    algorithm
        tolerance=1e-08 emtolerance=0.01 emiteration=500 nriterations=500;
    startvalues
        seed=1 sets=100 tolerance=1e-05 iterations=100;
    montecarlo
        seed=1 replicates=500 tolerance=1e-008;
    bayes
        categorical=1 variances=1 latent=1 poisson=1;
    missing includeall;
    output       
    	parameters=first standarderrors profile reorderclasses iterationdetails;\n")
  
  #Specify parameter estimates
  contract = c(0.8252, -0.3192)
  q = c(-0.1506, -0.1584) 
  SBIgroep = c(-2.5392, 2.2426)
  baanduur = c(-4.6105, -3.4073)
  data_param = gsub(",","",toString(c(contract,q,SBIgroep,baanduur)))
  
  #Parameters for 10% measurement error
  if(ME==1){
    outfile_name = paste0("simDat1_iteration",seed,".dat")
    a2=a3=-log(18)  #Coefficients for measurement error matrix
    b22=b33=log(324)
    b32=b23=log(18)
    ME_coefs = c(a2,a3,b22,b32,b23,b33,"\n")
    ME_coefs = gsub(",","",toString(rep(ME_coefs,3)))
    parameters = paste("\n",data_param,"\n",ME_coefs,"}\nend model")
  }
  
  #Parameters for 20% measurement error
  else if(ME==2){
    outfile_name = paste0("simDat2_iteration",seed,".dat")
    a2=a3=-3*log(2) #Coefficients for measurement error matrix
    b22=b33=6*log(2)
    b32=b23=3*log(2)
    ME_coefs = c(a2,a3,b22,b32,b23,b33,"\n")
    ME_coefs = gsub(",","",toString(rep(ME_coefs,3)))
    parameters = paste("\n",data_param,"\n",ME_coefs,"}\nend model")
  }
  
  #Parameters for 30% measurement error
  else if(ME==3){
    outfile_name = paste0("simDat3_iteration",seed,".dat")
    a2=a3=-1.54045 #Coefficients for measurement error matrix
    b22=b33=3.0809
    b32=b23=1.54045
    ME_coefs = c(a2,a3,b22,b32,b23,b33,"\n")
    ME_coefs = gsub(",","",toString(rep(ME_coefs,3)))
    parameters = paste("\n",data_param,"\n",ME_coefs,"}\nend model")
  }
  
  #Parameters for realistic amount of measurement error
  else if(ME==4){
    outfile_name = paste0("simDat4_iteration",seed,".dat")
    Y1_a2= -4.4917
    Y1_a3= -4.94368
    Y1_b22=7.64123
    Y1_b32=4.55064
    Y1_b23=3.09876
    Y1_b33=5.6678
    Y2_a2=-6.14311
    Y2_a3=-2.63157
    Y2_b22=11.9482
    Y2_b32=1.53296
    Y2_b23=1.72427
    Y2_b33=5.03275
    ME_coefs_Y1 = c(Y1_a2,Y1_a3,Y1_b22,Y1_b32,Y1_b23,Y1_b33,"\n")
    ME_coefs_Y2 = c(Y2_a2,Y2_a3,Y2_b22,Y2_b32,Y2_b23,Y2_b33,"\n")
    ME_coefs = gsub(",","",toString(rep(c(ME_coefs_Y1,ME_coefs_Y2,ME_coefs_Y1),1))) 
    parameters = paste(data_param,"\n",ME_coefs,"}\nend model")
  }
  
  script_part2 = paste0("\toutfile '",outfile_name, "' simulation=1 seed=",seed,";
    variables
         caseid id;
         caseweight w;
         dependent Y1 nominal 3, Y2 nominal 3, Y3 nominal 3;
         independent q nominal, SBIgroep nominal, baanduur nominal;
         latent cluster nominal 3;
     equations
         cluster <- 1 + q + SBIgroep + baanduur;			
         Y1      <- 1 + cluster;	
         Y2      <- 1 + cluster;
         Y3      <- 1 + cluster;
{ ")
  
  #Combine parts of script
  script = paste0(script_part1,script_part2,parameters)
  writeLines(script, paste0("simDat",ME,"_iteration",seed,"_script.lgs"))
  
  #Execute Latent Gold script
  shell(paste0('"C:/Program Files/LatentGOLDnet6.0/lg60.exe" ', paste0("simDat",ME,"_iteration",seed,"_script.lgs"), ' /b'))
  
  #Import simulated data set
  simDat = as.data.frame(fread(outfile_name,dec=","))
  
  #To add extra 'problem covariate' category, get id of observations with contract 'other' in Y1
  id_other_Y1 = simDat[simDat$Y1==2,]$id   
  
  #Select a random 90% of these observations (like in the real data)
  id_add_extra_cat = sample(id_other_Y1,(0.9*length(id_other_Y1))) 
  
  #Find out how many categories both covariates have
  ncat1 = length(levels(factor(simDat[,which(names(simDat)=="baanduur")]))) 
  ncat2 = length(levels(factor(simDat[,which(names(simDat)=="SBIgroep")]))) 
  
  #Assign the selected observations to a new covariate category
  simDat[id_add_extra_cat,which(names(simDat)=="baanduur")] = ncat1 + 1  
  simDat[id_add_extra_cat,which(names(simDat)=="SBIgroep")] = ncat2 + 1   
  
  return(simDat)
}

#####################################################################################################
## 1. Perform simulation study 2                                                                   ##
#####################################################################################################

## Specification of simulation conditions
ind = c(2,3)
N = c(10000)
ME = c(1:4)
iteration = c(1)
cov_problem = c("NULL","baanduur","baanduur-SBIgroep")

## Create lists to store results in
LC_models = list()
LCT_models = list()
treeMILC_models = list()

m = 0 

#Keep track of potential errors
errors = 0
error_vec = vector()

## Execute simulations
for(l in iteration){
  for(k in ME){
    for(i in ind){
      if(i==2){
        cov_ok = "q"
      } else {
        cov_ok = NULL
      }
      
      for(j in N){
        for(n in cov_problem){
          m = m + 1
          name = paste0(i,"-",ifelse(is.null(cov_ok),"NULL",cov_ok),"-",n,"-",j,"-",k)
          print(paste0("---------------", m,"/",length(iteration)*length(cov_problem)*length(ind)*length(N)*length(ME) ,"------------------"))

          #Get parameters in right format
          if(n=="NULL"){
            n=NULL
          } else if(n=="baanduur-SBIgroep"){
            n=c("baanduur","SBIgroep")
          }
          
          #Perform LC, LCT and tree-MILC
          LC = perform_lc(l,i,cov_ok,n,j,k,folder="F:\\Documents\\Thesis\\Simulatie\\Simulatie_8_met_twee_cov\\LC\\")
          LC_models = append(LC_models,list(LC))
          print(paste0("LC model ",name," complete."))
          
          LCT = perform_lct(l,i,cov_ok,n,j,k,folder="F:\\Documents\\Thesis\\Simulatie\\Simulatie_8_met_twee_cov\\LCT\\")
          LCT_models = append(LCT_models,list(LCT))
          print(paste0("LCT model ",name," complete."))
          
          treeMILC = perform_treeMILC(l,i,cov_ok,n,j,k,folder="F:\\Documents\\Thesis\\Simulatie\\Simulatie_8_met_twee_cov\\treeMILC\\")
          treeMILC_models  = append(treeMILC_models ,list(treeMILC))
          
          #Keep track of errors
          if((LC[[3]]=="Error")|(LCT[[3]]=="Error")|(treeMILC[[3]]=="Error")){
            errors = errors + 1
            error_vec = c(error_vec,LC[[1]]$id)
          }
          print(paste0("treeMILC model ",name," complete."))
          print(paste0("Number of errors caught: ",errors, " in ",error_vec))
        }
      }
    }
    
    #If a certain data set does not exist, create it
    if(exists(paste0("simDat",k,"_iteration",l))){
      rm(list=paste0("simDat",k,"_iteration",l))
    }
  }
}

##################################################################################################
## Get results of the simulation study                                                          ##
##################################################################################################

## True proportions in simulated data
true_proportions = c(0.6178, 0.2473, 0.1349)
names(true_proportions) = c("Permanent","Other","Flexible")

## Create (true) ME matrices
ME_matrix1 = create_ME_matrix(-log(18),-log(18),log(324),log(324),log(18),log(18)) #10% ME
ME_matrix2 = create_ME_matrix(-3*log(2),-3*log(2),6*log(2),6*log(2),3*log(2),3*log(2)) #20% ME
ME_matrix3 = create_ME_matrix(-1.54045,-1.54045,3.0809,3.0809,1.54045,1.54045) #30% ME
ME_matrix4a = create_ME_matrix(-4.4917,-4.94368,7.64123,5.6678,3.09876,4.55064) #Realistic 7% ME
ME_matrix4b = create_ME_matrix(-6.14311,-2.63157,11.9482,5.03275,1.72427,1.53296) #Realistic 7% ME
ME_matrix4 = list(ME_matrix4a,ME_matrix4b) #Realistic 7% ME

#Get results
LC_results = get_results(LC_models)
LCT_results = get_results(LCT_models)
treeMILC_results = get_results(treeMILC_models)

#Get summary
LC_summary = get_summary("LC", LC_results)
LCT_summary = get_summary("LCT", LCT_results)
treeMILC_summary = get_summary("tree-MILC", treeMILC_results)

#Merge summaries
all_joined = rbind(LC_summary,LCT_summary,treeMILC_summary)


