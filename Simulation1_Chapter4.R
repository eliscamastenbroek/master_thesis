
##################################################################################################
## Simulate a data set                                                                          ##
## @param seed (int): Seed                                                                      ##
## @param ME (int): Amount of measurement error (1=10%, 2=20%, 3=30%, 4=realistic)              ##
## @param folder (string): Folder to save files in                                              ##
## @returns (data.frame): A simulated data set (n=10,000)                                       ##
##################################################################################################

simulate_data = function(seed, ME, folder){
  
  #Create data structure file required by Latent GOLD
  filepath_input = paste0(folder,"exampleData.dat")
  exampleData = paste0("id q Y1 Y2 Y3 Y4 n
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
         caseweight n;
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
## Perform simulation study 1 (without missing covariates)                                      ##
##################################################################################################

## Specification of simulation conditions
ind = c(2:4)
N = c(1000,10000)
ME = c(1:4)
iteration = c(1:50)

## Create lists to store results in
LC_models = list()
LCT_models = list()
treeMILC_models = list()

## Execute simulations
for(l in iteration){
  for(k in ME){
    for(i in ind){
      
      #Add covariate if the number of indicators is two
      if(i==2){
        cov_ok="q"
      } else{
        cov_ok=NULL
      }
      
      for(j in N){
        #Create model name
        m = m + 1
        name = paste0(i,"-",cov_ok,"-",j,"-",k)
        row=c(l,i,cov,j,k,name)
        print(paste0("-----------------------------"))
        
        ## LC
        LC_models = append(LC_models,list(perform_lc(l,i,cov_ok,cov_problem=NULL,j,k,folder="F:\\Documents\\Thesis\\Simulatie\\Simulatie_1_zonder_cov\\LC\\")))
        print(paste0("LC model ",name," complete."))
        
        ## LCT
        LCT_models = append(LCT_models,list(perform_lct(l,i,cov_ok,cov_problem=NULL,j,k,folder="F:\\Documents\\Thesis\\Simulatie\\Simulatie_1_zonder_cov\\LCT\\")))
        print(paste0("LCT model ",name," complete."))
         
        ## tree-MILC
        treeMILC_models = append(treeMILC_models,list(perform_treeMILC(l,i,cov_ok,cov_problem=NULL,j,k,folder="F:\\Documents\\Thesis\\Simulatie\\Simulatie_1_zonder_cov\\treeMILC\\")))
        print(paste0("tree-MILC model ",name," complete."))
      }
    }
    
    ## Remove data set from global environment if not needed anymore   
    if(exists(paste0("simDat",k,"_iteration",l))){
      rm(list=paste0("simDat",k,"_iteration",l))
    }
  }
}

##################################################################################################
## Get results of the simulation study                                                          ##
##################################################################################################

## True proportions in simulated data
true_proportions = c(0.6061153, 0.2577143, 0.1361704)
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
LC_summary = get_summary(LC_results)
LCT_summary = get_summary(LCT_results)
treeMILC_summary = get_summary(treeMILC_results)

#Merge summaries
all_joined = rbind(LC_summary,LCT_summary,treeMILC_summary)


