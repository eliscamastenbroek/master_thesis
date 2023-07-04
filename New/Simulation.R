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

## Specification of simulation conditions
ind = c(2:4)
N = c(1000,10000)
ME = c(1:4)
iteration = c(150:200)

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
        cov=1
      } else{
        cov=0
      }
      
      for(j in N){
        #Create model name
        m = m + 1
        name = paste0(i,"-",cov,"-",j,"-",k)
        row=c(l,i,cov,j,k,name)
        print(paste0("-----------------------------"))
        
        ## LC
        LC_models = append(LC_models,list(perform_lc(l,i,cov,j,k,folder="F:\\Documents\\Thesis\\Simulatie\\Simulatie_1_zonder_cov\\LC\\")))
        print(paste0("LC model ",name," complete."))
        
        ## LCT
        LCT_models = append(LCT_models,list(perform_lct(l,i,cov,j,k,folder="F:\\Documents\\Thesis\\Simulatie\\Simulatie_1_zonder_cov\\LCT\\")))
        print(paste0("LCT model ",name," complete."))
         
        ## tree-MILC
        treeMILC_models = append(treeMILC_models,list(perform_treeMILC(l,i,cov,j,k,folder="F:\\Documents\\Thesis\\Simulatie\\Simulatie_1_zonder_cov\\treeMILC\\")))
        print(paste0("tree-MILC model ",name," complete."))
      }
    }
    
    ## Remove data set from global environment if not needed anymore   
    if(exists(paste0("simDat",k,"_iteration",l))){
      rm(list=paste0("simDat",k,"_iteration",l))
    }
  }
}
