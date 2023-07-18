tm_var <- function(results){  
  results = results[[2]]   #Ignore first data frame with model information
  
  n = nrow(results[[2]])
  prop_per_bootstrap = list()
  for(i in 1:length(results)){
    prop_per_bootstrap = append(prop_per_bootstrap,list(summary(factor(results[[i]]$cluster))/nrow(results[[i]])))
  }
  
  #Pool results
  pooled_proportions = bind_rows(prop_per_bootstrap) 
  pooled_proportions = as.data.frame(pooled_proportions)
  pooled_proportions = split(pooled_proportions, seq(5))
  
  res_var = lapply(pooled_proportions, fun_var, nsize = n) # within variance
  tvarmat = colMeans(bind_rows(res_var)) + apply(as.matrix(bind_rows(pooled_proportions)), 2, var) * (1 + 1/5)
  
  return(tvarmat)
  
}
