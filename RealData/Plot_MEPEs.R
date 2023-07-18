########################################## Plot_MEPEs.R ########################################## 
## This file obtains the results for the measurement error probability estimates (MEPEs) of the ##
## real data analyses and displays the results (see Section 6.3-6.4).                           ##
##################################################################################################

#Load required packages
library(ggplot2)
library(RColorBrewer)
library(scales)
library(dplyr)

##################################################################################################
## Define measurement error probability matrices as estimated by the HMMs                       ##
##################################################################################################

ER.2016 <- matrix(c(0.932, 0.001, 0.067, 0.004, 0.996, 0, 0.107, 0.001, 0.892), ncol = 3, byrow = TRUE)
LFS.2016 <- matrix(c(0.985, 0.01, 0.005, 0.04, 0.934, 0.026, 0.332, 0.071, 0.597), ncol = 3, byrow = TRUE)
HMM_ME.2016 <- list(ER.2016, LFS.2016)

ER.2017 <- matrix(c(0.929, 0.001, 0.07, 0.001, 0.998, 0.001, 0.089, 0, 0.911), ncol = 3, byrow = TRUE)
LFS.2017 <- matrix(c(0.982, 0.014, 0.005, 0.042, 0.928, 0.03, 0.296, 0.071, 0.633), ncol = 3, byrow = TRUE)
HMM_ME.2017 <- list(ER.2017, LFS.2017)

ER.2018 <- matrix(c(0.931, 0.002, 0.068, 0.004, 0.996, 0, 0.083, 0.001, 0.916), ncol = 3, byrow = TRUE)
LFS.2018 <- matrix(c(0.982, 0.011, 0.007, 0.04, 0.933, 0.027, 0.302, 0.074, 0.623), ncol = 3, byrow = TRUE)
HMM_ME.2018 <- list(ER.2018, LFS.2018)

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

###################################################################################################
## Obtain ME probability estimates for the LC and the tree-MILC models                           ##                                                          ##
###################################################################################################

## 1=Permanent, 2=Other, 3=Flexible

# ME probability matrices for (best) approach with direct effects and parameter restrictions (Section 6.1-6.3)
get_ME(LC_2016_restricted)
get_ME(LC_2017_restricted)
get_ME(LC_2018_restricted)
get_ME(treeMILC_2016_restricted)
get_ME(treeMILC_2017_restricted)
get_ME(treeMILC_2018_restricted)

# ME probability matrices for additional approach without missing covariates (Section 6.4)
get_ME(LC_2016_ok)
get_ME(LC_2017_ok)
get_ME(LC_2018_ok)
get_ME(treeMILC_2016_ok)
get_ME(treeMILC_2017_ok)
get_ME(treeMILC_2018_ok)

# ME probability matrices for additional approach without HMM recodings (Section 6.4)
get_ME(LC_2016_recoded)
get_ME(LC_2017_recoded)
get_ME(LC_2018_recoded)
get_ME(treeMILC_2016_recoded)
get_ME(treeMILC_2017_recoded)
get_ME(treeMILC_2018_recoded)

###################################################################################################
## Plot differences in ME probability estimates                                                  ##
###################################################################################################

# Differences in ME probability estimates for (best) approach with direct effects and parameter restrictions (Section 6.1-6.3)
HM_LC_2016_restricted <- get_ME_heatmap(get_ME(LC_2016_restricted), HMM_ME.2016) 
HM_LC_2017_restricted <- get_ME_heatmap(get_ME(LC_2017_restricted), HMM_ME.2017) 
HM_LC_2018_restricted <- get_ME_heatmap(get_ME(LC_2018_restricted), HMM_ME.2018)
HM_treeMILC_2016_restricted <- get_ME_heatmap(get_ME(treeMILC_2016_restricted), HMM_ME.2016)  
HM_treeMILC_2017_restricted <- get_ME_heatmap(get_ME(treeMILC_2017_restricted), HMM_ME.2017) 
HM_treeMILC_2018_restricted <- get_ME_heatmap(get_ME(treeMILC_2018_restricted), HMM_ME.2018) 

# Differences in ME probability estimates for additional approach without missing covariates (Section 6.4)
HM_LC_2016_ok <- get_ME_heatmap(get_ME(LC_2016_ok), HMM_ME.2016)
HM_LC_2017_ok <- get_ME_heatmap(get_ME(LC_2017_ok), HMM_ME.2017) 
HM_LC_2018_ok <- get_ME_heatmap(get_ME(LC_2018_ok), HMM_ME.2018) 
HM_treeMILC_2016_ok <- get_ME_heatmap(get_ME(treeMILC_2016_ok), HMM_ME.2016) 
HM_treeMILC_2017_ok <- get_ME_heatmap(get_ME(treeMILC_2017_ok), HMM_ME.2017) 
HM_treeMILC_2018_ok <- get_ME_heatmap(get_ME(treeMILC_2018_ok), HMM_ME.2018)

# Differences in ME probability estimates for additional approach without HMM recodings (Section 6.4)
HM_LC_2016_recoded <- get_ME_heatmap(get_ME(LC_2016_recoded), HMM_ME.2016) 
HM_LC_2017_recoded <- get_ME_heatmap(get_ME(LC_2017_recoded), HMM_ME.2017)
HM_LC_2018_recoded <- get_ME_heatmap(get_ME(LC_2018_recoded), HMM_ME.2018) 
HM_treeMILC_2016_recoded <- get_ME_heatmap(get_ME(treeMILC_2016_recoded), HMM_ME.2016) 
HM_treeMILC_2017_recoded <- get_ME_heatmap(get_ME(treeMILC_2017_recoded), HMM_ME.2017)
HM_treeMILC_2018_recoded <- get_ME_heatmap(get_ME(treeMILC_2018_recoded), HMM_ME.2018)
