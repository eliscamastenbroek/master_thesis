
################################## Plot_Simulation_1_and_2.R ##################################### 
## This file contains the code that is required to create the plots for the simulation study    ##
## in Chapter 4. The data is prepared for plotting in the second half. The plots are created    ##
## in the second half.                                                                          ##
##                                                                                              ##
## Note that the prepared data is already made available in the .RData files. This is because   ##
## not all objects that are required for preparing the data (i.e. the objects 'LC_models',      ##
## 'LCT_models' and 'treeMILC_models' could be uploaded due to storage limits in GitHub. As a   ## 
## result, creating the measurement error (ME) heatmaps is only possible when these objects are ##
## created manually using the 'Simulation.R' file.
##################################################################################################

# Load packages
library(ggplot2)
library(RColorBrewer)
library(scales)

# Set working directory 
setwd("add working directory here")

# Load data from either simulation study 1 (Chapter 4) or simulation study 2 (Chapter 5)
load("Simulation1.RData")
load("Simulation2.RData")

# Define which colours to use in the plots
colours <- c(LC = "#F8766D", LCT = "#00BA38", "tree-MILC" ="#619CFF") 

##################################################################################################
## Prepare data to plot population proportion estimates (PPEs) (i.e. expected value and RMSE)   ##                                                                    ##
##################################################################################################

# Create data frames for each contract type
perm <- all_results[,c(names(all_results)[1:6],"prop1","rmse_prop1","sd_prop1","type","ME1.1","ME2.1","sd_ME1.1","sd_ME2.1","diag","entropy")]
other <- all_results[,c(names(all_results)[1:6],"prop2","rmse_prop2","sd_prop2","type","ME1.2","ME2.2","sd_ME1.2","sd_ME2.2","diag","entropy")]
flex <- all_results[,c(names(all_results)[1:6],"prop3","rmse_prop3","sd_prop3","type","ME1.3","ME2.3","sd_ME1.3","sd_ME2.3","diag","entropy")]
names(perm) <- names(flex) <- names(other) <- c(names(all_results)[1:6],"prop","rmse","sd","type","ME1","ME2","sd_ME1","sd_ME2","diag","entropy")
perm$Contract <- "Permanent"
other$Contract <- "Other"
flex$Contract <- "Flexible"

# Add true proportions
perm$hline <- true_proportions[1]
other$hline <- true_proportions[2]
flex$hline <- true_proportions[3]

# Combine data frames 
plot_df <- rbind(perm, flex, other)
rm(perm, flex, other)
plot_df$Contract <- as.factor(plot_df$Contract)
plot_df$Contract <- factor(plot_df$Contract, levels = c("Permanent", "Flexible", "Other"))

# Get variables into the right format
plot_df$ME <- as.factor(plot_df$ME)
plot_df$ME <- relevel(plot_df$ME, 4)
levels(plot_df$ME) <- c("Realistic 7%", "10%", "20%", "30%")

# Add true ME probabilities (10%, 20%, and 30%)
plot_df$mline <- NA
plot_df[plot_df$Contract == "Permanent" & plot_df$ME == "10%", ]$mline <- ME_matrix1[1, 1] 
plot_df[plot_df$Contract == "Other" & plot_df$ME == "10%", ]$mline <- ME_matrix1[2, 2] 
plot_df[plot_df$Contract == "Flexible" & plot_df$ME == "10%", ]$mline <- ME_matrix1[3, 3] 
plot_df[plot_df$Contract == "Permanent" & plot_df$ME == "20%", ]$mline <- ME_matrix2[1, 1] 
plot_df[plot_df$Contract == "Other" & plot_df$ME == "20%", ]$mline <- ME_matrix2[2, 2] 
plot_df[plot_df$Contract == "Flexible" & plot_df$ME == "20%", ]$mline <- ME_matrix2[3, 3] 
plot_df[plot_df$Contract == "Permanent" & plot_df$ME == "30%", ]$mline <- ME_matrix3[1, 1] 
plot_df[plot_df$Contract == "Other" & plot_df$ME == "30%", ]$mline <- ME_matrix3[2, 2] 
plot_df[plot_df$Contract == "Flexible" & plot_df$ME == "30%", ]$mline <- ME_matrix3[3, 3]

##################################################################################################
## Prepare data to plot expected value of measurement error probability estimates (MEPEs)       ##
##################################################################################################

# Create data frame that contains the expected value of MEPEs for Y1 and Y3 in conditions with a realistic 7% ME
plot_rmse_ME2 <- plot_df[!is.na(plot_df$ME2), c("indicator","cov_ok","cov_problem","N","ME","id","type","ME2","sd_ME2","Contract","mline")]
levels(plot_rmse_ME2$ME)[1] <- "Real. 7% (ind. 2)"
plot_rmse_ME2[plot_rmse_ME2$Contract == "Permanent" & plot_rmse_ME2$ME == "Real. 7% (ind. 2)", ]$mline <- ME_matrix4b[1, 1] 
plot_rmse_ME2[plot_rmse_ME2$Contract == "Other" & plot_rmse_ME2$ME == "Real. 7% (ind. 2)", ]$mline <- ME_matrix4b[2, 2] 
plot_rmse_ME2[plot_rmse_ME2$Contract == "Flexible" & plot_rmse_ME2$ME == "Real. 7% (ind. 2)", ]$mline <- ME_matrix4b[3, 3] 
names(plot_rmse_ME2)[8:9] <- c("ME1","sd_ME1")

# Create data frame that contains the expected value of MEPEs for Y2 and Y4 in conditions with a realistic 7% ME
plot_rmse_ME1 <- plot_df[, c("indicator","cov_ok","cov_problem","N","ME","id","type","ME1","sd_ME1","Contract","mline")]
levels(plot_rmse_ME1$ME)[1] <- "Real. 7% (ind. 1)"
plot_rmse_ME1[plot_rmse_ME1$Contract == "Permanent" & plot_rmse_ME1$ME == "Real. 7% (ind. 1)", ]$mline <- ME_matrix4a[1, 1] 
plot_rmse_ME1[plot_rmse_ME1$Contract == "Other" & plot_rmse_ME1$ME == "Real. 7% (ind. 1)", ]$mline <- ME_matrix4a[2, 2] 
plot_rmse_ME1[plot_rmse_ME1$Contract == "Flexible" & plot_rmse_ME1$ME == "Real. 7% (ind. 1)", ]$mline <- ME_matrix4a[3, 3] 

# Combine data frames
plot_ME <- rbind(plot_rmse_ME2, plot_rmse_ME1)
rm(plot_rmse_ME1, plot_rmse_ME2)
plot_ME$ME <- as.factor(plot_ME$ME)
plot_ME$ME <- relevel(plot_ME$ME, "Real. 7% (ind. 2)")
plot_ME$ME <- relevel(plot_ME$ME, "Real. 7% (ind. 1)")

##################################################################################################
## Prepare data to plot RMSE of measurement error probability estimators (MEPEs)                ##
##################################################################################################

# Compute RMSE of ME probability estimators
LC_ME_rmse <- get_rmse_ME(LC_models, LC_results)
LCT_ME_rmse <- get_rmse_ME(LCT_models, LCT_results)
treeMILC_ME_rmse <- get_rmse_ME(treeMILC_models, treeMILC_results)
LC_ME_rmse$rmse <- as.numeric(LC_ME_rmse$rmse)
LCT_ME_rmse$rmse <- as.numeric(LCT_ME_rmse$rmse)
treeMILC_ME_rmse$rmse <- as.numeric(treeMILC_ME_rmse$rmse)
LC_ME_rmse$type <- "LC"
LCT_ME_rmse$type <- "LCT"
treeMILC_ME_rmse$type <- "tree-MILC"

# Combine results
plot_rmse_ME <- rbind(LC_ME_rmse, LCT_ME_rmse, treeMILC_ME_rmse)
plot_rmse_ME$ME <- as.numeric(plot_rmse_ME$ME)
plot_rmse_ME$ind <- as.numeric(plot_rmse_ME$ind)

# Combine indicators Y1/Y3 and Y2/Y4
plot_rmse_ME[plot_rmse_ME$ME == 4 & (plot_rmse_ME$ind == 1 | plot_rmse_ME$ind == 3), ]$ind <- 1
plot_rmse_ME[plot_rmse_ME$ME == 4 & (plot_rmse_ME$ind == 2 | plot_rmse_ME$ind == 4), ]$ind <- 2

# Compute mean RMSE over all indicators
plot_rmse_ME_not4 <- plot_rmse_ME[plot_rmse_ME$ME != 4, ] 
plot_rmse_ME_not4 <- as.data.frame(plot_rmse_ME_not4 %>% group_by(indicator, covariate, N, ME, id, Contract, type) %>% summarise_at(c("rmse"), mean))
plot_rmse_ME_not4$ind <- NA
plot_rmse_ME_4 <- plot_rmse_ME[plot_rmse_ME$ME == 4, ] 
plot_rmse_ME_4 <- as.data.frame(plot_rmse_ME_4 %>% group_by(indicator, covariate, N, ME, id, ind, Contract, type) %>% summarise_at(c("rmse"), mean))
plot_rmse_ME <- rbind(plot_rmse_ME_4, plot_rmse_ME_not4)
rm(plot_rmse_ME_4, plot_rmse_ME_not4)

# Get variables into the right format
plot_rmse_ME$Contract <- as.factor(plot_rmse_ME$Contract)
plot_rmse_ME$ME <- as.character(plot_rmse_ME$ME)
plot_rmse_ME[plot_rmse_ME$ME == 4 & plot_rmse_ME$ind == 1, ]$ME <- "Real. 7% (ind. 1)"
plot_rmse_ME[plot_rmse_ME$ME == 4 & plot_rmse_ME$ind == 2, ]$ME <- "Real. 7% (ind. 2)"
plot_rmse_ME$ME <- as.factor(plot_rmse_ME$ME)
plot_rmse_ME$ME <- relevel(plot_rmse_ME$ME, "Real. 7% (ind. 2)")
plot_rmse_ME$ME <- relevel(plot_rmse_ME$ME, "Real. 7% (ind. 1)")
levels(plot_rmse_ME$ME)[3:5] <- c("10%", "20%", "30%")
levels(plot_rmse_ME$Contract) <- c("Permanent", "Other", "Flexible")
plot_rmse_ME$Contract <- factor(plot_rmse_ME$Contract, levels = c("Permanent", "Flexible", "Other"))

##################################################################################################
## Create plots for PPEs (i.e. expected value and RMSE)                                         ##
##################################################################################################

# Expected value (all)
plot_df %>% ggplot(aes(x=indicator, y=prop, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="Model", values=colours) +
  labs(y= "Expected value", x = "Number of indicators") +
  facet_grid(Contract + N ~ ME, labeller=label_both) +
  geom_hline(aes(yintercept=hline), color="black", linewidth=0.3) +
  geom_errorbar(aes(ymin=prop-sd, ymax=prop+sd), width=.25, size=0.3, position=position_dodge(.9)) +
  theme(legend.position="top")
  
# RMSE (all)
plot_df %>% ggplot(aes(x=indicator, y=rmse, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="Model", values=colours) +
  labs(y= "RMSE", x = "Number of indicators") +
  facet_grid(Contract + N ~ ME, labeller=labeller(.rows=label_value, .cols=label_both)) +
  theme(legend.position="top")
  
# Expected value (permanent and n=1,000)
plot_df[plot_df$N==1000 & plot_df$Contract=="Permanent",] %>% ggplot(aes(x=indicator, y=prop, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="", values=colours) +
  labs(y= "Expected\nvalue", x = "Number of indicators") +
  facet_wrap(~ME, ncol=4, labeller=label_value) +
  geom_hline(aes(yintercept=hline), color="black", size=0.3) +
  geom_errorbar(aes(ymin=prop-sd, ymax=prop+sd), width=.25, size=0.3, position=position_dodge(.9)) +
  guides(fill = guide_legend(direction = "horizontal")) +
  theme(plot.title = element_text(size = 13), legend.justification="right")
  
# RMSE (permanent and n=1,000)
plot_df[plot_df$N==1000 & plot_df$Contract=="Permanent",] %>% ggplot(aes(x=indicator, y=rmse, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="", values=colours) +
  labs(y= "Expected\nvalue", x = "Number of indicators") +
  facet_wrap(~ME, ncol=4, labeller=label_value) +
  theme(legend.position="top") +
  guides(fill = guide_legend(direction = "horizontal")) +
  theme(plot.title = element_text(size = 13), legend.justification="right")

##################################################################################################
## Create plots for MEPEs (i.e. expected value and RMSE)                                        ##
##################################################################################################

# Expected value (all)
plot_ME %>% ggplot(aes(x=indicator, y=ME1, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="Model", values=colours) +
  labs(y= "(Mean) expected value", x = "Number of indicators") +
  geom_errorbar(aes(ymin=ME1-sd_ME1, ymax=ME1+sd_ME1), width=.25, size=0.3, position=position_dodge(.9)) +
  facet_grid(Contract + N ~ ME, labeller=label_both) +
  geom_hline(aes(yintercept=mline), color="black", size=0.3) +
  theme(legend.position="top")

# RMSE (all)
plot_rmse_ME %>% ggplot(aes(x=indicator, y=rmse, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="Model", values=colours) +
  labs(y= "RMSE", x = "Number of indicators") +
  facet_grid(Contract + N ~ ME, labeller=labeller(.rows=label_both, .cols=label_both)) +
  theme(legend.position="top")

# Expected value (permanent and n=1,000)
plot_ME[plot_ME$N==1000 & plot_ME$Contract=="Permanent",] %>% ggplot(aes(x=indicator, y=ME1, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="", values=colours) +
  labs(y= "(Mean) expected value", x = "Number of indicators") +
  theme(legend.position=c(1, 1.65), legend.justification="right") +
  geom_errorbar(aes(ymin=ME1-sd_ME1, ymax=ME1+sd_ME1), width=.25, size=0.3, position=position_dodge(.9)) +
  guides(fill = guide_legend(direction = "horizontal")) +
  facet_wrap(~ME, ncol=5, labeller=label_value) +
  geom_hline(aes(yintercept=mline), color="black", size=0.3)

# RMSE (permanent and n=1,000)
plot_rmse_ME[plot_rmse_ME$N==1000 & plot_rmse_ME$Contract=="Permanent",] %>% ggplot(aes(x=indicator, y=rmse, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="Model", values=colours) +
  labs(y= "(Mean) RMSE", x = "Number of indicators") +
  guides(fill = guide_legend(direction = "horizontal")) +
  theme(plot.title = element_text(size = 13), legend.position="top") +
  facet_wrap(~ ME, ncol=5, labeller=label_value)

##################################################################################################
## Plot entropy R2                                                                              ##
##################################################################################################

#Entropy R2 (all)
plot_df[plot_df$type!="tree-MILC",] %>% ggplot(aes(x=indicator, y=entropy, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="Model", values=colours) +
  labs(y= "Mean entropy R-Squared", x = "Number of indicators") +
  facet_grid(rows = vars(N), cols=vars(ME), labeller=label_both) +
  theme(legend.position="top")

#Entropy R2 (n=1,000)
plot_df[plot_df$type!="tree-MILC" & plot_df$N==1000,] %>% ggplot(aes(x=indicator, y=entropy, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="Model", values=colours) +
  labs(y= "Mean entropy R-Squared", x = "Number of indicators") +
  facet_wrap(~ME, ncol=4, labeller=label_both) +
  theme(legend.position="top")

##################################################################################################
## Plot mean summed bias                                                                        ##
##################################################################################################

plot_df %>% ggplot(aes(x=indicator, y=diag, fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(name="", values=colours) +
  labs(y= "Mean summed bias", x = "Number of indicators") +
  facet_wrap(~ ME, ncol=5, labeller=label_value) +
  theme(legend.position="top")

##################################################################################################
## Create plots for heatmaps                                                                    ##
##################################################################################################

# Simulation study 1
ind <- 2:4
cov_problem <- "null"

# Simulation study 2
ind <- 2:3
cov_problem <- c("baanduur", "SBIgroep", "")

for (j in c(1000, 10000)) {
  for (m in cov_problem) {
    for (k in ind) {
      if (k == 2) {
        cov_ok <- "q"
      } else {
        cov_ok <- NULL
      }
      # Plot bias and variance for 10%, 20%, and 30% ME
      bias1 <- get_ME_heatmap(LC_models, LCT_models, treeMILC_models, LC_results, 1, k, cov_ok, cov_problem, j)
      variance1 <- get_ME_heatmap(LC_models, LCT_models, treeMILC_models, LC_results, 2, k, cov_ok, cov_problem, j)
      
      # Plot bias and variance for realistic 7% ME
      bias2 <- get_ME_heatmap_realistic(LC_models, LCT_models, treeMILC_models, LC_results, 1, k, cov_ok, cov_problem, j)
      variance2 <- get_ME_heatmap_realistic(LC_models, LCT_models, treeMILC_models, LC_results, 1, k, cov_ok, cov_problem, j)

    }  
  }
}
