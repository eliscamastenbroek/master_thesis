# Introduction

This GitHub repository contains the R code that is used for the MSc thesis 'Comparing the performance of LC, LCT and tree-MILC analysis to correct for measurement error in the Employment Register and the Labour Force Survey' by Elisca Mastenbroek [(link to thesis)](https://github.com/eliscamastenbroek/master_thesis/blob/main/MSc_thesis.pdf)).

In this README.md file, instructions are provided to reproduce the results of the two simulation studies (see Chapters 4-5) and the analyses of the real data (see Chapter 6). Note that the real data from the ER and the LFS is not publicly available.

# Software requirements
For this project, the following software was used:
- RStudio 2022.02.01
- R version 4.1.3
- Latent GOLD 6.0

To perform the analyses, the following R packages were used:
- data.table 1.14.8
- dplyr 1.1.2

To plot the results, the following R packages were used:
- x


# 1. Instructions to reproduce Simulation study 1 (Chapter 4)
To reproduce the first simulation study in R, the following steps should be taken:
1. Load the functions in the files '[Helpfunctions_General.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_General.R)', '[Helpfunctions_Simulations.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_Simulations.R)', and '[Methods_BestApproach.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Methods_BestApproach.R)', and '[Helpfunctions_Performance_Measures_and_Plots.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_Performance_Measures_and_Plots.R)' in R.
2. Open the file '[Perform_Simulation1.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Analyses/Perform_Simulation1.R)'.
3. Change the working directory in this file.
4. Change the argument _folder = ""_ in lines 177, 182, and 187 to where you would like to store the model results.
5. Execute the rest of the code. Note that this takes approximately 5 days.

Note that the functions _perform_lc_, _perform_lct_ and _perform_treeMILC_ in the files above can also be used to estiate individual models. For example, the following lines can be used to estimate LC, LCT and tree-MILC models with two indicators and the (non-missing) covariate "q" using a data set of n=1000 with 10% ME that is generated with the seed 1:
```{r}
#LC
perform_lc(iteration = 1, ind = 2, cov_ok = "q", cov_problem = NULL, N = 1000, ME = 1, folder="your_folder_here")

#LCT
perform_lct(iteration = 1, ind = 2, cov_ok = "q", cov_problem = NULL, N = 1000, ME = 1, folder="your_folder_here")

#tree-MILC
perform_treeMILC(iteration = 1, ind = 2, cov_ok = "q", cov_problem = NULL, N = 1000, ME = 1, folder="your_folder_here")
```

To plot the results of the first simulation study in R, the following steps should be taken:
1. Open the file '[Plot_Simulation_1_and_2.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Plots/Plot_Simulation_1_and_2.R)'.
2
   a. Execute the code in the first half to prepare the data for plotting. Note that this part requires that the simulation in the file '[Perform_Simulation1.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Analyses/Perform_Simulation1.R)' has finished running.
   b. The part in 2a can be avoided by directly loading the content of the file '[Simulation1_Reduced.RData](https://github.com/eliscamastenbroek/master_thesis/blob/main/RData/Simulation1_Reduced.RData)' in R.
3. Execute the code in the second half to create the desired plots. 

# 2. Instructions to reproduce Simulation study 2 using the best approach (Chapter 5)
To reproduce the second simulation study, the following steps should be taken:
1. X

To plot the results of the first simulation study in R, the following additional steps should be taken:
1. X
   
# 3. Instructions to reproduce Simulation study 2 using the less optimal approach (Section 5.1.1)
To reproduce the second simulation study, the following steps should be taken:
1. Load the functions in the files X into the global environment.
2. Follow steps 2-X as described above in 2. Instructions to reproduce Simulation study 2 using the best approach (Chapter 5).

# 4. Instructions to reproduce the analyses on real data fomr the ER and the LFS (Chapter 6)
