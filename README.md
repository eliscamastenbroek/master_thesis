# Introduction

This GitHub repository contains the R code used in the MSc thesis 'Comparing the performance of LC, LCT and tree-MILC analysis to correct for measurement error in the Employment Register and the Labour Force Survey' by Elisca Mastenbroek [(link)](https://github.com/eliscamastenbroek/master_thesis/blob/main/MSc_thesis.pdf)).

In this README.md file, instructions are provided for reproducing the results of the two simulation studies (see Chapters 4-5) and the analyses of the real data (see Chapter 6). Note that the real data from the ER and the LFS is not publicly available.

# Software requirements
For this project, the following software was used:
- RStudio 2022.02.01
- R version 4.1.3
- Latent GOLD 6.0

To perform the analyses, the following R packages were used:
- x

To plot the results, the following R packages were used:
- x


# 1. Instructions to reproduce Simulation study 1 (Chapter 4)
To reproduce the first simulation study in R, the following steps should be taken:
1. Load the functions in the files '[Helpfunctions_General.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_General.R)', '[Helpfunctions_SimulatedData.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_SimulatedData.R)', and '[Methods_BestApproach.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Methods_BestApproach.R)', and '[Helpfunctions_Performance_Measures_and_Plots.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_Performance_Measures_and_Plots.R)' into your global environment.
2. Open the file '[Simulation1.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Simulation/Simulation1.R)'.
3. Change the working directory in this file.
4. Change the argument _folder = ""_ in lines 177, 182, and 187 to where you would like to store the model results.
5. Execute the rest of the code. Note that this takes approximately 5 days.

Note that the functions _perform_lc_, _perform_lct_ and _perform_treeMILC_ in the files above can also be used to create individual models. For example, the following lines create LC, LCT and tree-MILC models with two indicators and the (non-missing) covariate "q". Here, a data set is used of n=1000 with 10% ME that is generated using the seed 1:
```{r}
#LC
perform_lc(iteration = 1, ind = 2, cov_ok = "q", cov_problem = NULL, N = 1000, ME = 1, folder="your_folder_here")

#LCT
perform_lct(iteration = 1, ind = 2, cov_ok = "q", cov_problem = NULL, N = 1000, ME = 1, folder="your_folder_here")

#tree-MILC
perform_treeMILC(iteration = 1, ind = 2, cov_ok = "q", cov_problem = NULL, N = 1000, ME = 1, folder="your_folder_here")
```

To plot the results of the first simulation study in R, the following additional steps should be taken:
1. Open the file '[Plot_Simulation_1_and_2.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Simulation/Plot_Simulation_1_and_2.R)'.
2. Execute the code in the first half of this file to prepare the data for plotting.
3. Execute the code in the second half of this file to create the desired plots. 

# 2. Instructions to reproduce Simulation study 2 (Chapter 5)
To reproduce the second simulation study, the functions in the following files should be loaded into the global environment.



