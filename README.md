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
1. Load the functions in the file '[Helpfunctions_General.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_General.R)', '[Helpfunctions_SimulatedData.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_SimulatedData.R)', and '[Methods_BestApproach.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Methods_BestApproach.R)', and '[Helpfunctions_Performance_Measures_and_Plots.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_Performance_Measures_and_Plots.R)' into the global environment.
2. Open the file '[Simulation1.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Simulation/Simulation1.R)'.
3. Change the working directory in this file.
4. Change the folder = "" in lines 177, 182, and 187 to where you would like to store the model results.
5. Execute the rest of the code. Note that this takes approximately 5 days.

# 2. Instructions to reproduce Simulation study 2 (Chapter )
To reproduce the first simulation study, the functions in the following files should be loaded into the global environment.



