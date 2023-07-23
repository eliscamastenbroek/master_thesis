# Introduction

This GitHub repository contains the R code that is used for the MSc thesis 'Correcting for measurement error in the
Employment Register (ER) and the Labour Force Survey (LFS) using LC, LCT and tree-MILC analysis' by Elisca Mastenbroek ([link to thesis](https://github.com/eliscamastenbroek/master_thesis/blob/main/MSc_thesis.pdf)).

In this README.md file, instructions are provided to reproduce the results of the two simulation studies (see Chapters 4-5) and the analyses of the real data (see Chapter 6). Note that the real data from the ER and the LFS is not publicly available.

General remark: The sample size is referred to as _n_ throughout the thesis, but as _N_ throughout the code. Similarly, the covariates _Z<sub>1</sub>_, _Z<sub>2</sub>_, and _Q_ in the thesis are referred to as _baanduur_, _SBIgroep_, and _q_ throughout the code.

# Software requirements
For this project, the following software is required:
| **Software** | **Version** |
|--------------|-------------|
| RStudio      | 2022.02.01  |
| R            | 4.1.3       |
| Latent GOLD  | 6.0         |

To perform the analyses and plot the results, the following R packages are required:
| **Package**  | **Version** |
|--------------|-------------|
| data.table   | 1.14.8      |
| dplyr        | 1.1.2       |
| ggplot2      | 3.4.2       |
| reshape2     | 1.4.4       |
| RColorBrewer | 1.1-3       |
| scales       | 1.2.1       |
| tidyverse    | 2.0.0       |

# 1. Instructions to reproduce and/or plot Simulation study 1 (Chapter 4)

### Instructions to reproduce the simulation study:
1. Load the functions in the files '[Helpfunctions_General.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_General.R)', '[Helpfunctions_Simulations.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_Simulations.R)', and '[Methods_Best_Approach.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Methods_Best_Approach.R)', and '[Helpfunctions_Performance_Measures_and_Plots.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_Performance_Measures_and_Plots.R)' into R.
2. Open the file '[Perform_Simulation_1.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Analyses/Perform_Simulation_1.R)'.
3. Change the working directory in line 24.
4. Change the argument _folder = ""_ in lines 213, 218, and 223 to where you would like to store the model results. Make sure to end the folder name with a "/".
5. Execute the rest of the code. Note that this takes approximately 5 days.

The functions _perform_lc_, _perform_lct_ and _perform_treeMILC_ in the files above can also be used to estimate individual models. For example, the following lines can be used to estimate LC, LCT and tree-MILC models with two indicators and the (non-missing) covariate "q" using a data set of n=1,000 with 10% ME that is generated with the seed 1:
```{r}
#LC
perform_lc(iteration = 1, ind = 2, cov_ok = "q", cov_problem = NULL, N = 1000, ME = 1, folder="your_folder_here/")

#LCT
perform_lct(iteration = 1, ind = 2, cov_ok = "q", cov_problem = NULL, N = 1000, ME = 1, folder="your_folder_here/")

#tree-MILC
perform_treeMILC(iteration = 1, ind = 2, cov_ok = "q", cov_problem = NULL, N = 1000, ME = 1, folder="your_folder_here/")
```

### Instructions to plot the results (except for the heatmaps):
1. Follow either step 1a **or** 1b:
   
   a. Run the code in the file '[Prepare_Simulation_Results_for_Plotting.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Plots/Prepare_Simulation_Results_for_Plotting.R)' to prepare the data for plotting. Note that this part requires that the simulation in the file '[Perform_Simulation_1.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Analyses/Perform_Simulation_1.R)' has been executed.
   
   b. Load the contents of the file '[Simulation_1_Reduced.RData](https://github.com/eliscamastenbroek/master_thesis/blob/main/RData/Simulation_1_Reduced.RData)' into R.
3. Run the code in the file '[Plot_Simulation_1.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Plots/Plot_Simulation_1.R)' to create the desired plots.

### Instructions to plot the heatmaps:
1. (If not already loaded): Load the functions in the files '[Helpfunctions_General.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_General.R)', '[Helpfunctions_Simulations.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_Simulations.R)', and '[Helpfunctions_Performance_Measures_and_Plots.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_Performance_Measures_and_Plots.R)' into R.
2. (If the simulation in the file '[Perform_Simulation_1.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Analyses/Perform_Simulation_1.R)' has not been executed): Load the content of the file 'Simulation_1_Complete.RData' in R. Note that this file only available in the folder 'Stage_Elisca ' at CBS and not on GitHub due to storage limits.
3. Specify a working directory in line 36 to store the plots in.
4. Run the code in the file '[Plot_Simulation_1_Heatmaps.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Plots/Plot_Simulation_1_Heatmaps.R)' to create the desired plots.

# 2. Instructions to reproduce Simulation study 2 using the best approach (Chapter 5)

### Instructions to reproduce the simulation study:
1. Load the functions in the files '[Helpfunctions_General.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_General.R)', '[Helpfunctions_Simulations.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_Simulations.R)', and '[Methods_Best_Approach.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Methods_Best_Approach.R)', and '[Helpfunctions_Performance_Measures_and_Plots.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_Performance_Measures_and_Plots.R)' into R.
2. Open the file '[Perform_Simulation_2.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Analyses/Perform_Simulation_2.R)'.
3. Change the working directory in line 32.
4. Change the argument _folder = ""_ in lines 236, 240, 244 to where you would like to store the model results. Make sure to end the folder name with a "/".
5. Execute the rest of the code. Note that this takes approximately 5 days.

### Instructions to plot the results (except for the heatmaps):
1. Follow either step 1a **or** 1b:

   a. Run the code in the file '[Prepare_Simulation_Results_for_Plotting.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Plots/Prepare_Simulation_Results_for_Plotting.R)' to prepare the data for plotting. Note that this part requires that the simulation in the file '[Perform_Simulation_2.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Analyses/Perform_Simulation_2.R)' has been executed.

   b. Load the contents of the file '[Simulation_2_Reduced.RData](https://github.com/eliscamastenbroek/master_thesis/blob/main/RData/Simulation_2_Reduced.RData)' into R.
2. Run the code in the file '[Plot_Simulation_2_Best_Approach.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Plots/Plot_Simulation_2_Best_Approach.R)' to create the desired plots.

### Instructions to plot the heatmaps:
1. (If not already loaded): Load the functions in the files '[Helpfunctions_General.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_General.R)', '[Helpfunctions_Simulations.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_Simulations.R)', and '[Helpfunctions_Performance_Measures_and_Plots.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_Performance_Measures_and_Plots.R)' into R.
2. (If the simulation in the file '[Perform_Simulation_2.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Analyses/Perform_Simulation_2.R)' has not been executed): Load the content of the file 'Simulation2_Complete.RData' in R. Note that this file is not available on GitHub due to storage limits.
3. Specify a working directory in line 36 to store the plots in.
4. Run the code in the file '[Plot_Simulation_2_Best_Approach_Heatmaps.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Plots/Plot_Simulation_2_Best_Approach_Heatmaps.R)' to create the desired plots.

# 3. Instructions to reproduce Simulation study 2 using the less optimal approach (Section 5.1.1)

### Instructions to reproduce the simulation study:
1. Load the functions in the files '[Helpfunctions_General.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_General.R)', '[Helpfunctions_Simulations.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_Simulations.R)', and '[Methods_Less_Optimal_Approach.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Methods_Less_Optimal_Approach.R)', and '[Helpfunctions_Performance_Measures_and_Plots.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_Performance_Measures_and_Plots.R)' into R.
2. Open the file '[Perform_Simulation_2.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Analyses/Perform_Simulation_2_.R)'.
3. Change the working directory in line 32.
4. Change the argument _folder = ""_ in lines 236, 240, 244 to where you would like to store the model results. Make sure to end the folder name with a "/".
5. Execute the rest of the code. Note that this takes approximately 5 days.

### Instructions to plot the results (except for the heatmaps):
Nog afmaken.

# 4. Instructions to reproduce the analyses on real data fomr the ER and the LFS (Chapter 6)

Note that the data from the ER and the LFS is not publicly available.

### Instructions to reproduce the analyses:
1. Load the functions in the files '[Helpfunctions_General.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_General.R)', '[Helpfunctions_Real_Data.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/[Helpfunctions_Real_Data.R)', and '[Methods_Real_Data.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Methods_Real_Data.R)', and '[Helpfunctions_Performance_Measures_and_Plots.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_Performance_Measures_and_Plots.R)' into R.
2. Open the file '[Perform_Real_Data_Analyses.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Analyses/Perform_Real_Data_Analyses.R)'.
3. Change the working directory in line 19.
4. Change the argument _folder = ""_ in lines 82-103 to where you would like to store the model results. Make sure to end the folder name with a "/".
5. Execute the rest of the code. Note that this may take a few hours.

### Instructions to obtain and plot the results:
1. (If not already loaded): Load the functions in the files '[Helpfunctions_General.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_General.R)', '[Helpfunctions_Real_Data.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/[Helpfunctions_Real_Data.R)', and '[Helpfunctions_Performance_Measures_and_Plots.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Functions/Helpfunctions_Performance_Measures_and_Plots.R)' into R.
2. (If the models in the file '[Perform_Real_Data_Analyses.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Analyses/Perform_Real_Data_Analyses.R)' were not created): Load the content of the file 'Real_Data_Analyses.RData' in R. Note that this file is not publicly available.
3. Run the code in the file '[Plot_Real_Data_Analyses.R](https://github.com/eliscamastenbroek/master_thesis/blob/main/Plots/Plot_Real_Data_Analyses.R)' to obtain the results and to create the desired plots.
