## Description

This project implements the experiments section and provides the appendices of the paper "Improved Multi-Label Propagation for Small Data with Multi-Objective Optimization", K.Musayeva and M.Binois, 2023.

The repository is as follows:
    
* main.R
* src: contains the R code of the project
* data: contains the data used in the paper
* results: the results of the experiments are saved into this folder


## Running instructions

main.R : this is where the project should be run. This file contains the examples for running each method considered in the paper.

The following variables should be specified: 
* name of the method
* decision function / thresholding
* hyperparameter list
            
The evaluation measures used and the parameters of the cross-validation should be modified in the global_var.R file.

To run the harmonic function as a stacking method, you need to first run the method the predictions of which you want to 

use in the stacking (the results are saved automatically into results/name_data_*.obj file). Then the compromise solution, 

more precisely, the hyperparameter configuration corresponding to the compromise solution should be computed from this .obj file. 

This hyperparameter configuration is then used in the stacking. 

The reported results are obtained by applying the ks_solution to the results/name_data_*.obj file.

### Typos in the paper

* Section 3: $D^{'}_{jj}$ should be
```math
D^{'}_{jj}=\sum_{i=1}^{n} W_{ij}
```

* Section 4: $A_{ii}$ should be
```math
A_{ii}=\sum_{j=1}^C \tilde{\mathbf{Y}}_{ij}
```
and $B_{jj}$ should be
```math
B_{jj}=\sum_{i=1}^n \tilde{\mathbf{Y}}_{ij}
```
