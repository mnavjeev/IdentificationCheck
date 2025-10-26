# Documentation

This repository contains code to check whether a certain class of parameters considered in [Navjeevan, Pinto, and Santos (2025)](https://arxiv.org/abs/2310.05311) -- roughly, those considered in Corollary 4.4 of the paper. Conditional on being identified, the code also implements the estimator proposed in Section 5 of the paper, reporting both point estimates and standard errors. The two code files in the repository are:

1. `check_identified.R`: This file contains code to check whether the parameters of interest are identified based on the conditions outlined in Corollary 4.4 of the paper and estimate them if they are identified.

2. `check_identified_test.R`: This file contains some example code to illustrate how to use the functions in `check_identified.R`. 

## How to use the code 

The code is written in R. Each function takes some combination of the following inputs: 

1. "responseTypes": This is a matrix of shape $N_Z \times N_S$ whose entries are elements of $`\{0,1,\dots,N_T - 1\}`$ where $N_Z$ is the number of instruments, $N_S$ is the number of response types. The assignment of treatment values to numbers need not represent any particular ordering over the response types but is rather a convention adopted so that the code can easily check which parameters are identified.

2. "ellMat": This is a matrix of shape $N_S \times N_X$ whose entries are again elements of $`\{0,1,\dots,N_T - 1\}`$ where $N_X$ is the number of covariate values if there is a discrete covariate. This vector represents values of the function $\ell(T^\star, X)$ evaluated at each response type and covariate value. The ordering of the response types corresponds to the ordering in the "responseTypes" matrix, i.e the first column of responseTypes corresponds to the first row of ellMat, and so on. If the covariate $X$ does not enter the function $\ell(T^\star, X)$, this can be a vector of length $N_S$.

3. "formatted_data": This is a list containing the data in a specific format. The list should contain the following elements:
   - "Y": A vector of length $N$ containing the outcome variable with values in $\mathbb{R}$.
   - "T": A vector of length $N$ containing the observed treatment variable with values in $`\{0,1,\dots,N_T - 1\}`$.
   - "Z": A vector of length $N$ containing the instrument variable with values in $`\{0,1,\dots,N_Z - 1\}`$. As with the treatment variable, the numbering of the instrument values does not represent any particular ordering but is rather a convention adopted for ease of implementation.
   - "X":  A matrix of size $N \times N_X$ containing the covariate variable with values in $\mathbb{R}$. If there are no covariates set this to a vector of ones of length $N$.
