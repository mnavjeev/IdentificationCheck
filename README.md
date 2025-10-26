# Documentation

This repository contains code to check whether a certain class of parameters considered in [Navjeevan, Pinto, and Santos (2025)](https://arxiv.org/abs/2310.05311) -- roughly, those considered in Corollary 4.4 of the paper. Conditional on being identified, the code also implements the estimator proposed in Section 5 of the paper, reporting both point estimates and standard errors. The two code files in the repository are:

1. `check_identified.R`: This file contains code to check whether the parameters of interest are identified based on the conditions outlined in Corollary 4.4 of the paper and estimate them if they are identified.

2. `check_identified_test.R`: This file contains some example code to illustrate how to use the functions in `check_identified.R`. 

## How to use the code 

The code is written in R. Each function takes some combination of the following inputs: 

1. "responseTypes": This is a matrix of length $N_Z \times N_S$ whose entries are elements of $`\{0,1,\dots,N_T - 1\}`$ where $N_Z$ is the number of instruments, $N_S$ is the number of response types 
