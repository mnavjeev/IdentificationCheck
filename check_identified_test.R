#########################################################################
# Test Using LATE Model
#########################################################################
source("check_identified.R")
library(mvtnorm)

# Preliminaries
lateResponse <- cbind(c(1,1), c(0,0), c(0,1))
ellMat <- cbind(c(0,0,1))

# Some global hyper parameters
always_taker_share = 0.2
complier_share = 0.5

# Complier means
y1_complier_mean = 2
y0_complier_mean = 0

# Always taker mean
y1_always_taker_mean = 1

# Never taker mean
y0_never_taker_mean = -1

# Function to Create LATE Data
lateDGP <- function(n) {
  
  # Assign instrument conditional on covariates
  X <- rmvnorm(n, mean = c(0,0,0,0,0))
  X.init <- 2*(rbinom(n*5, 1, 0.5) - 0.5)
  X <- matrix(X.init, nrow = n, ncol = 5)
  Z <- 1*( ((c(0.25,0.25,0.25,0.25,0.25) %*% t(X)) + rnorm(n)) >= 0)
  # Z <- t(Z)
  
  # Divide into compliers, always takers, and never takers
  AT_upper <- always_taker_share
  C_upper <- always_taker_share + complier_share
 
  # Assign type conditional on the covariates 
  Type.init <- (c(0.25,0.25,0.25,0.25,0.25) %*% t(X)) + rnorm(n)
  F_hat <- ecdf(Type.init)
  Type <- F_hat(Type.init)
  
  D <- integer(n)
  D[Type <= AT_upper] <- 1 # Always Takers
  D[(Type > AT_upper & Type <= C_upper)] <- Z[(Type > AT_upper & Type <= C_upper)] # Compliers
  D[Type > C_upper] <- 0 # Never Takers
  
  # Means of Y
  Y <- integer(n)
  numAT <- length(D[Type <= AT_upper])
  numCompliers <- length(D[(Type > AT_upper & Type <= C_upper)])
  numNT <- length(D[Type > C_upper])
  Y[Type <= AT_upper] <- rnorm(numAT, y1_always_taker_mean) 
  Y[(Type > AT_upper & Type <= C_upper)] <- rnorm(numCompliers, y0_complier_mean) * 
    (1 - D[(Type > AT_upper & Type <= C_upper)]) + 
    rnorm(numCompliers, y1_complier_mean)*D[(Type > AT_upper & Type <= C_upper)] 
  Y[(Type > C_upper)] <- rnorm(numNT, y0_never_taker_mean) 
  
  # Return all in dataframe
  list("Y" = Y, "Z" = Z, "T" = D, "X" = X)
  
}

########################################
# Set Global Hyperparameters
########################################

# Type probabilities
always_taker_share = 0.2
complier_share = 0.5

# Complier means
y1_complier_mean = 2
y0_complier_mean = 0

# Always taker mean
y1_always_taker_mean = 1

# Never taker mean
y0_never_taker_mean = -1

# Sample Size
n = 1000

#######################################
# Generate Data and Perform Estimation
#######################################

lateData <- lateDGP(n)
typeEstimation <- estimateAllTypes(lateData, lateResponse)
outcomeEstimation <- estimateAllOutcomes(lateData, lateResponse)

#######################################
# Test using only Always Takers
#######################################
atResponse <- cbind(
  c(0,0,0),
  c(1,1,1),
  c(2,2,2)
)

atDGP <- function(n) {
  
  # Assign instrument conditional on covariates
  X <- rmvnorm(n, mean = c(0,0,0,0,0))
  Z.init <- ((c(0.25,0.25,0.25,0.25,0.25) %*% t(X)) + rnorm(n))
  F_hat <- ecdf(Z.init)
  Z.type <- F_hat(Z.init)
  Z <- integer(n)
  Z[Z.type <= 0.3] = 0
  Z[(Z.type > 0.3 & Z.type <= 0.65)] = 1
  Z[(Z.type > 0.65)] = 2
  Z <- matrix(Z, nrow = 1, ncol = n)
  Z <- t(Z)
  
  # Assign type conditional on the covariates 
  Type.init <- (c(0.25,0.25,0.25,0.25,0.25) %*% t(X)) + rnorm(n)
  F_hat <- ecdf(Type.init)
  Type <- F_hat(Type.init)
  
  D <- integer(n)
  D[Type <= 0.3] <- 0 # Always Takers
  D[(Type > 0.3 & Type <= 0.65)] <- 1 # Compliers
  D[Type > 0.65] <- 2 # Never Takers
  
  # Means of Y
  Y <- integer(n)
  numAT.1 <- length(D[Type <= 0.3])
  numAT.2 <- length(D[(Type > 0.3 & Type <= 0.65)])
  numAT.3 <- length(D[Type > 0.65])
  Y[Type <= 0.3] <- rnorm(numAT.1, 1) 
  Y[(Type > 0.3 & Type <= 0.65)] <- rnorm(numAT.2, 2)
  Y[(Type > 0.65)] <- rnorm(numAT.3, 3) 
  
  # Return all in dataframe
  list("Y" = Y, "Z" = Z, "T" = D, "X" = X)
}

atData <- atDGP(n)

typeEstimation.AT <- estimateAllTypes(atData, atResponse)
outcomeEstimation.AT <- estimateAllOutcomes(atData, atResponse)



