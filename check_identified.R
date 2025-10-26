# Goal is to build code to check what is identified given a set of restrictions
# For now, assume that the instruments and covariates are discrete. 

# We will assume that the researcher hands off a list of response types
# Since instruments are discrete, a response type will look like a vector in R^{N_Z}
# where N_Z is the number of instruments 

# We will assume that the parameter we are interested in is of the form E[g(Y,X)l(T*,X)]
# The function l(T*,X) will be handed to us in a {N_T* x N_X} matrix
# where N_T* is the number of response types and N_X is the number of covariate values

#########################################################################
# General code to check if l(T*,X) is identified
#########################################################################

isIdentified <- function(ellMat, responseTypes) {
  
  # locate the number of instruments
  NZ <- dim(responseTypes)[1]
  
  # locate the number of response types
  NTypes <- dim(responseTypes)[2]
  
  # locate the number of covariates
  NX <- dim(ellMat)[2]
  
  # collect the set of treatments
  responseMat <- matrix(responseTypes, nrow = NZ)
  treatments <- c()
  for (i in 1:NZ) {
    for (j in 1:NTypes) {
      t = responseMat[i,j]
      if (!t %in% treatments) {
        treatments <- c(treatments, t)
      }
    }
  }
  
  # Locate the number of treatments
  NT <- length(treatments)
  
  # Build the total design matrix
  designMat <- cbind()
  for (t in 0:(NT - 1)) {
    Bt <- 1*(responseMat == t)
    designMat <- cbind(designMat, t(Bt))
  }
  
  # For each row of ellMat, check if it is in the range of linear operator Upsilon
  identified = 1
  nu = cbind()
  for (xval in 1:NX) {
    ellVec <- ellMat[, xval]
    temp <- lm(ellVec ~ designMat + 0)
    temp.sum <- summary(temp)
    if (temp.sum$r.squared < 0.999) {
      identified = 0
      # print("Functional is not Identified")
    }
    else {
      beta <- temp$coefficients
      beta[is.na(beta)] <- 0
      nu = cbind(nu, beta)
    }
  }
  # if (identified == 1) {
  #   print("Functional is Identified")
  # }
  list("isIdentified" = identified, nu = nu)
}

#########################################################################
# Building on this code, check which type probabilities are identified
#########################################################################

typesIdentified <- function(responseTypes) {
  
  # Generate all the 'ell' vectors
  NTypes <- dim(responseTypes)[2]
  idMat <- diag(nrow = NTypes)
  
  # Check if each response type is identified
  idVec <- c()
  for (type in 1:NTypes) {
    ellMat <- matrix(idMat[,type])
    idCheck <- isIdentified(ellMat, responseTypes)
    if (idCheck$isIdentified == 1) {
      print(paste("Pr(T* = t", type, "*) is identified", sep = ""))
    }
    idVec <- c(idVec, idCheck$isIdentified)
  }
  idVec
}

#########################################################################
# Check which outcome functionals are identified
#########################################################################

outcomesIdentified <- function(responseTypes) {
  
  # locate the number of instruments
  NZ <- dim(responseTypes)[1]
  
  # locate the number of response types
  NTypes <- dim(responseTypes)[2]
  idMat <- diag(nrow = NTypes)
  
  # locate the number of covariates
  NX <- dim(ellMat)[2]
  
  # collect the set of treatments
  responseMat <- matrix(responseTypes, nrow = NZ)
  treatments <- c()
  for (i in 1:NZ) {
    for (j in 1:NTypes) {
      t = responseMat[i,j]
      if (!t %in% treatments) {
        treatments <- c(treatments, t)
      }
    }
  }
  
  # Locate the number of treatments
  NT <- length(treatments)
  
  
  # loop through the treatments
  for (t in treatments) {
    Bt <- t(1*(responseMat == t))
  
    # Loop through the response types
    for (type in 1:NTypes) {
      ellMat <- matrix(idMat[, type])
      temp <- lm(ellMat ~ Bt + 0)
      temp.sum <- summary(temp)
      if (temp.sum$r.squared > 0.999) {
        print(paste(
          "E[g(Y",t,") | T* = t",type,"*] is identified",
          sep = ""
        ))
      }
    }
  }
}

#########################################################################
# Estimate parameters of interest given data
#########################################################################
library(glmnet)
library(pracma)

estimateTypeProb <- function(formatted_data, ellMat, responseTypes) {
  
  # Check if the parameter is identified
  test <- isIdentified(ellMat, responseTypes)
  if (test$isIdentified == 0) {
    print("Parameter is not Identified")
    break
  }
  
  # Collect needed attributed from data
  Y <- formatted_data$Y
  D <- formatted_data$T
  Z <- formatted_data$Z
  X <- formatted_data$X
  n <- dim(X)[1]
  if (dim(Z)[2] != n) {
    Z <- t(Z)
  }
 
  # locate the number of instruments
  NZ <- dim(responseTypes)[1]
  
  # locate the number of response types
  NTypes <- dim(responseTypes)[2]
  idMat <- diag(nrow = NTypes)
  
  # locate the number of covariates
  NX <- dim(ellMat)[2]
  
  # collect the set of treatments
  responseMat <- matrix(responseTypes, nrow = NZ)
  treatments <- c()
  for (i in 1:NZ) {
    for (j in 1:NTypes) {
      t = responseMat[i,j]
      if (!t %in% treatments) {
        treatments <- c(treatments, t)
      }
    }
  }
  
  # Locate the number of treatments
  NT <- length(treatments)
  
  # Collect nu from test 
  nuVec <- test$nu
  
  # Go from the nu vector to a function
  nuFunction <- function(t, z, xIndex = 1) {
    
    # Collect nu from test again
    nuVec2 <- test$nu[,xIndex]
    nuMat <- t(matrix(nuVec2, nrow = NZ))
    
    nu = 0
    for (j in 0:(NT-1)) {
      for (k in 0:(NZ-1)) {
        nu = nu + nuMat[j + 1,k + 1]*(t == j)*(z == k)
      }
    }
   
    nu 
  }
  
  # Setup for Estimation Procedure
  basis <- createBasis(Z, X)
  muBasis <- createMuBasis(Z, X)
  nu <- 0 
  for (j in 1:NX) {
    nu <- nu + nuFunction(D, Z, j)
  }
  
  # Estimate Nuisance Models
  cm.model <- cv.glmnet(basis, nu, intercept = T)
  cm.fitted <- predict(cm.model, basis)
  cm.mufitted <- predict(cm.model, muBasis)
  
  rr.model <- estimateRR(formatted_data)
  rr.fitted <- predict(rr.model, basis)
  
  # Construct orthogonal score
  psi <- rr.fitted * (t(nu) - cm.fitted) + cm.mufitted
  psi <- NZ*psi
  
  # Return psi, mean, and standard deviation
  list("scores" = psi, "point.estimate" = mean(psi), "se" = sd(psi)/sqrt(n), "nu" = unlist(as.list(nu)))
  
}

estimateAllTypes <- function(formatted_data, responseTypes) {
  
  # Generate all the 'ell' vectors
  NTypes <- dim(responseTypes)[2]
  idMat <- diag(nrow = NTypes)
  
  # Set up vectors to store identified parameters, point estimates, and scores
  labels.identified <- c()
  point.estimates <- c()
  standard.errors <- c()
  scores <- cbind()
  
  # Check if each response type is identified and estimate if so
  for (type in 1:NTypes) {
    ellMat <- matrix(idMat[,type])
    idCheck <- isIdentified(ellMat, responseTypes)
    if (idCheck$isIdentified == 1) {
      
      # Conduct estimation
      estimation <- estimateTypeProb(formatted_data, ellMat, responseTypes)
      label = paste("Pr(T* = t", type, "*)", sep = "") 
      print(paste(label, "is identified."))
      print(paste("Point Estimate: ", round(estimation$point.estimate, 5), 
                  ", Standard Error: ", round(estimation$se, 5), 
                  sep = ""))
      
      # Store relevant quantities
      labels.identified <- c(labels.identified, label)
      point.estimates <- c(point.estimates, estimation$point.estimate)
      standard.errors <- c(standard.errors, estimation$se)
      scores <- cbind(scores, estimation$scores)
    }
  }
  
  names(point.estimates) = labels.identified
  names(standard.errors) = labels.identified 
  vcov = genVCOV.types(scores, labels.identified)
  
  list("identified.params" = labels.identified, 
       "point.estimates" = point.estimates, 
       "standard.errors" = standard.errors,
       "vcov" = vcov,
       "scores" = scores)
}

estimateOutcome <- function(formatted_data, ellMat, responseTypes, tVal = 0) {
  
  # Collect attributes from the data
  Y <- formatted_data$Y
  D <- formatted_data$T
  Z <- formatted_data$Z
  X <- formatted_data$X
  n <- dim(X)[1]
  if (dim(Z)[2] != n) {
    Z <- t(Z)
  }
 
  # locate the number of instruments
  NZ <- dim(responseTypes)[1]
  
  # locate the number of response types
  NTypes <- dim(responseTypes)[2]
  idMat <- diag(nrow = NTypes)
  
  # locate the number of covariates
  NX <- dim(ellMat)[2]
  
  # Check if parameter is identified
  responseMat <- matrix(responseTypes, nrow = NZ)
  Bt <- t(1*(responseMat == tVal))
  temp <- lm(ellMat ~ Bt + 0)
  # temp.sum <- summary(temp) 
  if (var(temp$residuals) > 0.001) {
    print("Parameter is not identified")
    break
  }
  
  temp$coefficients[is.na(temp$coefficients)] = 0
  nuVec <- as.numeric(temp$coefficients)
  # Turn the nu vector into a function
  nuFunction <- function(t, z) {
    nu <- 0
    for (j in 0:(NZ - 1)) {
      nu <- nu + nuVec[j + 1] * ( t == tVal) * (z == j)
    }
    
    nu
  }
  
  # Setup for estimation procedure
  basis <- createBasis(Z, X)
  muBasis <- createMuBasis(Z, X)
  nu <- Y * nuFunction(D, Z)
  
  # Estimate Nuisance Models
  cm.model <- cv.glmnet(basis, nu, intercept = T)
  cm.fitted <- predict(cm.model, basis)
  cm.mufitted <- predict(cm.model, muBasis)
  
  rr.model <- estimateRR(formatted_data)
  rr.fitted <- predict(rr.model, basis)
  
  # Construct orthogonal score
  psi <- rr.fitted * (t(nu) - cm.fitted) + cm.mufitted
  psi <- NZ*psi

  # Report means and standard deviations
  list("scores" = psi, 
       "point.estimate" = mean(psi), 
       "se" = sd(psi)/sqrt(n))
  
}

estimateAllOutcomes <- function(formatted_data, responseTypes) {
 
  # locate the number of instruments
  NZ <- dim(responseTypes)[1]
  n <- dim(formatted_data$X)[1]
  
  # locate the number of response types
  NTypes <- dim(responseTypes)[2]
  idMat <- diag(nrow = NTypes)
  
  # locate the number of covariates
  NX <- dim(ellMat)[2]
  
  # collect the set of treatments
  responseMat <- matrix(responseTypes, nrow = NZ)
  treatments <- c()
  for (i in 1:NZ) {
    for (j in 1:NTypes) {
      t = responseMat[i,j]
      if (!t %in% treatments) {
        treatments <- c(treatments, t)
      }
    }
  }
  
  # Locate the number of treatments
  NT <- length(treatments)
  
  # Set up vectors to store identified parameters, point estimates, and scores
  labels.identified <- c()
  point.estimates <- c()
  standard.errors <- c()
  outcome.scores <- cbind()
  type.scores <- cbind()
  
  # loop through the treatments
  for (type in 1:NTypes) {
    ellMat <- matrix(idMat[, type])
    
    # Loop through the response types
    for (t in treatments) {
      
      # Check if parameter is identified
      # First, setup to only use observations for which T = t
      Bt <- t(1*(responseMat == t))
      
      # Check if there is a kappa that maps to \ell under Upsilon
      temp <- lm(ellMat ~ Bt + 0)
      
      # If parameter is identified, residuals should be zero
      if (var(temp$residuals) < 0.001) {
        
        # Conduct Estimation
        # We are interested in conditional means, so need E[Y(t)*1{T* = t*}]
        # as well as Pr(T* = t*). Point estimate is a ratio of these two
        typeEstimation <- estimateTypeProb(formatted_data, ellMat, responseTypes) 
        outcomeEstimation <- estimateOutcome(formatted_data, ellMat, responseTypes, t)
        point.estimate <- outcomeEstimation$point.estimate / typeEstimation$point.estimate
        
        # Calculate the asymptotic variance
        var.outcome <- (outcomeEstimation$se) ** 2
        var.type <- (typeEstimation$se) ** 2 
        pe.outcome <- outcomeEstimation$point.estimate 
        pe.type <- typeEstimation$point.estimate
        
        # Use Delta Method to get variance for the conditional mean
        cov <- cov(outcomeEstimation$scores, typeEstimation$scores) / n
        vcov <- matrix(c(var.outcome, cov, cov, var.type), nrow = 2)
        g <- c(1/pe.type, -1*pe.outcome/ ( pe.type **2)) 
        avar <- t(g) %*% vcov %*% g
        
        # Print relevant quantities and store the name of identified parameter
        label = paste("E[Y(",t,") | T* = t", type,"*]", sep = "")
        print(paste(label, "is identified."))
        print(paste("Point Estimate: ", round(point.estimate, 5),
                    ", Standard Error: ", round(sqrt(avar), 5),
                    sep = ""))
        
        # Store all needed quantities
        labels.identified <- c(labels.identified, label)
        point.estimates <- c(point.estimates, point.estimate) 
        standard.errors <- c(standard.errors, sqrt(avar)) 
        outcome.scores <- cbind(outcome.scores, outcomeEstimation$scores) 
        type.scores <- cbind(type.scores, typeEstimation$scores)
      }
    }
  }
  
  # Generate variance covariance matrix of estimated quantities
  vcov = genVCOV.outcomes(outcome.scores, type.scores, labels.identified)
  
  # Apply labels
  names(point.estimates) = labels.identified
  names(standard.errors) = labels.identified 
  colnames(type.scores) = labels.identified
  colnames(outcome.scores) = labels.identified
  
  # Return all relevant parameters
  list("identified.params" = labels.identified, 
       "point.estimates" = point.estimates,
       "standard.errors" = standard.errors, 
       "outcome.scores" = outcome.scores, 
       "vcov" = vcov,
       "type.scores" = type.scores)
}

#########################################################################
# Helper Functions for Estimation
#########################################################################
library(matrixcalc)

createBasis <- function(dataZ, dataX) {
  
  # Collect information from the data
  NX <- dim(dataX)[2]
  instrumentValues <- as.numeric(names(table(dataZ)))
  
  # Create basis by interacting each Xj, Xj^2 with instrument indicators
  basis <- cbind()
  for (j in 1:NX) {
    for (inst in instrumentValues) {
      Z.Xj <- dataX[,j] * (dataZ == inst)
      Z.Xj2 <- (dataX[,j]**2) * (dataZ == inst)
      basis <- cbind(basis, t(Z.Xj), t(Z.Xj2))
    }
  }
  
  # Need the basis to be non-singular for the solver
  if (is.singular.matrix(t(basis) %*% basis)) {
    basis <- cbind()
    for (j in 1:NX) {
      for (inst in instrumentValues) {
        Z.Xj <- dataX[,j] * (dataZ == inst) 
        basis <- cbind(basis, t(Z.Xj))
      }
    }
  }
  
  basis
}

createMuBasis <- function(dataZ, dataX) {
    
  # Collect information from the data
  NX <- dim(dataX)[2]
  instrumentValues <- as.numeric(names(table(dataZ)))
  
  # Create basis by interacting each Xj, Xj^2 with instrument indicators
  basis <- cbind()
  for (j in 1:NX){
    for (inst in instrumentValues) {
      Z.Xj <- dataX[,j] / length(instrumentValues)
      Z.Xj2 <- (dataX[,j]**2) / length(instrumentValues)
      basis <- cbind(basis, Z.Xj, Z.Xj2)
    }
  }
  
  # Check if initial basis would be singular
  basis.init <- cbind()
  for (j in 1:NX) {
    for (inst in instrumentValues) {
      Z.Xj <- dataX[,j] * (dataZ == inst)
      Z.Xj2 <- (dataX[,j]**2) * (dataZ == inst)
      basis.init <- cbind(basis.init, t(Z.Xj), t(Z.Xj2))
    }
  }
  
  # If initial basis would be singular, do not include interaction terms
  if (is.singular.matrix(t(basis.init) %*% basis.init)) {
    basis <- cbind()
    for (j in 1:NX) {
      for (inst in instrumentValues) {
        Z.Xj <- dataX[,j] / length(instrumentValues)
        basis <- cbind(basis, Z.Xj)
      }
    }
  }
  
  basis
}

estimateRR <- function(formatted_data) {
  
  # Get information from the data
  Z <- formatted_data$Z 
  X <- formatted_data$X
  n <- dim(X)[1]
  if (dim(Z)[2] != n) {
    Z <- t(Z)
  }
  instrumentValues <- as.numeric(names(table(Z)))
  NZ <- length(instrumentValues)
  
  # Create mu 
  muBasis <- createMuBasis(Z, X)
  
  # Create the regular basis
  basis <- createBasis(Z, X)
  
  # Complete the square by solving for the auxiliary Y
  auxY <- basis %*% solve( t(basis) %*% basis) %*% colSums(muBasis)
  
  # Estimate the reisz representer
  cv.init  <- cv.glmnet(basis, auxY, intercept = T)
  glmnet(basis, auxY, lambda = cv.init$lambda.1se, intercept = T)
}

genVCOV.types <- function(scores, labels) {
  num.types <- dim(scores)[2]
  n <- dim(scores)[1]
  vcov <- matrix(0, nrow = num.types, ncol = num.types)
  for (j in 1:num.types) {
    for (k in 1:num.types) {
      vcov[j,k] = cov(scores[,j], scores[,k])/n
    }
  }
  
  rownames(vcov) = labels 
  colnames(vcov) = labels
  vcov
}

genVCOV.outcomes <- function(outcome.scores, type.scores, 
                             labels) {
  
  # Get some information from the data
  num.outcomes <- dim(outcome.scores)[2]
  n <- dim(outcome.scores)[1]
  outcome.estimates = colMeans(outcome.scores)
  type.estimates = colMeans(type.scores)
  
  # Set up the original variance covariance matrix of outcomes and type probabilities
  vcov.init <- matrix(0, nrow = 2*num.outcomes, ncol = 2*num.outcomes)
  for (j in 1:num.outcomes) {
    for (k in 1:num.outcomes) {
      vcov.init[(2*j - 1), (2*k - 1)] = cov(outcome.scores[,j], outcome.scores[,k])/n 
    }
  }
  for (j in 1:num.outcomes) {
    for (k in 1:num.outcomes) {
      vcov.init[(2*j), (2*k)] = cov(type.scores[,j], type.scores[,k])/n
    }
  }
  for (j in 1:num.outcomes) {
    for (k in 1:num.outcomes) {
      vcov.init[2*j - 1, 2*k] = cov(outcome.scores[,j], type.scores[,k])/n
      vcov.init[2*k, 2*j - 1] = cov(outcome.scores[,j], type.scores[,k])/n
    }
  }
  
  # Set up delta method derivative
  g.init <- matrix(0, nrow = num.outcomes, ncol = 2*num.outcomes)
  for (j in 1:num.outcomes) {
    g.init[j, (2*j - 1)] = 1/type.estimates[j]
    g.init[j, (2*j)] = -1*outcome.estimates[j]/ (type.estimates[j])**2
  }
  
  # Get final vcov matrix from delta method
  vcov.final = g.init %*% vcov.init %*% t(g.init)
  
  # Apply labels 
  rownames(vcov.final) = labels 
  colnames(vcov.final) = labels
  
  vcov.final
}


