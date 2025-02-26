#' Source functions for the paper "Improved order selection method for hidden Markov models: a case study with movement data."

#' Computing the logarithm of a sum given the logarithms of the summands
#' 
#' @param logx Logarithm of x
#' @param logy Logarithm of y
#' @return The extended logarithm of the sum of x and y given as inputs the extended logarithm of x and y based on http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf.

logspace_add <- function(logx,logy){
  if(logx==-Inf || logy==-Inf){
    if(logx==-Inf){
      return(logy)
    }else{
      return(logx)
    }
  }else{
    
    if(logx>logy){
      return(logx+log(1+exp(logy-logx)))
    }else{
      return(logy+log(1+exp(logx-logy)))
    }
  }
}

#' Ordering the multidimensional parameter
#' 
#' @param N Upper bound (for the number of states) for double penalized method. 
#' @param par Vector of parameters to sort with the group-fuse method from Manole, T. and Khalili, A. (2021).
#' @return Sorted parameters.

ordering <- function(N, par) {
  # Initialize the sequence of indices
  alpha <- numeric(N) 
  # Choose the index with the minimum value as the first index
  alpha[1] <- which.min(apply(par,1,FUN=function(x)sqrt(sum(x^2))))  # done with max in  Manole, T. and Khalili, A. (2021).
  for (k in 2:N) {
    # Exclude indices already selected
    available_indices <- setdiff(1:N, alpha[1:(k-1)]) 
    # Calculate distances
    distances <- sapply(available_indices, function(i)sqrt(sum((par[i]-par[alpha[k-1]])^2)))
    # Choose index with minimum distance
    alpha[k] <- available_indices[which.min(distances)] 
  }
  return(alpha)
}

#' Numerically stable forward algorithm for non-stationary HMM based on the code provided in McClintock, B. T. (2021).
#' 
#' @param pi Vector of initial distribution. It should be a vector of size nbStates.
#' @param trMat Transition Matrix. dim(trMat) =  c(nbSteps,nbStates, nbStates).
#' @param lnProbs Diagonal matrix of log emission distributions for all observations (for all individuals).
#' @param nbSteps Length of the time-series.
#' @param animal Number of the individual (from 1 to nbAnimals)
#' @param aInd First observation for each individual
#' @return Joint negative log-likelihood of the model.

forward_algNH <- function(pi,trMat,lnProbs,nbSteps,aInd)
{
  # Define the number of states
  nbStates<- ncol(trMat)
  # Initialize matrices, arrays and vectors
  ltrMat <- array(NA,dim=c(nbSteps,nbStates,nbStates))
  lpiG <- lalpha <- lnewalpha <- matrix(NA,nrow=nbStates,ncol=1)
  jnll <- 0
  nbAnimals <-  length(aInd)
  animal <- 1
  for(t in (1:nbSteps)){
    if(animal<=nbAnimals && t==aInd[animal]){
      # Start likelihood at 0, and log(0) = -Inf
      sumalpha <-  -Inf
      # Formula for joint log-likelihood at t=1 of individual "animal" 
      for(j in (1:nbStates)){
        lpiG[j] <-  -Inf
        for(i in (1:nbStates)){
          # Take the log of the matrix
          ltrMat[t,i,j] <-  log(trMat[t,i,j])
          # Likelihood of S_1 = i and S_2 = j.
          lpiG[j] <-  logspace_add(lpiG[j],log(pi[i])+ltrMat[t,i,j])
        }
        # Likelihood of S_1 = i and S_2 = j and emission probabilities at time 2.
        lalpha[j] <-  lpiG[j]+lnProbs[t,j]
        sumalpha  <-  logspace_add(sumalpha,lalpha[j])
      }
      # Negative log-likelihood.
      jnll <- jnll- sumalpha
      lalpha <- lalpha -sumalpha
      animal <-  animal + 1
    } else {
      sumalpha <-  -Inf
      for(j in (1:nbStates)){
        logalpha <-  -Inf
        for(i in (1:nbStates)){
          ltrMat[t,i,j] <-  log(trMat[t,i,j])
          # Likelihood of S_{t-1} = i and S_t = j + likelihood of everything before.
          logalpha <-  logspace_add(logalpha,lalpha[i]+ltrMat[t,i,j])
        }
        # Likelihood up to time t + emission probabilities at time t.
        lnewalpha[j] <-  logalpha + lnProbs[t,j]
        sumalpha <-  logspace_add(sumalpha,lnewalpha[j])
      }
      # Negative log-likelihood.
      jnll <-  jnll - sumalpha
      lalpha = lnewalpha - sumalpha
    }
  }
  return(jnll)
}


#' Numerically stable forward algorithm for stationary HMM based on the code provided in McClintock, B. T. (2021).
#' 
#' @param pi Stationary distribution. It should be a vector of size nbStates.
#' @param trMat Transition Matrix. dim(trMat) =  c(nbStates, nbStates) for stationary.
#' @param lnProbs Diagonal matrix of log emission distributions for all observations (for all individuals).
#' @param nbSteps Length of the time-series.
#' @param animal Number of the individual (from 1 to nbAnimals)
#' @param aInd First observation for each individual
#' @return Joint negative log-likelihood of the model.

forward_alg <- function(pi,trMat,lnProbs,nbSteps,aInd){
  #Same comments as for forward_algNH. Dimensions for ltrMat are different.
  nbStates<- ncol(trMat)
  ltrMat <- array(NA,dim=c(nbStates,nbStates))
  lpiG <- lalpha <- lnewalpha <- matrix(NA,nrow=nbStates,ncol=1)
  jnll <- 0
  nbAnimals <-  length(aInd)
  animal <- 1
  for(t in (1:nbSteps)){
    if(animal<= nbAnimals && t==aInd[animal]){
      sumalpha <-  -Inf
      for(j in (1:nbStates)){
        lpiG[j] <-  -Inf
        for(i in (1:nbStates)){
          ltrMat[i,j] <-  log(trMat[i,j])
          lpiG[j] <-  logspace_add(lpiG[j],log(pi[i])+ltrMat[i,j])
        }
        lalpha[j] <-  lpiG[j]+lnProbs[t,j]
        sumalpha  <-  logspace_add(sumalpha,lalpha[j])
      }
      jnll <- jnll- sumalpha
      lalpha <- lalpha -sumalpha
      animal <-  animal + 1
    } else {
      sumalpha <-  -Inf
      for(j in (1:nbStates)){
        logalpha <-  -Inf
        for(i in (1:nbStates)){
          logalpha <-  logspace_add(logalpha,lalpha[i]+ltrMat[i,j])
        }
        lnewalpha[j] <-  logalpha + lnProbs[t,j]
        sumalpha <-  logspace_add(sumalpha,lnewalpha[j])
      }
      jnll <-  jnll - sumalpha
      lalpha = lnewalpha - sumalpha
    }
  }
  return(jnll)
}

#' Numerically stable forward algorithm based on the code provided in McClintoCk, B. T. (2021). Only used for non-stationary HMM.
#' 
#' @param pi Stationary distribution (or its estimate if non-stationary HMM). It should be a vector of size nbStates.
#' @param trMat Transition Matrix.
#' @param lnProbs Diagonal matrix of log emission distributions. Different from the one provided for forward_alg. lnProbs only corresponds to the observation from individual "animal".
#' @param nbSteps Length of the time-series.
#' @param animal Number of the animal (from 1 to nbAnimals)
#' @param nbAnimals Total number of animals
#' @return Negative log-likelihood for individual "animal".

forward_alg_ind <- function(pi,trMat,lnProbs,nbSteps,animal,nbAnimals)
{
  nbStates<- ncol(trMat)
  ltrMat <- array(dim=dim(trMat))
  lpiG <- lalpha <- lnewalpha <- matrix(NA,nrow=nbStates,ncol=1)
  jnll <- 0
  
  if(animal!=nbAnimals){
    TT = nrow(lnProbs)+aInd[animal]-1
  }else{
    TT = nrow(lnProbs)+aInd[animal]-1
  }
  
  LL = TT-aInd[animal]+1
  sumalpha <-  -Inf
  for(j in (1:nbStates)){
    lpiG[j] <-  -Inf
    for(i in (1:nbStates)){
      ltrMat[1,i,j] <-  log(trMat[1,i,j])
      lpiG[j] <-  logspace_add(lpiG[j],log(pi[i])+ltrMat[1,i,j])
    }
    lalpha[j] <-  lpiG[j]+lnProbs[1,j]
    sumalpha  <-  logspace_add(sumalpha,lalpha[j])
  }
  jnll <- jnll- sumalpha
  lalpha <- lalpha -sumalpha
  
  for(t in (2:LL)){
    sumalpha <-  -Inf
    for(j in (1:nbStates)){
      logalpha <-  -Inf
      for(i in (1:nbStates)){
        ltrMat[t,i,j] <-  log(trMat[t,i,j])
        logalpha <-  logspace_add(logalpha,lalpha[i]+ltrMat[t,i,j])
      }
      lnewalpha[j] <-  logalpha + lnProbs[t,j]
      sumalpha <-  logspace_add(sumalpha,lnewalpha[j])
    }
    jnll <-  jnll - sumalpha
    lalpha = lnewalpha - sumalpha
  }
  return(jnll)
}

#' SCAD penalty function
#' 
#' @param eta Vector of sorted parameters 
#' @param lambda Hyperparameter of penalty function
#' @param a Second hyperparameter of penalty function, usually set to a = 3.7 as recommended by Fan, J. and Li, R. (2001). 
#' @return Value of the SCAD penalty function

pen <- function(eta,lambda,a){
  if(abs(eta)<=lambda){
    return(lambda*abs(eta))
  }else{
    if(lambda < abs(eta) && abs(eta) <=  a*lambda){
      return((2*a*lambda*abs(eta)-eta^2-lambda^2)/(2*(a-1)))
    }else{
      return(lambda^2*(a+1)/2)
    }
  }
}

#' local approximation of the SCAD penalty
#' @param eta Vector of sorted parameters from previous iteration
#' @param lambda Hyperparameter of penalty function
#' @param a Second hyperparameter of penalty functio, usually set to a = 3.7 as recommended by Fan, J. and Li, R. (2001). 
#' @param etaP Vector of sorted parameters from previous iteration.
#' @return Value of the  local approximation of the SCAD penalty, defined in Eq (16).

dscad_pen <- function(eta,lambda,a,etaP, nbAnimals){
  if(length(eta)!=length(etaP)){
    stop('eta and etaP are not the same length')
  }
  
  tildeP <- array(NA,dim=length(eta))
  Dpen <- function(x){
    if(abs(x)<=lambda){
      DD <- lambda
    }else{
      DD <- max(a*lambda-abs(eta),0)/(a-1)
    }
    return(DD)
  }
  
  for(k in (1:length(eta))){
    tildeP[k] <- nbAnimals*(pen(etaP[k],lambda,a)+Dpen(etaP[k])*(eta[k]-etaP[k]))
  }
  return(sum(tildeP))
}

#' Numerically stable forward log-probabilities for non-stationary HMM based on McClintock, B. T. (2021).
#' @param pi Stationary distribution (or its estimate if non-stationary HMM). It should be a vector of size nbStates.
#' @param trMat Transition Matrix. dim(trMat) = c(nbSteps,nbStates, nbStates).
#' @param lnProbs Diagonal matrix of log emission distributions for all observations (for all individuals).
#' @param nbSteps Length of the time-series.
#' @param aInd First observation for each individual
#' @return Value of forward log-probabilities. 

logAlphaNH <- function(pi, trMat, lnProbs, nbSteps, aInd)
{
  # Initialize parameters
  nbStates <- ncol(trMat)
  elnalpha <- array(NA,dim=c(nbSteps,nbStates))
  ltrMat  <- array(NA,dim=c(nbSteps,nbStates,nbStates))
  lpiG <- array(NA,dim=c(1,nbStates))
  nbAnimals <- length(aInd)
  animal <- 1
  for(t in (1:nbSteps)){
    # Formula for forward probability at t=1 of individual "animal" 
    if( animal <= nbAnimals && t==aInd[animal]){
      for(j in (1:nbStates)){
        # Start the likelihood at 0, log(0) = -Inf
        lpiG[j] = -Inf
        for(i in (1:nbStates)){
          # Take the log of the matrix
          ltrMat[t,i,j] = log(trMat[t,i,j])
          # Likelihood of S_1 = i and S_2 = j.
          lpiG[j] = logspace_add(lpiG[j],log(pi[i])+ltrMat[t,i,j])
        }
        elnalpha[t,j] = lpiG[j]+lnProbs[t,j]
      }
      animal <- animal+1
    } else {
      for(j in (1:nbStates)){
        logalpha = -Inf
        for(i in (1:nbStates)){
          ltrMat[t,i,j] = log(trMat[t,i,j])
          logalpha <-  logspace_add(logalpha,elnalpha[t-1,i]+ltrMat[t,i,j])
        }
        # Forward probability at time t and emission probabilities at time t.
        elnalpha[t,j] <-  logalpha+lnProbs[t,j]
      }
    }
  }
  return(elnalpha)
}


#' Numerically stable forward log-probabilities for stationary HMM based on McClintock, B. T. (2021) and http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf.
#' @param pi Stationary distribution (or its estimate if non-stationary HMM). It should be a vector of size nbStates.
#' @param trMat Transition Matrix. dim(trMat) =  c(nbStates, nbStates).
#' @param lnProbs Diagonal matrix of log emission distributions for all observations (for all individuals).
#' @param nbSteps Length of the time-series.
#' @param aInd First observation for each individual
#' @return Value of forward log-probabilities. 

logAlpha <- function(pi, trMat, lnProbs, nbSteps, aInd){
  # Initialize parameters
  nbStates <- ncol(trMat)
  elnalpha <- array(NA,dim=c(nbSteps,nbStates))
  ltrMat  <- array(NA,dim=c(nbStates,nbStates))
  lpiG <- array(NA,dim=c(1,nbStates))
  nbAnimals <- length(aInd)
  animal <- 1
  for(t in (1:nbSteps)){
    if( animal <= nbAnimals && t==aInd[animal]){
      for(j in (1:nbStates)){
        lpiG[j] = -Inf
        for(i in (1:nbStates)){
          ltrMat[i,j] = log(trMat[i,j])
          lpiG[j] = logspace_add(lpiG[j],log(pi[i])+ltrMat[i,j])
        }
        elnalpha[t,j] = lpiG[j]+lnProbs[t,j]
      }
      animal <- animal+1
    } else {
      for(j in (1:nbStates)){
        logalpha = -Inf
        for(i in (1:nbStates)){
          logalpha <-  logspace_add(logalpha,elnalpha[t-1,i]+ltrMat[i,j])
        }
        elnalpha[t,j] <-  logalpha+lnProbs[t,j]
      }
    }
  }
  return(elnalpha)
}


#' Numerically stable backward log-probabilities based on McClintock, B. T. (2021) and and http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf.
#' @param pi Stationary distribution (or its estimate if non-stationary HMM). It should be a vector of size nbStates.
#' @param trMat Transition Matrix. dim(trMat) =  c(nbStates, nbStates) for stationary, and c(nbSteps,nbStates, nbStates) for non-stationary HMM.
#' @param lnProbs Diagonal matrix of log emission distributions for all observations (for all individuals).
#' @param nbSteps Length of the time-series.
#' @param aInd First observation for each individual
#' @return Value of backward log-probabilities. 

logBetaNH <- function(pi, trMat,lnProbs,nbSteps, aInd)
{ 
  # Initialize parameters
  nbStates <- ncol(trMat)
  elnbeta <- array(NA,dim=c(nbSteps,nbStates))
  ltrMat  <- array(NA,dim=c(nbSteps,nbStates,nbStates))
  lpiG <- array(NA,dim=c(1,nbStates))
  nbAnimals <- length(aInd)
  animal <- nbAnimals
  for(j in (1:nbStates)){
    elnbeta[nbSteps,j] <- 0
    for(t in 1:nbSteps){
      for(i in (1:nbStates)){
        ltrMat[t,i,j] <-  log(trMat[t,i,j])
      }
    }
  }
  for(t in (nbSteps-1):1){
    if(animal>1 && t==(aInd[animal]-1)){
      animal <- animal-1
      for(j in (1:nbStates)){
        elnbeta[t,j] <- 0
      }
    } else {
      for(i in (1:nbStates)){
        logbeta = -Inf
        for(j in (1:nbStates)){
          logbeta <-  logspace_add(logbeta,ltrMat[t+1,i,j]+lnProbs[t+1,j]+elnbeta[t+1,j])
        }
        elnbeta[t,i] <-  logbeta
      }
    }
  }
  return(elnbeta)
}

#' Numerically stable backward log-probabilities for stationary HMM based on McClintock, B. T. (2021).
#' @param pi Stationary distribution (or its estimate if non-stationary HMM). It should be a vector of size nbStates.
#' @param trMat Transition Matrix. dim(trMat) =  c(nbStates, nbStates).
#' @param lnProbs Diagonal matrix of log emission distributions for all observations (for all individuals).
#' @param nbSteps Length of the time-series.
#' @param aInd First observation for each individual
#' @return Value of backward log-probabilities. 

logBeta <- function(pi, trMat,lnProbs,nbSteps, aInd){
  # Initialize parameters
  nbStates <- ncol(trMat)
  elnbeta <- array(NA,dim=c(nbSteps,nbStates))
  ltrMat  <- array(NA,dim=c(nbStates,nbStates))
  lpiG <- array(NA,dim=c(1,nbStates))
  nbAnimals <- length(aInd)
  animal <- nbAnimals
  for(j in (1:nbStates)){
    elnbeta[nbSteps,j] <- 0
    for(i in (1:nbStates)){
      ltrMat[i,j] <-  log(trMat[i,j])
    }
  }
  for(t in (nbSteps-1):1){
    if(animal>1 && t == (aInd[animal]-1)){
      animal <- animal-1
      for(j in (1:nbStates)){
        elnbeta[t,j] <- 0
      }
    } else {
      for(i in (1:nbStates)){
        logbeta = -Inf
        for(j in (1:nbStates)){
          logbeta <-  logspace_add(logbeta,ltrMat[i,j]+lnProbs[t+1,j]+elnbeta[t+1,j])
        }
        elnbeta[t,i] <-  logbeta
      }
    }
  }
  return(elnbeta)
}


#' Numerically stable derivation of the expectation step of u in the EM algorithm for non-stationary HMM.
#' @param pi Vector of initial distributions. It should be a vector of size nbStates.
#' @param trMat Transition Matrix. dim(trMat) = c(nbSteps,nbStates, nbStates) for non-stationary HMM.
#' @param lnProbs Diagonal matrix of log emission distributions for all observations (for all individuals).
#' @param nbSteps Length of the time-series.
#' @param aInd First observation for each individual
#' @return Expectation of vector u, conditionally to the parameters from Eq (13)

logUNH <- function(pi, trMat,lnProbs,nbSteps, aInd){
  # Initialize parameters
  nbStates <- ncol(trMat)
  elnU <- array(NA,dim=c(nbSteps,nbStates))
  elnalpha <- logAlphaNH(pi, trMat, lnProbs, nbSteps, aInd) 
  elnbeta <- logBetaNH(pi, trMat, lnProbs, nbSteps, aInd) 
  
  for(t in (1:nbSteps)){
    nn <- -Inf
    for(i in (1:nbStates)){
      elnU[t,i] <- elnalpha[t,i]+elnbeta[t,i] #because it is already the jnllk
      nn <- logspace_add(nn, elnU[t,i]) #nn was to check that normlizer = sum(alpha*beta), thus elnU is good.
    }
    for(i in (1:nbStates)){
      elnU[t,i] = elnU[t,i]  - nn
    }
  }
  return(elnU)
}


#' Numerically stable derivation of the expectation step of v in the EM algorithm for non-stationary HMM.
#' @param pi Vector of initial distributions. It should be a vector of size nbStates.
#' @param trMat Transition Matrix. dim(trMat) = c(nbSteps,nbStates, nbStates) for non-stationary HMM.
#' @param lnProbs Diagonal matrix of log emission distributions for all observations (for all individuals).
#' @param nbSteps Length of the time-series.
#' @param aInd First observation for each individual
#' @return Expectation of vector u, conditionally to the parameters from Eq (13) for non-stationary HMM.

logVNH <-  function(pi, trMat,lnProbs,nbSteps, aInd){
  #Initialize parameters
  nbStates <- ncol(trMat)
  nbAnimals <- length(aInd)
  elnV <- array(NA,dim=c((nbSteps-1),nbStates,nbStates)) #because it is define for t=2,...,T (in Sanchez, for one individual)
  elnalpha <- logAlphaNH(pi, trMat, lnProbs, nbSteps, aInd) 
  elnbeta <- logBetaNH(pi, trMat, lnProbs, nbSteps, aInd) 
  ltrMat <- array(dim=c(nbSteps,nbStates,nbStates))
  
  animal <-  1
  for(t in (1:(nbSteps-1))){
    nn <- -Inf
    if(animal < nbAnimals && t==(aInd[animal+1]-1)){
      for(i in (1:nbStates)){
        for(j in (1:nbStates)){
          elnV[t,i,j] <- 0
          nn <- logspace_add( nn, elnV[t,i,j])
        }
      }
      animal <- animal+1
    }else{
      for(i in (1:nbStates)){
        for(j in (1:nbStates)){
          ltrMat[t+1,i,j] <- log(trMat[t+1,i,j])
          elnV[t,i,j] <- elnalpha[t,i]+ltrMat[t+1,i,j]+lnProbs[t+1,j]+elnbeta[t+1,j]#+normalizer
          nn <- logspace_add( nn, elnV[t,i,j])
        }
      }
    }
    for(i in (1:nbStates)){
      for(j in (1:nbStates)){
        elnV[t,i,j] = elnV[t,i,j]  - nn
      }
    }
  }
  return(elnV)
}

#' Numerically stable derivation of the expectation step of u in the EM algorithm for stationary HMM.
#' @param pi Stationary distribution. It should be a vector of size nbStates.
#' @param trMat Transition Matrix. dim(trMat) =  c(nbStates, nbStates).
#' @param lnProbs Diagonal matrix of log emission distributions for all observations (for all individuals).
#' @param nbSteps Length of the time-series.
#' @param aInd First observation for each individual
#' @return Expectation of vector u, conditionally to the parameters from Eq (13) non-stationary HMM.

logU <- function(pi, trMat,lnProbs,nbSteps, aInd){
  # Initialize parameters
  nbStates <- ncol(trMat)
  elnU <- array(NA,dim=c(nbSteps,nbStates))
  elnalpha <- logAlpha(pi, trMat, lnProbs, nbSteps, aInd) 
  elnbeta <- logBeta(pi, trMat, lnProbs, nbSteps, aInd) 
  # normalizer <- forward_alg(pi,trMat,lnProbs,nbSteps,aInd)
  
  for(t in (1:nbSteps)){
    nn <- -Inf
    for(i in (1:nbStates)){
      elnU[t,i] <- elnalpha[t,i]+elnbeta[t,i]
      nn <- logspace_add(nn, elnU[t,i]) 
    }
    for(i in (1:nbStates)){
      elnU[t,i] = elnU[t,i]  - nn
    }
  }
  return(elnU)
}

#' Numerically stable derivation of the expectation step of v in the EM algorithm for stationary HMM.
#' @param pi Stationary distribution. It should be a vector of size nbStates.
#' @param trMat Transition Matrix. dim(trMat) =  c(nbStates, nbStates).
#' @param lnProbs Diagonal matrix of log emission distributions for all observations (for all individuals).
#' @param nbSteps Length of the time-series.
#' @param aInd First observation for each individual
#' @return Expectation of vector u, conditionally to the parameters from Eq (13).

logV <-  function(pi, trMat,lnProbs,nbSteps, aInd){
  # Initialize parameters
  nbStates <- ncol(trMat)
  nbAnimals <- length(aInd)
  elnV <- array(NA,dim=c((nbSteps-1),nbStates,nbStates)) #because it is define for t=2,...,T (in Sanchez, for one individual)
  elnalpha <- logAlpha(pi, trMat, lnProbs, nbSteps, aInd) 
  elnbeta <- logBeta(pi, trMat, lnProbs, nbSteps, aInd) 
  # normalizer <- forward_alg(pi,trMat,lnProbs,nbSteps,aInd)
  ltrMat <- array(dim=c(nbStates,nbStates))
  animal <-  1
  for(t in (1:(nbSteps-1))){
    nn <- -Inf
    if(animal < nbAnimals && t==(aInd[animal+1]-1)){
      for(i in (1:nbStates)){
        for(j in (1:nbStates)){
          elnV[t,i,j] <- 0
          nn <- logspace_add( nn, elnV[t,i,j])
        }
      }
      animal <- animal+1
    }else{
      for(i in (1:nbStates)){
        for(j in (1:nbStates)){
          ltrMat[i,j] = log(trMat[i,j])
          elnV[t,i,j] <- elnalpha[t,i]+ltrMat[i,j]+lnProbs[t+1,j]+elnbeta[t+1,j]#+normalizer
          nn <- logspace_add( nn, elnV[t,i,j])
        }
      }
    }
    for(i in (1:nbStates)){
      for(j in (1:nbStates)){
        elnV[t,i,j] = elnV[t,i,j]  - nn
      }
    }
  }
  return(elnV)
}


#' Compute the part of (II) of Eq (14) that is not penalized by the SCAD. It is the same for stationary and non-stationary HMMs.
#' @param lsd Log of the standard deviation parameter of the gamma distribution.
#' @param lmu Log of the mean of the gamma distribution.
#' @param pi Stationary distribution (or its estimate if non-stationary HMM). It should be a vector of size nbStates.
#' @param trMat Diagonal Matrix of the HMM.
#' @param nbSteps Length of time-series
#' @param aInd First observation for each individual
#' @param elnU Expectation of vector u, conditionally to the parameters from Eq (13), obtained from logU function.
#' @return (II) without the SCAD penalty function 

ComputeIISd <- function(lsd,lmu,pi, trMat,nbSteps, aInd,elnU){
  
  mu <- exp(lmu)
  sd <- exp(lsd)
  shape = (mu * mu)/ (sd * sd)
  scale = (sd * sd/ mu)
  
  nbStates <- ncol(trMat)
  # to prevent for overflow max_x(f(x)) = max_x(g=f(x)*a), for a positive constant not depending on x. Here, a is exp(-max(elnU)).
  b <- max(elnU) 
  Theta <- t(matrix(c(shape,scale),ncol=2,nrow=nbStates))
  
  lnProbs <- array(0,dim=c(nbSteps,nbStates))
  lnProbs[1,] <- apply(Theta,2,FUN=function(x)dgamma(Obs[1,1],shape=x[1],scale=x[2],log=TRUE))
  for(t in (2:nbSteps)){
    if(!is.na(Obs[t,1])){lnProbs[t,] <- apply(Theta,2,FUN=function(x)dgamma(Obs[t,1],shape=x[1],scale=x[2],log=TRUE))}
  } 
  lll <- array(NA,dim=c(nbSteps,nbStates))
  for(t in (1:nbSteps)){
    for(i in (1:nbStates)){
      lll[t,i] <- lnProbs[t,i]*exp(elnU[t,i]-b)
    }
  }
  
  # we take negative bc optim does minimization by default
  lll <- -sum(lll) 
  return(lll)
}



#' Optimizes Eq (14) with respect to the standard deviation of the gamma distribution.
#' @param lmu Current estimate of the mean of the gamma distributions.
#' @param lsd0 Initial value of the gamma standard deviation parameters
#' @param pi Stationary distribution (or its estimate if non-stationary HMM). It should be a vector of size nbStates.
#' @param trMat Diagonal Matric of the HMM.
#' @param nbSteps Length of time-series
#' @param aInd First observation for each individual
#' @param elnU Expectation of vector u, conditionally to the parameters from Eq (13), obtained from logU function.
#' @return Values for standard deviation of gamma distributions that maximizes Eq (14).

SdUpdate <- function(lmu,lsd0,pi, trMat,nbSteps, aInd,elnU)
{
  res <- optim(lsd0,ComputeIISd,lmu=lmu,pi=pi, trMat=trMat,nbSteps=nbSteps, aInd=aInd,elnU=elnU,method=c("BFGS"))
  return(res$par)
} 


#' Beta coefficient to transition probability matrix
#' @param nbStates Number of states of the HMM
#' @param beta Regression coefficient for the covariates in the t.p.m
#' @param X Vector of covariates, the first column should be intercept. Hence ncol(X) = 1+ "Nb Covariates".
#' @param nbSteps Length of time-series
#' @return Transition Matrix.

HMM.beta2tpm <- function(nbStates,beta,X,nbSteps)
{
  beta <- matrix(beta,ncol=ncol(X))
  tpm <- array(NA,dim=c(nbSteps,nbStates,nbStates))
  pi <- array(dim=c(nbSteps,nbStates))
  # compute the sum of the product of the covariates and the regression coefficient
  prod <- X%*%t(beta)
  # Take the exponential as in Eq (5) minus the maximum of the sum of the product to prevent from overflow
  qq <- exp(prod-max(prod))
  # Set off-diagonal elements of matrix, starting from 1st col and moving down
  for(t in 1:nbSteps){
    int <- tpm[t,,]
    diag(int) <- exp(-max(prod)) # since qq is Eq (5) times exp(-max(prod)) was computed line 661, it has to be multiplied by exp(-max(prod)) everywhere.
    # Take only off-diagonal elements of resulting matrix
    int[col(int)!=row(int)] <- qq[t,]
    int = t(int)
    # to prevent value to be exactly equal to zero, since the logarithm of trMat is used later (log(0) does not exist).
    tpm[t,,] <- int+1e-15 
    # Divide each element by the row sum, so the rows will sum to 1
    tpm[t,,] <- tpm[t,,]/apply(tpm[t,,],1,sum)
  }
  return(list(tpm=tpm))
}


#' Transform natural parameters for trMat to working parameters for stationary HMMs.
#' @param trMat t.p.m of dim = c(nbStates,nbStates)
#' @param aInd First observation for each individual
#' @return Working parameters for trMat for stationary HMMs.

pn2pw <- function(trMat,aInd){
  nbAnimals <- length(aInd)
  nbStates <- ncol(trMat)
  #Transform gamma by log(gamma_ij/gamma_ii)
  foo <- log(trMat/diag(trMat))
  #Take only off-diagonal elements of resulting matrix
  ttrMat <- as.vector(foo[!diag(nbStates)])
  return(ttrMat)
}


#' Transform working parameters to trMat for stationary HMMs.
#' @param ttrMat Working parameters for trMat.
#' @param aInd First observation for each individual
#' @param nbStates Number of states, i.e., upper bound of the MPLE
#' @return trMat for stationary HMMs.

pw2pn <-  function(ttrMat,nbStates){
  trMat <- array(dim=c(nbStates,nbStates))
  pi <- array(dim=c(nbStates))
  gamma<- diag(nbStates)
  # Set off-diagonal elements of matrix, starting from 1st col and moving down
  gamma[!gamma] <- exp(ttrMat)
  # Divide each element by the row sum, so the rows will sum to 1
  gamma <- gamma/apply(gamma,1,sum)
  trMat <- gamma
  return(trMat)
}


#' Compute the part of (I) of Eq (14) for non-stationary HMM (i.e., with an estimate of pi).
#' @param beta Regression coefficient for the covariates in the t.p.m
#' @param X Vector of covariates, the first column should be intercept. Hence ncol(X) = 1+ "Nb Covariates".
#' @param trMat Current estimate of the t.p.m.
#' @param nbSteps Length of time-series
#' @param aInd First observation for each individual
#' @param pi Stationary distribution (or its estimate if non-stationary HMM). It should be a vector of size nbStates.
#' @param lnProbs Diagonal matrix of log emission distributions for all observations (for all individuals).
#' @param Cn Hyper parameter of the MPLE.
#' @param elnV Expectation of vector u, conditionally to the parameters from Eq (13), obtained from logV function.
#' @return Value of (I) for non-stationary HMMs.

ComputeItpmNH<- function(beta,X,trMat,nbSteps, aInd,pi,lnProbs,Cn,elnV)
{
  # Initialize parameters
  nbStates <- ncol(trMat)
  # The maximum of the log(v) vector is used to prevent from overflow.
  b <- max(elnV) 
  nbAnimals = length(aInd)
  tpm <- HMM.beta2tpm(nbStates,beta,X,nbSteps)$tpm
  
  lProbaState = array(dim=c(nbSteps,nbStates)) 
  lalpha = logAlphaNH(pi,tpm,lnProbs,nbSteps,aInd)
  lbeta = logBetaNH(pi,tpm,lnProbs,nbSteps,aInd)
  
  # Compute HatPi, based on forward and backward probabilities.
  II <- array(0,dim=c(nbStates,nbStates))
  animal=1
  for(t in (1:(nbSteps-1))){
    if(animal<= nbAnimals && t==aInd[animal]){
      if(animal==nbAnimals){negllk = forward_alg_ind(pi,tpm[aInd[animal]:nbSteps,,],lnProbs[aInd[animal]:nbSteps,],nrow(lnProbs[aInd[animal]:nbSteps,]),animal,nbAnimals)}
      else{negllk = forward_alg_ind(pi,tpm[aInd[animal]:(aInd[animal+1]-1),,],lnProbs[aInd[animal]:(aInd[animal+1]-1),],nrow(lnProbs[aInd[animal]:(aInd[animal+1]-1),]),animal,nbAnimals)}
      animal=animal+1
    }
    # Compute I from Eq (14)
    for(i in (1:nbStates)){
      lProbaState[t,i] = lalpha[t,i] + lbeta[t,i] +negllk
      for(j in (1:nbStates)){
        # Multiply Eq (14) by exp(-b) to prevent from overflow
        II[i,j] <- II[i,j]+ exp(elnV[t,i,j]-b)*log(tpm[t,i,j])
      }
    }
  }
  for(j in (1:nbStates)){
    lProbaState[nbSteps,j] = lalpha[nbSteps,j] + lbeta[nbSteps,j]  +negllk
  }
  # Hatpi = exp(lProbaState)
  HatPi = apply(exp(lProbaState),2,mean)
  # Multiply the other side of Eq (14) by exp(-b), to obtain equivalent problem of maximization than maximizing Eq (14) w respect to beta 
  return(-sum(II)-exp(-b)*Cn*sum(log(HatPi)))
}

#' Maximize (I) with respect to beta
#' @param beta0 Initial value for regression coefficient for the covariates in the t.p.m
#' @param X Vector of covariates, the first column should be intercept. Hence ncol(X) = 1+ "Nb Covariates".
#' @param trMat Current estimate of the t.p.m.
#' @param nbSteps Length of time-series
#' @param aInd First observation for each individual
#' @param pi Stationary distribution (or its estimate if non-stationary HMM). It should be a vector of size nbStates.
#' @param lnProbs Diagonal matrix of log emission distributions for all observations (for all individuals).
#' @param Cn Hyper parameter of the MPLE.
#' @param elnV Expectation of vector v, conditionally to the parameters from Eq (13), obtained from logV function.
#' @return Values of beta that maximizes (I).


UpdateBeta <- function(beta0,X,trMat,nbSteps, aInd,pi,lnProbs,Cn,elnV)
{
  res = optim(beta0,ComputeItpmNH,X=X,trMat=trMat,nbSteps=nbSteps, aInd=aInd,pi=pi,lnProbs=lnProbs,Cn=Cn,elnV=elnV,method=c("BFGS"))
  return(res$par)
}  

#' Compute the part of (I) of Eq (14), stationary HMM.
#' @param ttrMat Working parameters for the t.p.m
#' @param nbSteps Length of time-series
#' @param aInd First observation for each individual
#' @param elnU Expectation of vector u, conditionally to the parameters from Eq (13), obtained from logU function.
#' @param elnV Expectation of vector v, conditionally to the parameters from Eq (13), obtained from logV function.
#' @param nbStates Number of states.
#' @param Cn Hyperparameter of MPLE.
#' @return Value of (I) 

ComputetrMat <- function(ttrMat,nbSteps, aInd,elnU,elnV,nbStates,Cn)
{
  
  nbAnimals <- length(aInd)
  trMat <- pw2pn(ttrMat,nbStates)
  ltrMat <- array(NA,dim=c(nbStates,nbStates))
  ll <- lll <- 0
  pi<-base::solve(t(diag(nbStates)-trMat+1),rep(1,nbStates))
  # The maximum of the log(u) and log(v) vectors is used to prevent from overflow.
  b <- max(elnU,elnV)
  for(animal in (1:nbAnimals)){
    if(animal<nbAnimals){
      TT <- aInd[animal+1]-1
    }else{
      TT <- nbSteps-1
    }
    for(i in (1:nbStates)){
      # Multiply Eq (14) by exp(-b) to prevent from overflow
      lll <- lll+exp(elnU[aInd[animal],i]-b)*log(pi[i])
      # Compute I from Eq (14)
      for(j in (1:nbStates)){
        ltrMat[i,j] = log(trMat[i,j])
        for(t in (aInd[animal]:TT)){
          ll <- ll+exp(elnV[t,i,j]-b)*ltrMat[i,j]
        }
      }
    }
  }
  # Multiply the other side of Eq (14) by exp(-b), to obtain equivalent problem of maximization than maximizing eq (14) w respect to ttrMat 
  nllk <-  -exp(-b)*Cn*sum(log(pi))-ll-lll
  
  return(nllk)
}


#' Update t.p.m (stationary HMM)
#' @param ttrMat Initival value for working parameters for the t.p.m
#' @param nbSteps Length of time-series
#' @param aInd First observation for each individual
#' @param elnU Expectation of vector u, conditionally to the parameters from Eq (13), obtained from logU function.
#' @param elnV Expectation of vector v, conditionally to the parameters from Eq (13), obtained from logV function.
#' @param nbStates Number of states.
#' @param Cn Hyperparameter of MPLE.
#' @return Updated working parameters.

trMatUpdate <- function(ttrMat,nbSteps, aInd,elnU,elnV,nbStates,Cn){
  res <- optim(ttrMat,ComputetrMat,nbSteps=nbSteps,aInd=aInd,elnU=elnU,elnV=elnV,nbStates=nbStates,Cn=Cn)
  return(res$par)
}


#' Update the vector of initial probabilities from Zucchini et al.
#' @param pi Stationary distribution (or its estimate if non-stationary HMM). It should be a vector of size nbStates.
#' @param trMat Transition Matrix. dim(trMat) =  c(nbStates, nbStates) for stationary, and c(nbSteps,nbStates, nbStates) for non-stationary HMM.
#' @param lnProbs Diagonal matrix of log emission distributions for all observations (for all individuals).
#' @param nbSteps Length of the time-series.
#' @param aInd First observation for each individual
#' @return Updated vector of initial probabilities.

DeltaUpdate <- function(pi, trMat,lnProbs,nbSteps, aInd)
{
  nbStates <- ncol(trMat)
  nbAnimals <- length(aInd)
  elnU <- logUNH(pi, trMat,lnProbs,nbSteps, aInd)
  pi <- array(0,dim=c(nbStates))
  
  for(i in (1:nbStates)){
    # page 72 of Zucchini et al.
    pi[i] <- sum(sapply(c(1:nbAnimals),FUN=function(x)exp(elnU[aInd[x],i])))
  }
  pi <- pi/nbAnimals
  return(pi)
}


#' Compute the part of (II) of Eq (14) that is penalized by the SCAD. Here: mean of step-length (mu) and concentration of turning angle (kappa)
#' @param lpar Log of the vector of parameter estimates that are penalized by the SCAD function. ncol(lpar) = number of parameters to penalize.
#' @param lds Log of the standard deviation of gamma distribution.
#' @param etaP Vector of sorted parameters (mu,kappa) from current estimates
#' @param pi Stationary distribution (or its estimate if non-stationary HMM). It should be a vector of size nbStates.
#' @param lambda Hyperparameter of penalty function
#' @param a Hyperparameter of penalty function, a = 3.7.
#' @param Obs Vector of observations. 
#' @param trMat Transition Matrix. dim(trMat) =  c(nbStates, nbStates) for stationary, and c(nbSteps,nbStates, nbStates) for non-stationary HMM.
#' @param nbSteps Length of time-series
#' @param aInd First observation for each individual
#' @param elnU Expectation of vector u, conditionally to the parameters from Eq (13), obtained from logU function.
#' @return (II) from Eq (14) without the SCAD penalty function 

ComputeIIMuKappaNH <- function(lpar,lsd,etaP,pi,lambda,a,Obs, trMat,nbSteps, aInd,elnU){
  # Initialize parameters
  nbStates <- ncol(trMat)
  nbAnimals <- length(aInd)
  par <- matrix(exp(lpar),ncol=2,nrow=nbStates)
  mu <- par[,1]
  kappa <- par[,2]
  sd <- exp(lsd)
  # from mean and sd, obtain shape and scale of gamma distribution
  shape = (mu * mu)/ (sd * sd)
  scale = (sd * sd/ mu)
  # The maximum of the log(u) vector is used to prevent from overflow.
  b <- max(elnU)
  
  Theta <- t(matrix(c(shape,scale),ncol=2,nrow=nbStates))
  # Diagonal matrix of log emission distributions for all observations (for all individuals).
  lnProbs <- array(0,dim=c(nbSteps,nbStates))
  lnProbs[1,] <- apply(Theta,2,FUN=function(x)dgamma(Obs[1,1],shape=x[1],scale=x[2],log=TRUE))
  for(t in (2:nbSteps)){
    if(!is.na(Obs[t,1])){lnProbs[t,] <- apply(Theta,2,FUN=function(x)dgamma(Obs[t,1],shape=x[1],scale=x[2],log=TRUE))}
    if(!is.na(Obs[t,2])){lnProbs[t,] <- lnProbs[t,]+sapply(kappa,FUN=function(x)dvon_mises(Obs[t,2],mu=0,kappa=x,log=TRUE))}
  }
  # Compute II
  lll <- array(NA,dim=c(nbSteps,nbStates))
  for(t in (1:nbSteps)){
    for(i in (1:nbStates)){
      # Multiply by exp(-b) to prevent from overflow
      lll[t,i] <- lnProbs[t,i]*exp(elnU[t,i]-b)
    }
  }
  # Sort the parameters
  idx <- ordering(nbStates,par)
  par <- par[idx, ]
  # Vector of differences of sorted parameters
  eta <- numeric(nrow(par)-1)
  
  # Compute L2 norm of differences between successive rows
  for (i in 1:(nrow(par)-1)) {
    eta[i] <- sqrt(sum((par[i+1, ] - par[i, ])^2))
  }
  # Compute II & Multiply the other side of Eq (14) by exp(-b)
  lll <- -sum(lll)+exp(-b)*dscad_pen(eta,lambda,a,etaP,nbAnimals)
  
  return(lll)
}


#' Maximize II with respect to the mean of step-length (mu) and concentration of turning angle (kappa): multi-dimension version
#' @param lpar Log of the vector of parameter estimates that are penalized by the SCAD function. ncol(lpar) = number of parameters to penalize.
#' @param lds Log of the standard deviation of gamma distribution.
#' @param etaP Vector of sorted parameters (mu,kappa) from current estimates
#' @param pi Stationary distribution (or its estimate if non-stationary HMM). It should be a vector of size nbStates.
#' @param lambda Hyperparameter of penalty function
#' @param a Hyperparameter of penalty function, a = 3.7.
#' @param Obs Vector of observations. 
#' @param trMat Transition Matrix. dim(trMat) =  c(nbStates, nbStates) for stationary, and c(nbSteps,nbStates, nbStates) for non-stationary HMM.
#' @param nbSteps Length of time-series
#' @param aInd First observation for each individual
#' @param elnU Expectation of vector u, conditionally to the parameters from Eq (13), obtained from logU function.
#' Updated mean and concentration.


MuKappaUpdate <-  function(lpar,lsd,etaP,pi,lambda,a,Obs, trMat,nbSteps, aInd,elnU){
  res <- optim(lpar,ComputeIIMuKappaNH,lsd=lsd,etaP=etaP,pi=pi,lambda=lambda,a=a,Obs=Obs, trMat=trMat,
               nbSteps=nbSteps, aInd=aInd,elnU=elnU)
  return(res$par)
}


#' Compute II with the SCAD penalty
#' @param lmu  mean of the gamma distributions.
#' @param lsd  standard deviation of the gamma distributions.
#' @param lmuP  current estimate (sorted) of the mean parameter of the gamma distributions.
#' @param pi Stationary distribution (or its estimate if non-stationary HMM). It should be a vector of size nbStates.
#' @param lambda Hyperparameter of penalty function
#' @param a Hyperparameter of penalty function, a = 3.7.
#' @param trMat Transition Matrix. dim(trMat) =  c(nbStates, nbStates) for stationary, and c(nbSteps,nbStates, nbStates) for non-stationary HMM.
#' @param nbSteps Length of time-series
#' @param aInd first observation for each individual
#' @param elnU Expectation of vector u, conditionally to the parameters from Eq (13), obtained from logU function.
#' @param Obs Vector of observations. 
#' Value of II.

ComputeIIMu <- function(lmu,lsd,lmuP,pi,lambda,a, trMat,nbSteps, aInd,elnU,Obs)
{ 
  # Initialize parameters
  nbAnimals <- length(aInd)
  mu <- exp(lmu)
  sd <- exp(lsd)
  # from mean and sd, obtain shape and scale of gamma distribution
  shape = (mu * mu)/ (sd * sd)
  scale = (sd * sd/ mu)
  # The maximum of the log(u) vector is used to prevent from overflow.
  b <- max(elnU)
  NbStates <- ncol(trMat)
  
  
  Theta <- t(matrix(c(shape,scale),ncol=2,nrow=NbStates))
  # Diagonal matrix of log emission distributions for all observations (for all individuals).
  lnProbs <- apply(Theta,2,FUN=function(x)dgamma(Obs,shape=x[1],scale=x[2],log=TRUE))  
  lll <- array(NA,dim=c(nbSteps,NbStates))
  for(t in (1:nbSteps)){
    for(i in (1:NbStates)){      
     # Multiply by exp(-b) to prevent from overflow
      lll[t,i] <- lnProbs[t,i]*exp(elnU[t,i]-b)
    }
  }
  # Sort parameters estimates from previous iteration
  muP <- exp(lmuP)
  etaP <- sort(muP,decreasing=TRUE)
  # Sort current parameters estimates
  eta <- sort(mu,decreasing=TRUE)
  #Difference of sorted parameters
  eta <- c(rev(abs(diff(eta))))
  etaP <- c(rev(abs(diff(etaP))))
  # Multiply the other side of Eq (14) by exp(-b)
  # Compute II.
  lll <- -sum(lll)+exp(-b)*dscad_pen(eta,lambda,a,etaP, nbAnimals)
  
  return(lll)
}


#' Maximize II with respect to the mean of step-length (mu): univariate version
#' @param lmu0  Initial value for mean of the gamma distributions.
#' @param lsd  Standard deviation of the gamma distributions.
#' @param lmuP  Current estimate (sorted) of the mean parameter of the gamma distributions.
#' @param pi Stationary distribution (or its estimate if non-stationary HMM). It should be a vector of size nbStates.
#' @param lambda Hyperparameter of penalty function
#' @param a Hyperparameter of penalty function, a = 3.7.
#' @param trMat Transition Matrix. dim(trMat) =  c(nbStates, nbStates) for stationary, and c(nbSteps,nbStates, nbStates) for non-stationary HMM.
#' @param nbSteps Length of time-series
#' @param aInd First observation for each individual
#' @param elnU Expectation of vector u, conditionally to the parameters from Eq (13), obtained from logU function.
#' @param Obs Vector of observations. 
#' Updated mean parameter.


MuUpdate <-  function(lmu0,lsd,lmuP,pi,lambda,a, trMat,nbSteps, aInd,elnU,Obs)
{
  res <- optim(lmu0,ComputeIIMu,lsd=lsd,lmuP=lmuP,pi=pi,lambda=lambda,a=a, trMat=trMat,
               nbSteps=nbSteps, aInd=aInd,elnU=elnU,Obs=Obs)
  return(res$par)
}

