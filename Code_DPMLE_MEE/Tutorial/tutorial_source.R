#' Function to fit both DPMLE methods from the paper "Improved order selection method for hidden Markov models: a case study with movement data"
source("sourcefunctions.R")

#' DPMLE function: EM algorithm to perform double penalized maximum likelihood with gamma or gamma + von Mises. "function.NH" means that the non-stationary DPMLE is applied. No ".MH" means the stationary DPMLE is applied.
#' @param NbIter Number maximum of iterations for the EM algorithm
#' @param epsilon Stopping rule.
#' @param Pi0 Initial distribution for stationary (stationary HMM) or initial distributions (non-stationary HMM). It should be a vector of size nbStates.
#' @param trMat0 Initial t.p.m of dimension nbSteps x nbStates x nbStates if non-stationary HMM and nbStates x nbStates otherwise.
#' @param beta0 Initial regression vector for the t.p.m, only for non-stationary DPMLE.
#' @param lmu0 Logarithm of the initial mean values for gamma distributions
#' @param lsd0 Logarithm of the initial standard deviation values for gamma distributions
#' @param lkappa0 (Optional) Logarithm of the concentration parameters of von Mises
#' @param nbSteps Time-length of the time-series.
#' @param aInd Vector of time of first observation for each individual (Vector of size NbAnimals)
#' @param Obs Vector of observations, with ncol = number of datastreams.
#' @param lambda Hyperparameter for double penalized likelihood, for the SCAD function.
#' @param a a = 3.7
#' @param Cn Hyperparameter for double penalized likelihood, for the penalty applied on the stationary probabilities.
#' @param type list of two elements to describe which method to be applied: "multi-uni" dimensional and "NS-S" for non-stationary or stationary.
#' @return List of estimated parameters: emission distribution, transition matrix, stationary probability + initial probability if relevant
#' vector of negative joint log likelihood, vector of double penalized likelihood, estimated number of states, vector of diagonal matrix of log emissions,
#' values of hyperparameters for double penalized likelihood: lambda and Cn, and convergence: 0 if it converged, 1 otherwise.

gamma.EMpen <- function(NbIter,epsilon,Pi0,trMat0,lmu0,lsd0,nbSteps,aInd,Obs,lambda,a,Cn){
  
  convergence = 1 #1 if did not converged, 0 o/w.
  #Initialize parameters
  nllk <- nllpen <- array()
  pi <- Pi0
  nbAnimals <- length(aInd)
  trMat <- trMat0
  nbStates <- ncol(trMat)
  
  mu0 <- exp(lmu0)
  sd0 <- exp(lsd0)
  
  if(length(lmu0) != nbStates){
    stop(paste0("The number of parameters don't match, it should be", nbStates, "parameters."))
  }
  
  if(length(lsd0) != nbStates){
    stop(paste0("The number of parameters don't match, it should be", nbStates, "parameters."))
  }
  
  #from mean and sd, obtain shape and scale of gamma distribution
  lshape = log((mu0 * mu0 )/ (sd0 * sd0))
  lscale = log(sd0 * sd0 / mu0)
  
  #Vector of log parameters of the gamma distributions mean and sd
  lXi <- t(matrix(c(lmu0,lsd0),ncol=2,nrow=nbStates))
  
  #Vector of log parameters of the gamma distributions shape and scale
  lTheta <- t(matrix(c(lshape,lscale),ncol=2,nrow=nbStates))
  
  #Diagonal matrix of log emission distributions for all observations (for all individuals).
  lnProbs <- array(0,dim=c(nbSteps,nbStates))
  lnProbs[1,] <- apply(exp(lTheta),2,FUN=function(x)dgamma(Obs[1,1],shape=x[1],scale=x[2],log=TRUE))
  for(t in (2:nbSteps)){
    if(!is.na(Obs[t,1])){lnProbs[t,] <- apply(exp(lTheta),2,FUN=function(x)dgamma(Obs[t,1],shape=x[1],scale=x[2],log=TRUE))}
  }

  for(i in (1:NbIter)){
    # get the working parameter for t.p.m
    ttrMat <- pn2pw(trMat,aInd)
    
    # Update u and v with current parameters
    elnU <- logU(pi, trMat,lnProbs,nbSteps, aInd)
    elnV <- logV(pi, trMat,lnProbs,nbSteps, aInd)
    
    # Update working parameter for t.p.m
    ttrMat_new <-  trMatUpdate(ttrMat,nbSteps, aInd,elnU,elnV,nbStates,Cn)
    # Update t.p.m
    trMat_new <- pw2pn(ttrMat_new,nbStates)
    
    # Update stationary distribution from updated t.p.m
    pi_new<- base::solve(t(diag(nbStates)-trMat_new+1),rep(1,nbStates))
    
    # Update sd
    lsd_new <-  SdUpdate(lXi[1,],lXi[2,],pi, trMat,nbSteps, aInd,elnU)
    
    # Update all parameters before updating mu
    
    # Update the logshape
    lTheta[1,] <-lXi[1,]+lXi[1,] - (lsd_new+lsd_new)
    # Update the logscale
    lTheta[2,] <-lsd_new+lsd_new-lXi[1,]
    # Update standard deviation
    lXi[2,] <- lsd_new
    # Update the t.p.m and stationary probability
    trMat <- trMat_new
    pi  <- pi_new
    
    # Diagonal matrix of log emission distributions for all observations (for all individuals).
    
    lnProbs <- array(0,dim=c(nbSteps,nbStates))
    lnProbs[1,] <- apply(exp(lTheta),2,FUN=function(x)dgamma(Obs[1,1],shape=x[1],scale=x[2],log=TRUE))
    for(t in (2:nbSteps)){
      if(!is.na(Obs[t,1])){lnProbs[t,] <- apply(exp(lTheta),2,FUN=function(x)dgamma(Obs[t,1],shape=x[1],scale=x[2],log=TRUE))}
    }
    
    # Update logU
    elnU <- logU(pi, trMat,lnProbs,nbSteps, aInd)

    # Current log estimate of mu
    lmuP <- lXi[1,] 
    
    # Updated log estimate of mu
    lmu_new <- MuUpdate(lXi[1,],lXi[2,],lmuP,pi,lambda,a, trMat,nbSteps, aInd,elnU)

    # Finish updating parameters
    lTheta[1,] <- (lmu_new+lmu_new)-(lsd_new+lsd_new)
    lTheta[2,] <-lsd_new+lsd_new-lXi[1,]
  
    lXi[1,] <- lmu_new
    
    # Diagonal matrix of log emission distributions for all observations (for all individuals).
    lnProbs <- array(0,dim=c(nbSteps,nbStates))
    lnProbs[1,] <- apply(exp(lTheta),2,FUN=function(x)dgamma(Obs[1,1],shape=x[1],scale=x[2],log=TRUE))
    for(t in (2:nbSteps)){
      if(!is.na(Obs[t,1])){lnProbs[t,] <- apply(exp(lTheta),2,FUN=function(x)dgamma(Obs[t,1],shape=x[1],scale=x[2],log=TRUE))}
    }    
    mu <- exp(lmu_new)
    sd <- exp(lsd_new)
    
    # Updated joint negative log likelihood
    nllk[i] <- forward_alg(pi,trMat,lnProbs,nbSteps,aInd)
    
    
    # Updated double penalized likelihood and estimated number of states K (hatN in the paper).
    eta <- sort(mu,decreasing=TRUE)
    eta <- c(rev(abs(diff(eta))))
    
    # threshold for overlaps is 0.1, it is arbitrary, does not affect the algorithm and can be changed after fitting model.
    idx <- which(abs(eta) <= 0.1)
    eta[idx]=0
    K <- length(which(eta!=0))+1
    eta <- array(eta)
    
    nllpen[i] <- forward_alg(pi,trMat,lnProbs,nbSteps,aInd)-Cn*sum(log(pi))+nbAnimals*sum(sapply(eta,FUN=function(x)pen(x,lambda,a)))

    # stopping rule
    if(i>1 && abs(nllpen[i]-nllpen[i-1])<epsilon){
      convergence = 0
      break()
    }
  }
  
  return(list(theta=c(mu=mu,sd=sd),pi=pi,trMat=trMat,nllk=nllk,nllpen=nllpen,K=K,lnProbs=lnProbs,lambda=lambda,Cn=Cn,convergence=convergence))
}

gamma.EMpen.NH <- function(NbIter,epsilon,Pi0,beta0,lmu0,lsd0,nbSteps,aInd,Obs,lambda,a,Cn,X){
  convergence <- 1 # 1 if did not converged, 0 o/w.
  

  if(ncol(X) != ncol(beta0)){
    stop("Covariates and coefficient dimensions do not match.")
  }
  if(nrow(beta0)!=(nbStates*(nbStates-1))|| length(mu0)!= nbStates || length(sd0)!= nbStates){
    stop("Parameter has the wrong number of states.")
  }
  
  # Initialize parameters
  nbStates <- length(Pi0)
  nbAnimals <- length(aInd)
  
  nllk <- nllpen <- array()
  pi <- Pi0
  
  mu0 <- exp(lmu0)
  sd0 <- exp(lsd0)
  
  # from mean and sd, obtain shape and scale of gamma distribution
  lshape = log((mu0 * mu0 )/ (sd0 * sd0))
  lscale = log(sd0 * sd0 / mu0)
  
  # Vector of log parameters of the gamma distributions shape and scale
  lTheta <- t(matrix(c(lshape,lscale),ncol=2,nrow=nbStates))
  # Diagonal matrix of log emission distributions for all observations (for all individuals).
  lnProbs <- array(0,dim=c(nbSteps,nbStates))
  lnProbs[1,] <- apply(exp(lTheta),2,FUN=function(x)dgamma(Obs[1,1],shape=x[1],scale=x[2],log=TRUE))
  for(t in (2:nbSteps)){
    if(!is.na(Obs[t,1])){lnProbs[t,] <- apply(exp(lTheta),2,FUN=function(x)dgamma(Obs[t,1],shape=x[1],scale=x[2],log=TRUE))}
  }
  
  pi_new <- pi
  beta <- c(beta0)
  
  # Get t.p.m from beta
  pp <- HMM.beta2tpm(nbStates,beta0,X,nbSteps)
  trMat <- pp$tpm
  
  #Vector of log parameters of the gamma distributions: logmean and logsd
  lXi <- t(matrix(c(lmu0,lsd0),ncol=2,nrow=nbStates))
  
  for(i in (1:NbIter)){
    # print(i)
    
    # Update u and v with current parameters
    elnU <- logUNH(pi, trMat,lnProbs,nbSteps, aInd)
    elnV <- logVNH(pi, trMat,lnProbs,nbSteps, aInd)
    
    # Update beta and initial distribution
    beta_new <- UpdateBeta(beta,X,trMat,nbSteps, aInd,pi,lnProbs,Cn,elnV)
    pi_new <- DeltaUpdate(pi, trMat,lnProbs,nbSteps, aInd)
    
    # Update standard deviation of gamma distribution
    lsd_new <-  SdUpdate_NH(lXi[1,],lXi[2,],pi, trMat,nbSteps, aInd,elnU)
    
    # Update parameters before updating mu
    lTheta[2,] <-lsd_new+lsd_new-lXi[1,]
    lTheta[1,] <-lXi[1,]+lXi[1,]-(lsd_new+lsd_new)
    lXi[2,] <- lsd_new
    # get updated t.p.m from updated beta
    trMat <-  HMM.beta2tpm(nbStates,beta_new,X,nbSteps)$tpm
    pi  <- pi_new
    beta <- beta_new
    
    # Diagonal matrix of log emission distributions for all observations (for all individuals).
    lnProbs <- array(0,dim=c(nbSteps,nbStates))
    lnProbs[1,] <- apply(exp(lTheta),2,FUN=function(x)dgamma(Obs[1,1],shape=x[1],scale=x[2],log=TRUE))
    for(t in (2:nbSteps)){
      if(!is.na(Obs[t,1])){lnProbs[t,] <- apply(exp(lTheta),2,FUN=function(x)dgamma(Obs[t,1],shape=x[1],scale=x[2],log=TRUE))}
    }
    
    elnU <- logUNH(pi_new, trMat,lnProbs,nbSteps, aInd)
    # Current log estimate of mu
    lmuP <- lXi[1,] 
    # Update mu
    lmu_new <- MuUpdate(lXi[1,],lXi[2,],lmuP,pi,lambda,a, trMat,nbSteps, aInd,elnU,Obs)
    
    
    # Finish updating parameters
    lTheta[1,] <- (lmu_new+lmu_new)-(lsd_new+lsd_new)
    lXi[1,] <- lmu_new
    lTheta[2,] <-lsd_new+lsd_new-lXi[1,]
    
    #Updated matrix of log emission distributions for all observations (for all individuals).
    lnProbs <- array(0,dim=c(nbSteps,nbStates))
    lnProbs[1,] <- apply(exp(lTheta),2,FUN=function(x)dgamma(Obs[1,1],shape=x[1],scale=x[2],log=TRUE))
    for(t in (2:nbSteps)){
      if(!is.na(Obs[t,1])){lnProbs[t,] <- apply(exp(lTheta),2,FUN=function(x)dgamma(Obs[t,1],shape=x[1],scale=x[2],log=TRUE))}
    }
    
    #Update mu and sd
    mu <- exp(lmu_new)
    sd <- exp(lsd_new)
    
    eta <- sort(mu,decreasing=TRUE)
    eta <- c(rev(abs(diff(eta))))
    
    #threshold for overlaps is 0.1, it is arbitrary, does not affect the algorithm and can be changed after fitting model.
    idx <- which(abs(eta) <= 0.1)
    eta[idx]=0
    #Estimate of the number of states (hatN in the paper)
    K <- length(which(eta!=0))+1
    
    #Joint log-likelihood
    alpha <- forward_alg(pi,trMat,lnProbs,nbSteps,aInd)
    
    lProbaState = array(dim=c(nbSteps,nbStates)) 
    # Updated forward probability
    lalpha = logAlphaNH(pi,trMat,lnProbs,nbSteps,aInd)
    # Updated backward probability
    lbeta = logBetaNH(pi,trMat,lnProbs,nbSteps,aInd)
    
    # compute estimate of stationary distribution denoted HatPi
    animal=1
    for(t in (1:(nbSteps-1))){
      if(animal<=nbAnimals && t==aInd[animal]){
        # Compute likelihood of individual "animal" (i.e., the one observed at time t).
        if(animal==nbAnimals){negllk = forward_alg_ind_NH(pi,trMat[aInd[animal]:nbSteps,,],lnProbs[aInd[animal]:nbSteps,],nrow(lnProbs[aInd[animal]:nbSteps,]),animal,nbAnimals)}
        else{negllk = forward_alg_ind_NH(pi,trMat[aInd[animal]:(aInd[animal+1]-1),,],lnProbs[aInd[animal]:(aInd[animal+1]-1),],nrow(lnProbs[aInd[animal]:(aInd[animal+1]-1),]),animal,nbAnimals)}
        animal=animal+1
      }
      for(j in (1:nbStates)){
        # Estimated log-probability of being  at state j at time t, from Zucchini et al. 
        lProbaState[t,j] = lalpha[t,j] + lbeta[t,j] + negllk
      }
    }
    
    for(j in (1:nbStates)){
      # Estimated log-probability of being  at state j at time nbSteps, from Zucchini et al. 
      lProbaState[nbSteps,j] = lalpha[nbSteps,j] + lbeta[nbSteps,j] + negllk
    }
    
    HatPi = apply(exp(lProbaState),2,mean)
    
    
    Beta = matrix(beta,  nrow=(nbStates*(nbStates-1)),ncol=ncol(X))
    from = 0
    rnames = array(dim=c(nbStates*(nbStates-1)))
    for(pp in 1:(nbStates)){
      from = from+1
      states = seq(1,nbStates)[-from]
      for(idx in 1:(nbStates-1)){
        rnames[(pp-1)*(nbStates-1)+idx] = paste0(from," -> ",states[idx])
      }
    }
    rownames(Beta) = rnames
    # Updated double penalized likelihood
    nllpen[i] <- alpha-Cn*sum(log(HatPi))+nbAnimals*sum(sapply(eta,FUN=function(x)pen(x,lambda,a)))
    
    # Updated joint negative log likelihood
    nllk[i] <- alpha
    
    # stopping rule
    if(i>1 && abs(nllpen[i]-nllpen[i-1])<epsilon){
      convergence <- 0
      break()
    }
  }
  return(list(theta=c(mu=mu,sd=sd),pi=pi,HatStationary=HatPi,beta=Beta,nllk=nllk,nllpen=nllpen,trMat=trMat,K=K,lnProbs=lnProbs,convergence=convergence))
}

gamma.vonMises.EMpen.NH <- function(NbIter,epsilon,Pi0,beta0,lmu0,lsd0,lkappa0,nbSteps,aInd,Obs,lambda,a,Cn,X){
  nbStates <- length(Pi0)
  if(ncol(X) != ncol(beta0)){
    stop("Covariates and coefficient dimensions do not match.")
  }
  if(nrow(beta0)!=(nbStates*(nbStates-1))|| length(lmu0)!= nbStates || length(lsd0)!= nbStates || length(lkappa0)!= nbStates){
    stop("Parameter has the wrong number of states.")
  }
  
  if(length(ncol(X)) < 1 || ncol(X) != ncol(beta0)){
    stop("Vector of covariates not well specificied. For no covariate, a column of 1 has to be specified")
  }
  
  convergence <- 1  # 1 if did not converged, 0 o/w.
  # Initialize parameters
  nbAnimals <- length(aInd)
  nllk <- nllpen <- array()
  pi <- Pi0
  
  lkappa <- (lkappa0)
  kappa <- exp(lkappa0)
  lsd <- lsd0
  lmu <- lmu0
  mu0 <- exp(lmu0)
  sd0 <- exp(lsd0)
  
  # from mean and sd, obtain shape and scale of gamma distribution
  lshape <- log((mu0 * mu0 )/ (sd0 * sd0))
  lscale <- log(sd0 * sd0 / mu0)

  # Vector of log parameters of the gamma distributions shape and scale
  lTheta <- t(matrix(c(lshape,lscale),ncol=2,nrow=nbStates))
  
  # Vector of log parameters of the gamma distributions mean and sd
  lXi <- t(matrix(c(lmu,lsd),ncol=2,nrow=nbStates))
  
  # Vector of mean of the gamma distributions and von Mises: parameters that are constrained by the SCAD.
  lpar <-  c(matrix(cbind(lXi[1,],lkappa),ncol=2,nrow=nbStates))
  # Vector of parameters of the gamma distributions and von Mises.
  par <- exp(lpar)
  
  # Diagonal matrix of log emission distributions for all observations (for all individuals).
  lnProbs <- array(0,dim=c(nbSteps,nbStates))
  lnProbs[1,] <- apply(exp(lTheta),2,FUN=function(x)dgamma(Obs[1,1],shape=x[1],scale=x[2],log=TRUE))
  
  for(t in (2:nbSteps)){
    if(!is.na(Obs[t,1])){lnProbs[t,] <- apply(exp(lTheta),2,FUN=function(x)dgamma(Obs[t,1],shape=x[1],scale=x[2],log=TRUE))}
    if(!is.na(Obs[t,2])){lnProbs[t,] <- lnProbs[t,]+sapply(kappa,FUN=function(x)dvon_mises(Obs[t,2],mu=0,kappa=x,log=TRUE))}
  }
  
  
  beta <- c(beta0)
  # Compute trMat from beta
  pp <- HMM.beta2tpm(nbStates,beta0,X,nbSteps)
  trMat <- pp$tpm
  etaP <- array()
  
  for(i in (1:NbIter)){
    # print(i)
    
    # Update u and v with current parameters
    elnU <- logUNH(pi, trMat,lnProbs,nbSteps, aInd) # NH refers to non-homogeneous, i.e., t.p.m has nbSteps rows.
    elnV <- logVNH(pi, trMat,lnProbs,nbSteps, aInd) # NH refers to non-homogeneous, i.e., t.p.m has nbSteps rows.
    
    # Update beta and initial distribution
    beta_new <- UpdateBeta(beta,X,trMat,nbSteps, aInd,pi,lnProbs,Cn,elnV) 
    pi_new <- DeltaUpdate(pi, trMat,lnProbs,nbSteps, aInd)
    
    # Update standard deviation of gamma distribution
    lsd_new <-  SdUpdate(lXi[1,],lXi[2,],pi, trMat,nbSteps, aInd,elnU)
    
    # Update parameters before updating mu
    lTheta[2,] <-lsd_new+lsd_new-lXi[1,]
    lTheta[1,] <-lXi[1,]+lXi[1,]-(lsd_new+lsd_new)
    lXi[2,] <- lsd_new
    
    trMat <-  HMM.beta2tpm(nbStates,beta_new,X,nbSteps)$tpm
    pi  <- pi_new
    beta <- beta_new
    
    lnProbs <- array(0,dim=c(nbSteps,nbStates))
    lnProbs[1,] <- apply(exp(lTheta),2,FUN=function(x)dgamma(Obs[1,1],shape=x[1],scale=x[2],log=TRUE))
    for(t in (2:nbSteps)){
      if(!is.na(Obs[t,1])){lnProbs[t,] <- apply(exp(lTheta),2,FUN=function(x)dgamma(Obs[t,1],shape=x[1],scale=x[2],log=TRUE))}
      if(!is.na(Obs[t,2])){lnProbs[t,] <- lnProbs[t,]+sapply(exp(lkappa),FUN=function(x)dvon_mises(Obs[t,2],mu=0,kappa=x,log=TRUE))}
    }
    
    elnU <- logUNH(pi, trMat,lnProbs,nbSteps, aInd)
    
    # Sort the current vector of parameter estimates
    t = matrix(par,ncol=2,nrow=nbStates)
    idx <- ordering(nbStates,t)
    t <- t[idx, ]
    for (j in 1:(nrow(t)-1)){
      etaP[j] <- sqrt(sum((t[j+1, ] - t[j, ])^2))
    }
    
    # Update Mu and Kappa
    lpar_new <- MuKappaUpdate(lpar,lXi[2,],etaP,pi,lambda,a,Obs, trMat,nbSteps, aInd,elnU)
    lmu_new <- lpar_new[1:nbStates]
    lkappa_new <- lpar_new[(nbStates+1):(2*nbStates)]
    
    
    # Finish updating parameters
    lTheta[1,] <- (lmu_new+lmu_new)-(lsd_new+lsd_new)
    lTheta[2,] <-lsd_new+lsd_new-lXi[1,]
    lXi[1,] <- lmu_new
    lkappa <- lkappa_new
    lpar <- lpar_new
    par <- exp(lpar)
    
    lnProbs <- array(0,dim=c(nbSteps,nbStates))
    lnProbs[1,] <- apply(exp(lTheta),2,FUN=function(x)dgamma(Obs[1,1],shape=x[1],scale=x[2],log=TRUE))
    for(t in (2:nbSteps)){
      if(!is.na(Obs[t,1])){lnProbs[t,] <- apply(exp(lTheta),2,FUN=function(x)dgamma(Obs[t,1],shape=x[1],scale=x[2],log=TRUE))}
      if(!is.na(Obs[t,2])){lnProbs[t,] <- lnProbs[t,]+sapply(exp(lkappa),FUN=function(x)dvon_mises(Obs[t,2],mu=0,kappa=x,log=TRUE))}
    }
    
    
    # Update mu and sd
    mu <- exp(lmu_new)
    sd <- exp(lsd_new)
    
    # Sort the current vector of parameter estimates
    t = matrix(par,ncol=2,nrow=nbStates)
    idx <- ordering(nbStates,t)
    t <- t[idx, ]
    eta = array()
    for (j in 1:(nrow(t)-1)){
      eta[j] <- sqrt(sum((t[j+1, ] - t[j, ])^2))
    }
    
    # threshold for overlaps is 0.1, it is arbitrary, does not affect the algorithm and can be changed after fitting model.
    idx <- which(eta <= 0.1)
    eta[idx]=0
    # Number of states estimate (hatN in the paper).
    K <- length(which(eta!=0))+1
    eta <- array(eta)
    
    lProbaState = array(dim=c(nbSteps,nbStates)) 
    # Updated forward probabilities
    lalpha = logAlphaNH(pi,trMat,lnProbs,nbSteps,aInd)
    # Updated backward probabilities
    lbeta = logBetaNH(pi,trMat,lnProbs,nbSteps,aInd)
    # Updated negative joint log-likelihood
    alpha <- forward_algNH(pi,trMat,lnProbs,nbSteps,aInd)
    
    #Derive Beta from beta, with s_i->s_j
    Beta = matrix(beta,  nrow=(nbStates*(nbStates-1)),ncol=ncol(X))
    from = 0
    rnames = array(dim=c(nbStates*(nbStates-1)))
    for(pp in 1:(nbStates)){
      from = from+1
      states = seq(1,nbStates)[-from]
      for(idx in 1:(nbStates-1)){
          rnames[(pp-1)*(nbStates-1)+idx] = paste0(from," -> ",states[idx])
      }
    }
    rownames(Beta) = rnames
    # compute lHatPi
    animal=1
    for(t in (1:(nbSteps-1))){
      if(animal<= nbAnimals && t==aInd[animal]){
        # Compute likelihood of individual "animal" (i.e., the one observed at time t).
        if(animal==nbAnimals){negllk = forward_alg_ind(pi,trMat[aInd[animal]:nbSteps,,],lnProbs[aInd[animal]:nbSteps,],nrow(lnProbs[aInd[animal]:nbSteps,]),animal,nbAnimals)}
        else{negllk = forward_alg_ind(pi,trMat[aInd[animal]:(aInd[animal+1]-1),,],lnProbs[aInd[animal]:(aInd[animal+1]-1),],nrow(lnProbs[aInd[animal]:(aInd[animal+1]-1),]),animal,nbAnimals)}
        animal=animal+1
      }
      for(j in (1:nbStates)){
        # Estimated log-probability of being  at state j at time t, from Zucchini et al. 
        lProbaState[t,j] = lalpha[t,j] + lbeta[t,j] + negllk
      }
    }
    
    for(j in (1:nbStates)){
      # Estimated log-probability of being  at state j at time nbSteps, from Zucchini et al. 
      lProbaState[nbSteps,j] = lalpha[nbSteps,j] + lbeta[nbSteps,j] + negllk
    }
    
    HatPi = apply(exp(lProbaState),2,mean)
    
    # Updated double penalized likelihood
    nllpen[i] <- alpha-Cn*sum(log(HatPi))+nbAnimals*sum(sapply(eta,FUN=function(x)pen(x,lambda,a)))
    # Updated joint negative log likelihood
    nllk[i] <- alpha
  
    # stopping rule
    if(i>1 && abs(nllpen[i]-nllpen[i-1])<epsilon){
      convergence <- 0
      break()
    }
  }
  
  return(list(theta=c(mu=mu,sd=sd,Kappa=exp(lkappa)),pi=pi,HatStationary=HatPi,beta=Beta,nllk=nllk,nllpen=nllpen,trMat=trMat,K=K,lnProbs=lnProbs,convergence=convergence))
}

gamma.vonMIses.EMpenMulti <- function(NbIter,epsilon,Pi0,trMat0,lmu0,lsd0,lkappa0,nbSteps,aInd,Obs,lambda,a,Cn){
  
  convergence = 1 # 1 if did not converged, 0 o/w.
  
  # Initialize parameters
  nllk <- nllpen <- array()
  pi <- Pi0
  nbAnimals <- length(aInd)
  
  trMat <- trMat0
  nbStates <- ncol(trMat)
  lkappa <- lkappa0
  
  mu0 <- exp(lmu0)
  sd0 <- exp(lsd0)
  
  kappa = exp(lkappa0)
  lsd  <- lsd0
  lmu <- lmu0
  # from mean and sd, obtain shape and scale of gamma distribution
  lshape <- log((mu0 * mu0 )/ (sd0 * sd0))
  lscale <- log(sd0 * sd0 / mu0)
  # Vector of log-parameters of gamma distribution: shape and scale.
  lTheta <- t(matrix(c(lshape,lscale),ncol=2,nrow=nbStates))
  # Vector of log-parameters of gamma distribution: mean and standard deviation
  lXi <- t(matrix(c(lmu0,lsd0),ncol=2,nrow=nbStates))
  # Vector of log-parameters that are penalized by SCAD.
  lpar <-  c(matrix(cbind(lXi[1,],lkappa),ncol=2,nrow=nbStates))
  # Vector of parameters that are penalized by SCAD
  par <- exp(lpar)
  
  theta <- exp(lTheta)
  # Diagonal matrix of log emission distributions for all observations (for all individuals).
  lnProbs <- array(0,dim=c(nbSteps,nbStates))
  lnProbs[1,] <- apply(theta,2,FUN=function(x)dgamma(Obs[1,1],shape=x[1],scale=x[2],log=TRUE))
  for(t in (2:nbSteps)){
    if(!is.na(Obs[t,1])){lnProbs[t,] <- apply(theta,2,FUN=function(x)dgamma(Obs[t,1],shape=x[1],scale=x[2],log=TRUE))}
    if(!is.na(Obs[t,2])){lnProbs[t,] <- lnProbs[t,]+sapply(kappa,FUN=function(x)dvon_mises(Obs[t,2],mu=0,kappa=x,log=TRUE))}
  }
  
  pi_new <- pi
  etaP <- eta<- array()
  
  for(i in (1:NbIter)){
    # print(i)
    # Obtain trMat from working parameters
    ttrMat <- pn2pw(trMat,aInd)
    
    # Update u and v with current parameters
    elnU <- logU(pi, trMat,lnProbs,nbSteps, aInd)
    elnV <- logV(pi, trMat,lnProbs,nbSteps, aInd)
   
    # Update working parameter for t.p.m
    ttrMat_new <-  trMatUpdate(ttrMat,nbSteps, aInd,elnU,elnV,nbStates,Cn)
    
    # Update t.p.m
    trMat_new <- pw2pn(ttrMat_new,nbStates)
    # Update stationary probability
    pi_new<- base::solve(t(diag(nbStates)-trMat_new+1),rep(1,nbStates))
    # Update standard deviation of gamma distribution
    lsd_new <- SdUpdate(lXi[1,],lsd,pi, trMat,nbSteps, aInd,elnU)
    
    # Update parameters before updating mu
    lTheta[2,] <-lsd_new+lsd_new-lXi[1,]
    lTheta[1,] <- (lXi[1,]+lXi[1,])-(lsd_new+lsd_new)
    lXi[2,] <- lsd_new
    
    trMat <- trMat_new
    pi  <- pi_new
    
    #Update diagonal matrix of log emission distributions for all observations (for all individuals).
    lnProbs <- array(0,dim=c(nbSteps,nbStates))
    lnProbs[1,] <- apply(exp(lTheta),2,FUN=function(x)dgamma(Obs[1,1],shape=x[1],scale=x[2],log=TRUE))
    for(t in (2:nbSteps)){
      if(!is.na(Obs[t,1])){lnProbs[t,] <- apply(exp(lTheta),2,FUN=function(x)dgamma(Obs[t,1],shape=x[1],scale=x[2],log=TRUE))}
      if(!is.na(Obs[t,2])){lnProbs[t,] <- lnProbs[t,]+sapply(exp(lkappa),FUN=function(x)dvon_mises(Obs[t,2],mu=0,kappa=x,log=TRUE))}
    }
    
    elnU <- logU(pi, trMat,lnProbs,nbSteps, aInd)
    
    
    t = matrix(par,ncol=2,nrow=nbStates)
    idx <- ordering(nbStates,t)
    t <- t[idx, ]
    for (j in 1:(nrow(t)-1)){
      etaP[j] <- sqrt(sum((t[j+1, ] - t[j, ])^2))
    }
    
    lpar_new <- MuKappaUpdate(lpar,lXi[2,],etaP,pi,lambda,a,Obs, trMat,nbSteps, aInd,elnU)
    
    lmu_new <- lpar_new[1:nbStates]
    lkappa_new <- lpar_new[(nbStates+1):(2*nbStates)]
    
    
    #Finish updating parameters
    lTheta[1,] <- (lmu_new+lmu_new)-(lsd_new+lsd_new)
    lTheta[2,] <-lsd_new+lsd_new-lmu_new
    lXi[1,] <- lmu_new
    lkappa <- lkappa_new
    lpar <- lpar_new
    par <- exp(lpar)
    
    lnProbs <- array(0,dim=c(nbSteps,nbStates))
    lnProbs[1,] <- apply(exp(lTheta),2,FUN=function(x)dgamma(Obs[1,1],shape=x[1],scale=x[2],log=TRUE))
    for(t in (2:nbSteps)){
      if(!is.na(Obs[t,1])){lnProbs[t,] <- apply(exp(lTheta),2,FUN=function(x)dgamma(Obs[t,1],shape=x[1],scale=x[2],log=TRUE))}
      if(!is.na(Obs[t,2])){lnProbs[t,] <- lnProbs[t,]+sapply(exp(lkappa),FUN=function(x)dvon_mises(Obs[t,2],mu=0,kappa=x,log=TRUE))}
    }
    
    
    # Update mu and sd
    mu <- exp(lmu_new)
    sd <- exp(lsd_new)
    
    # Sort the current vector of parameter estimates
    t = matrix(par,ncol=2,nrow=nbStates)
    idx <- ordering(nbStates,t)
    t <- t[idx, ]
    for (j in 1:(nrow(t)-1)){
      eta[j] <- sqrt(sum((t[j+1, ] - t[j, ])^2))
    }
    
    # threshold for overlaps is 0.1, it is arbitrary, does not affect the algorithm and can be changed after fitting model.
    idx <- which(eta <= 0.1)
    eta[idx]=0
    # Estimate of the number of states.
    K <- length(which(eta!=0))+1
    eta <- array(eta)
    
    # Updated joint negative log likelihood
    nllk[i] <- forward_alg(pi,trMat,lnProbs,nbSteps,aInd)
    # Updated double penalized likelihood
    nllpen[i] <- forward_alg(pi,trMat,lnProbs,nbSteps,aInd)-Cn*sum(log(pi))+nbAnimals*sum(sapply(eta,FUN=function(x)pen(x,lambda,a)))

    #Stopping rule
    if(i>1 && abs(nllpen[i]-nllpen[i-1])<epsilon){
      convergence = 0
      break()
    }
  }
  return(list(theta=c(mu=mu,sd=sd,Kappa=exp(lkappa)),pi=pi,trMat=trMat,nllk=nllk,nllpen=nllpen,K=K,lnProbs=lnProbs,lambda=lambda,Cn=Cn,convergence = convergence))
}


DPMLE <- function(NbIter,epsilon,Pi0,lmu0,lsd0,nbSteps,aInd,Obs,lambda,a,Cn,type, lkappa0=NULL, X = NULL,trMat0 = NULL, beta0 = NULL,...)
{
 if(is.null(lambda) || is.null(Cn)){
   stop("Specify hyperpamarameters.")
 }
 if (length(type) != 2){
   stop("Wrong specification of method")
 } 
  if (is.null(trMat0) == is.null(beta0)){
    stop("One type of  t.p.m has to be formulated: trMat0 or beta0")
  } 
  if(length(lmu0) != length(lsd0)){
    stop("The mean and sd vectors need to be the same length")
  }
  if(type[[1]]=="multi"){
    if(length(lmu0) != length(lkappa0)){
      stop("The mean and concentration vectors need to be the same length")
    }
    if(type[[2]]=="NS"){
      return(gamma.vonMises.EMpen.NH(NbIter,epsilon,Pi0,beta0,lmu0,lsd0,lkappa0,nbSteps,aInd,Obs,lambda,a,Cn,X))
    }else{
      return(gamma.vonMIses.EMpenMulti(NbIter,epsilon,Pi0,trMat0,lmu0,lsd0,lkappa0,nbSteps,aInd,Obs,lambda,a,Cn))
    }
  }else{
    if(type[[2]]=="NS"){
      return(gamma.EMpen.NH(NbIter,epsilon,Pi0,beta0,lmu0,lsd0,nbSteps,aInd,Obs,lambda,a,Cn,X))
    }else{
      return(gamma.EMpen(NbIter,epsilon,Pi0,trMat0,lmu0,lsd0,nbSteps,aInd,Obs,lambda,a,Cn))
    }
  }
}




gamma.forecast_getcumdistr <- function(xt,tt,lalpha,shape,scale,trMat){
  # to prevent from overflow
  lafact <- apply (lalpha ,1, max)
  # Probability that observation at time tt is lower or equal than xt
  Px = apply(cbind(shape,scale),1,FUN=function(x)pgamma(xt,shape=x[1],scale=x[2]))
  # Conditional probability that observation at time t is lower or equal than xt, given all the previous observations, based on Zucchini et al.
  Proba = sum(exp(lalpha[tt-1,]-lafact[tt-1]) %*% trMat%*%diag(Px))/sum(exp(lalpha[tt-1,]-lafact[tt-1]))
  return(Proba)
}


gamma.conditional_distribution <- function(xt,tt,lalpha,lbeta,shape,scale,trMat){
  # to prevent from overflow
  lafact <- apply (lalpha ,1, max)
  lbfact <- apply (lbeta ,1, max)
  # Probability that observation at time tt is lower or equal than xt
  Px = apply(cbind(shape,scale),1,FUN=function(x)pgamma(xt,shape=x[1],scale=x[2])) 
  # Probability that observation at time tt is lower or equal than xt given all the other observations other than at time t.
  foo <- (exp (lalpha[tt-1,]- lafact[tt]) %*% trMat) * exp (lbeta[tt,]- lbfact [tt])
  foo <- foo / sum (foo)
  Proba = sum(foo %*% diag(Px))
  return(Proba)
}

vmises.forecast_getcumdistr <- function(xt,tt,lalpha,kappa,trMat){
  # to prevent from overflow
  lafact <- apply (lalpha ,1, max)
  if(!is.na(xt)){ #pvon_mises provides warnings if xt == NA.
  # Probability that observation at time tt is lower or equal than xt
  Px = sapply(kappa,FUN=function(x)pvon_mises(xt,mu=0,kappa=x)) 
  # Conditional probability that observation at time t is lower or equal than xt, given all the previous observations, based on Zucchini et al.
  Proba = sum(exp(lalpha[tt-1,]-lafact[tt-1]) %*% trMat%*%diag(Px))/sum(exp(lalpha[tt-1,]-lafact[tt-1]))
  }else{
    Proba = 0
  }
  return(Proba)
}


vmises.conditional_distribution <- function(xt,tt,lalpha,lbeta,kappa,trMat){
  # to prevent from overflow
  lafact <- apply (lalpha ,1, max)
  lbfact <- apply (lbeta ,1, max)
  if(!is.na(xt)){#pvon_mises provides warnings if xt == NA.
  # Probability that observation at time tt is lower or equal than xt
  Px = sapply(kappa,FUN=function(x)pvon_mises(xt,mu=0,kappa=x)) 
  # Conditional probability that observation at time t is lower or equal than xt, given all the previous observations, based on Zucchini et al.
  foo <- (exp (lalpha[tt-1,]- lafact[tt]) %*% trMat) * exp (lbeta[tt,]- lbfact [tt])
  foo <- foo / sum (foo)
  Proba = sum(foo %*% diag(Px))
  }else{
  Proba = 0
  }
  return(Proba)
}


# Same as without "_NH" except that the trMat now depends on time.
gamma.forecast_getcumdistr_NH <- function(xt,tt,lalpha,shape,scale,trMat){
  lafact <- apply (lalpha ,1, max)
  Px = apply(cbind(shape,scale),1,FUN=function(x)pgamma(xt,shape=x[1],scale=x[2]))
  Proba = sum(exp(lalpha[tt-1,]-lafact[tt-1]) %*% trMat[tt,,]%*%diag(Px))/sum(exp(lalpha[tt-1,]-lafact[tt-1]))
  return(Proba)
}

# Same as "gamma.conditional_distribution" except that the trMat now depends on time.
gamma.conditional_distribution_NH <- function(xt,tt,lalpha,lbeta,shape,scale,trMat){
  lafact <- apply (lalpha ,1, max)
  lbfact <- apply (lbeta ,1, max)
  Px = apply(cbind(shape,scale),1,FUN=function(x)pgamma(xt,shape=x[1],scale=x[2])) 
  foo <- (exp (lalpha[tt-1,]- lafact [tt]) %*% trMat[tt,,]) * exp (lbeta[tt,]- lbfact [tt])
  foo <- foo / sum (foo)
  proba = sum(foo %*% diag(Px))
  return(proba)
}

# Same as "vmises.forecast_getcumdistr" except that the trMat now depends on time.
vmises.forecast_getcumdistr_NH <- function(xt,tt,lalpha,kappa,trMat){
  lafact <- apply (lalpha ,1, max)
  if(!is.na(xt)){
    Px = sapply(kappa,FUN=function(x)pvon_mises(xt,mu=0,kappa=x)) 
    Proba = sum(exp(lalpha[tt-1,]-lafact[tt-1]) %*% trMat[tt,,]%*%diag(Px))/sum(exp(lalpha[tt-1,]-lafact[tt-1]))
    
  }else{
    Proba = 0
  }
  return(Proba)
}

# Same as "vmises.conditional_distribution" except that the trMat now depends on time.
vmises.conditional_distribution_NH <- function(xt,tt,lalpha,lbeta,kappa,trMat){
  lafact <- apply (lalpha ,1, max)
  lbfact <- apply (lbeta ,1, max)
  if(!is.na(xt)){
    
    Px = sapply(kappa,FUN=function(x)pvon_mises(xt,mu=0,kappa=x)) 
    foo <- (exp (lalpha[tt-1,]- lafact [tt]) %*% trMat[tt,,]) * exp (lbeta[tt,]- lbfact [tt])
    foo <- foo / sum (foo)
    Proba = sum(foo %*% diag(Px))
  }else{
    Proba = 0
  }
  return(Proba)
}


pseudo_residuals <- function(Obs,type=NULL,dist=NULL,lalpha,lbeta=NULL,delta,shape=NULL,scale=NULL,trMat,kappa=NULL){
n = length(Obs)
  if(length(type[[1]])==0){
  type ="ordinary"
}
if(type[[2]]=="H"){  
if(dist=="gamma"){
  if(type[[1]] == "forecast"){
    dists <- sapply(1:length(Obs),FUN=function(x)gamma.forecast_getcumdistr.step(Obs[x],x,lalpha,shape,scale,trMat))
  }else{
    dists <- sapply(1:length(Obs),FUN=function(x)gamma.conditional_distribution(Obs[x],x,lalpha,lbeta,shape,scale,trMat))
  }
  npsr <- rep(NA, n)
  npsr <- qnorm(dists)
  P = apply(cbind(shape,scale),1,FUN=function(x)pgamma(Obs[1],shape=x[1],scale=x[2]))
  npsr<-c(qnorm( delta%*% P),npsr)
}else{
  if(type[[1]] == "forecast"){
    dists <- sapply(2:length(Obs),FUN=function(x)vmises.forecast_getcumdistr.step(Obs[x],x,lalpha,kappa,trMat))
  }else{
    dists <- sapply(2:length(Obs),FUN=function(x)vmises.conditional_distribution(Obs[x],x,lalpha,lbeta,kappa,trMat))
  }
  npsr <- rep(NA, n)
  npsr <- qnorm(dists)
}
}else{
  
  
  if(dist=="gamma"){
    if(type[[1]] == "forecast"){
      dists <- sapply(1:length(Obs),FUN=function(x)gamma.forecast_getcumdistr_NH(Obs[x],x,lalpha,shape,scale,trMat))
    }else{
      dists <- sapply(1:length(Obs),FUN=function(x)gamma.conditional_distribution_NH(Obs[x],x,lalpha,lbeta,shape,scale,trMat))
    }
    npsr <- rep(NA, n)
    npsr <- qnorm(dists)
    P = apply(cbind(shape,scale),1,FUN=function(x)pgamma(Obs[1],shape=x[1],scale=x[2]))
    npsr<-c(qnorm( delta%*% P),npsr)
  }else{
    if(type[[1]] == "forecast"){
      dists <- sapply(2:length(Obs),FUN=function(x)vmises.forecast_getcumdistr_NH(Obs[x],x,lalpha,kappa,trMat))
    }else{
      dists <- sapply(2:length(Obs),FUN=function(x)vmises.conditional_distribution_NH(Obs[x],x,lalpha,lbeta,kappa,trMat))
    }
    npsr <- rep(NA, n)
    npsr <- qnorm(dists)
  }
  
}
  return (npsr)
}




























