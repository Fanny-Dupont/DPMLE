#################################################################
###                Simulation Setting                         ###
 #               three simulated states                        #
 #      code to fit stationary and non-stationary HMMs         #
###                 with momentuHMM                           ###
#################################################################

require(momentuHMM)
require(boot)
require(tidyverse)
require(pracma)
require(dirmult)

nbStates <- N <- 3             # number of states
M <- 1                         # number of individuals
nobs <- n <- 5000              # number of observation per individuals
obsVect <- rep(list(nobs),M)   # number of observations per animal (length of obsVect must be a factor of M)
T_m <- rep(obsVect,M/length(obsVect))


# observation distribution parameters, only simulated gamma
dist <- list(step = "gamma")

#set parameters, same a in the Pohle et al. (2017).
mu0 <- c(5.5,3,1)
shape0<- c(12,4,1.5)
sd0<- mu0/sqrt(shape0)
beta0<- matrix(c(rep(log(0.1/(1-0.2)),6)),nrow=1,ncol=6)
stateNames <- c("travelling","foraging","resting")

#Scenario 1: Benchmark
#Simulate the data
id = 1 # for simulation id 1.
set.seed(id)  # id of the job
exData <- simData(nbAnimals = M, nbStates,dist=dist,Par=list(step = c(mu0,sd0)),beta=beta0,delta=rep(1/N, N),obsPerAnimal=obsVect)
exData <- exData[-which(is.na(exData$step)),]

#Scenario 2: Outliers
# Adding uniform [10,20] to 0.5% of the step-length of the dataset?
#   data <- exData 
#   idx <- seq(1,(nobs-1))
#   idx1 <- sample(idx,ceiling(0.005*(nobs)),replace=FALSE,prob=rep(1/(nobs-1),nobs-1))
#   data$step[idx1] <-  data$step[idx1]+runif(ceiling(0.005*nobs),10,20)
#   exData <- data


#Scenario 3: Heterogeneity in t.p.m
# M = 10
# nobs <- n <- 500               # number of observation per individuals
# obsVect <- rep(list(nobs),M)   # number of observations per animal (length of obsVect must be a factor of M)
# T_m <- rep(obsVect,M/length(obsVect))
# pie <- matrix(c(0.5,0.5),nrow=1,ncol=2)         # mixture probabilities 
# # mixture probabilities 
# beta0 <- list(beta=rbind(matrix(c(rep(log(0.1/(1-0.2)),2),rep(log(0.1/(1-0.2)),2),log(0.8/0.1),log(0.1/0.1)),nrow=1,ncol=6),matrix(c(rep(log(0.1/(1-0.2)),6)),nrow=1,ncol=6)), pi=pie) 
# # 1 -> 2 et 2 <- 1 for mix 1 is row1
# delta0 <- matrix(c(0.00001,0.9991,0.00001,0.0001,0.999,0.00001),nrow=2,ncol=3)#1 and 2 for mix 1 is first row
# exData <- momentuHMM::simData(nbAnimals=M ,nbStates=nbStates,dist=dist,Par=list(step = c(mu0,sd0)),delta=delta0,beta=beta0,obsPerAnimal=obsVect,mixtures=2,states=TRUE)
# exData <- exData[-which(is.na(exData$step)),]

#Scenario 4: Heterogeneity in observation process
# mean0s1 <- rshifted_lnorm(M, meanlog = log(5.5), sdlog =  0.15, shift = 0)
#
# ##generate data
# ##set parameters 
# DataInd <- list()
# #parameters shared by every individuals
# for(i in (1:M)){#simulate each individual according to their mean for state 2, and then grouping them
#   #in one dataset
#   mu <- c(mean0s1[i],3,1)
#   sigma <- c(mean0s1[i]/sqrt(12), sd0[2],sd0[3]) 
#   DataInd[[i]] <- simData(nbAnimals = 1, nbStates,dist=dist,Par=list(step = c(mu,sigma)),beta=beta0,obsPerAnimal=obsVect[[i]])
#   levels(DataInd[[i]]$ID) <-  as.numeric(factor(DataInd[[i]]$ID)) +i-1
#   if(i==1){exData <- DataInd[[i]]}else{
#     exData <- rbind(exData,DataInd[[i]])
#   }
# }
# exData <- exData[-which(is.na(exData$step)),]

#Scenario 5: Available upon request, custom coded.
#time-varying mean parameter of the state-dependent gamma distributions for slow movement 
#generated using autoregressive processes of order 1, persistence coefficient 0.85 and noise variance de 0.25.

#Scenario 6: Temporal variation in hidden process
# M = 1
# nobs <- n <- 5000               # number of observation per individuals
# obsVect <- rep(list(nobs),M)   # number of observations per animal (length of obsVect must be a factor of M)
# T_m <- rep(obsVect,M/length(obsVect))
# beta <- matrix(c(logit(0.1),0,0,logit(0.1),0,0,logit(0.2),0,0,logit(0.3),logit(0.2),0,logit(0.2),0,0,logit(0.187),-logit(0.2),0),nrow=3,ncol=6)#intercept is: 1 -> 2, 1 -> 3, 2 -> 1, 2 -> 3, 3 -> 1, 3 -> 2
# 
# covs<-data.frame(hour96=0:95)
# formula= ~cosinor(hour96,96)
# rownames(beta)<-c("(Intercept)",
#                   "cosinorCos(hour96, 96)",
#                   "cosinorSin(hour96, 96)")
# 
# 
# exData <- simData(nbAnimals = M, nbStates,dist=dist,Par=list(step = c(mu0,sd0)),beta=beta,formula=formula,covs=covs,obsPerAnimal=obsVect)
# exData <- exData[-which(is.na(exData$step)),]



Nrep <- 150                 # number of random initial values explored.
data <- exData

# Fit standard stationary HMM with 2 states

nbStates <- 2 # number of states
stateNames <- c("directed","wandering")
dist <- list(step = "gamma")
Par0 <- list()
# Set vectors 
inter <- list()
loglike <- array(dim = c(Nrep,M))

for (i in (1:Nrep)){
  skip_to_next <- FALSE
  #Initialization of parameters, with two states, from Pohle et al. (2017)
  mu0<-c(runif(2,min=min(exData$step),max=max(exData$step)))
  shape0<-c(rgamma(2,shape=1,scale=2.5))
  sigma0<- mu0/sqrt(shape0)
  beta<-NULL
  j <- 2
  beta<-array(NA,dim=c(j,j))
  for(l in 1:j){
    alpha<-rep(0.2/(j-1),j)
    alpha[l]<-0.8
    alpha<-alpha*10
    beta[l,]<-rdirichlet(1,alpha)
  }
  beta <- logit(beta)
  beta <- t(beta)
  beta <- matrix(beta[col(beta) != row(beta)],ncol=2,nrow=1)
  beta <- rbind(beta)
  
  Par0 <- list(step=c(mu0,sigma0))
  
  tryCatch(a <- fitHMM(data, nbState = nbStates, 
                       dist = dist, Par0 = Par0, beta0=beta,stateNames = stateNames,stationary=TRUE), error = function(e) { skip_to_next <- TRUE}) # Trycatch function required for code submitted in Compute Canada
  
  if(skip_to_next){ next }
  inter[[i]] <- a
  loglike[i] <- inter[[i]]$mod$minimum
}
idx <- which.min(loglike)
null2 <- inter[[idx]]


# Fit standard stationary HMM with 3 states
nbStates <- 3 # number of states
stateNames <- c("directed","wandering","resting")
# re-initialize vectors 
inter <- list()
loglike <- array(dim = c(Nrep))
Par0 <- list()

for (i in (1:Nrep)){
  skip_to_next <- FALSE
  #Initialization of parameters, with three states, from Pohle et al. (2017)
  mu0<-c(runif(3,min=min(exData$step),max=max(exData$step)))
  shape0<-c(rgamma(3,shape=1,scale=2.5))
  sigma0<- mu0/sqrt(shape0)
  beta<-NULL
  j <- 3
  beta<-array(NA,dim=c(j,j))
  for(l in 1:j){
    alpha<-rep(0.2/(j-1),j)
    alpha[l]<-0.8
    alpha<-alpha*10
    beta[l,]<-rdirichlet(1,alpha)
  }
  beta <- logit(beta)
  beta <- t(beta)
  beta <- matrix(beta[col(beta) != row(beta)],ncol=6,nrow=1)
  beta <- rbind(beta)
  Par0 <- list(step=c(mu0,sigma0))
  
  tryCatch(a <- fitHMM(data, nbState = nbStates, 
                       dist = dist, Par0 = Par0, beta0=beta,stateNames = stateNames,stationary=TRUE), error = function(e) { skip_to_next <- TRUE})
  
  if(skip_to_next){ next }
  inter[[i]] <- a
  loglike[i] <- inter[[i]]$mod$minimum
}
idx <- which.min(loglike)
null3 <- inter[[idx]]


# Fit standard HMM with 4 states
nbStates <- 4 # number of states
stateNames <- c("directed","wandering","resting","other")

# Setting up the starting values
# Re-initialize set vectors 
inter <- list()
loglike <- array(dim = c(Nrep))
Par0 <- list()
for (i in (1:Nrep)){
  skip_to_next <- FALSE
  #Initialization of parameters, with four states, from Pohle et al. (2017)
  mu0<-c(runif(4,min=min(exData$step),max=max(exData$step)))
  shape0<-c(rgamma(4,shape=1,scale=2.5))
  sigma0<- mu0/sqrt(shape0)
  beta<-NULL
  j <- 4
  beta<-array(NA,dim=c(j,j))
  for(l in 1:j){
    alpha<-rep(0.2/(j-1),j)
    alpha[l]<-0.8
    alpha<-alpha*10
    beta[l,]<-rdirichlet(1,alpha)
  }
  beta <- logit(beta)
  beta <- t(beta)
  beta <- matrix(beta[col(beta) != row(beta)],ncol=12,nrow=1)
  beta <- rbind(beta)
  Par0 <- list(step=c(mu0,sigma0))
  
  tryCatch(a <- fitHMM(data, nbState = nbStates, 
                       dist = dist, Par0 = Par0, beta0=beta,stateNames = stateNames,stationary=TRUE), error = function(e) { skip_to_next <- TRUE})
  
  if(skip_to_next){ next }
  inter[[i]] <- a
  loglike[i] <- inter[[i]]$mod$minimum
}
idx <- which.min(loglike)
null4 <- inter[[idx]]


#Save Results.
saveRDS(null2, file=paste("1e_scenario_3_5000_s2_",id, '.RData', sep = ""))
saveRDS(null3, file=paste("1e_scenario_3_5000_s3_",id, '.RData', sep = ""))
saveRDS(null4, file=paste("1e_scenario_3_5000_s4_",id, '.RData', sep = ""))




# Code to fit non-stationary HMM with time of day a linear covariate (data every 15 minutes)

original_matrix <- matrix(c(rep(1:96, length.out = (n-1))),ncol=1)
# Repeat the matrix M times using rbind
tod <- do.call(rbind, replicate(M, original_matrix, simplify = FALSE))

# fit model
formula <- ~ tod
data$tod = tod


#Fit standard non-stationary HMM with tod as covariate

# Fit standard HMM with 2 states
nbStates <- 2 # number of states
stateNames <- c("directed","wandering")
dist <- list(step = "gamma")
Par0 <- list()
# set vectors 
inter <- list()
loglike <- array(dim = c(Nrep))

for (i in (1:Nrep)){
  skip_to_next <- FALSE
  mu0<-c(runif(2,min=min(exData$step),max=max(exData$step)))
  shape0<-c(rgamma(2,shape=1,scale=2.5))
  sigma0<- mu0/sqrt(shape0)
  beta<-NULL
  j <- 2
  beta<-array(NA,dim=c(j,j))
  for(l in 1:j){
    alpha<-rep(0.2/(j-1),j)
    alpha[l]<-0.8
    alpha<-alpha*10
    beta[l,]<-rdirichlet(1,alpha)
  }
  beta <- logit(beta)
  beta <- t(beta)
  beta <- matrix(beta[col(beta) != row(beta)],ncol=2,nrow=1)
  beta <- rbind(beta,c(0,0))
  Par0 <- list(step=c(mu0,sigma0))
  
  tryCatch(a <- fitHMM(data, nbState = nbStates, 
                       dist = dist, Par0 = Par0, beta0=beta,formula=formula,stateNames = stateNames,stationary=FALSE), error = function(e) { skip_to_next <- TRUE})
  
  if(skip_to_next){ next }
  inter[[i]] <- a
  loglike[i] <- inter[[i]]$mod$minimum
}
idx <- which.min(loglike)
null2NS <- inter[[idx]]


# Fit standard HMM with 3 states
nbStates <- 3           #number of states
stateNames <- c("directed","wandering","resting")
nbAnimals <- length(unique(exData$ID))

# set vectors 
minsim <- array(dim = c(1,Nrep))
inter <- list()
loglike <- array(dim = c(Nrep))
Par0 <- list()
for (i in (1:Nrep)){
  skip_to_next <- FALSE
  mu0<-c(runif(3,min=min(exData$step),max=max(exData$step)))
  shape0<-c(rgamma(3,shape=1,scale=2.5))
  sigma0<- mu0/sqrt(shape0)
  beta<-NULL
  j <- 3
  beta<-array(NA,dim=c(j,j))
  for(l in 1:j){
    alpha<-rep(0.2/(j-1),j)
    alpha[l]<-0.8
    alpha<-alpha*10
    beta[l,]<-rdirichlet(1,alpha)
  }
  beta <- logit(beta)
  beta <- t(beta)
  beta <- matrix(beta[col(beta) != row(beta)],ncol=6,nrow=1)
  beta <- rbind(beta,c(0,0,0,0,0,0))
  Par0 <- list(step=c(mu0,sigma0))
  
  tryCatch(a <- fitHMM(data, nbState = nbStates, 
                       dist = dist, Par0 = Par0, beta0=beta,formula=formula,stateNames = stateNames,stationary=FALSE), error = function(e) { skip_to_next <- TRUE})
  
  if(skip_to_next){ next }
  inter[[i]] <- a
  loglike[i] <- inter[[i]]$mod$minimum
}
idx <- which.min(loglike)
null3NS <- inter[[idx]]

# Fit standard HMM with 4 states
nbStates <- 4 # number of states
stateNames <- c("directed","wandering","resting","other")

# Setting up the starting values
# set vectors 
minsim <- array(dim = c(1,Nrep))
inter <- list()
loglike <- array(dim = c(Nrep))
Par0 <- list()
for (i in (1:Nrep)){
  skip_to_next <- FALSE
  mu0<-c(runif(4,min=min(exData$step),max=max(exData$step)))
  shape0<-c(rgamma(4,shape=1,scale=2.5))
  sigma0<- mu0/sqrt(shape0)
  beta<-NULL
  j <- 4
  beta<-array(NA,dim=c(j,j))
  for(l in 1:j){
    alpha<-rep(0.2/(j-1),j)
    alpha[l]<-0.8
    alpha<-alpha*10
    beta[l,]<-rdirichlet(1,alpha)
  }
  beta <- logit(beta)
  beta <- t(beta)
  beta <- matrix(beta[col(beta) != row(beta)],ncol=12,nrow=1)
  beta <- rbind(beta,rep(0,12))
  Par0 <- list(step=c(mu0,sigma0))
  
  tryCatch(a <- fitHMM(data, nbState = nbStates, 
                       dist = dist, Par0 = Par0, beta0=beta,formula=formula,stateNames = stateNames,stationary=FALSE), error = function(e) { skip_to_next <- TRUE})
  
  if(skip_to_next){ next }
  inter[[i]] <- a
  loglike[i] <- inter[[i]]$mod$minimum
}
idx <- which.min(loglike)
null4NS <- inter[[idx]]

saveRDS(null2NS, file=paste("1e_scenario_3_5000_cov_s2_",id, '.RData', sep = ""))
saveRDS(null3NS, file=paste("1e_scenario_3_5000_cov_s3_",id, '.RData', sep = ""))
saveRDS(null4NS, file=paste("1e_scenario_3_5000_cov_s4_",id, '.RData', sep = ""))


















