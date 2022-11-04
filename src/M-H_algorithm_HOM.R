############################################################################
# Function for Metropolis-Hasting under HOM model
############################################################################
#' @title Simulation of MCMC samples using Metropolis-Hasting under HOM model
#'
#' @description This function simulates MCMC samples using Metropolis-Hasting algorithm.
#' @param theta.k Initial values of parameters
#' @param data data sets contain the original data
#' @param sigma.error measurement error variation
#' @param mu_f mean of X given W 
#' @param var_f variance of X given W 
#' @param alpha0 percentage of missing values
#' @param alpha missing data mechanism parameters
#' @param M Number of MC size
#' @return M simulated values for U, X and missing Y
#' @example 



# Metropolis-Hasting Algorithm
metropolis.HOM <- function (n, q, data, theta.k, sigma.error, mu_f, var_f, alpha0, alpha, M, burnin) 
{ 
  Ustar<-matrix(0,nrow= M,ncol=n) # Mxn matrix to be filled with simulated u (random effects).
  xstar<-matrix(0,nrow= M,ncol=n*q) # Mxn matrix to be filled with simulated x (measurement error covariate).
  ystar<-matrix(NA,nrow= M,ncol=n*q) # Mxn matrix to be filled with simulated y
  ap.x<-rep(0,n*q) # to store the acceptance probability for candidate for x.
  ap.u<-rep(0,n) # to store the acceptance probability for candidate for u.
  ap.y<-rep(0,n*q) # to store the acceptance probability for candidate for y.
  ones <- matrix(rep(1,q),q,1) # qx1 vector pf 1's.
  z<-matrix(0,nrow=n*q,ncol=n) # nqxn, Design matrix
  for(j in 1:n)
  {
    z[(1+(j-1)*q):(q*j),j]=1
  }
  
  #generate initial samples
  Ustar[1,]<-rnorm(n,0, sqrt(theta.k[5])) # initial simulation of random effects.
  ## Generate initial error covariate x
  xstar[1,] <- rnorm(n*q, 1, sqrt(2)) # X ~ N(mu_f,sqrt(var_f)) ## initial simulation of X
  # Generate initial sample of missing responses
  y.obs <- data$Y.obs
  id.na <- which(is.na(y.obs))
  fixed.y <- theta.k[1]+theta.k[2]*xstar[1,][id.na]+theta.k[3]*data$Zij[id.na]+theta.k[4]*data$Zistar[id.na]
  eta.y <- fixed.y+(z%*%Ustar[1,])[id.na]
  pij.y <- exp(eta.y)/(1+exp(eta.y)) # probability of Bernoulli trail.
  y.obs[id.na] <- rbinom(length(id.na),1,pij.y)
  ystar[1,] <- y.obs # nqx1; copies of the response.
  
  ## Metropolis steps
  for (m in 2:M)
  { 
    ## Generate candidate samples for y, u and x.
    u.candidate<- rnorm(n,0, sqrt(theta.k[5])) 
    ##x.candidate<- matrix(rnorm(n*q, mu_f, sqrt(var_f)),n,q)
    x.candidate<-matrix(rnorm(n*q,1,sqrt(2)),n,q)
    y.obs[id.na] <- rbinom(length(id.na),1,pij.y)
    y.candidate <- y.obs
    
    ## Calculate Denominator of Acceptance criterion vectors
    xmat <-matrix(xstar[m-1,],n,q) # nxq; copies of previously simulated x values.
    umat <-matrix(rep(Ustar[m-1,],1),n,1)# nx1; copies of previously simulated random effects.
    ymat <-matrix(ystar[m-1,],n,q)
    
    fixed.denom<- theta.k[1]+theta.k[2]*xmat+theta.k[3]*data$Zij+theta.k[4]*data$Zistar # fixed effects
    Fixed.denom<- matrix(rep(fixed.denom,1),n*q,1,byrow=T) # nqx1; copies of fixed effects.
    eta.denom<- Fixed.denom+z%*%umat # nqx1, copies of linear predictors based on previously simulated u's.
    pij.denom<- exp(eta.denom)/(1+exp(eta.denom)) # probability of Bernoulli trail.
    # For U
    density.y.denom <- dbinom(ymat,1,pij.denom) # nqx1; Binomial density at the linear predictor.
    # For X:
    ym <- matrix(rep(ymat,1),n*q,1,byrow=T)
    # missing data process
    # Probability of observing
    if (alpha[1]==!0) {p.miss <- 1/(1+ exp(- alpha0 - alpha[1]*ym[1:(n*q-n),])) #MAR
    } else if (alpha[2]==!0 | alpha[3]==!0) { p.miss <- 1/(1+ exp(- alpha0 - alpha[2]*ym[1:(n*q),]-alpha[3]*xmat))#MNAR
    } else { p.miss <- 1/(1+ exp(- alpha0)) } #MCAR
    R.denom <- rbinom(n*q,1,p.miss)
    density.R.denom <- dbinom(R.denom,1,p.miss)
    density.X.denom <- matrix(rep(density.y.denom,1),n*q,1,byrow=T)*density.R.denom 

    density.U.denom<- matrix(density.y.denom,n,q) # nxq; q copies of Binomial density for U
    density.X.denom<- matrix(density.X.denom,n,q) # nxq; q copies of product of two Binomial density for X
    density.Y.denom<- matrix(density.R.denom,n,q) # nxq; q copies of Binomial density for Y

    apdenominator.u <- density.U.denom%*%ones # nx1; jth element is the product of Binomial pmfs at previously simulated values of uj.
    apdenominator.x <- density.X.denom
    apdenominator.y <- density.Y.denom
    
    ## Calculate Numerator of Acceptance criterion vectors
    
    ## For U
    umat2 <- u.candidate # nx1; elements with new candidates.  
    xmat2 <- x.candidate # nqx1; elements with new candidates.
    ymat2 <- y.candidate # nqx1; elements with new candidates.
    
    fixed.numer <-theta.k[1]+theta.k[2]*xmat2+theta.k[3]*data$Zij+theta.k[4]*data$Zistar # nqx1 fixed effects in linear predictor.
    Fixed.numer <- matrix(rep(fixed.numer,1),n*q,1,byrow=T) # nqx1; copies of fixed effects.
    eta.numer <- Fixed.numer+z%*%umat2 # nqx1; copies of linear predictors based on previously simulated u's.
    pij.numer <-exp(eta.numer)/(1+exp(eta.numer))
    
    # For U
    density.y.numer <- dbinom(ymat2,1,pij.numer) # nqx1; jth column has Binomial pmfs with all u's equal 
    # to previously simulated values, except with uj replaced by new candidate.
    density.y.numer <- matrix(density.y.numer,n,q) # nxq; q copies of Binomial density.
    # For X:
    ym <- matrix(rep(ymat2,1),n*q,1,byrow=T)
    # missing data process
    # Probability of observing
    if (alpha[1]==!0) {p.miss <- 1/(1+ exp(- alpha0 - alpha[1]*ym[1:(n*q-n),])) #MAR
    } else if (alpha[2]==!0 | alpha[3]==!0) { p.miss <- 1/(1+ exp(- alpha0 - alpha[2]*ym[1:(n*q),]-alpha[3]*xmat2))#MNAR
    } else { p.miss <- 1/(1+ exp(- alpha0)) } #MCAR
    R.numer <- rbinom(n*q,1,p.miss)
    density.R.numer <- dbinom(R.numer,1,p.miss)
    density.X.numer <- matrix(rep(density.y.numer,1),n*q,1,byrow=T)*density.R.numer
    
    density.U.numer<- matrix(density.y.numer,n,q) # nxq; q copies of Binomial density for U
    density.X.numer<- matrix(density.X.numer,n,q) # nxq; q copiesof product of two Binomial density for X
    density.Y.numer<- matrix(density.R.numer,n,q) # nxq; q copies of Binomial density for Y
    
    apnumerator.u <- density.U.numer%*%ones # nx1; jth element is the product of Binomial pmfs at previously simulated values of uj.
    apnumerator.x <- density.X.numer
    apnumerator.y <- density.Y.numer
    
    ## Acceptance criterion for u and x.      
    ap.u <-ifelse(apnumerator.u < apdenominator.u,apnumerator.u/apdenominator.u,1)
    #acceptance probability of jth element for candidate uj.
  
    ap.x <-ifelse(apnumerator.x < apdenominator.x,apnumerator.x/apdenominator.x,1) 
    #acceptance probability of jth element for candidate xj.
    
    ap.y <-ifelse(apnumerator.y < apdenominator.y,apnumerator.y/apdenominator.y,1) 
    #acceptance probability of jth element for candidate xj.
    
    ## Acceptance rule.
    # Generate a bernoulli random variable with the candidate u acceptance probabilities.
    ber.u<-rbinom(n,1,ap.u) 
    ber.x<-rbinom(n*q,1,ap.x)
    ber.y<-rbinom(n*q,1,ap.y)
    Ustar[m,]<-ber.u*u.candidate+(1-ber.u)*Ustar[m-1,]  # nx1; new u is either previous u or new candidate.
    xstar[m,]<-ber.x*x.candidate+(1-ber.x)*xstar[m-1,] # nx1; new x is either previous x or new candidate.
    ystar[m,]<-ber.y*y.candidate+(1-ber.y)*ystar[m-1,] # nx1; new x is either previous y or new candidate.
  }
  burnin <- burnin
  b <- burnin+1
  index <- b:M
  Ustar <- Ustar[index,]
  xstar <- xstar[index,]
  ystar <- ystar[index,]
  list(xstar=xstar,Ustar=Ustar, ystar=ystar) # Return M simulated values
}





