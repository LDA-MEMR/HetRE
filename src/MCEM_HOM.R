############################################################################
# Function for MCEM under HOM model
############################################################################
#' @title Monte Carlo Expected-Maximization algorithm steps.
#' @description This function evaluate the MCEM algorithm.
#' @param n sample size
#' @param q number for follow-up for each individual
#' @param data data sets contain the original data
#' @param theta.old initial values of parameters
#' @param alpha.old initial parameters of missing data mechanism
#' @param mu_f mean of X given W 
#' @param var_f variance of X given W 
#' @param M Number of MC size
#' @param tol given tolerance value for EM algorithm
#' @return MLEs of model parameter of interest
#' @example 

### Function for calculating the E and M-steps of MCEM algorithm.
MCEM.HOM <- function(n, q, data, M=1000, tol = 0.0005, theta.old, burnin = 500, alpha.old, mu_f=0, var_f= 1)
{ 
  burnin <- burnin
  Mm <- M-burnin
  k <- 0
  theta.k <- c(theta.old,alpha.old)
  beta <- c(theta.k[1],theta.k[2],theta.k[3],theta.k[4])
  varu <- theta.k[5]
  alphaA <- c(theta.k[6], theta.k[7],theta.k[8],theta.k[9])
  alpha0 <- theta.k[6]
  alpha <- c(theta.k[7],theta.k[8],theta.k[9])

  Zij <- data$Zij
  Zistar <- data$Zistar
  eij <- data$eij
  
  ## E-M step of EM algorithm. 
  convergence <- F
  while(1-convergence){
    k <- k+1
    # Implementation of M-H algorithm
    sample <- metropolis.HOM(n, q, data=data, theta.k = theta.k, sigma.error=sigma.error, mu_f=mu_f, 
                             var_f=var_f, alpha0=alpha0, alpha=alpha, M=M, burnin = burnin) 
    #sample <- metropolis(theta.k, M)
    U <- sample$Ustar # Mxn
    Xij <- t(sample$xstar) # nqxM
    # Generate surroagte using classical additive measurement error model
    Wij <- Xij+eij
    Yij <- t(sample$ystar) # nqxM
    
    yvector <- matrix(Yij,ncol=1,byrow=T) 
    
  #########################################################
  # Likelihood function for fixed effects estimate    #
  #########################################################   
    ll.beta <- function(beta)
    {
      beta0<-beta[1]
      beta1<-beta[2]
      beta2<-beta[3]
      beta3<-beta[4]
      fixed<-beta0 + beta1*Wij+beta2*Zij+beta3*Zistar
      Fixed<-matrix(rep(fixed,1),n*q*Mm,1,byrow=T)
      z<-matrix(0,nrow=n*q,ncol=n) # nqxn, Design matrix.
      for(j in 1:n)
      {
        z[(1+(j-1)*q):(q*j),j]=1
      }
      uvector <- matrix(z%*%t(U),ncol=1,byrow=T)
      eta <- Fixed+uvector
      p <- exp(eta)/(1+exp(eta))
      d <- dbinom(yvector,1,p, log=TRUE)
      sum(d)/Mm
    }
  #########################################################
  # Likelihood function for random effects estimate    #
  #########################################################
    ll.random<-function(varu)
    {
      varu2 <- varu
      uvector2 <- matrix(U,ncol=1,byrow=T)
      z<-matrix(0,nrow=n*q,ncol=n) # nqxn, Design matrix.
      for(j in 1:n)
      {
        z[(1+(j-1)*q):(q*j),j]=1
      }
      uvector<-matrix(z%*%t(U),ncol=1,byrow=T)
      sum(log(dnorm(uvector,0,sqrt(varu2))))/Mm
    }
    
  #########################################################
  # Likelihood function for missing data mechanism    #
  #########################################################
    ll.alpha <- function(alphaA){
      alpha0 <-alphaA[1]
      alpha1 <-alphaA[2]
      alpha2 <-alphaA[3]
      alpha3 <-alphaA[4]
      yc <- yvector 
      # Probability of observing
      if (alpha[1]==!0) { p.miss <- 1/(1+ exp(- alpha0 - alpha1*yc[1:(n*q*Mm-n*Mm),])) #MAR
      } else if (alpha[2]==!0 | alpha[3]==!0) { p.miss <- 1/(1+ exp(- alpha0 - alpha2*yc[1:(n*q*Mm),]-alpha3*Wij)) #MNAR
      } else { p.miss <- 1/(1+ exp(- alpha0)) } #MCAR
      R.numer <- rbinom(n*q,1,p.miss)
      density.R <- log(dbinom(R.numer,1,p.miss))
      sum(density.R)/Mm
    }
    
    # updated estimates for parameters   
    beta.update <- optim(fn = ll.beta, par = c(theta.k[1],theta.k[2],theta.k[3],theta.k[4]),control=list(fnscale=-1),
                                                        method="BFGS",hessian = TRUE)
    beta.new <- beta.update$par
    # note that this may give some values that are too large or not finite:
    ifelse(beta.new<0 | beta.new>5, next, beta.new)
    # updated estimates for random parameter  
    varu.new <- optimize(ll.random,interval=c(0,5),maximum=TRUE)$maximum
    
    # updated estimates for missing data parameter 
    alpha.update <- optim(fn = ll.alpha, par = c(theta.k[6],theta.k[7],theta.k[8],theta.k[9]),control=list(fnscale=-1),
                          method="L-BFGS-B",hessian = TRUE)
    alpha.new <- alpha.update$par
    ifelse(alpha.new<0 | alpha.new>5, next, alpha.new)
    
    theta.new <-list(beta.new, varu.new, alpha.new)
    theta.new<-unlist(theta.new)
    # Checking the convergence   
    if (max(abs(c(theta.new[1]-theta.k[1],theta.new[2]-theta.k[2],theta.new[3]-theta.k[3],theta.new[4]-theta.k[4],
                  theta.new[5]-theta.k[5],theta.new[6]-theta.k[6],theta.new[7]-theta.k[7],theta.new[8]-theta.k[8],
                  theta.new[9]-theta.k[9])))>tol){
      
      theta.k <-theta.new      
      beta <-c(theta.new[1],theta.new[2],theta.new[3],theta.new[4])
      varu <-c(theta.new[5])
      alpha0<-theta.new[6]
      alpha <-c(theta.new[7],theta.new[8],theta.new[9])
      alphaA <- c(theta.new[6],theta.new[7],theta.new[8],theta.new[9])
    }
    else convergence <- T
  }
  if(convergence){
    print(paste('Convergence was reached after',k,'iteration(s)'))
    mle <- theta.new  
    ##### Calculation of standard errors for parameters.
    # Calculation of standard errors for parameters from "optim" function.
    # NEGATIVE of the hessian is the "fishers observed information"
    I.beta <- solve(-beta.update$hessian)
    se.beta <- sqrt(diag(I.beta))
    
    se.rand <- sqrt(2*mle[5]^2/(q*Mm)) # sd of RE variance
    
    I.alpha <- solve(-alpha.update$hessian)
    se.alpha <- sqrt(abs(diag(I.alpha)))
    
    se<-c(se.beta,se.rand,se.alpha)  
    
    estimates <- c(mle,se)
    return(estimates)
  }
  else print(paste('Convergence was not reached'))
}

