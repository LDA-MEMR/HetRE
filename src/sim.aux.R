############################################################################
## AUXILLARY FUNCTIONS FOR SIMULATION OF LONGITUDINAL DATA ##     
############################################################################

# Required libraries
#library(mvtnorm)

############################################################################
# Function for generating random effects by proposed approach (using MCD)
############################################################################
#' @title Function for generating random effects using MCD
#'
#' @description This function generate random effects.
#'
#' @param n sample size
#' @param q number for follow-up for each individual
#' @param delta GARP parameters
#' @param lambda IV parameters
#' @param Zi subject specific variable
#' @return A list
#' @examples
#' n=10 
#' q=5
#' delta=c(0.5, 0.3)
#' lambda=c(0.1, 0.2)
#' Zi = matrix(rep(sample(c(0,1), n, replace = TRUE),q),n,q)
#' library(mvtnorm)
#' uu.sim = u.sim.HET(n,q,delta = delta, Zi=Zi, lambda = lambda)
#' head(uu.sim)
#'@export

u.sim.HET <- function(n,q, delta, Zi, lambda){
  u.all <- matrix(NA,n,q)
  # Indicator variable to generate GAR and IV parameters
  I <- matrix(NA,q,q)  
  for (i in 1:n){
    for (j in 1:q){
      for (t in (1:q)){
        #I[j,t]=ifelse((abs(j-t)==1),1,0)
        I[j,t]=if(abs(j-t)==1) 1 else 0
      }
    }
    phi <- delta[1]*I+delta[2]*I*Zi[i,]
    T.hat<-phi+diag(1,q,q)
    T.hat[upper.tri(T.hat)] <-0
    T.hat[col(T.hat) < row(T.hat)] <- -1*T.hat[col(T.hat)<row(T.hat)]
    sigma_sq <- sqrt(exp(lambda[1]+lambda[2]*Zi[i,]))
    D.hat<-diag(sigma_sq,q,q)
    Sigma<-solve(T.hat)%*%D.hat%*%t(solve(T.hat))
    mean<-rep(0,q)
    ## Load relevant R library
    #require(mvtnorm)
    u.all[i,]<-rmvnorm(1,mean,Sigma)
  }
  U.all=matrix(t(u.all),ncol=1,byrow=T)
  return(list(U.all=U.all,u.all=u.all,I=I))
}


############################################################################ 
# Function to Simulate longitudinal data using heterogeneous RE model
############################################################################
#' @title Simulation of longitudinal data using heterogeneous RE model
#'
#' @description This function simulates longitudinal copula data using heterogeneous RE model.
#'
#' @param n sample size
#' @param q number for follow-up for each individual
#' @param sigma.error measurement error variation
#' @param beta regression parameters 
#' @param alpha0 percentage of missing values
#' @param alpha missing data mechanism parameters
#' @param delta GARP parameters
#' @param lambda IV parameters
#' @return A data frame
#' @examples
#'
#' library(mvtnorm)
#' n = 50
#' q = 5
#' sigma.error = 0.5
#' beta = c(-5.7,1.5,-0.75,1)
#' alpha0 = 1.4 #generate 20% missing values
#' alpha = c(0,0,0)
#' delta = c(0.5, 0.3)
#' lambda = c(0.1, 0.2)
#' data = sim.data.HET(n=n,q=q,sigma.error=sigma.error,beta=beta,alpha0=alpha0,alpha=alpha,delta=delta,lambda=lambda)
#' head(data)
#' @export

# Data generation
sim.data.HET <- function(n=100, q=5, sigma.error=0, beta, alpha0=1.4, alpha, delta, lambda){
  
## Simulate the true model
id <- rep(c(1:n), rep(q, n))
## Generate covariates (Z and Zstar) based on information given in Section 4.
vi<-rnorm(n*q,0,0.5)
Eij<-rnorm(n*q,0,0.5)
Zij<-vi+Eij
Zi <- matrix(rep(sample(c(0,1), n, replace = TRUE),q),n,q)
Zistar<-as.vector(matrix(Zi,ncol=1,byrow=T))
## Generate error covariate X based on information given in Section 4.
mu_x<-1
ai<-rnorm(n*q,0,1)
eps_ij<-rnorm(n*q,0,1)
Xij<- mu_x + ai + eps_ij
# Generate surrogate using classical additive measurement error model
sigma <- sigma.error
eij<-rnorm(n*q,0,sigma)
Wij<-Xij+eij

# response model using surrogate variable
fixed <- beta[1]+beta[2]*Wij+beta[3]*Zij+beta[4]*Zistar 
Fixed <- matrix(rep(fixed,1),n*q,1,byrow=T)
# Generate random effects
uvector <- u.sim.HET(n,q,delta = delta, Zi = Zi, lambda = lambda)
uvector <- uvector$U.all
eta <- Fixed+uvector
p.y <- exp(eta)/(1+exp(eta))
y.true <- rbinom(n*q,1,p.y) 
ym <- matrix(y.true,ncol=1,byrow=T) 

#### simulate missing data process
alpha <- alpha 
# Probability of observing
if (alpha[1]==!0) {p.miss <- 1/(1+ exp(- alpha0 - alpha[1]*ym[1:(n*q-n),])) #MAR
} else if (alpha[2]==!0 | alpha[3]==!0) { p.miss <- 1/(1+ exp(- alpha0 - alpha[2]*ym[1:(n*q),]-alpha[3]*Wij))#MNAR
} else { p.miss <- 1/(1+ exp(- alpha0)) } #MCAR

# =====================================
# R = 1 for observed; R = 0 for missing
# =====================================
R <- matrix(rbinom(n*q,1,p.miss),n*q,1,byrow=T)
for (i in 1:(n*q)){if (R[i,]==0) R[i,]=NA}
Y.obs <- R*ym  # observed Y including missing values

data.sim <- data.frame(id, y.true, p.y, Xij, eij, Wij, Zij, Zistar, R, Y.obs)
return(data.sim)
}



############################################################################ 
# Function to Simulate longitudinal data using homogeneous RE model
############################################################################

#' @title Simulation of longitudinal data using homogeneous RE model
#'
#' @description This function simulates longitudinal copula data using homogeneous RE model.
#'
#' @param n sample size
#' @param q number for follow-up for each individual
#' @param sigma.error measurement error variation
#' @param beta regression parameters 
#' @param alpha0 percentage of missing values
#' @param alpha missing data mechanism parameters
#' @param delta GARP parameters
#' @param lambda IV parameters
#' @return A data frame
#' @examples
#'
#' library(mvtnorm)
#' n = 50
#' q = 5
#' sigma.error = 0.4
#' beta = c(-5.7,1.5,-0.75,1)
#' alpha0 = 1.4 #generate 20% missing values
#' alpha = c(0,1,-1)
#' lambda = c(0.1, 0.2)
#' varu=1.2
#' data = sim.data.HOM(n=n,q=q,sigma.error=sigma.error,beta=beta,varu=varu,alpha0=alpha0,alpha=alpha,lambda=lambda)
#' head(data)
#' @export

# Data generation
sim.data.HOM <- function(n=100, q=5, sigma.error=0, beta, varu=1.2, alpha0=1.4, alpha, lambda){
  
  ## Simulate the true model
  id <- rep(c(1:n), rep(q, n))
  ## Generate covariates (Z and Zstar) based on information given in Section 4.
  vi<-rnorm(n*q,0,0.5)
  Eij<-rnorm(n*q,0,0.5)
  Zij<-vi+Eij
  Zi<-matrix(rep(sample(c(0,1), n, replace = TRUE),q),n,q)
  Zistar<-as.vector(matrix(Zi,ncol=1,byrow=T))
  ## Generate error covariate X based on information given in Section 4.
  mu_x<-1
  ai<-rnorm(n*q,0,1)
  eps_ij<-rnorm(n*q,0,1)
  Xij<- mu_x + ai + eps_ij
  # Generate surrogate using classical additive measurement error model
  sigma <- sigma.error
  eij<-rnorm(n*q,0,sigma) #sigma=0.4
  Wij<-Xij+eij
  # response model using surrogate variable
  fixed <- beta[1]+beta[2]*Wij+beta[3]*Zij+beta[4]*Zistar 
  Fixed <- matrix(rep(fixed,1),n*q,1,byrow=T)
  # Generate random effects
  #sigma_sq.u <- exp(lambda[1]+lambda[2]*Zi)  # average of IV parameterization
  #varu <- round(mean(sigma_sq.u[,1]),2)
  uvector <- rnorm(n, 0, sqrt(varu))
  z <- matrix(0,nrow=n*q,ncol=n) # nqxn, Design matrix.
  for(j in 1:n)
  {
    z[(1+(j-1)*q):(q*j),j]=1
  }
  uvector <- matrix(z%*%uvector, ncol=1, byrow=T)
  
  eta <- Fixed+uvector
  p.y <- exp(eta)/(1+exp(eta))
  y.true <- rbinom(n*q,1,p.y) 
  ym <- matrix(y.true,ncol=1,byrow=T) 
  
  #### simulate missing data process
  alpha <- alpha 
  # Probability of observing
  if (alpha[1]==!0) {p.miss <- 1/(1+ exp(- alpha0 - alpha[1]*ym[1:(n*q-n),])) #MAR
  } else if (alpha[2]==!0 | alpha[3]==!0) { p.miss <- 1/(1+ exp(- alpha0 - alpha[2]*ym[1:(n*q),]-alpha[3]*Wij))#MNAR
  } else { p.miss <- 1/(1+ exp(- alpha0)) } #MCAR
  
  # =====================================
  # R = 1 for observed; R = 0 for missing
  # =====================================
  R <- matrix(rbinom(n*q,1,p.miss),n*q,1,byrow=T)
  for (i in 1:(n*q)){if (R[i,]==0) R[i,]=NA}
  Y.obs <- R*ym  # observed Y including missing values
  
  data.sim <- data.frame(id, y.true, p.y, Xij, eij, Wij, Zij, Zistar, R, Y.obs)
  return(data.sim)
}




