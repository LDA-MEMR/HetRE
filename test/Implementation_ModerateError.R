
############################################################################
# Sample script to demonstrate HET and HOM models
############################################################################
# load Required Packages
library(mvtnorm)

# Load source files
source("sim.aux.R")
source("M-H_algorithm_HET.R")
source("MCEM_HET.R")
source("HET.R")

source("M-H_algorithm_HOM.R")
source("MCEM_HOM.R")
source("HOM.R")

# Load the auxiliary functions located in 'R' folder which is
#source("src/sim.aux.R")
#source("src/M-H_algorithm_HET.R")
#source("src/MCEM_HET.R")
#source("src/HET.R")
#source("src/M-H_algorithm_HOM.R")
#source("src/MCEM_HOM.R")
#source("src/HOM.R")

###########################################################################
# An Illustration of analysis with simulated data
###########################################################################

# set model parameters

n <- 100
q <- 5
sigma.error <- 0.5 # moderate measurement error
beta <- c(1,2,2,2)
delta <- c(0.5, 0.3)
lambda <- c(0.1, 0.2)
alpha0 <- 1.4 # generate 20% missing values
alpha <- c(1,0,0) # missing data scenario 2 (MAR)

set.seed(1234)
sigma_sq.u = exp(lambda[1]+lambda[2]*(matrix(rep(sample(c(0,1), n, replace = TRUE),q),n,q)))
varu = round(mean(sigma_sq.u[,1]),2)

###########################################################################
# generate data 
###########################################################################

set.seed(1234)
data.sim <- sim.data.HET(n=n,q=q,sigma.error=sigma.error,beta=beta,
                         alpha0=alpha0,alpha=alpha,delta=delta,lambda=lambda)

###########################################################################
## Implementation of HET model
###########################################################################
start <- Sys.time()
estimate  <- tryCatch(HET(data = data.sim, n=n, q=q, sigma.error = sigma.error, 
                            beta = beta, alpha0 = alpha0, alpha = alpha, 
                            delta = delta, lambda = lambda),
                        error=function(e){return(NA)}
                        )
end<- Sys.time() 
com.time <- difftime(end, start, units = "secs")
print(com.time)
print(estimate)


###########################################################################
## Implementation of HOM model
###########################################################################

start <- Sys.time()
estimate  <- tryCatch(HOM(data = data.sim, n, q, sigma.error = sigma.error, 
                          beta = beta, varu= varu, alpha0 = alpha0, alpha = alpha, 
                          delta = delta, lambda = lambda),
                      error=function(e){return(NA)}
)
end<- Sys.time() 
com.time <- difftime(end, start, units = "secs")
print(com.time)
print(estimate)



  



  
  
  
  
