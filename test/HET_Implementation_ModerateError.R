
############################################################################
# Sample script to demonstrate HET model 
############################################################################
# load Required Packages
library(mvtnorm)

# Load source files
source("sim.aux.R")
source("M-H_algorithm_HET.R")
source("MCEM_HET.R")
source("HET.R")

##################################################################################################
# An Illustration of HET model with n=100, q=5 under missing data scenario 2 (alpha = c(1,0,0))
##################################################################################################

# set model parameters
n <- 100
q <- 5
sigma.error <- 0.5 # moderate measurement error
beta <- c(1,2,2,2)
delta <- c(0.5, 0.3)
lambda <- c(0.1, 0.2)
alpha0 <- 1.4 # generate 20% missing values
alpha <- c(1,0,0) # missing data scenario 1 (MCAR)

#set.seed(2022)
#sigma_sq.u = exp(lambda[1]+lambda[2]*(matrix(rep(sample(c(0,1), n, replace = TRUE),q),n,q)))
#varu = round(mean(sigma_sq.u[,1]),2)

# generate data
set.seed(1234)
data.sim <- sim.data.HET(n=n,q=q,sigma.error=sigma.error,beta=beta,alpha0=alpha0,alpha=alpha,delta=delta,lambda=lambda)

## Implementation of HET model
# pars <- c(beta, delta,lambda,alpha0,alpha)

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



  
# Load the auxiliary functions located in 'R' folder which is
# needed to generate data
#source(here::here("SCRIPT", "R", "LDvine.aux.R"))
#source(/src/simulation_functions.R"))

#source("src/sim.aux.R")
#source("src/M-H_algorithm_HET.R")
#source("src/MCEM_HET.R")
#source("src/HET.R")
  
  
  
  
