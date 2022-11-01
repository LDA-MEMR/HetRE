############################################################################
# SIMULATION EXAMPLE: n=100, q=5 under MAR
############################################################################

# Required libraries
library(mvtnorm)
source("Sim.aux.R")

data.gen <- function (n=100,q=5,sigma.error = 0, beta = c(1,2,2,2),alpha0 = 1.4,alpha = c(0,0,0),
                      delta = c(0.5, 0.3),lambda = c(0.1, 0.2), B=200){
set.seed(2022)
sigma_sq.u = exp(lambda[1]+lambda[2]*(matrix(rep(sample(c(0,1), n, replace = TRUE),q),n,q)))  # average of IV parameterization
varu = round(mean(sigma_sq.u[,1]),2)

Full.data <- vector("list", B)

for(k in 1:B){
  set.seed(2022+k)
  Full.data[[k]] <- sim.data.HET(n=n,q=q,sigma.error=sigma.error,beta=beta,alpha0=alpha0,alpha=alpha,
                                 delta=delta,lambda=lambda)
}

return(list(Full.data = Full.data, 
            beta=beta, varu = varu, alpha0 = alpha0, alpha = alpha, sigma.error = sigma.error, 
            delta = delta, lambda = lambda))
}


# Let's generate data on n = 100 subjects each with 5 visits

# No error
data <- data.gen (n=100,q=5,sigma.error = 0.0,beta = c(1,1,2,2), alpha0 = 1.4, alpha = c(1,0,0),
                      delta = c(0.5, 0.3),lambda = c(0.1, 0.2), B=5)

save(data, file="SimHET_MAR_mildMEq5.Rdata")

# Moderate error
data <- data.gen (n=100,q=5,sigma.error = 0.5,beta = c(1,2,2,2),alpha0 = 1.4,alpha = c(1,0,0),
                  delta = c(0.5, 0.3),lambda = c(0.1, 0.2), B=200)

save(data, file="SimHET_MAR_moderateMEq5.Rdata")

# Severe error
data <- data.gen (n=100,q=5,sigma.error = 0.9,beta = c(1,2,2,2),alpha0 = 1.4,alpha = c(1,0,0),
                  delta = c(0.5, 0.3),lambda = c(0.1, 0.2), B=5)

save(data, file="SimHET_MAR_severeMEq5.Rdata")



## Different initial value:

# Moderate error
data <- data.gen (n=100,q=5,sigma.error = 0.5,beta = c(1,1,2,2),alpha0 = 1.4,alpha = c(1,0,0),
                  delta = c(0.5, 0.3),lambda = c(0.1, 0.2), B=200)

save(data, file="SimHET_MAR_moderateMEq5.Rdata")

