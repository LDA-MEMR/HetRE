############################################################################
# Implementation of HOM model
############################################################################

HOM <- function(data, n, q, sigma.error = 0, beta, varu=1.2, alpha0 =1.4, alpha = c(0,0,0), delta, lambda){
  n<-n
  q<-q
  data <- data
  # Initial parameters for M-H algorithm
  mean.x <- 1
  var.x <- 2
  mu_f <- (mean.x*sigma.error^2 + data$Wij*var.x)/(sigma.error^2 + var.x)
  var_f <- (sigma.error^2 * var.x) /(sigma.error^2 + var.x)
  theta.old <- c(beta, varu)
  alpha.old <- c(alpha0,alpha)
  
  estimate <- MCEM.HOM(n, q, data=data, M=1000, tol = 0.05, theta.old=c(beta,varu), burnin = 500, alpha.old=c(alpha0,alpha),
                        mu_f = mu_f, var_f = var_f)
  
  return(estimate = estimate)
}

