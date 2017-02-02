

################################################################################
################################################################################
## For the generic Metropolis sampler
## Samples from gamma and beta from the H-O-T scale mixing distribution.
## params is the vector (gamma, beta)', in the notation in the manuscript.
##
## data............................... a n.t-vector of scaling factors
## params............................. a 2-vector of parameters (theta in dhuser.thibaud):
##                                     theta[1] = gamma (from H-O-T (2017))
##                                     theta[2] = beta (from H-O-T(2017))
##

HOT.mix.likelihood <- function(data, params) {
  # The params vector is the same as the theta vector in the likelihood function
  return(sum(dhuser.thibaud(data, params, log=TRUE)))
}

#                                                                              #
################################################################################


################################################################################
################################################################################
## For the generic Metropolis sampler
## Samples from gamma and beta from the H-O-T scale mixing distribution.
## params is the vector (gamma, beta)', in the notation in the manuscript.
##
## data............................... a (n.s x n.t) matrix of normal r.v.s, with
##                                     mean zero and mu[i] having sd sigma.s[i]
## params............................. a 2-vector of correlation parameters:
##                                     theta[1] = rho (range parameter)
##                                     theta[2] = nu (smoothness parameter)
## h ................................. a (n.s x n.s) matrix of distances
## sigma.s ........................... a n.t-vector of standard deviations
##

normal.mu.likelihood <- function(data, params, h, sigma.s) {

  mu <- data
  theta.corr <- params
  
  mu.standardized <- sweep(mu, 2, sigma.s, "/")     # divide each column by its standard dev.
  Sigma <- corr.fn(h, log(theta.corr))                   # now each column has unit standard dev.
  #**********   Eigendecomposition of the correlation matrix   **************#
  eig.Sigma <- eigen(Sigma, symmetric=TRUE)
  V <- eig.Sigma$vectors
  d <- eig.Sigma$values
  
  return(dmvn.eig(mu.standardized, V, 1/d)) 

}

#                                                                              #
################################################################################


################################################################################
## For the generic Metropolis sampler
## Samples from the scale factors in a scale mixture of Gaussians, where the   
## mixing distribution comes from from the H-O-T scale mixing distribution.
## 
##
## data............................... a n.s-vector of N(0, Sigma) r.v.s, 
##                                     Sigma = VDV'
## params............................. a single scaling parameter
## V ................................. a (n.s x n.s) matrix whose columns are
##                                     the eigenvectors of Sigma
## d ................................. a n.s-vector eigenvalues of Sigma
##

normal.scale.likelihood <- function(data, params, V, d) {
  
  mu <- data
  sigma.s <- params
  
  return(dmvn.eig(mu, V, 1/sigma.s^2 * 1/d)) 
}

#                                                                              #
################################################################################


















