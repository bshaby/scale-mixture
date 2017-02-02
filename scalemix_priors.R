

################################################################################
#  Computes two independent half Cauchy densities with scale hyper.params.     #

half.cauchy.half.cauchy <- function(params, hyper.params) {
  theta.1 <- params[1]
  theta.2 <- params[2]
  scale.theta.1 <- hyper.params[1]
  scale.theta.2 <- hyper.params[2]
  
  if (min(params)<0) return(-Inf)
  return(log(2) + dcauchy(theta.1, scale=scale.theta.1, log=TRUE) + 
         log(2) + dcauchy(theta.2, scale=scale.theta.2, log=TRUE))
}
#                                                                              #
################################################################################



################################################################################
#  Computes two independent half Cauchy densities with scale hyper.params.     #

huser.thibaud.prior <- function(params, hyper.params) {
  r <- params[1]
  return(dhuser.thibaud(r, theta=hyper.params, log=TRUE))
}
#                                                                              #
################################################################################

