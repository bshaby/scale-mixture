



################################################################################
################################################################################
## Compute the density of the mixing distribution in Huser-Opitz-Thibaud (2017).
## Theta is the vector (gamma, beta)', in the notation in the manuscript.
## The support of r is [1, infinity)
##
## r ................................. value at which to evaluate the density
## theta ............................. a 2-vector of parameters:
##                                     theta[1] = gamma (from H-O-T (2017))
##                                     theta[2] = beta (from H-O-T(2017))
##
dhuser.thibaud <- function(r, theta, log=TRUE) {
  
  gamma <- theta[1]
  beta <- theta[2]
  
  eps <- 1e-06
  
  in.support <- r > 1
  r <- r[in.support]

  if (log) {
    dens <- rep(-Inf, max(1, length(in.support)))
    if (beta > eps) {
      dens[in.support] <- log(gamma) + (beta-1)*log(r) - (gamma/beta * (r^beta-1))
    } else {
      dens[in.support] <- log(gamma) - (gamma+1)*log(r)
    }
  } else {
    dens <- rep(0, max(1, length(in.support)))
    if (beta > eps) {
      dens[in.support] <- exp(-gamma/beta * (r^beta-1)) * gamma * r^(beta-1)
    } else {
      dens[in.support] <- gamma * r^(-(gamma+1))
    }    
  }
  return(dens)
}
##
################################################################################
################################################################################

























################################################################################
################################################################################
## Draws from the mixing distribution in Huser-Opitz-Thibaud (2017).
## Theta is the vector (gamma, beta)', in the notation in the manuscript.
## The support of r is [1, infinity)
##
## n ................................. the number of draws to make
## theta ............................. a 2-vector of parameters:
##                                     theta[1] = gamma (from H-O-T (2017))
##                                     theta[2] = beta (from H-O-T(2017))
##
rhuser.thibaud <- function(n, theta) {
  require(geoR)
  
  
  gamma <- theta[1]
  beta <- theta[2]
  
  eps <- 1e-06

  u <- runif(n)
  
  if (beta > eps) {
    return((-beta/gamma * log(1-u) + 1)^(1/beta))
  } else {
    return((1-u)^(-1/gamma))
  }

}
##
################################################################################
################################################################################





################################################################################
################################################################################
## Compute the Matern correlation function from a matrix of pairwise distances
## and a vector of parameters
##
## d ................................. a square matrix of distances
## theta ............................. a 2-vector of correlation parameters:
##                                     theta[1] = log(phi) (from the geoR 'matern' function)
##                                                **usually rho**
##                                     theta[2] = log(kappa) (from the geoR 'matern' function)
##                                                **usually nu**
##
corr.fn <- function(d, theta) {
  require(geoR)
  
  phi <- exp(theta[1])
  kappa <- exp(theta[2])
  
  Sigma <- matern(d, phi, kappa)
  diag(Sigma) <- 1
  return(Sigma)
}
##
################################################################################
################################################################################

################################################################################
################################################################################
## Multivariate normal log density of R, where each column of 
## R iid N(0,Q), where Q is the precision (inverse covariance matrix
## It essentially computes the log density of each column of R, then takes
## the sum.  Faster than looping over the columns, but not as transparent.
dmvn <- function(R,Q){
  
  0.5*ncol(R)*determinant(Q)$modulus - 0.5*sum(R*(Q%*%R))
} 
##
################################################################################
################################################################################


################################################################################
################################################################################
## Multivariate normal log density of R, where each column of 
## R iid N(0,C'C), where C'C is the covariance matrix
## It essentially computes the log density of each column of R, then takes
## the sum.  Faster than looping over the columns, but not as transparent.
##
## Is there a faster way to to this? i.e. by using the quadform function?
##
## R ................................. a (s x t) matrix, where each column is
##                                     multivariage normal, independently from
##                                     the other columns
## C ................................. a (s x s) upper triangular matrix 
##                                     such that C'C=Sigma, as from a call to
##                                     chol(Sigma)
##
dmvn.chol <- function(R, C){
  
  return(-0.5*ncol(R)*chol2logdet(C) - 0.5 * sum(R * chol2inv.times.vector(C, R)))
} 
##
################################################################################
################################################################################



################################################################################
################################################################################
## Multivariate normal log density of R, where each column of 
## R iid N(0,V'DV), where V'DV is the covariance matrix
## It essentially computes the log density of each column of R, then takes
## the sum.  Faster than looping over the columns, but not as transparent.
##
## Is there a faster way to to this? i.e. by using the quadform function?
##
## R ................................. a (s x t) matrix, where each column is
##                                     multivariage normal, independently from
##                                     the other columns
## V ................................. a (s x s) matrix of eigenvectors of the
##                                     covariance matrix of the columns of R
## d.inv ............................. a (s) vector of the inverted eigenvalues,
##                                     d.inv_i = 1/D_ii
##
dmvn.eig <- function(R, V, d.inv){
  if (is.vector(R)) n.rep <- 1 else n.rep <- ncol(R)
  return(-0.5*n.rep*eig2logdet(1/d.inv) - 0.5 * sum(R * eig2inv.times.vector(V, d.inv, R)))
} 
##
################################################################################
################################################################################



################################################################################
################################################################################
## This is slow
normal.mean.conjugate.update.from.chol <- function(y, prior.mean, C.Sigma, 
                                                   sigma.y, sigma.prior) {
  ns <- length(y)
  
  regularized.precision <- 1/sigma.y^2 * chol2inv(C.Sigma)
  diag(regularized.precision) <- diag(regularized.precision) + 1/sigma.prior^2
  
  C.regularized.precision <- chol(regularized.precision)
  
  post.mean <- chol2inv.times.vector(C.regularized.precision, y) / sigma.prior^2
  # return(post.mean)
  return(post.mean + backsolve(C.regularized.precision, rnorm(ns), transpose=TRUE))
}
##
################################################################################
################################################################################


################################################################################
################################################################################
## This is slow
normal.mean.conjugate.update.from.inverse <- function(y, prior.mean, Sigma.inv, 
                                                      sigma.y, sigma.prior) {
  ns <- length(y)
  
  regularized.precision <- 1/sigma.y^2 * Sigma.inv
  diag(regularized.precision) <- diag(regularized.precision) + 1/sigma.prior^2
  
  C.regularized.precision <- chol(regularized.precision)
  
  post.mean <- chol2inv.times.vector(C.regularized.precision, y) / sigma.prior^2
  # return(post.mean)
  return(post.mean + backsolve(C.regularized.precision, rnorm(ns), transpose=TRUE))
}
##
################################################################################
################################################################################


################################################################################
################################################################################
## This is faster
normal.mean.conjugate.update.from.eigen <- function(y, prior.mean, V, d, 
                                                    sigma.y, sigma.prior) {
  ns <- length(y)
  
  regularized.d.inv <- 1/d * 1/sigma.y^2 + 1/sigma.prior^2
  post.mean <- eig2inv.times.vector(V, 1/regularized.d.inv, y) / sigma.prior^2
  
  return(post.mean + eig2inv.times.vector(V, 1/sqrt(regularized.d.inv), rnorm(ns)))
}
##
################################################################################
################################################################################


################################################################################
################################################################################
## Assumes that A = C'C, where C is a upper triangular matrix (e.g. returned 
## by a call to chol(A))
##
## Computes A^{-1}x
##
## C ................................. an upper triangular Cholesky factor, as above
## x ................................. a vector x, as above
##
chol2inv.times.vector <- function(C, x) {
  return(backsolve(C, forwardsolve(C, x, upper.tri=TRUE, transpose=TRUE)))
}
##
################################################################################
##################################################################################


################################################################################
################################################################################
## Assumes that A = C'C, where C is a upper triangular matrix (e.g. returned 
## by a call to chol(A))
##
## Computes x'A^{-1}x
##
## C ................................. an upper triangular Cholesky factor, as above
## x ................................. a vector x, as above
##
chol2inv.quadform.vector <- function(C, x) {
  bs <- backsolve(C, x, transpose=TRUE)
  return(drop(crossprod(bs)))
}
##
################################################################################
################################################################################


################################################################################
################################################################################
## Assumes that A = C'C, where C is a upper triangular matrix (e.g. returned 
## by a call to chol(A))
##
## Computes log(|A|)
##
## C ................................. an upper triangular Cholesky factor, as above
## vector ............................ a vector x, as above
##
chol2logdet <- function(C) {
  return(2*sum(log(diag(C))))
}
##
################################################################################
################################################################################



################################################################################
################################################################################
## Assumes that A = VDV', where D is a diagonal vector of eigenvectors of A, and
## V is a matrix of normalized eigenvectors of A.
## Computes A^{-1}x
##
## V ................................. a matrix of eigenvectors, as above
## d.inv ............................. a vector of the inverted eigenvalues, with
##                                     d.inv_i = 1/D_ii
## x ................................. a vector x, as above
##
eig2inv.times.vector <- function(V, d.inv, x) {
  return(V %*% (d.inv * crossprod(V, x)))
}
##
################################################################################
################################################################################

################################################################################
################################################################################
## Assumes that A = VDV', where D is a diagonal vector of eigenvectors of A, and
## V is a matrix of normalized eigenvectors of A.
##
## Computes x'A^{-1}x
##
## V ................................. a matrix of eigenvectors, as above
## d.inv ............................. a vector of the inverted eigenvalues, with
##                                     d.inv_i = 1/D_ii
## x ................................. a vector x, as above
##
eig2inv.quadform.vector <- function(V, d.inv, x) {
  cp <- crossprod(V, x)
  return(crossprod(cp, d.inv*cp))
}
##
################################################################################
################################################################################

################################################################################
################################################################################
## Just computes x'Ax
##
## A ................................. a square matris
## x ................................. a vector
##
quadform.vector <- function(A, x) {
  return(crossprod(x, A %*% x))
}
##
################################################################################
################################################################################



################################################################################
################################################################################
## Assumes that A = VDV', where D is a diagonal vector of eigenvectors of A, and
## V is a matrix of normalized eigenvectors of A.
##
## log(|A|)
##
## d ................................. a vector of the eigenvalues, with
##                                     d_i = D_ii
##
eig2logdet <- function(d) {
  return(sum(log(d)))
}
##
################################################################################
################################################################################


################################################################################
################################################################################
##                     FOR MARGINAL TRANSFORMATIONS                           ##
################################################################################
################################################################################


################################################################################
################################################################################
## The thing that gets integrated dr to result in the marginal CDF of X
##
## 
##
## theta ............................. A vector, (beta, gamma), of parameters
##                                     that get passed to dhuser.thibaud
##
mix.distn.integrand <- function(r, x, theta) {
  prod <- rep(0, length(r))
  in.support <- r > 0
  
  prod[in.support] <- pnorm(x, sd=r[in.support]) * dhuser.thibaud(r[in.support], theta, log=FALSE)
  return(prod)
}
##
################################################################################
################################################################################

################################################################################
################################################################################
## Integrates mix.distn.integrand dr to result in the marginal CDF of X
##
## theta ............................. A vector, (beta, gamma), of parameters
##                                     that get passed to dhuser.thibaud
##
pmixture.uni <- function(x, theta) {
  integ <- integrate(mix.distn.integrand, 0, Inf, x=x, theta=theta)
  return(integ$value)
}
pmixture <- Vectorize(pmixture.uni, "x")

##
################################################################################
################################################################################

################################################################################
################################################################################
## Numerically inverts the marginal CDF of X using uniroot to get the maringal
## quantile function at probability p.  Only searches for a root on the
## interval (lower, upper).
##
## theta ............................. A vector, (beta, gamma), of parameters
##                                     that get passed to dhuser.thibaud
##
qmixture.uni <- function(p, theta, lower=-1000, upper=1000) {
  if (p == 0) return(NA)
  if (p == 1) return(NA)
  
  return(uniroot((function (x) pmixture(x, theta) - p), lower=lower, upper=upper)$root)
}
qmixture <- Vectorize(qmixture.uni, "p")
##
################################################################################
################################################################################

################################################################################
################################################################################
## Approximates the marginal quantile function by taking values of the
## marginal CDF of X and doing linear interpolation.  If no values of the CDF
## are supplied, it computes n.x of them, for x in (lower, upper).
##
## theta ............................. A vector, (beta, gamma), of parameters
##                                     that get passed to dhuser.thibaud
##
qmixture.uni.interp <- function(p, cdf.vals) {
  return(NLSstClosestX(cdf.vals, p))
}
qmixture.vec.interp <- Vectorize(qmixture.uni.interp, "p")

qmixture.interp <- function(p, cdf.vals = NULL, x.vals=NULL,
                            n.x=100, lower=-100, upper=100, theta=NULL) {
  
  if (is.null(x.vals)) {
    x.vals <- seq(lower, upper, length=n.x)
    cdf.vals <- pmixture(x.vals, theta)
  } else {
    if (is.null(cdf.vals)) {
      cdf.vals <- pmixture(x.vals, theta)
    }
  }
  cdf.vals <- sortedXyData(x.vals, cdf.vals)

  return(qmixture.vec.interp(p, cdf.vals))
  
}
##
################################################################################
################################################################################




