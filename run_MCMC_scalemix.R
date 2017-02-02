
#Test it out:
library(fields)
library(geoR)

source("scalemix_aux_fns.R")

# Simulation settings

n.s <- 200        # Number of sites
n.t <- 50         # Number of time points
sigma.e <- 0.1    # Nugget SD
beta <- 0.5       # Spatial variance ~ H-O-T(gamma, beta)
gamma <- 1
rho <- 0.1          # Matern range
nu <- 1.0         # Matern smoothness
beta.0 <- 10      # Overall mean


# Threshold for fitting
thresh <- 13


experiment.name <- "test_HOT"

# Generate fake data

theta  <- c(rho,nu)   

# S     <- cbind(runif(ns),runif(ns))
S     <- cbind(seq(0, 1, length=n.s), rep(1, n.s))
Cor   <- corr.fn(rdist(S), log(theta))
C.Cor <- chol(Cor)

# set.seed(1111)
sigma.s <- rhuser.thibaud(n.t, c(gamma, beta))

Y <- matrix(NA, n.s, n.t)
mu <- matrix(NA, n.s, n.t)
for(t in 1:n.t) {
  mu[ ,t] <- sigma.s[t]*crossprod(C.Cor, rnorm(n.s))
  Y[ ,t] <- mu[ ,t] + sigma.e*rnorm(n.s)
  plot(Y[ ,t], type="l", lwd=3)
  lines(mu[ ,t], col=2, lwd=2)
  abline(h=thresh, lty=2, lwd=3, col="gray80")
  # readline()
}
Y <- Y + beta.0

true.sigma.s <- sigma.s
sum(Y > thresh) / length(Y)

# Fit the model
source("scalemix_sampler_eigen.R")

starting.prop.Sigma.theta.corr <- matrix(c( 9.6e-06, -1.4e-05,
                                            -1.4e-05,  7.0e-05), 2, 2, byrow=TRUE)
starting.prop.Sigma.theta.scale <- matrix(c( 0.04, -0.04,
                                            -0.04,  0.06), 2, 2, byrow=TRUE)

fit.2 <- MCMC_scalemix_HOT(Y, S, thresh=thresh, 
                           n.iters=10000, thin=10, echo.interval=100,
                           experiment.name=experiment.name,
                           true.sigma.s=true.sigma.s, 
                           true.sigma.e=sigma.e,
                           true.beta.0=beta.0,
                           true.gamma=gamma, true.beta=beta,
                           true.rho=rho, true.nu=nu,
                           true.mu=mu,
                           prop.Sigma=list(theta.corr=starting.prop.Sigma.theta.corr,
                                           theta.scale=starting.prop.Sigma.theta.scale))






################################################################################
# Test scale param sampling
sigma.s <- rhuser.thibaud(50, c(gamma, beta))
sigma.m.theta.scale <- 0.25
n.updates.test <- 10000
metr.out <- static.metr(sigma.s, c(0.5, 0.5), HOT.mix.likelihood,
                        half.cauchy.half.cauchy, c(1,1),
                        n.updates=n.updates.test,
                        prop.Sigma=diag(2),
                        sigma.m=sigma.m.theta.scale)
cat(metr.out$acc.prob, "\n")
par(mfrow=c(2,1))
plot(metr.out$trace[ ,1], type="l")
abline(h=gamma, col=2, lwd=3)
plot(metr.out$trace[ ,2], type="l")
abline(h=beta, col=2, lwd=3)






