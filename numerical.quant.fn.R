
source("scalemix_aux_fns.R")


mix.distn.integrand <- function(r, x, theta) {
  prod <- rep(0, length(r))
  in.support <- r > 0
  
  prod[in.support] <- pnorm(x, sd=r[in.support]) * dhuser.thibaud(r[in.support], theta, log=FALSE)
  return(prod)
}

beta <- 0.5       # Spatial variance ~ H-O-T(gamma, beta)
gamma <- 1
theta <- c(beta, gamma)

xx <- rhuser.thibaud(10000, theta)
hist(xx, probability = TRUE, breaks=50)
x.vals<- seq(0, 50, length=500)
lines(x.vals, dhuser.thibaud(x.vals, theta, log=FALSE))

n.samps <- 1e04
r <- rhuser.thibaud(n.samps, theta)
z <- rnorm(n.samps)

y <- r*z

emp.cdf <- ecdf(y)
plot(emp.cdf)
p.vals <- seq(0, 1, length=500)
emp.quantile <- quantile(y, probs=p.vals)
plot(p.vals, emp.quantile, type="l")

r.vals <- seq(0, 100, length=200)
plot(r.vals, mix.distn.integrand(r.vals, 2, theta), type="l")

integrate(mix.distn.integrand, 0, 100, x=2, theta=theta)


pmixture.uni <- function(x, theta) {
  integ <- integrate(mix.distn.integrand, 0, Inf, x=x, theta=theta)
  return(integ$value)
}
pmixture <- Vectorize(pmixture.uni, "x")

qmixture.uni <- function(p, theta) {
  if (p == 0) return(NA)
  if (p == 1) return(NA)
  
  lower <- -100
  upper <- 100
  return(uniroot((function (x) pmixture(x, theta) - p), lower=lower, upper=upper)$root)
}
qmixture <- Vectorize(qmixture.uni, "p")


qmixture.uni(0.5, theta)
pmixture.uni(1, theta)

qmixture.uni(pmixture.uni(2, theta), theta)


## Test CDF
n <- 1000
x.vals <- seq(-100, 100, length=n)
system.time(cdf.vals <- pmixture(x.vals, theta))
plot(emp.cdf)
lines(x.vals, cdf.vals, col="red", lwd=3, lty=2)

## Test quantile function
n.samps <- 1e06
n <- 1000
p.vals <- seq(0, 1, length=n)
system.time({r <- rhuser.thibaud(n.samps, theta); 
             z <- rnorm(n.samps);
             y <- r*z;
             emp.quantile <- quantile(y, probs=p.vals)})
system.time(quant.vals <- qmixture(p.vals, theta))
plot(p.vals, emp.quantile, type="l")
lines(p.vals, quant.vals, col="red", lwd=3, lty=2)

## Try logspline
library(logspline)
n.samps <- 1e04
n <- 1000
x.vals <- seq(-100, 100, length=n)
# CDF
system.time({r <- rhuser.thibaud(n.samps, theta); 
             z <- rnorm(n.samps);
             y <- r*z;
             logspline.fit <- logspline(y);
             plogspline.vals <- plogspline(x.vals, logspline.fit)})

system.time(cdf.vals <- pmixture(x.vals, theta))
plot(x.vals, plogspline.vals, type="l")
lines(x.vals, cdf.vals, col="red", lwd=3, lty=2)

# Quantile
p.vals <- seq(0, 1, length=n)
system.time(qlogspline.vals <- qlogspline(p.vals, logspline.fit))
system.time(quant.vals <- qmixture(p.vals, theta))
plot(p.vals, quantile(y, probs=p.vals), type="l")
plot(p.vals, quant.vals, type="l")
lines(p.vals, qlogspline.vals, col="red", lwd=3, lty=2)

# Use interpolation to find the quantile function
qinterp.uni <- function(p, data) {
  return(NLSstClosestX(data, p))
}
qinterp <- Vectorize(qinterp.uni, "p")

system.time(qinterp.vals <- qinterp(p.vals, data=sortedXyData(x.vals, cdf.vals)))
plot(p.vals, quant.vals, type="l")
lines(p.vals, qinterp.vals, col="red", lwd=3, lty=2)

plot(quant.vals - qlogspline.vals)
points(quant.vals - qinterp.vals, col="red")
points(quant.vals - quantile(y, probs=p.vals), col="green")


################################################################################
## Try the sourced functions
rm(list=ls())
source("scalemix_aux_fns.R")

beta <- 0.5       # Spatial variance ~ H-O-T(gamma, beta)
gamma <- 1
theta <- c(beta, gamma)


n <- 1000
p.vals <- seq(0, 1, length=n)
system.time(quant.vals.interp <- qmixture.interp(p.vals, theta=theta, n.x=500))
system.time(quant.vals.numer <- qmixture(p.vals, theta))

plot(p.vals, quant.vals.numer, type="l")
lines(p.vals, quant.vals.interp, col="red", lwd=3, lty=2)




