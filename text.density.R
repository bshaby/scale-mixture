

source("scalemix_aux_fns.R")

n <- 500
beta.grid <- seq(0, 10, length=n)
r <- 5
gamma <- 1

dens <- rep(NA, n)

for (i in 1:n) {
  dens[i] <- dhuser.thibaud(r, c(gamma, beta.grid[i]), log=FALSE)
}

plot(beta.grid, dens, type="l")



r.grid <- seq(0, 10, length=n)
beta <- 0.5
gamma <- 1

plot(r.grid, dhuser.thibaud(r.grid, c(gamma, beta), log=FALSE), type="l")

r.dens <- dhuser.thibaud(r.grid, c(gamma, beta), log=FALSE)
r.samps <- rhuser.thibaud(1000000, c(gamma, beta))
plot(density(r.samps, bw=0.1), xlim=c(0, 10), ylim=c(0, max(r.dens)))
lines(r.grid, dhuser.thibaud(r.grid, c(gamma, beta), log=FALSE), col=2)

