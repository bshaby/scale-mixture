
source("scalemix_aux_fns.R")
source("scalemix_likelihoods.R")
source("scalemix_priors.R")
source("~/Dropbox/code_bank/generic.samplers.R")


#######################################################
#
# MODEL
#
# Y[,t]    = beta.0 + mu[,t] + rnorm(n.s, 0, sig2e)
# mu[,t]   = N(0,sig2s[t]*Matern)
# sig2s[2] ~ InvGamma(a,b)
#
#
# INPUTS
#
# Y        := ns x nt matrix of observations
# S        := ns x  2 matrix of spatial coordinates
# thresh   := censoring threshold
#
# mean_nu  := log(nu) ~ N(mean_nu,sd_nu)
# sd_nu
# range_nu := log(rangee) ~ N(mean_range,sd_range)
# sd_range
#
# n.iters    := number of MCMC iteration
# burn     := length of burn-in
#
#######################################################

MCMC_scalemix_HOT <- function(Y, S, thresh,
                              mean_nu=0, sd_nu=1,
                              mean_range=-1, sd_range=1,
                              n.iters=1000, thin=10,
                              experiment.name="scalemix", echo.interval=50,
                              sigma.m.theta.scale=1, sigma.m.theta.corr=1,
                              sigma.m.scale=NULL,
                              prop.Sigma=NULL, adapt.cov=TRUE,
                              true.sigma.s=NULL, true.mu=NULL,
                              true.sigma.e=NULL, true.beta.0=NULL,
                              true.gamma=NULL, true.beta=NULL,
                              true.rho=NULL, true.nu=NULL) {
  

  n.metr.updates.theta.scale <- 5
  n.metr.updates.theta.corr <- 2
  n.metr.updates.scale <- 2
  hyper.params.scale <- c(1,1)
  hyper.params.corr <- c(1,1)
  
  eps <- 1e-06
  
  # Bookkeeping
  n.s <- nrow(Y)
  n.t <- ncol(Y)
  cen <- Y < thresh
  h <- rdist(S)
  diag(h)  <- 0
  theta.mn <- c(mean_range,mean_nu)
  theta.sd <- c(sd_range,sd_nu)
  
  # Initial values
  beta.0 <- mean(Y) # Intercept
  mu      <- 0.9*(Y-beta.0)
  tau.s    <- 1/apply(mu, 2, var)
  tau.e    <- 25*mean(tau.s)
  sigma.e <- 1/sqrt(tau.e)
  sigma.s <- pmax(1.001, 1/sqrt(tau.s))

  theta.scale <- c(1, 1)  # (gamma, beta)
  theta.corr <- c(0.1, 1.0)
  Sigma <- corr.fn(h, log(theta.corr))
  
  #**********   Eigendecomposition of the correlation matrix   ****************#
  eig.Sigma <- eigen(Sigma, symmetric=TRUE)
  V <- eig.Sigma$vectors
  d <- eig.Sigma$values
  
  # Allocate trace objects
  theta.corr.trace <- matrix(NA, n.iters, 2)
  theta.scale.trace <- matrix(NA, n.iters, 2)
  params.trace <-  matrix(NA, n.iters, 2)
  sigma.s.trace <- matrix(NA, n.iters, n.t)
  mu.trace <- array(NA, c(n.iters, n.s, n.t))
  Y.trace <- array(NA, c(n.iters, n.s, n.t))
  
  colnames(theta.corr.trace) <- c("range", "nu")
  colnames(theta.scale.trace) <- c("gamma", "beta")
  colnames(params.trace) <- c("nugget", "beta.0")
  
  # Fill in trace objects with initial values
  theta.corr.trace[1, ] <- theta.corr
  theta.scale.trace[1, ] <- theta.scale
  sigma.s.trace[1, ] <- sigma.s
  params.trace[1, ] <- c(sigma.e, beta.0)

  # Temperature parameters for adaptation
  c.0 <- 10
  c.1 <- 0.8
  k <- 3
  metr.opt.theta.scale <- 0.35
  metr.opt.theta.corr <- 0.35
  metr.opt.scale <-0.41
  
  # For Metropolis updates of theta
  if (is.null(sigma.m.theta.scale)) sigma.m.theta.scale <- 1
  if (is.null(sigma.m.theta.corr)) sigma.m.theta.corr <- 1
  if (is.null(sigma.m.scale)) sigma.m.scale <- rep(1, n.t)
  r.hat.theta.scale <- NA
  r.hat.theta.corr <- NA
  r.hat.scale <- rep(NA, n.t)
  
  if (is.null(prop.Sigma$theta.scale)) {
    prop.Sigma.theta.scale <- diag(2)
  } else {
    prop.Sigma.theta.scale <- prop.Sigma$theta.scale
  }
  if (is.null(prop.Sigma$theta.corr)) {
    prop.Sigma.theta.corr <- diag(2)
  } else {
    prop.Sigma.theta.corr <- prop.Sigma$theta.corr
  }
  
  # Start timer
  time.0 <- proc.time()[3]
  
    
  ######################################################################
  # Begin main loop
  
  for(i in 1:n.iters){
    
    for (thin.iter in 1:thin) {
      ##############################################:
      #####      IMPUTE LATENT Y (Gibbs)     #######:
      ##############################################:
  
      # Draws a collection of univariate truncated normals
      B      <- pnorm(thresh, beta.0 + mu[cen], sigma.e)
      U      <- B * runif(sum(cen))
      Y[cen] <- qnorm(U, beta.0 + mu[cen], sigma.e)
      
      # Plot the imputations
      n.pl <- which.min(apply(cen, 2, sum))  # Choose a year to plot impuations
      par(mfrow=c(1,1))
      plot(S[ ,1], Y[ ,n.pl], type="l", ylim=c(0, 20), lwd=2)
      clip(min(S[ ,1]), max(S[ ,1]), -10, thresh)
      lines(S[ ,1], Y[ ,n.pl], col="gray80", lwd=3)
      abline(h=thresh)
          
      
  
      ##############################################:
      #####      UDPATE MU (Gibbs)           #######:
      ##############################################:
      
      for(t in 1:n.t){
        mu[ ,t] <- normal.mean.conjugate.update.from.eigen(Y[ ,t]-beta.0, 0, V, d, sigma.s[t], sigma.e)
      }
      
      # FOR DEBUGGING
      # mu <- true.mu
      
      ##############################################:
      #####       MEAN PARAMETERS (Gibbs)    #######:
      ##############################################:
      
      post.prec.beta.0 <- tau.e*n.s*n.t + 0.001
      post.mean.beta.0 <- tau.e*sum(Y-mu)
      beta.0 <- rnorm(1, post.mean.beta.0/post.prec.beta.0, 1/sqrt(post.prec.beta.0))
      
      # FOR DEBUGGING
      # beta.0 <- true.beta.0
      
      ##############################################:
      #####     SCALING FACTORS              #######:
      ##############################################:
      
      for(t in 1:n.t){
        metr.out <- static.metr(mu[ ,t], sigma.s[t], normal.scale.likelihood,
                                huser.thibaud.prior, hyper.params=theta.scale,
                                n.updates=n.metr.updates.scale,
                                prop.Sigma=1,
                                sigma.m=sigma.m.scale[t],
                                V=V, d=d)
        r.hat.scale[t] <- metr.out$acc.prob
        sigma.s[t] <- metr.out$trace[n.metr.updates.scale, ]
      }
      
      # FOR DEBUGGING
       # sigma.s <- true.sigma.s
      
      ##############################################:
      #####     NUGGET VARIANCE (Gibbs)      #######:
      ##############################################:
      
      tau.e <- rgamma(1, n.s*n.t/2+0.1, sum((Y-mu-beta.0)^2)/2+0.1) 
      sigma.e <- 1/sqrt(tau.e)
      # FOR DEBUGGING
      # tau.e <- 1/true.sigma.e^2
      # sigma.e <- true.sigma.e
      
      ##############################################:
      #### CORRELATION PARAMETERS (Metropolis) #####:
      ##############################################:
      metr.out.theta.corr <- static.metr(mu, theta.corr, normal.mu.likelihood,
                                         half.cauchy.half.cauchy, hyper.params=hyper.params.corr,
                                         n.updates=n.metr.updates.theta.corr,
                                         prop.Sigma=prop.Sigma.theta.corr,
                                         sigma.m=sigma.m.theta.corr,
                                         sigma.s=sigma.s, h=h)
      r.hat.theta.corr <- metr.out.theta.corr$acc.prob
      theta.corr <- metr.out.theta.corr$trace[n.metr.updates.theta.corr, ]

      # FOR DEBUGGING
      # theta.corr <- c(true.rho, true.nu)
      
            
      # This is inefficient since we've already done this decomposition in the
      # above likelihood, but doing it again makes the code more modular.
      Sigma <- corr.fn(h, log(theta.corr))                   # now each column has unit standard dev.
      eig.Sigma <- eigen(Sigma, symmetric=TRUE)
      V <- eig.Sigma$vectors
      d <- eig.Sigma$values
      
          
      
      ##############################################:
      ####    UPDATE SCALING PARAMETERS (METR)  ####:
      ##############################################:
      metr.out.theta.scale <- static.metr(sigma.s, theta.scale, HOT.mix.likelihood,
                                          half.cauchy.half.cauchy, hyper.params.scale,
                                          n.updates=n.metr.updates.theta.scale,
                                          prop.Sigma=prop.Sigma.theta.scale,
                                          sigma.m=sigma.m.theta.scale)
      r.hat.theta.scale <- metr.out.theta.scale$acc.prob
      theta.scale <- metr.out.theta.scale$trace[n.metr.updates.theta.scale, ]
      # FOR DEBUGGING
      # theta.scale <- c(true.gamma, true.beta)
      
  
      #################################################################
      #  Adapt proposal for theta                                     #
      gamma1 <- c.0 / (i + k)^(c.1)
      gamma2 <- 1 / (i + k)^(c.1)
      sigma.m.theta.corr <- exp(log(sigma.m.theta.corr) +
                                gamma1*(r.hat.theta.corr - metr.opt.theta.corr))
      sigma.m.theta.scale <- exp(log(sigma.m.theta.scale) +
                                 gamma1*(r.hat.theta.scale - metr.opt.theta.scale))
      for (t in 1:n.t) {
        sigma.m.scale[t] <- exp(log(sigma.m.scale[t]) +
                                    gamma1*(r.hat.scale[t] - metr.opt.scale))
        
      }
      
      if (adapt.cov) {
        prop.Sigma.theta.corr <- prop.Sigma.theta.corr +
          gamma2*(cov(metr.out.theta.corr$trace) - prop.Sigma.theta.corr)
        prop.Sigma.theta.scale <- prop.Sigma.theta.scale +
          gamma2*(cov(metr.out.theta.scale$trace) - prop.Sigma.theta.scale)
        
        while(is(try(chol(prop.Sigma.theta.corr), silent=TRUE), "try-error")) {
          prop.Sigma.theta.corr <- prop.Sigma.theta.corr + eps*diag(p)
          cat("Oops. Proposal covariance matrix is now:\n")
          print(prop.Sigma.theta.corr)
        }
        while(is(try(chol(prop.Sigma.theta.scale), silent=TRUE), "try-error")) {
          prop.Sigma.theta.scale <- prop.Sigma.theta.scale + eps*diag(p)
          cat("Oops. Proposal covariance matrix is now:\n")
          print(prop.Sigma.theta.scale)
        }
      }
      #                                                               #
      #################################################################

    } # End thin loop
    
    ##############################################:
    #####        KEEP TRACK OF STUFF       #######:
    ##############################################:
    
    theta.corr.trace[i, ] <- theta.corr
    theta.scale.trace[i, ] <- theta.scale
    sigma.s.trace[i, ] <- sigma.s
    params.trace[i, ] <- c(sigma.e, beta.0)
    Y.trace[i, , ] <- Y
    mu.trace[i, , ] <- mu
    
    
    
    
    ##############################################:
    #####       PLOT RESULTS SO FAR        #######:
    ##############################################:
    
    if ((i %% echo.interval) == 0) {
      gibbs.out <- list(theta.corr.trace=theta.corr.trace, 
                        theta.scale.trace=theta.scale.trace, 
                        sigma.s.trace=sigma.s.trace, 
                        params.trace=params.trace,
                        mu.trace=mu.trace, Y.trace=Y.trace,
                        Y=Y, S=S, thresh=thresh,
                        sigma.m.theta.scale=sigma.m.theta.scale,
                        sigma.m.theta.corr=sigma.m.theta.corr,
                        i=i)
      save(gibbs.out, file=paste(experiment.name, "_progress.RData", sep=""))
      pdf(file=paste(experiment.name, "_progress.pdf", sep=""))
      par(mar=c(2, 3, 3, 1))
      par(mfrow=c(3,2))
      plot(theta.corr.trace[1:i, 1], type="l",
           xlab="MCMC iteration",ylab="Sample",
           main=expression(rho))
      if (!is.null(true.rho)) abline(h=true.rho, col="red", lty=2, lwd=3)
      plot(theta.corr.trace[1:i, 2], type="l",
           xlab="MCMC iteration",ylab="Sample",
           main=expression(nu))
      if (!is.null(true.nu)) abline(h=true.nu, col="red", lty=2, lwd=3)
      plot(theta.scale.trace[1:i, 1], type="l",
           xlab="MCMC iteration",ylab="Sample",
           main=expression(gamma))
      if (!is.null(true.gamma)) abline(h=true.gamma, col="red", lty=2, lwd=3)
      plot(theta.scale.trace[1:i, 2], type="l",
           xlab="MCMC iteration",ylab="Sample",
           main=expression(beta))
      if (!is.null(true.beta)) abline(h=true.beta, col="red", lty=2, lwd=3)
      plot(params.trace[1:i, 1], type="l",
           xlab="MCMC iteration",ylab="Sample",
           main=expression(sigma[e]))
      if (!is.null(true.sigma.e)) abline(h=true.sigma.e, col="red", lty=2, lwd=3)
      plot(params.trace[1:i, 2], type="l",
           xlab="MCMC iteration",ylab="Sample",
           main=expression(beta[0]))
      if (!is.null(true.beta.0)) abline(h=true.beta.0, col="red", lty=2, lwd=3)
      par(mfrow=c(5, 3))
      for (t in 1:min(n.t, 30)) {
        plot(sigma.s.trace[1:i, t], type="l", ylim=range(sigma.s.trace, na.rm=TRUE),
             xlab="MCMC iteration",ylab="Sample",
             main=substitute(sigma[s][index], list(index=t)))
        if (!is.null(true.sigma.s)) abline(h=true.sigma.s[t], col="red", lty=2, lwd=3)
      }
      dev.off()
      elapsed.time <- round((proc.time()[3] - time.0) / 60 / 60, 2)
      remaining.time <- round((elapsed.time/i) * (n.iters-i), 2)
      cat("Iteration", i, "out of", n.iters, "complete.\n")
      cat("  Elapsed time:", elapsed.time, "hours\n")
      cat("  Remaining time:", remaining.time, "hours\n\n")
    }
    # End main loop
    ######################################################################
    
    
  }
  
  return(gibbs.out)
  
}


