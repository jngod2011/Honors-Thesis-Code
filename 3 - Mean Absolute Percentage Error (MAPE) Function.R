########################################################################
##### MEAN ABSOLUTE PERCENTAGE ERROR (MAPE) FUNCTION
########################################################################

# Calculates the MAPE of a Log ACD1(p,q) model using EE parameter estimates and AR-method initial values
# x is the full-length original observations vector
# n is the number of points used to fit the Log ACD model
# L is the number of steps to predict from forecast origin n
# This function only works for orders (1,1) and (2,1)
# momentsdist is the distribution used in the EE parameter estimation procedure
# For exponential moments, par1=lambda
# For weibull moments, par1=alpha and par2=beta
# For gamma moments, par1=k (shape) and par2=theta (scale)
# For none (no distribution specified, look at user input moments

fmape.logacd1 <- function(x, n, L, p, q, momentsdist, par1=1, par2=1) {
  # make sure that n and L are valid
  if (n < 0 | L < 0) { cat("Error: n and L must be positive") } else {
    if (n > length(x)) { cat("Error: Please specify an n <= the length of the observed data") } else {
      if (n+L > length(x)) { cat("Error: The sum of the number of points used to fit the model and the number of points used for prediction must not exceed the length of the observations vector") } else {
        # define fit and prediction parts of x
        fit <- x[1:n]
        pred.actual <- x[(n+1):(n+L)]
        # psi values associated with fit part of x
        psi.fit <- log(x[1:n])
        # Log ACD1(1,1) Case   
        if (p==1 & q==1) {
          # EE Parameter estimation
          if (momentsdist=="exponential") {
            param <- festeq.logacd1(fit, 1,1, momentsdist="exponential", par1, par2, arinitval="TRUE")
          } else if (momentsdist=="gamma") {
            param <- festeq.logacd1(fit, 1,1, momentsdist="gamma", par1, par2, arinitval="TRUE")
          } else if (momentsdist=="weibull") {
            param <- festeq.logacd1(fit, 1,1, momentsdist="weibull", par1, par2, arinitval="TRUE")
          } else { cat("Error: Unrecognized moments distribution") }
          # define omega, alpha, beta estimates
          omega <- param[1]; alpha <- param[2]; beta <- param[3]
          # Prediction (j=1:L steps from forecast origin n) given n values of x
          pred.psi <- rep(NA, L)
          pred.xhat <- rep(NA, L)
          # For j=1
          pred.psi[1] <- omega + alpha*log(x[(n+1)-1]) + beta*psi.fit[(n+1)-1]
          pred.xhat[1] <- exp(pred.psi[1])
          # For j=2:L steps from n	
          for (j in 2:L) {
            pred.psi[j] <- omega + alpha*log(x[(n+j)-1]) + beta*pred.psi[j-1]
            pred.xhat[j] <- exp(pred.psi[j])
          }
          # MAPE Calculation			
          mape <- 100*(mean(abs((pred.actual-pred.xhat)/pred.actual)))
          cat("The mean absolute percentage error (MAPE) is:", mape, "\n")
          mape
        } else if (p==2 & q==1) {
          # EE Parameter estimation
          if (momentsdist=="exponential") {
            param <- festeq.logacd1(fit, 2,1, momentsdist="exponential", par1, par2, arinitval="TRUE")
          } else if (momentsdist=="gamma") {
            param <- festeq.logacd1(fit, 2,1, momentsdist="gamma", par1, par2, arinitval="TRUE")
          } else if (momentsdist=="weibull") {
            param <- festeq.logacd1(fit, 2,1, momentsdist="weibull", par1, par2, arinitval="TRUE")
          } else { cat("Error: Unrecognized moments distribution") }
          # define omega, alpha1, alpha2, beta estimates
          omega <- param[1]; alpha1 <- param[2]; alpha2 <- param[3]; beta <- param[4]
          # Prediction (j=1:L steps from forecast origin n) given n values of x
          pred.psi <- rep(NA, L)
          pred.xhat <- rep(NA, L)
          # For j=1
          pred.psi[1] <- omega + alpha1*log(x[(n+1)-1]) + alpha2*log(x[(n+1)-2]) + beta*psi.fit[(n+1)-1]
          pred.xhat[1] <- exp(pred.psi[1])
          # For j=2:L steps from n	
          for (j in 2:L) {
            pred.psi[j] <- omega + alpha1*log(x[(n+j)-1]) + alpha2*log(x[(n+j)-2]) + beta*pred.psi[j-1]
            pred.xhat[j] <- exp(pred.psi[j])
          }
          # MAPE Calculation			
          mape <- 100*(mean(abs((pred.actual-pred.xhat)/pred.actual)))
          cat("The mean absolute percentage error (MAPE) is:", mape, "\n")
          mape
        } else { cat("This function is only defined for the Log ACD1(1,1) and Log ACD1(2,1) cases") }	
      } 
    }
  }
}