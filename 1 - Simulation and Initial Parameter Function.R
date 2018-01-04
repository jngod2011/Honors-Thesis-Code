########################################################################
##### SIMULATION FUNCTION
########################################################################

# Log ACD1(m,q) Simulation function: Simulate from the Log ACD1(m,q) model
	# n is the final sample size (after burn-in of nb time points) 
	# parameters are omega, alpha, beta
		# alpha and beta can be vectors, depending on the number of lags to be included
		# alpha = (alpha_1, alpha_2, ..., alpha_m )
		# beta = (beta_1, beta_2, ..., beta_q )
	# feps is the distribution of the errors epsilon
		# par1 and par2 are parameters of feps 
		# par2 defaults to 1
		# function is defined for the exponential, weibull, and gamma distributions
			# in the exponential function, par1=lambda
			# in the weibull function, par1=alpha and par2=beta
			# in the gamma function, par1=k (shape) and par2=theta (scale)

fsim.logacd1 <- function (n, nb, omega, alpha, beta, feps, par1, par2=1) {
# Check for stationarity constraint
if (sum(alpha)+sum(beta)>=1) { cat("Error: The simulated process is not weakly stationary")
} else {
	# total number of simulated datapoints
		nt <- n + nb		
	# Initialize psis, xs, x, and psi
		psis <- rep(1,nt)
		xs <-rep(1,nt)
		x <-rep(NA,n)		# we'll save the last n entries of xs into x 
		psi <-rep(NA,n)		# we'll save the last n entries of psis into psi
	if (feps=="exponential") {
		# randomly generate errors from exponential distribution
			eps <- rexp(nt,par1) 
		# generate xs and psis
			# Initial psis and xs values (depends on maximum number of lags)
			for (t in (1:max(length(alpha), length(beta))) ) {
				# define temporary lagged.logxs and lagged.psis vectors
					# initialization
						lagged.logxs <- rep(NA, length(alpha))
						lagged.psis <- rep(NA, length(beta))
					# calculation 
						for (i in 1:length(alpha)) { 
							# if the index t-i is less than 0, then define lagged.logxs[i]=0. Else, define lagged.logxs[i] = log(xs[t-i])
								index <- t-i
								if (index <= 0) {lagged.logxs[i] <- 0} else { lagged.logxs[i] <- log(xs[index]) }
						}		
						for (j in 1:length(beta)) { 
							# if the index t-j is less than 0, then define lagged.psis[j]=0. Else, define lagged.psis[j] = log(xs[t-j])
								index <- t-j
								if (index <= 0) {lagged.psis[j] <- 0} else { lagged.psis[j] <- psis[index] }
						}		
				psis[t] = omega + alpha%*%lagged.logxs + beta%*%lagged.psis 
				xs[t] = exp(psis[t])*eps[t]
			}	
			# For t = max(length(alpha), length(beta))+1 : nt
			for (t in (max(length(alpha), length(beta))+1):nt) {
				# define temporary lagged.logxs and lagged.psis vectors
					# initialization
						lagged.logxs <- rep(NA, length(alpha))
						lagged.psis <- rep(NA, length(beta))
					# calculation 
						for (i in 1:length(alpha)) { lagged.logxs[i] <- log(xs[t-i]) }
						for (j in 1:length(beta)) { lagged.psis[j] <- psis[t-j] }
				# calculate psis and xs
				psis[t] = omega + alpha%*%lagged.logxs + beta%*%lagged.psis 
				xs[t] = exp(psis[t])*eps[t]
			}
	# save the simulation data: keep last n data and psi after burn-in save into x and psi
		x <- ts(xs[(nb+1):nt])
		psi <- ts(psis[(nb+1):nt])	
			# x is the simulated durations series (after the first nb points are deleted)	
	# output as dataframe
		output <- as.data.frame(cbind(x,psi))		
	} else if (feps=="weibull") {
		# randomly generate errors from weibull distribution
			eps <- rweibull(nt, par1, par2) 
		# generate xs and psis
			# Initial psis and xs value (depends on maximum number of lags)
			for (t in (1:max(length(alpha), length(beta))) ) {
				# define temporary lagged.logxs and lagged.psis vectors
					# initialization
						lagged.logxs <- rep(NA, length(alpha))
						lagged.psis <- rep(NA, length(beta))
					# calculation 
						for (i in 1:length(alpha)) { 
							# if the index t-i is less than 0, then define lagged.logxs[i]=0. Else, define lagged.logxs[i] = log(xs[t-i])
								index <- t-i
								if (index <= 0) {lagged.logxs[i] <- 0} else { lagged.logxs[i] <- log(xs[index]) }
						}		
						for (j in 1:length(beta)) { 
							# if the index t-j is less than 0, then define lagged.psis[j]=0. Else, define lagged.psis[j] = log(xs[t-j])
								index <- t-j
								if (index <= 0) {lagged.psis[j] <- 0} else { lagged.psis[j] <- psis[index] }
						}		
				psis[t] = omega + alpha%*%lagged.logxs + beta%*%lagged.psis 
				xs[t] = exp(psis[t])*eps[t]
			}	
			# For t = max(length(alpha), length(beta))+1 : nt
			for (t in (max(length(alpha), length(beta))+1):nt) {
				# define temporary lagged.logxs and lagged.psis vectors
					# initialization
						lagged.logxs <- rep(NA, length(alpha))
						lagged.psis <- rep(NA, length(beta))
					# calculation 
						for (i in 1:length(alpha)) { lagged.logxs[i] <- log(xs[t-i]) }
						for (j in 1:length(beta)) { lagged.psis[j] <- psis[t-j] }
				# calculate psis and xs
				psis[t] = omega + alpha%*%lagged.logxs + beta%*%lagged.psis 
				xs[t] = exp(psis[t])*eps[t]
			}
	# save the simulation data: keep last n data and psi after burn-in save into x and psi
		x <- ts(xs[(nb+1):nt])
		psi <- ts(psis[(nb+1):nt])	
			# x is the simulated durations series (after the first nb points are deleted)	
	# output as dataframe
		output <- as.data.frame(cbind(x,psi))		
	} else if (feps=="gamma") {
		# randomly generate errors from gamma distribution
			eps <- rgamma(nt,par1,par2) 
		# generate xs and psis
			# Initial psis and xs value (depends on maximum number of lags)
			for (t in (1:max(length(alpha), length(beta))) ) {
				# define temporary lagged.logxs and lagged.psis vectors
					# initialization
						lagged.logxs <- rep(NA, length(alpha))
						lagged.psis <- rep(NA, length(beta))
					# calculation 
						for (i in 1:length(alpha)) { 
							# if the index t-i is less than 0, then define lagged.logxs[i]=0. Else, define lagged.logxs[i] = log(xs[t-i])
								index <- t-i
								if (index <= 0) {lagged.logxs[i] <- 0} else { lagged.logxs[i] <- log(xs[index]) }
						}		
						for (j in 1:length(beta)) { 
							# if the index t-j is less than 0, then define lagged.psis[j]=0. Else, define lagged.psis[j] = log(xs[t-j])
								index <- t-j
								if (index <= 0) {lagged.psis[j] <- 0} else { lagged.psis[j] <- psis[index] }
						}		
				psis[t] = omega + alpha%*%lagged.logxs + beta%*%lagged.psis 
				xs[t] = exp(psis[t])*eps[t]
			}	
			# For t = max(length(alpha), length(beta))+1 : nt
			for (t in (max(length(alpha), length(beta))+1):nt) {
				# define temporary lagged.logxs and lagged.psis vectors
					# initialization
						lagged.logxs <- rep(NA, length(alpha))
						lagged.psis <- rep(NA, length(beta))
					# calculation 
						for (i in 1:length(alpha)) { lagged.logxs[i] <- log(xs[t-i]) }
						for (j in 1:length(beta)) { lagged.psis[j] <- psis[t-j] }
				# calculate psis and xs
				psis[t] = omega + alpha%*%lagged.logxs + beta%*%lagged.psis 
				xs[t] = exp(psis[t])*eps[t]
			}
	# save the simulation data: keep last n data and psi after burn-in save into x and psi
		x <- ts(xs[(nb+1):nt])
		psi <- ts(psis[(nb+1):nt])	
			# x is the simulated durations series (after the first nb points are deleted)			
	# output as dataframe
		output <- as.data.frame(cbind(x,psi))
	} else { cat("Error: The error distribution must be specified as exponential, gamma, or weibull") }
}
}	

# EXAMPLE:
# test <- fsim.logacd1(100, 1000, 1, .5, .3, feps="exponential", 1); x <- test$x

########################################################################
##### AR INITIAL VALUES FUNCTION
########################################################################

# Log ACD1 initial parameters function: Calculate initial values for parameters in the log ACD(1,1) and log ACD(2,1) model by fitting an AR(m,n) model to the observed (simulated) data
	# Given observed/simulated data x 
	# Order of the Log ACD1 model is given by (p,q)
		# Fit an ARMA(m,n) model to the data, where m=max(p,q) and n=q
		
finitval.logacd1 <- function (x, p=1, q=1) 
	{
		# Take the log transform of the observed (simulated) data and define it as y
			y <- log(x)
		# Fit an ARMA(m,n) model to the log transformed data y, using a 'high enough' order m
			ar.model <- arima(y, order=c(max(p,q),0,q), include.mean=TRUE)
		# OMEGA: omega vector is the intercept
			#omega.hat <- ar.model$coef[length(ar.model$coef)]
				MEAN <- ar.model$coef[length(ar.model$coef)]
				intercept <- MEAN*(1-sum(ar.model$coef[1:p]))
				omega.hat <- intercept
		# BETA: beta vector is equivalent to MA coefficients
				beta.hat <- -1*ar.model$coef[(length(ar.model$coef)-q):(length(ar.model$coef)-1)]
		# ALPHA:
			# phi vector equals AR coefficients
				phi.hat <- ar.model$coef[1:max(p,q)]
			# theta vector has length equal to that of the phi vector, and values equal to the MA coefficients (betas) for the first q entries
				theta.hat <- rep(0, max(p,q))
				for(i in 1:length(beta.hat)) { theta.hat[i] <- beta.hat[i] }
			# alpha vector equals phi-theta
				alpha.hat <- phi.hat-theta.hat		
					
		# Results
			results <- c(omega.hat, alpha.hat, beta.hat)
		}


