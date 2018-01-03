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

########################################################################
##### ESTIMATING EQUATIONS FUNCTION
########################################################################

# Note: This code was modified from code provided by Dr. Nalini Ravishanker to include initial parameter estimates (#INITIAL VALUES section) 
# Log ACD1 estimating equations function: Calculate parameter estimates for the log ACD1(1,1) function using the estimating equations approach
	# Given observed/simulated data x 
		# n is defined as the length of the observation vector 
	# Order of the Log ACD1 model is given by (p,q)
		# This function will only work for orders (1,1) and (2,1)
		# p is the length of alpha and q is the length of beta
	# moments denotes the distribution from which we derive the first four central moments of epsilon
		# For exponential moments, par1=lambda
		# For weibull moments, par1=alpha and par2=beta
		# For gamma moments, par1=k (shape) and par2=theta (scale)
		# For none (no distribution specified, look at user input moments
		# by default, the moments are drawn from an exponential(1) distribution
		# par1 and par2 (optional) are parameters of the distribution of epsilon 
			# the default values are par1=1 and par2=1
	# arinitval tells function whether or not to use AR(m) fit to get initial values
		# by default, arinitval is TRUE and m is 20
		# if false, function will look at user-input initial values in initval
	# CHANGE_mod: change for different Duration Models

festeq.logacd1 <- function (x, p=1, q=1, momentsdist="exponential", usermoments, par1=1, par2=1, arinitval="TRUE", initval, m=20) {
	n <- length(x)
	if (p==1 & q==1) {
			pdim <- 1+p+q 			# number of parameters to be estimated
		# MOMENTS	
			if (momentsdist=="exponential") {
				lamda=par1
				# central moments, exponential errors
					mue=1/lamda
					vare=1/lamda**2
					skewe=2/lamda**3   
					kurte=9/lamda**4    
			} else if (momentsdist=="weibull") {
				walpha=par1
				wbeta=par2
				#raw moments, weibull errors
					fimom=wbeta*gamma(1+walpha**-1)
					smom=(wbeta**2)*gamma(1+2*walpha**-1)
					tmom=(wbeta**3)*gamma(1+3*walpha**-1)
					fomom=(wbeta**4)*gamma(1+4*walpha**-1)
				#central moments, weibull errors
					mue=fimom
					vare=smom-fimom**2
					skewe=tmom-3*vare*fimom-fimom**3
					kurte=fomom-4*skewe*fimom-6*vare*fimom**2-fimom**4
			} else if (momentsdist=="gamma") {
				k=par1
				gtheta=par2
				#raw moments, gamma errors
					fimom=gtheta*gamma(1+k)/gamma(k)
					smom=(gtheta**2)*gamma(2+k)/gamma(k)
					tmom=(gtheta**3)*gamma(3+k)/gamma(k)
					fomom=(gtheta**4)*gamma(4+k)/gamma(k)
				#central moments, gamma errors
					mue=fimom
					vare=smom-fimom**2
					skewe=tmom-3*vare*fimom-fimom**3
					kurte=fomom-4*skewe*fimom-6*vare*fimom**2-fimom**4
			} else if (momentsdist=="none") {
				# check to see if user-input moments (usermoments) are valid
				if (is.numeric(usermoments)==TRUE & is.vector(usermoments)==TRUE & length(usermoments)==4) {
				# user-input central moments
					mue=usermoments[1]
					vare=usermoments[2]
					skewe=usermoments[3]
					kurte=usermoments[4]
				} else {cat("Error: Unrecognized moments specification", "\n")}
			} else { cat("Error: Unrecognized moments specification", "\n") }
		# INITIALIZATION
			# identity matrix
				iden = diag(pdim)	
			# moments of x
				mu <- rep(1,n) 	   # mu(i)
				sigsq <- rep(1,n)  # sigsq(i)
				gamma <- rep(1,n)  # third central moment of x (not skewness)
				kappa <- rep(1,n)  # fourth central moment of x (not kurtosis)
			# psi hat
				psih <- rep(1,n)   
			# k matrix (variance-covariance) and k inverse (observed information)
				kmat = array(NA, dim = c(pdim, pdim, n)) 
				kinv = array(NA, dim = c(pdim, pdim, n))
			# parameter estimates for each iteration
				thehat = array(NA, dim = c(pdim, 1, n))
			# derivative of psi and second derivative of psi
				derpsi<-matrix(rep(0),pdim,1)
				der2psi<-matrix(rep(0),pdim,pdim)
			# derivative of mu, sigsq; second derivates of mu, sigsq
				dermu<-matrix(rep(0),pdim,1)
				dersigsq<-matrix(rep(0),pdim,1)
				der2mu<-matrix(rep(0),pdim,pdim)
				der2sigsq<-matrix(rep(0),pdim,pdim)
			# derivative of m, M, quadratic variation, quadratic covariation, eta
				derm<-matrix(rep(0),pdim,1)
				derqm<-matrix(rep(0),pdim,1)
				dervm<-matrix(rep(0),pdim,1)
				dervqm<-matrix(rep(0),pdim,1)
				dereta<-matrix(rep(0),pdim,1)
			# optimal a and b
				astr<-matrix(rep(0),pdim,1)
				bstr<-matrix(rep(0),pdim,1)
		# INITIAL VALUES
			# CHANGE_mod
			if (arinitval=="TRUE") {
				initial<-finitval.logacd1(x, 1, 1)
				initial<- as.numeric(c(initial[1], initial[2], initial[3]))
				cat("The initial values from an AR(m) fit are:", initial, "\n")
					# omega, alpha, beta initial values
			} else if (arinitval=="FALSE") {
				# check to see if user-input initial values (initval) are valid 
				if (is.numeric(initval)==TRUE & length(initval)==pdim) {
					initial <- c(initval[1], initval[2], initval[3])
					cat("The user-input initial values are:",initial, "\n")
				} else { cat("Error: invalid initial values input", "\n")}
			}
		# ESTIMATING EQUATIONS
			# Put initial values into initial positions of arrays
				thehat[,,1]=initial
				# Initial observed information and var-cov matrices
				# CHANGE_mod
					kinv[,,1]=diag(c(1/(initial[1]/2)**2,1/(initial[2]/2)**2,1/(initial[3]/2)**2))	
					kmat[,,1]=solve(kinv[,,1])
			# Recursive Estimation
				# t=1
					psih[1]=thehat[1,1,1]    # omega 		#CHANGE_mod
				# t=2:n
					for (t in 2:n){
						#CHANGE_mod (change psi, derivative of psi and second derivative of psi)
							# Define psi
								psih[t]=thehat[1,1,t-1]+thehat[2,1,t-1]*log(x[t-1])+thehat[3,1,t-1]*psih[t-1]
							# Define derivatives of psih(t) wrt theta: pdim*1 vector
								# First derivative
									derpsi=matrix(c(1,log(x[t-1]),psih[t-1]),pdim,1)   
								# Second derivative: 
									der2psi=matrix(rep(0),pdim,pdim)    
						# mu_t, sigsq_t,gamma_t, kappa_t
							mu=mue*exp(psih[t]) 
							sigsq=vare*exp(2*psih[t]) 
							gamma=skewe*exp(3*psih[t])  # recall skewe is third central moment
							kappa=kurte*exp(4*psih[t])  # recall kurte is fourth central moment
						# Compute m(t) and M(t)
							m=x[t]-mu
							qm=m**2-sigsq
						# Compute Quadratic variations and covariance
							vm=sigsq*exp(2*psih[t]) 
							vqm=(kurte-vare**2)*exp(4*psih[t]) 
							vmqm=skewe*exp(3*psih[t])
						# Define rho^2(t) and eta(t)
							#rho = vare*(kurte-vare**2)/(vare*(kurte-vare**2)-skewe**2)
								termr= 1-(vmqm**2/(vm*vqm))
								rho=1/termr
							#eta = skewe/(vare*(kurte-vare**2)*exp(3*psih[t]))
								eta=vmqm/(vm*vqm)
						# First Derivatives of mu(t) and sigsq(t)
							dermu=mue*derpsi*exp(psih[t])
							dersigsq=2*vare*derpsi*exp(2*psih[t])
						# Second Derivatives of mu(t) and sigsq(t)
							der2mu=mue*exp(psih[t])*(der2psi+derpsi%*%t(derpsi))
							der2sigsq=2*vare*exp(2*psih[t])*(der2psi+2*derpsi%*%t(derpsi))
						# Define vectors astr and bstr
							astr=rho*(-dermu/vm +dersigsq*eta)
							bstr=rho*(dermu*eta - dersigsq/vqm)
						# Define Derivatives of m(i) and M(i) 
							derm=-mue*exp(psih[t])*derpsi
							derqm=2*m*derm - dersigsq
						# Derivatives of <m>(i) and <M>(i)
							dervm=2*vare*exp(2*psih[t])*derpsi
							dervqm=4*(kurte-vare**2)*exp(4*psih[t])*derpsi
						# Derivative of eta(i)
							dereta=-3*skewe*derpsi/(vare*(kurte-vare**2)*exp(3*psih[t]))
						# Derivatives of astr and bstr
							# astr
							terma1=(vm*der2mu -dermu%*%t(derm))/vm**2
							terma2=der2sigsq*eta+dersigsq%*%t(dereta)
							derastr=-rho*terma1 + rho*terma2
							# bstr
							termb1=der2mu*eta + dermu%*%t(dermu)
							termb2=(vqm*der2sigsq - dersigsq%*%t(derqm))/vqm**2
							derbstr=rho*termb1-rho*termb2
						# Compute Kinv(t) 
							termk1=astr%*%t(derm) + m*derastr 
							termk2=bstr%*%t(derqm)+ qm*derbstr
							kinv[,,t] = kinv[,,t-1] - (termk1+termk2)
						# Invert to get K(t)    
							kmat[,,t]=solve(kinv[,,t])
						# compute thehat[t]
							termt=astr*m + bstr*qm
							thehat[,,t]=thehat[,,t-1]+kmat[,,t]%*%termt
						}
			# ESTIMATES
				# print(thehat[,,1:n])
				cat("The parameter estimates from the EE method are:", thehat[,,n], "\n")
				finalest <- thehat[,,n]
	} else if (p==2 & q==1) {
				pdim <- 1+p+q 			# number of parameters to be estimated
			# MOMENTS	
				if (momentsdist=="exponential") {
					lamda=par1
					# central moments, exponential errors
						mue=1/lamda
						vare=1/lamda**2
						skewe=2/lamda**3   
						kurte=9/lamda**4    
				} else if (momentsdist=="weibull") {
					walpha=par1
					wbeta=par2
					#raw moments, weibull errors
						fimom=wbeta*gamma(1+walpha**-1)
						smom=(wbeta**2)*gamma(1+2*walpha**-1)
						tmom=(wbeta**3)*gamma(1+3*walpha**-1)
						fomom=(wbeta**4)*gamma(1+4*walpha**-1)
					#central moments, weibull errors
						mue=fimom
						vare=smom-fimom**2
						skewe=tmom-3*vare*fimom-fimom**3
						kurte=fomom-4*skewe*fimom-6*vare*fimom**2-fimom**4
				} else if (momentsdist=="gamma") {
					k=par1
					gtheta=par2
					#raw moments, gamma errors
						fimom=gtheta*gamma(1+k)/gamma(k)
						smom=(gtheta**2)*gamma(2+k)/gamma(k)
						tmom=(gtheta**3)*gamma(3+k)/gamma(k)
						fomom=(gtheta**4)*gamma(4+k)/gamma(k)
					#central moments, gamma errors
						mue=fimom
						vare=smom-fimom**2
						skewe=tmom-3*vare*fimom-fimom**3
						kurte=fomom-4*skewe*fimom-6*vare*fimom**2-fimom**4
				} else if (momentsdist=="none") {
					# check to see if user-input moments (usermoments) are valid
					if (is.numeric(usermoments)==TRUE & is.vector(usermoments)==TRUE & length(usermoments)==4) {
					# user-input central moments
						mue=usermoments[1]
						vare=usermoments[2]
						skewe=usermoments[3]
						kurte=usermoments[4]
					} else {cat("Error: Unrecognized moments specification", "\n")}
				} else { cat("Error: Unrecognized moments specification", "\n") }
			# INITIALIZATION
				# identity matrix
					iden = diag(pdim)	
				# moments of x
					mu <- rep(1,n) 	   # mu(i)
					sigsq <- rep(1,n)  # sigsq(i)
					gamma <- rep(1,n)  # third central moment of x (not skewness)
					kappa <- rep(1,n)  # fourth central moment of x (not kurtosis)
				# psi hat
					psih <- rep(1,n)   
				# k matrix (variance-covariance) and k inverse (observed information)
					kmat = array(NA, dim = c(pdim, pdim, n)) 
					kinv = array(NA, dim = c(pdim, pdim, n))
				# parameter estimates for each iteration
					thehat = array(NA, dim = c(pdim, 1, n))
				# derivative of psi and second derivative of psi
					derpsi<-matrix(rep(0),pdim,1)
					der2psi<-matrix(rep(0),pdim,pdim)
				# derivative of mu, sigsq; second derivates of mu, sigsq
					dermu<-matrix(rep(0),pdim,1)
					dersigsq<-matrix(rep(0),pdim,1)
					der2mu<-matrix(rep(0),pdim,pdim)
					der2sigsq<-matrix(rep(0),pdim,pdim)
				# derivative of m, M, quadratic variation, quadratic covariation, eta
					derm<-matrix(rep(0),pdim,1)
					derqm<-matrix(rep(0),pdim,1)
					dervm<-matrix(rep(0),pdim,1)
					dervqm<-matrix(rep(0),pdim,1)
					dereta<-matrix(rep(0),pdim,1)
				# optimal a and b
					astr<-matrix(rep(0),pdim,1)
					bstr<-matrix(rep(0),pdim,1)
			# INITIAL VALUES
			# CHANGE_mod
				if (arinitval=="TRUE") {
					initial<-finitval.logacd1(x, 2, 1)
					initial<- as.numeric(c(initial[1], initial[2], initial[3], initial[4]))
					cat("The initial values from an AR(m) fit are:", initial, "\n")
						# omega, alpha1, alpha2, beta initial values
				} else if (arinitval=="FALSE") {
					# check to see if user-input initial values (initval) are valid 
					if (is.numeric(initval)==TRUE & length(initval)==pdim) {
						initial <- c(initval[1], initval[2], initval[3], initval[4])
						cat("The user-input initial values are:",initial, "\n")
					} else { cat("Error: invalid initial values input", "\n")}
				}
			# ESTIMATING EQUATIONS
				# Put initial values into initial positions of arrays
					thehat[,,1]=initial
					thehat[,,2]=thehat[,,1]
					# Initial observed information and var-cov matrices
					# CHANGE_mod
						kinv[,,1]=diag(c(1/(initial[1]/2)**2,1/(initial[2]/2)**2,1/(initial[3]/2)**2,1/(initial[4]/2)**2))
						kmat[,,1]=solve(kinv[,,1])
						kinv[,,2]=kinv[,,1]
						kmat[,,2]=kmat[,,1]
				# Recursive Estimation
					# t=1 and t=2
					#CHANGE_mod
						psih[1]=thehat[1,1,1]    # omega 		
						psih[2]=thehat[1,1,1]+thehat[2,1,1]*log(x[1])+thehat[4,1,1]*psih[1]		# omega + alpha1*logx1 + beta*psih1
					# t=3:n
						for (t in 3:n){
							#CHANGE_mod (change psi, derivative of psi and second derivative of psi)
								# Define psi
									psih[t]=thehat[1,1,t-1]+thehat[2,1,t-1]*log(x[t-1])+thehat[3,1,t-1]*log(x[t-2])+thehat[4,1,t-1]*psih[t-1]
								# Define derivatives of psih(t) wrt theta: pdim*1 vector
									# First derivative
										derpsi=matrix(c(1,log(x[t-1]),log(x[t-2]),psih[t-1]),pdim,1)   
									# Second derivative: 
										der2psi=matrix(rep(0),pdim,pdim)    
							# mu_t, sigsq_t,gamma_t, kappa_t
								mu=mue*exp(psih[t]) 
								sigsq=vare*exp(2*psih[t]) 
								gamma=skewe*exp(3*psih[t])  # recall skewe is third central moment
								kappa=kurte*exp(4*psih[t])  # recall kurte is fourth central moment
							# Compute m(t) and M(t)
								m=x[t]-mu
								qm=m**2-sigsq
							# Compute Quadratic variations and covariance
								vm=sigsq*exp(2*psih[t]) 
								vqm=(kurte-vare**2)*exp(4*psih[t]) 
								vmqm=skewe*exp(3*psih[t])
							# Define rho^2(t) and eta(t)
								#rho = vare*(kurte-vare**2)/(vare*(kurte-vare**2)-skewe**2)
									termr= 1-(vmqm**2/(vm*vqm))
									rho=1/termr
								#eta = skewe/(vare*(kurte-vare**2)*exp(3*psih[t]))
									eta=vmqm/(vm*vqm)
							# First Derivatives of mu(t) and sigsq(t)
								dermu=mue*derpsi*exp(psih[t])
								dersigsq=2*vare*derpsi*exp(2*psih[t])
							# Second Derivatives of mu(t) and sigsq(t)
								der2mu=mue*exp(psih[t])*(der2psi+derpsi%*%t(derpsi))
								der2sigsq=2*vare*exp(2*psih[t])*(der2psi+2*derpsi%*%t(derpsi))
							# Define vectors astr and bstr
								astr=rho*(-dermu/vm +dersigsq*eta)
								bstr=rho*(dermu*eta - dersigsq/vqm)
							# Define Derivatives of m(i) and M(i) 
								derm=-mue*exp(psih[t])*derpsi
								derqm=2*m*derm - dersigsq
							# Derivatives of <m>(i) and <M>(i)
								dervm=2*vare*exp(2*psih[t])*derpsi
								dervqm=4*(kurte-vare**2)*exp(4*psih[t])*derpsi
							# Derivative of eta(i)
								dereta=-3*skewe*derpsi/(vare*(kurte-vare**2)*exp(3*psih[t]))
							# Derivatives of astr and bstr
								# astr
								terma1=(vm*der2mu -dermu%*%t(derm))/vm**2
								terma2=der2sigsq*eta+dersigsq%*%t(dereta)
								derastr=-rho*terma1 + rho*terma2
								# bstr
								termb1=der2mu*eta + dermu%*%t(dermu)
								termb2=(vqm*der2sigsq - dersigsq%*%t(derqm))/vqm**2
								derbstr=rho*termb1-rho*termb2
							# Compute Kinv(t) 
								termk1=astr%*%t(derm) + m*derastr 
								termk2=bstr%*%t(derqm)+ qm*derbstr
								kinv[,,t] = kinv[,,t-1] - (termk1+termk2)
							# Invert to get K(t)    
								kmat[,,t]=solve(kinv[,,t])
							# compute thehat[t]
								termt=astr*m + bstr*qm
								thehat[,,t]=thehat[,,t-1]+kmat[,,t]%*%termt
							}
				# ESTIMATES
					# print(thehat[,,1:n])
					cat("The parameter estimates from the EE method are:", thehat[,,n], "\n")
					finalest <- thehat[,,n]
	} else { cat("Error: This function is only defined for the log ACD1 (1,1) and log ACD1 (2,1) cases") }				
}

# EXAMPLES:
# test <- fsim.logacd1(100, 1000, 1, .5, .3, feps="exponential", 1); x <- test$x
# ee <- festeq.logacd1(x, 1,1, momentsdist="exponential", arinitval="TRUE"); ee
# ee <- festeq.logacd1(x, 2,1, momentsdist="exponential", arinitval="TRUE"); ee
# festeq.logacd1(durations, 1,1, momentsdist="exponential", arinitval="TRUE")
	# The initial values from an AR(m) fit are: 0.6171753 0.1169692 0.5546814 
	# The parameter estimates from the EE method are: 0.607835 0.1172964 0.5670814 
# festeq.logacd1(durations[1:(length(durations)-4)], 1,1, momentsdist="exponential", arinitval="TRUE")
	# The initial values from an AR(m) fit are: 0.6217596 0.1172356 0.5520966 
	# The parameter estimates from the EE method are: 0.6133423 0.1175344 0.5637614 

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
