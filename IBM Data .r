########################################################################
##### IBM DATA PROCESSING 
########################################################################
# READ IN RAW IBM DATA (collected 9/13/2012)
	rawIBM <- read.table("C://Users//Lilian//Documents//UCONN Documents//2014 Spring//Stat Research 2014//GLARMA//IBM data//IBM raw data.txt",header=T)
	# Change the class of $time from a factor variable to date/time (note: date will be the current date)
		rawIBM$time <- strptime(rawIBM$time, format="%H:%M:%S")
	# rawIBM is in reverse chronological order by time. So invert it.
		rawIBM <- rawIBM[order(nrow(rawIBM):1),]
		rownames(rawIBM) <- as.character(seq(1, length(rawIBM$time),1))
	
# EXTRACT DURATIONS
	# An event occurs if the price difference between two consecutive trades exceeds some critical value
		crit <- 0.0125
	# keep only the time and price columns
		IBMdata <- rawIBM[,1:2]
		attach(IBMdata)
	# INITIALIZATION
	  # initialize a new variable event={1 if an event occurred, 0 otherwise} that is indexed by i
		event <- rep(NA, nrow(IBMdata))	
	  # initialize t, which will store the position at which an event occurs
		t <- 1
	  # initialize i, which will start at 2 since that is the time at which the first event can occur	
		i <- 2	
	# FOR LOOP
		for (i in 2:nrow(IBMdata)) {
		  # if an event occurs (i.e., |price at time i - price at previous event time| >= critical value)
			if (abs(price[i]-price[t])>=crit) { 
				# set event[i]=1
					event[i] <- 1
				# store the position at which the event occurred as t (so it can be compared to the next observed price in the next loop)	
					t <- i 
		  # if an event does not occur	
			} else { 
				# set event[i]=0
					event[i] <- 0 
				}
		}
	# bind event with IBMdata 
		IBMdata <- cbind(IBMdata, event)
	# define arrival times, which are the times associated with an event occurring (times such that event=1)
		arr.times <- time[event==1]
	# define durations, which are the non-zero times between events
		durations <- as.numeric(diff(arr.times)[diff(arr.times)!=0]	)
		durations <- durations[is.na(durations)==FALSE]
		
detach(IBMdata)


########################################################################
##### MODEL 1
########################################################################

x <- durations
L <- 5  					# observations withheld 
n <- length(durations)-L	# observations used to fit model

# EE parameter estimation
# param <- festeq.logacd1(durations[1:n], 1,1, momentsdist="weibull", 9, 1, arinitval="TRUE")
	# The initial values from an AR(m) fit are: 0.6171753 0.1169692 0.5546814 
	# The parameter estimates from the EE method are: 0.607835 0.1172964 0.5670814  

# MAPE	
#1
fmape.logacd1(durations, length(durations)-L, L, 2, 1, momentsdist="gamma", .5, 8)
#2
fmape.logacd1(durations, length(durations)-L, L, 1, 1, momentsdist="gamma", .5, 8)

#####################
# Graphs
#####################

# param <- festeq.logacd1(durations[1:n], 1,1, momentsdist="gamma", 0.5, 8, arinitval="TRUE")
# omega <- param[1]
# alpha <- param[2]
# beta <- param[3]

# # define fit and prediction parts of x
   # fit <- x[1:n]
   # pred.actual <- x[(n+1):(n+L)]
  # # psi values associated with fit part of x (known values of x)
   # psi.fit <- log(x[1:n]) 
   
   
# # Prediction (j=1:L steps from forecast origin n) given n values of x
		# pred.psi <- rep(NA, L)
		# pred.xhat <- rep(NA, L)
		# # For j=1
			# pred.psi[1] <- omega + alpha*log(x[(n+1)-1]) + beta*psi.fit[(n+1)-1]
			# pred.xhat[1] <- exp(pred.psi[1])
		# # For j=2:L steps from n	
			# for (j in 2:L) {
				# pred.psi[j] <- omega + alpha*log(x[(n+j)-1]) + beta*pred.psi[j-1]
				# pred.xhat[j] <- exp(pred.psi[j])
			# }
	  # # MAPE Calculation			
		# mape <- 100*(mean(abs((pred.actual-pred.xhat)/pred.actual)))
		# mape
			
# # Fit
	# # Initialize psi.fit and x.fit
	# psi.fit <- rep(1,length(x))
	# x.fit <-rep(1,length(x))

	# # For t=1
		# psi.fit[1] <- log(x[1])
		# x.fit[1] <- exp(psi.fit[1])
	# # For t=2:n
		# for (t in 2:length(x)) {
				# psi.fit[t] = omega + alpha*log(x[t-1]) + beta*psi.fit[t-1] 
				# x.fit[t] = exp(psi.fit[t] )
		# }

# # Graphs
# windows()
# plot(x, type="l", main="Fitted Log ACD model against IBM durations")
# lines(x.fit, type="l", col="gold")

# windows()
# par(mfrow=c(1,2))
# hist(x[1:n]/exp(psi.fit[1:n]), freq=FALSE, xlab="residuals", main="Histogram of Log ACD residuals")
# acf(x[1:n]/exp(psi.fit[1:n]), lag.max=100, main="ACF of Log ACD residuals")

param <- festeq.logacd1(durations[1:n], 2,1, momentsdist="gamma", 0.5, 8, arinitval="TRUE")
omega <- param[1]
alpha1 <- param[2]
alpha2 <- param[3]
beta <- param[4]

# define fit and prediction parts of x
   fit <- x[1:n]
   pred.actual <- x[(n+1):(n+L)]
  # psi values associated with fit part of x (known values of x)
   psi.fit <- log(x[1:n]) 
 
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

# Fit
	# Initialize psi.fit and x.fit
	psi.fit <- rep(1,length(x))
	x.fit <-rep(1,length(x))

	# For t=1
		psi.fit[1] <- log(x[1])
		x.fit[1] <- exp(psi.fit[1])
	# For t=2
		psi.fit[2] <- omega + alpha1*log(x[1]) + beta*psi.fit[1]
		x.fit[2] <- exp(psi.fit[2])
	# For t>2	
		for (t in 3:length(x)) {
				psi.fit[t] = omega + alpha1*log(x[t-1]) + alpha2*log(x[t-2]) + beta*psi.fit[t-1] 
				x.fit[t] = exp(psi.fit[t] )
		}

# Graphs
windows()
par(mfrow=c(2,2))
plot(x, type="l", main="Fitted Log ACD model against IBM durations")
lines(x.fit, type="l", col="gold")
acf(x, lag.max=100, main= "ACF of IBM durations")

hist(x[1:n]/exp(psi.fit[1:n]), freq=FALSE, xlab="residuals", main="Histogram of Log ACD residuals")
acf(x[1:n]/exp(psi.fit[1:n]), lag.max=100, main="ACF of Log ACD residuals")


