## Run after Functions Definition.r ##

#############################################################################
##### (1) MONTE CARLO SIMULATION - LOG ACD1(1,1) , nsim=100
#############################################################################
# Prior to using this code: Load custom functions

# number of simulations
nsim = 100
# initialize simulation number
isim = 1 
# initialize omgi, alphi, and beti to store initial values across nsim simulations:
omgi = rep(NA, nsim)
alpi = rep(NA, nsim)
beti = rep(NA, nsim)
# initialize omgh, alph, and beth to store parameter estimates across nsim simulations:
omgh = rep(NA, nsim)
alph = rep(NA, nsim)
beth = rep(NA, nsim)

# simulation setup
n <- 4000			# length of sample after burn-in
nb <- 1000			# burn-in
# parameters
# setup 1
# omg <- 1
# alp <- 0.5
# bet <- 0.3	
# setup 2 
# omg <- 3
# alp <- 0.2
# bet <- -0.4	
# setup 3
# omg <- 2
# alp <- -0.1
# bet <- .23	
# setup 4 
 omg <- -1.5
 alp <- -0.2
 bet <- 0.65
# distribution
# setups 1 and 2
# distr <- "exponential"
# par1 = 1
# par2 = 1
# setups 3 and 4
 distr <- "gamma"
 par1 = 2
 par2 = 0.5

###################################
## ARMA initial parameter estimates
###################################

set.seed(123457)   # random seed for generation

for (isim in 1:nsim) {
# simulation
sim <- fsim.logacd1(n, nb, omg, alp, bet, feps=distr, par1, par2)
x <- sim$x
# initial values
init <- finitval.logacd1(x, p=1, q=1)
# store initial parameter estimates into omgi, alphi, beti
omgi[isim] <- init[1]
alpi[isim] <- init[2]
beti[isim] <- init[3]
omgi <- as.numeric(omgi); alpi <- as.numeric(alpi); beti <- as.numeric(beti)
}
quantile(omgi,probs=c(0.05,.25,.5,.75,.95))
quantile(alpi,probs=c(0.05,.25,.5,.75,.95))
quantile(beti,probs=c(0.05,.25,.5,.75,.95))

# sample graphs
# ts.plot(x, main="Data generated from the Log ACD1(1,1) model with \n theta=(1,0.5,0.3) and exponential(1) errors")

###################################
## EE with ARMA initial parameter estimates
###################################

set.seed(123457)   # random seed for generation

for (isim in 1:nsim) {
# simulation
sim <- fsim.logacd1(n, nb, omg, alp, bet, feps=distr, par1, par2)
x <- sim$x
# estimation with AR initial values
try(EE <- festeq.logacd1(x, 1,1, momentsdist=distr, par1, par2, arinitval="TRUE", m=10))
# store estimated parameters into omgh, alph, beth
omgh[isim] <- EE[1]
alph[isim] <- EE[2]
beth[isim] <- EE[3]
}

quantile(omgh,probs=c(0.05,.25,.5,.75,.95), na.rm=TRUE)
quantile(alph,probs=c(0.05,.25,.5,.75,.95), na.rm=TRUE)
quantile(beth,probs=c(0.05,.25,.5,.75,.95), na.rm=TRUE)

# Save these values as omg1, alph1, bet1 since we'll be moving on to more estimation
omg1 <- omgh
alp1 <- alph
bet1 <- beth

###################################
## EE without ARMA initial parameter estimates
###################################

set.seed(123457)   # random seed for generation

# re-initialize omgh, alph, and beth to store parameter estimates across nsim simulations:
omgh = rep(NA, nsim)
alph = rep(NA, nsim)
beth = rep(NA, nsim)

for (isim in 1:nsim) {
# simulation
sim <- fsim.logacd1(n, nb, omg, alp, bet, feps=distr, par1, par2)
x <- sim$x
# estimation without AR initial values
initval1 <- runif(1, -5, 5)
initval2 <- runif(1, -1, 1)
initval3 <- runif(1, -1+abs(initval2), 1-abs(initval2))
try(EE <- festeq.logacd1(x, 1,1, momentsdist=distr, par1, par2, arinitval="FALSE", initval=c(initval1, initval2, initval3)))
# store estimated parameters into omgh, alph, beth
omgh[isim] <- EE[1]
alph[isim] <- EE[2]
beth[isim] <- EE[3]
}

quantile(omgh,probs=c(0.05,.25,.5,.75,.95), na.rm=TRUE)
quantile(alph,probs=c(0.05,.25,.5,.75,.95), na.rm=TRUE)
quantile(beth,probs=c(0.05,.25,.5,.75,.95), na.rm=TRUE)

omg2 <- omgh
alp2 <- alph
bet2 <- beth

### GRAPHS ###
windows()
par(mfcol=c(3,2))
hist(omg1, 12, xlim=c(min(omg1[is.na(omg1)==FALSE])-0.1, max(omg1[is.na(omg1)==FALSE])+0.1), main="Results using ARMA initial parameter method", xlab="Final estimates of omega")
hist(alp1, 12, xlim=c(min(alp1[is.na(alp1)==FALSE])-0.1, max(alp1[is.na(alp1)==FALSE])+0.1),  main="Results using ARMA initial parameter method", xlab="Final estimates of alpha")
hist(bet1,12, xlim=c(min(bet1[is.na(bet1)==FALSE])-0.1, max(bet1[is.na(bet1)==FALSE])+0.1),  main="Results using ARMA initial parameter method", xlab="Final estimates of beta")
hist(omg2, 12, xlim=c(min(omg2[is.na(omg2)==FALSE], omg1[is.na(omg1)==FALSE])-0.25, max(omg1[is.na(omg1)==FALSE], omg2[is.na(omg2)==FALSE])+0.25),  main="Results using uniform generation method", xlab="Final estimates of omega")
hist(alp2, 12, xlim=c(min(alp2[is.na(alp2)==FALSE], alp1[is.na(alp1)==FALSE])-0.25, max(alp1[is.na(alp1)==FALSE], alp2[is.na(alp2)==FALSE])+0.25),  main="Results using uniform generation method", xlab="Final estimates of alpha")
hist(bet2,12, xlim=c(min(bet2[is.na(bet2)==FALSE], bet1[is.na(bet1)==FALSE])-0.25, max(bet1[is.na(bet1)==FALSE], bet2[is.na(bet2)==FALSE])+0.25),  main="Results using uniform generation method", xlab="Final estimates of beta")




















#############################################################################
##### (2) MONTE CARLO SIMULATION - LOG ACD1(2,1) , nsim=100
#############################################################################


# number of simulations
	nsim = 100
# initialize simulation number
	isim = 1 
# initialize omgi, alphi, and beti to store initial values across nsim simulations:
	omgi = rep(NA, nsim)
	alpi1 = rep(NA, nsim)
	alpi2 = rep(NA, nsim)
	beti = rep(NA, nsim)
# initialize omgh, alph, and beth to store parameter estimates across nsim simulations:
	omgh = rep(NA, nsim)
	alph1 = rep(NA, nsim)
	alph2 = rep(NA, nsim)
	beth = rep(NA, nsim)

# simulation setup
	n <- 4000			# length of sample after burn-in
	nb <- 1000			# burn-in
	# parameters
		# setup 1
			omg <- 10
			alp1 <- 0.1
			alp2 <- -0.5
			bet <- 0.06	
		# setup 2
			# omg <- 5
			# alp1 <- -0.11
			# alp2 <- -0.6
			# bet <- -0.2	
	# distribution
		# setup 1
			 distr <- "exponential"
			 par1 = 1
			 par2 = 1
		# setup 2
			# distr <- "gamma"
			# par1 = 2
			# par2 = 0.5

###################################
## ARMA initial parameter estimates and EE using user-input initial parameter estimates from the AR(m) estimation
###################################

set.seed(123457)   # random seed for generation
		
for (isim in 1:nsim) {
	# simulation
		sim <- fsim.logacd1(n, nb, omg, c(alp1, alp2), bet, feps=distr, par1, par2)
		x <- sim$x
	# AR initial values
		try(init <- finitval.logacd1(x, p=2, q=1))
		# store initial parameter estimates into omgi, alphi, beti
			omgi[isim] <- init[1]
			alpi1[isim] <- init[2]
			alpi2[isim] <- init[3]
			beti[isim] <- init[4]
	# estimation with user-input AR initial values
		try(EE <- festeq.logacd1(x, p=2,q=1, momentsdist=distr, par1, par2, arinitval="FALSE", initval=c(as.numeric(omgi[isim]), as.numeric(alpi1[isim]), as.numeric(alpi2[isim]), as.numeric(beti[isim]))))
	# store EE estimated parameters into omgh, alph, beth
		omgh[isim] <- EE[1]
		alph1[isim] <- EE[2]
		alph2[isim] <- EE[3]
		beth[isim] <- EE[4]
}

omgi <- as.numeric(omgi)
alpi1 <- as.numeric(alpi1)
alpi2 <- as.numeric(alpi2)
beti <- as.numeric(beti)

quantile(omgi,probs=c(0.05,.25,.5,.75,.95), na.rm=TRUE)
quantile(alpi1,probs=c(0.05,.25,.5,.75,.95), na.rm=TRUE)
quantile(alpi2,probs=c(0.05,.25,.5,.75,.95), na.rm=TRUE)
quantile(beti,probs=c(0.05,.25,.5,.75,.95), na.rm=TRUE)

quantile(omgh,probs=c(0.05,.25,.5,.75,.95), na.rm=TRUE)
quantile(alph1,probs=c(0.05,.25,.5,.75,.95), na.rm=TRUE)
quantile(alph2,probs=c(0.05,.25,.5,.75,.95), na.rm=TRUE)
quantile(beth,probs=c(0.05,.25,.5,.75,.95), na.rm=TRUE)

omg1 <- omgh
alp11 <- alph1
alp12 <- alph2
bet1 <- beth

###################################
## EE without ARMA initial parameter estimates
###################################

set.seed(123457)   # random seed for generation

	# re-initialize parameter storage vectors:
		omgh = rep(NA, nsim)
		alph1 = rep(NA, nsim)
		alph2 = rep(NA, nsim)
		beth = rep(NA, nsim)
		
for (isim in 1:nsim) {
	# simulation
		sim <- fsim.logacd1(n, nb, omg, c(alp1, alp2), bet, feps=distr, par1, par2)
		x <- sim$x
	# estimation without AR initial values
		## setup 1 uniform dist
			initval1 <- runif(1, 0, 12)
			initval2 <- runif(1, 0, 0.5)
			initval3 <- runif(1, -1, 0)
			initval4 <- runif(1, 0, 0.2)
		## setup 2 uniform dist
			#initval1 <- runif(1, 4, 8)
			#initval2 <- runif(1, -0.3, 0)
			#initval3 <- runif(1, -0.9, -0.5)
			#initval4 <- runif(1, -0.3, 0)
		try(EE <- festeq.logacd1(x, p=2,q=1, momentsdist=distr, par1, par2, arinitval="FALSE", initval=c(initval1, initval2, initval3, initval4)))
	# store estimated parameters into omgh, alph, beth
		omgh[isim] <- EE[1]
		alph1[isim] <- EE[2]
		alph2[isim] <- EE[3]
		beth[isim] <- EE[4]
}

quantile(omgh,probs=c(0.05,.25,.5,.75,.95), na.rm=TRUE)
quantile(alph1,probs=c(0.05,.25,.5,.75,.95), na.rm=TRUE)
quantile(alph2,probs=c(0.05,.25,.5,.75,.95), na.rm=TRUE)
quantile(beth,probs=c(0.05,.25,.5,.75,.95), na.rm=TRUE)

omg2 <- omgh
alp21 <- alph1
alp22 <- alph2
bet2 <- beth

### GRAPHS ###
windows()
par(mfcol=c(4,2))
hist(omg1, 12, xlim=c(min(omg1)-0.1, max(omg1)+0.1), main="Results using ARMA initial parameter method", xlab="Final estimates of omega")
hist(alp11, 12, xlim=c(min(alp11)-0.1, max(alp11)+0.1),  main="Results using ARMA initial parameter method", xlab="Final estimates of alpha 1")
hist(alp12, 12, xlim=c(min(alp12)-0.1, max(alp12)+0.1),  main="Results using ARMA initial parameter method", xlab="Final estimates of alpha 2")
hist(bet1,12, xlim=c(min(bet1)-0.1, max(bet1)+0.1),  main="Results using ARMA initial parameter method", xlab="Final estimates of beta")
hist(omg2, 12, xlim=c(min(omg2[is.na(omg2)==FALSE], omg1)-0.25, max(omg1[is.na(omg1)==FALSE], omg2[is.na(omg2)==FALSE])+0.25),  main="Results using uniform generation method", xlab="Final estimates of omega")
hist(alp21, 12, xlim=c(min(alp21[is.na(alp21)==FALSE], alp11)-0.25, max(alp11[is.na(alp11)==FALSE], alp21[is.na(alp21)==FALSE])+0.25),  main="Results using uniform generation method", xlab="Final estimates of alpha 1")
hist(alp22, 12, xlim=c(min(alp22[is.na(alp22)==FALSE], alp12)-0.25, max(alp12[is.na(alp12)==FALSE], alp22[is.na(alp22)==FALSE])+0.25),  main="Results using uniform generation method", xlab="Final estimates of alpha 2")
hist(bet2,12, xlim=c(min(bet2[is.na(bet2)==FALSE], bet1)-0.25, max(bet1[is.na(bet1)==FALSE], bet2[is.na(bet2)==FALSE])+0.25),  main="Results using uniform generation method", xlab="Final estimates of beta")
