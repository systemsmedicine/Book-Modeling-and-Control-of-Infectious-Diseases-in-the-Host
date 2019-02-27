# R Codes for Chapter 3
# This code was developed by Kinh Nguyen
# You can also find at https://figshare.com/s/9f0c50984e470693839e

# Materials
# -------------------------------------------------------------------------
install.packages("deSolve"); library("deSolve")  	# see Note 2
install.packages("DEoptim"); library("DEoptim")
# Anything that follows the character ‘#’ is a comment and is not processed in R

# 3.1.	Preparing data
# -------------------------------------------------------------------------
# Read data from the downloaded CSV file
		
myData <- read.csv(“/path/to/myData.csv”)  		# see Notes 2, 5
myData  					# View data

# or input the data directly
# Write time in days
time <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 5, 5, 5, 5, 5, 7, 7, 7, 7, 7, 9, 9, 9, 9, 9)

# Write viral titres in log 10
V <- c(2.36, 3.23, 2.69, 2.57, 2.37, 4.87, 4.41, 5.04, 3.76, 4.16, 5.20, 5.43, 5.55, 5.60, 6.12, 4.96, 4.54, 4.65, 5.06, 4.80, 3.67, 3.22, 3.35, 2.84, 3.76, 1.61, 1.80, 2.00, 1.69, 2.22)
myData <- data.frame(time = time, V = V)
myData  					# View data

# 3.2.	Writing the model and the cost function
# -------------------------------------------------------------------------
myModel <- function(t,state,parameters) { 
  with(as.list(c(state,parameters)),{
	  # Fill in the equations			# see Note 6
    dU = - myBeta*U*V  				 
    dI = myBeta*U*V - myDelta*I 
    dV = myP*I - myC*V
    # Outputs of the model			# see Note 7
    list(c(dU, dI, dV))  			
  })
}

myStates <- c(U = 10^6, I = 0, V = 10)  		# see Note 8
myParams <- c("myBeta","myDelta","myP","myC") 	# list of parameters
myV      <- "V"  				# the component in model observed in data

modelTime  <- seq(from = 0, to = 10, by = 0.01) 	# see Note 9

myCostFn <- function(x) {
    parms        <- x[1:length(myParams)]
    names(parms) <- myParams
    yhat      <- ode(myStates, modelTime, myModel, 10^parms)
    yMatch    <- yhat[yhat[,1] %in% myData$time, ]  	# see Note 9
    nm        <- rle(myData$time)$lengths
    x         <- myData[, myV]- rep(log10(yMatch[, myV]), times = nm)
    rmse      <- sqrt(mean(x^2))
    return(rmse)
}  					# see Note 15

lower = log10(c(1e-7, 1e-2, 1e+0, 1e-1))  		# see Notes 10, 11
upper = log10(c(1e-3, 1e+2, 1e+2, 1e+2))

myOptions <- DEoptim.control(itermax = 10000, steptol = 50, reltol = 1e-8)  # see Note 12

fit <- do.call("DEoptim", list(myCostFn, lower, upper, myOptions))  # see Notes 13, 14, 15

(bestPar <- fit$optim$bestmem)  			# show parameter estimates (in log 10 scale)
names(bestPar) <- myParams			# assign names to the parameters
out <- ode(myStates, modelTime, myModel, 10^bestPar)	# solve the ODEs with the estimated values
plot(out[, "time"], log10(out[, "V"]), type = "l")		# plot fitted viral load in log 10
points(myData$time, myData$V) 			# superimpose observed data points

# 3.3.	Model comparison with Akaike information criteria (AIC).
# -------------------------------------------------------------------------
myAIC <- function(fit, np=NULL, rms=NULL, n=NULL) {
    if (is.null(n)) stop ("How many observations were used? n=#")
    if (is.null(np)) np <- length(fit$optim$bestmem)
    if (is.null(rms)) rms <- fit$optim$bestval
    return(2*np + n*log(rms))
}

myAIC(fit, n = 30)

# 3.4.	Likelihood profile of the model parameters
# -------------------------------------------------------------------------
myProfile <- function(lower, upper, bestPar) {
	pro.ll <-  NULL
	for (v in 1:length(bestPar)) {
	    # Create parameter sequence
	    tmpl <- seq(lower[v], bestPar[[v]], length.out = 100)
	    tmpl <- tmpl[order(tmpl, decreasing = TRUE)[cumsum(1:13)]]
	    tmpr <- seq(bestPar[[v]], upper[v], length.out = 100)
	    tmpr <- tmpr[cumsum(1:13)]
	    pars <- sort(unique(c(lower[v], tmpl, bestPar[[v]], tmpr, upper[v])))
	    ppl  <- NULL
	    # Run optimization for each, and record the parameters and RMSE
	    for (p in pars) {
	        DEargs  <- list(myCostFn, replace(lower, v, p), replace(upper, v, p), myOptions) 
	        fit     <- do.call("DEoptim", DEargs)
	        ppl     <- c(ppl, fit$optim$bestval)
	    }
	    pro.ll[[v]] <- cbind(pars, ppl)
	}
  return(pro.ll)
}

outProfiles <- myProfile(lower, upper, bestPar)

par(mfrow = c(2, 2))
sapply(1:4, function(x) plot(outProfiles[[x]], xlab = myParams[x], ylab = 'RMSE'))

# 3.5.	Bootstrapping parameters
# -------------------------------------------------------------------------
myBoot <- function(numboot = 1000, numpar = 4) {
	# numpar: number of parameters in the model
	# numboot: number of bootstrap samples
	results <- matrix(NA, numboot, numpar)
	original <- myData
	sampling <- function(x) sample(original$V[original$time==x], length(original$V[original$time==x]), replace = 1)
	for (i in 1:numboot) {
		message("Bootstrapping sample ", i)
		tmp    <- sapply(unique(original$time), sampling)
		myData <- cbind(original$time, as.vector(tmp))
		DEarguments <- list(myCostFn, lower, upper, myOptions)
		fit <- do.call("DEoptim", DEarguments)
		results[i, ] <- fit$optim$bestmem
	}
	results <- as.data.frame(results)
	colnames(results) <- myParams
	return(results)
}

bootResults <- myBoot()  			# see Note 17
par(mfrow = c(2,2))
sapply(1:4, function(x) hist(bootResults[, x]) )
apply(bootResults, 2, quantile, probs = c(.025,.975))


