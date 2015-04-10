"diffusionCoefficient" <- function (tr, byburst = FALSE, nsteps = 1000) {
	tr <- na.omit(tr)

	resultNames <- unique(id(tr))
	if (byburst) {
		resultNames <- burst(tr)
	}
	
	inputData <- list()
	for (n in resultNames) { inputData[[n]] <- matrix(0, nrow=3, ncol=0) }
	maxDiffCoeff <- rep(0, length(resultNames)) # For each element of the result, find a proper range to search
	names(maxDiffCoeff) <- resultNames

	# For each even numbered measurement in a burst (except the last if burst has even length):
	#  - compute the relevant parameters for the estimation of the diffusion coefficient:
	#      - the distance between the observation and the mean derived from the neighbouring observations
	#      - parameters used to compute the location variance from the diffusion coefficient
	#  - compute the ML value for the diffusion coefficient looking only at one bridge.
	#      When the diffusion coefficient gets larger than the maximum of these,
	#      the likelihood of the observations becomes a decreasing function of the diffusion coefficient.
	for (b in 1:length(tr)) {
		burst <- tr[[b]]

		if (nrow(burst) >= 3) { # We can't estimate the likelihood for shorter bursts
			burst$date <- as.double(burst$date) - min(as.double(burst$date))
		
			bdata <- matrix(NA, nrow=3, ncol=floor((nrow(burst)-1)/2))
		
			fac <- attr(burst, ifelse(byburst, "burst", "id"))
			# i runs over all even numbers that have a measurement before and after them
			for (i in seq(2, nrow(burst)-1, by=2)) {
				alpha <- (burst$date[i] - burst$date[i-1]) / (burst$date[i+1] - burst$date[i-1])
			
				# Squared distance from measured loc to estimated mean
				bdata[1, i/2] <- 
						(burst$x[i] - (1-alpha) * burst$x[i-1] - alpha * burst$x[i+1])^2 +
						(burst$y[i] - (1-alpha) * burst$y[i-1] - alpha * burst$y[i+1])^2
				# The coefficients for a linear function mapping diffusion coefficient to variance
				bdata[2, i/2] <- (burst$date[i+1]-burst$date[i-1]) * alpha * (1-alpha)
				bdata[3, i/2] <- (1-alpha)^2*burst$loc.var[i-1] + alpha^2*burst$loc.var[i+1]
			
				# Find the variance at which this bridge reaches its maximum likelihood
				maxVar <- 0.5 * bdata[1, i/2]
				maxDiffCoeff[fac] <- max(maxDiffCoeff[fac], (maxVar-bdata[3, i/2])/bdata[2, i/2])
			}
			
			inputData[[fac]] <- cbind(inputData[[fac]], bdata)
		}
	}
	
	result <- rep(NA, length(resultNames))
	names(result) <- resultNames
	for(fac in names(result)) {
		if (maxDiffCoeff[fac] <= 0) {
			# We already know that the ML is achieved for diff.coeff 0
			result[fac] <- 0
		} else {
			candidates <- seq(0,maxDiffCoeff[fac], length.out=nsteps)
		
			cResult <- .C("diffusion_static", double(1),
					c(inputData[[fac]]), ncol(inputData[[fac]]),
					candidates, length(candidates), PACKAGE="movementAnalysis")
			result[fac] <- cResult[[1]]
		}
	}
	
	return(result)
}

# New method for estimating diffusion coefficient, based on displacement over
# each bridge.
"diffusionCoefficient.new" <- function (tr, byburst = FALSE, nsteps = 1000) {
	tr <- na.omit(tr)

	resultNames <- unique(id(tr))
	if (byburst) {
		resultNames <- burst(tr)
	}
	
	inputData <- list()
	for (n in resultNames) { inputData[[n]] <- matrix(0, nrow=3, ncol=0) }
	maxDiffCoeff <- rep(0, length(resultNames)) # For each element of the result, find a proper range to search
	names(maxDiffCoeff) <- resultNames

	# For each relocation:
	#  - compute the relevant parameters for the estimation of the diffusion coefficient:
	#      - the displacement of this relocation
	#      - parameters used to compute the location variance from the diffusion coefficient
	#  - compute the ML value for the diffusion coefficient looking only at one bridge.
	#      When the diffusion coefficient gets larger than the maximum of these,
	#      the likelihood of the observations becomes a decreasing function of the diffusion coefficient.
	for (b in 1:length(tr)) {
		burst <- tr[[b]]

		if (nrow(burst) >= 2) { # We can't estimate the likelihood for shorter bursts
			burst$date <- as.double(burst$date) - min(as.double(burst$date))
		
			bdata <- matrix(NA, nrow=3, ncol=nrow(burst)-1)
		
			fac <- attr(burst, ifelse(byburst, "burst", "id"))
			# i runs over all even numbers that have a measurement before and after them
			for (i in 1:(nrow(burst)-1)) {			
				# Squared distance from measured loc to estimated mean
				bdata[1, i] <- 
						(burst$x[i+1] - burst$x[i])^2 + (burst$y[i+1] - burst$y[i])^2
				# The coefficients for a linear function mapping diffusion coefficient to variance
				bdata[2, i] <- (burst$date[i+1]-burst$date[i])
				bdata[3, i] <- burst$loc.var[i] + burst$loc.var[i+1]
			
				# Find the variance at which this bridge reaches its maximum likelihood
				maxVar <- 0.5 * bdata[1, i]
				maxDiffCoeff[fac] <- max(maxDiffCoeff[fac], (maxVar-bdata[3, i])/bdata[2, i])
			}
			
			inputData[[fac]] <- cbind(inputData[[fac]], bdata)
		}
	}
	
	result <- rep(NA, length(resultNames))
	names(result) <- resultNames
	for(fac in names(result)) {
		if (maxDiffCoeff[fac] <= 0) {
			# We already know that the ML is achieved for diff.coeff 0
			result[fac] <- 0
		} else {
			candidates <- seq(0,maxDiffCoeff[fac], length.out=nsteps)
		
			cResult <- .C("diffusion_static", double(1),
					c(inputData[[fac]]), ncol(inputData[[fac]]),
					candidates, length(candidates), PACKAGE="movementAnalysis")
			result[fac] <- cResult[[1]]
		}
	}
	
	return(result)
}
