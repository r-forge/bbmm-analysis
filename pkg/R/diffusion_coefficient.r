setGeneric("diffusion", function(object) standardGeneric("diffusion"))
setMethod(f = "diffusion",
	signature = c(object=".BBInfo"),
	definition = function (object) {
		object@diffusion
})

setGeneric("diffusion<-", function(object, value) standardGeneric("diffusion<-"))
setMethod(f = "diffusion<-",
		signature = c(object=".BBInfo"),
		definition = function (object, value) {
	## Convert the given values to the proper format
	diffusion <- sapply(value, diffusion.convert)

	if (!is.null(attr(value, 'grouping'))) {
		# This is the result of a call to diffusionCoefficient
		# with groupBy != NULL, expand the result
		diffusion <- diffusion[as.character(attr(value, 'grouping'))]
		names(diffusion) <- names(attr(value, 'grouping'))
	} else if (is.null(names(diffusion)) && inherits(object, "MoveBBStack")) {
		names(diffusion) <- levels(object@trackId)
	}
	
	object@diffusion <- diffusion
	object
})

setGeneric("diffusion.convert", function(dc.obj) standardGeneric("diffusion.convert"))

setMethod(f="diffusion.convert", signature=c(dc.obj="numeric"),
		definition = function(dc.obj) {
	# Turn the diffusion coefficients into (constant) functions of time
	sapply(dc.obj, function(d) { 
			dc.val <- d
			function(ts) { rep(dc.val, length(ts)) } 
	})
})

setMethod(f="diffusion.convert", signature=c(dc.obj="function"),
		definition = function(dc.obj) {
	dc.obj
})

setMethod(f="diffusion.convert", signature=c(dc.obj="matrix"),
		definition = function(dc.obj) {
	## Matrix should have two columns: time and diffusion coefficient
	## Each row specifies the diffusion coefficient from the given time
	## up to the timestamp in the next row
	function(ts) {
		.interpolate(dc.obj, ts, interpolation="constant")
	}
})

setMethod(f="diffusion.convert", signature=c(dc.obj="dBMvariance"),
		definition = function(dc.obj) {
	## Convert to matrix and redirect
	## brownian.motion.variance.dyn uses minutes by default, how to fix this?
	dc <- dc.obj@means / 60
	dc[is.na(dc)] <- 0
	diffusion.convert(cbind(dc.obj@timestamps, dc))
})

setGeneric("diffusionCoefficient", function(tr, groupBy=NULL, nsteps=1000,
	method=c("horne","new")) standardGeneric("diffusionCoefficient"))
setMethod(f = "diffusionCoefficient",
	signature = c(tr="MoveBB"),
	definition = function (tr, nsteps=1000, method=c("horne","new")) {
		.diffusionCoefficient(list(tr), NULL, nsteps, method)
})

setMethod(f = "diffusionCoefficient",
	signature = c(tr="MoveBBStack"),
	definition = function (tr, groupBy=NULL, nsteps=1000, method=c("horne","new")) {
		.diffusionCoefficient(split(tr), groupBy, nsteps, method)
})

".diffusionCoefficient" <- function(trs, groupBy, nsteps, method=c("horne","new")) {
	method <- match.arg(method)
	burstNames <- unlist(unname(lapply(trs, IDs, groupBy)))
	switch(method,
		horne=.diffusionCoefficient.horne(trs, burstNames, nsteps),
		new  =.diffusionCoefficient.new  (trs, burstNames, nsteps)
	)
}

"diffusion.LL" <- function(trs, dc.values) {
	# For each even numbered measurement in a burst (except the last if burst has even length):
	#  - compute the relevant parameters for the estimation of the diffusion coefficient:
	#      - the distance between the observation and the mean derived from the neighbouring observations
	#      - parameters used to compute the location variance from the diffusion coefficient
	#  - compute the ML value for the diffusion coefficient looking only at one bridge.
	#      When the diffusion coefficient gets larger than the maximum of these,
	#      the likelihood of the observations becomes a decreasing function of the diffusion coefficient.
	#	for (b in 1:length(trs)) {
	res <- sapply(trs, function(burst) {
		if (nrow(burst) >= 3) { # We can't estimate the likelihood for shorter bursts
			burst$date <- as.double(burst@timestamps) - min(as.double(burst@timestamps))
				
		#	bdata <- matrix(NA, nrow=3, ncol=floor((nrow(burst)-1)/2))
		
			## i runs over all even numbers that have a measurement before and after them
			ll <- rep(0, length(dc.values))
			for (i in seq(2, nrow(burst)-1, by=2)) {
				alpha <- (burst$date[i] - burst$date[i-1]) / (burst$date[i+1] - burst$date[i-1])
			
				## Squared distance from measured loc to estimated mean
				d2 <- sum((
						              burst@coords[i,]
						- (1-alpha) * burst@coords[i-1,]
						- alpha     * burst@coords[i+1,]
					)^2)
				## The coefficients for a linear function mapping diffusion coefficient to variance
				c1 <- (burst$date[i+1]-burst$date[i-1]) * alpha * (1-alpha)
				c0 <- (1-alpha)^2*burst@variance[i-1] + alpha^2*burst@variance[i+1]
				
				var <- c1 * dc.values + c0
				ll <- ll - log(2*pi) - log(var) - (0.5*d2 / var)
			}
			ll
		} else {
			rep(NA, length(dc.values))
		}
	})
	res <- matrix(res, nrow=length(dc.values))
	attr(res, "dc.values") <- dc.values
	rownames(res) <- as.character(dc.values)
	res
}

".diffusionCoefficient.horne" <- function(trs, burstNames, nsteps) {
		resultNames <- unique(burstNames)
		
		inputData <- list()
		for (n in resultNames) { inputData[[n]] <- matrix(0, nrow=3, ncol=0) }

		# For each element of the result, find a proper range to search
		maxDiffCoeff <- rep(0, length(resultNames))
		names(maxDiffCoeff) <- resultNames

		# For each even numbered measurement in a burst (except the last if burst has even length):
		#  - compute the relevant parameters for the estimation of the diffusion coefficient:
		#      - the distance between the observation and the mean derived from the neighbouring observations
		#      - parameters used to compute the location variance from the diffusion coefficient
		#  - compute the ML value for the diffusion coefficient looking only at one bridge.
		#      When the diffusion coefficient gets larger than the maximum of these,
		#      the likelihood of the observations becomes a decreasing function of the diffusion coefficient.
		for (b in 1:length(trs)) {
			burst <- trs[[b]]
			if (nrow(burst) >= 3) { # We can't estimate the likelihood for shorter bursts
				burst$date <- as.double(burst@timestamps) - min(as.double(burst@timestamps))
				
				bdata <- matrix(NA, nrow=3, ncol=floor((nrow(burst)-1)/2))
		
				fac <- as.character(burstNames[b])
				# i runs over all even numbers that have a measurement before and after them
				for (i in seq(2, nrow(burst)-1, by=2)) {
					alpha <- (burst$date[i] - burst$date[i-1]) / (burst$date[i+1] - burst$date[i-1])
			
					# Squared distance from measured loc to estimated mean
					bdata[1, i/2] <- sum((
							              burst@coords[i,]
							- (1-alpha) * burst@coords[i-1,]
							- alpha     * burst@coords[i+1,]
						)^2)
					# The coefficients for a linear function mapping diffusion coefficient to variance
					bdata[2, i/2] <- (burst$date[i+1]-burst$date[i-1]) * alpha * (1-alpha)
					bdata[3, i/2] <- (1-alpha)^2*burst@variance[i-1] + alpha^2*burst@variance[i+1]
			
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
				# We already know that the ML is achieved for diff.coeff == 0
				result[fac] <- 0
			} else {
				candidates <- seq(0,maxDiffCoeff[fac], length.out=nsteps)
		
				cResult <- .C("diffusion_static", double(1),
						c(inputData[[fac]]), ncol(inputData[[fac]]),
						candidates, length(candidates), PACKAGE="movementAnalysis")
				result[fac] <- cResult[[1]]
			}
		}
		
		if (any(names(burstNames) != burstNames)) { # groupBy != NULL
			names(burstNames) <- names(trs)
			attr(result, 'grouping') <- burstNames
		}
		
		LL <- rep(NA, length(trs))
		names(LL) <- names(trs)
		for (tn in 1:length(trs)) {
			LL[tn] <- diffusion.LL(trs[tn], result[burstNames[tn]])
		}
		
		attr(result, "LL") <- LL
		return(result)
}

# New method for estimating diffusion coefficient, based on displacement over
# each bridge.
".diffusionCoefficient.new" <- function(trs, burstNames, nsteps) {
	resultNames <- unique(burstNames)
		
	inputData <- list()
	for (n in resultNames) { inputData[[n]] <- matrix(0, nrow=3, ncol=0) }

	# For each element of the result, find a proper range to search
	maxDiffCoeff <- rep(0, length(resultNames))
	names(maxDiffCoeff) <- resultNames
	
	# For each relocation:
	#  - compute the relevant parameters for the estimation of the diffusion coefficient:
	#      - the displacement of this relocation
	#      - parameters used to compute the location variance from the diffusion coefficient
	#  - compute the ML value for the diffusion coefficient looking only at one bridge.
	#      When the diffusion coefficient gets larger than the maximum of these,
	#      the likelihood of the observations becomes a decreasing function of the diffusion coefficient.
	for (b in 1:length(trs)) {
		burst <- trs[[b]]

		if (nrow(burst) >= 2) { # We can't estimate the likelihood for shorter bursts
			burst$date <- as.double(burst@timestamps) - min(as.double(burst@timestamps))
		
			bdata <- matrix(NA, nrow=3, ncol=nrow(burst)-1)
		
			fac <- burstNames[b]
			# i runs over all even numbers that have a measurement before and after them
			for (i in 1:(nrow(burst)-1)) {
				# Squared distance between relocations
				bdata[1, i] <- sum((burst@coords[i+1,] - burst@coords[i,])^2)
				# The coefficients for a linear function mapping diffusion coefficient to variance
				bdata[2, i] <- (burst$date[i+1]-burst$date[i])
				bdata[3, i] <- burst@variance[i] + burst@variance[i+1]
			
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

## Interpolates values between x-coordinates of xy at the coordinates given by q
".interpolate" <- function(xy, q, interpolation=c("linear","constant")) {
	interpolation <- match.arg(interpolation)
	
	xy <- rbind(xy, 0)
	xy[nrow(xy),1] <- Inf
	
	## Find the interval in which each q value occurs
	intervals <- rep(NA, length(q))
	iv <- 0
	for (i in seq_len(length(q))) {
		iv.last <- iv
		d.iv <- 1
		## Exponential search to find upper bound for iv
		while (iv < nrow(xy) && q[i] >= xy[iv+1,1]) {
			iv <- min(iv+d.iv,nrow(xy))
			
			print(paste(iv.last, iv, d.iv, xy[iv,1], q[i]))
			
			d.iv <- d.iv*2
		}
		## The real iv is between iv.last and iv, find by binary search
		while (iv-iv.last > 0) {
			m <- ceiling((iv+iv.last)/2)
			if (q[i] >= xy[m,1]) {
				iv.last <- m
			} else {
				iv <- m-1
			}
		}
		print(iv)
		
		intervals[i] <- iv
	}
	
	res <- mapply(function(qe, iv) {
		if (iv == 0 || iv == nrow(xy)) {
			return(xy[nrow(xy),-1])
		}
		y1 <- xy[iv,-1]
		if (interpolation=="constant") {
			return(y1)
		} else {
			y2 <- xy[iv+1,-1]
			alpha <- (qe- xy[iv,1]) / (xy[iv+1,1]-xy[iv,1])
			
			return((1-alpha)*y1 + alpha*y2)
		}		
	}, q, intervals)
	names(res) <- q
	res
}
interpolate <- .interpolate
