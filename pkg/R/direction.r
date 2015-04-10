"direction" <- function(tr, time, time.scale=NA) {
	.direction_statistic(tr, time, time.scale, 0, .fn_mean_direction)[,"0",]
}

".fn_mean_direction" <- function(alpha, nu, theta) {
	theta # parameters are noncentrality (i.e. speed of the mean) and theta (angle of the mean)
}

"ddirection" <- function(d, tr, time, time.scale=NA) {
	.direction_statistic(tr, time, time.scale, d, .direction_pdf)
}

".direction_pdf" <- function(alpha, nu, theta) {
	beta <- theta - alpha
	
	exp(-0.5*nu^2) / (2*pi) +
	nu * cos(beta) / sqrt(2 * pi) * exp(-0.5*(nu * sin(beta))^2) * pnorm(nu * cos(beta))
}

"pdirection" <- function(d, tr, time, time.scale=NA, lower=0) {
	.direction_statistic(tr, time, time.scale, d,
			function(alpha, nu, theta) { .direction_cdf(alpha, nu, theta, lower) })
}

".direction_cdf" <- function(alpha, nu, theta, lower=0) {
	res <- sapply(alpha, function(a) {
		integrate(function(x) { .direction_pdf(x, nu, theta) }, lower, a)
	})
	
	unlist(res["value",])
}

#"qdirection" <- function(p, tr, time, time.scale=NA) {
#	.direction_statistic(tr, time, time.scale, p, lmomco::quarice)
#}

".direction_statistic" <- function(tr, time, time.scale, value, fn) {
	if (inherits(time, "POSIXct")) {
		time <- as.double(time)
	}
	if (is.na(time.scale)) {
		if (length(time) > 1) {
			# By default, use the average time between two elements of time argument
			time.scale <- (time[length(time)]-time[1])/(length(time)-1)
		} else {
			stop('Cannot determine appropriate time.scale')
		}
	}

	ids <- unique(adehabitatLT::id(tr))
	
	tr <- na.omit(tr) # This function does not like missing values
		
	# Construct an array indexed by time stamp and an ID.
	# Each element holds the requested statistic for this group at that time
	res <- sapply(time, function(t) {
		params <- velocity(tr, t, time.scale)[,,1]
		# If there is only one ID, its dimension was just dropped; restore it
		dim(params) <- c(length(ids), 3)
		dimnames(params) <- list(ids, c('x','y','var'))
		
		res <- array(NA, dim=c(length(ids),length(value)))
		for (i in 1:length(ids)) {
			nu <- sqrt(sum((params[i,c('x','y')])^2)) / sqrt(params[i, 'var'])
			theta <- atan2(params[i,'y'], params[i,'x'])
			
			if (all(!is.na(c(nu, theta)))) {
			#print(value)
			#print(do.call(fn, list(value, nu, theta)))
				res[i,] <- do.call(fn, list(value, nu, theta))
			} else {
				res[i,] <- NA
			}
		}
		return(res)
	})
	
	# Res is now a list; it really is a 3 dimensional array
	dim(res) <- c(length(ids), length(value), length(time))
	dimnames(res) <- list(ids, value, time)
	return(res)
}
