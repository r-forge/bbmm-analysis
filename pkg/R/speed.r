"speed" <- function(tr, time, time.scale=NA) {
	.speed_statistic(tr, time, time.scale, 0, .fn_mean_speed)[,"0",]
}

".fn_mean_speed" <- function(value, distr) {
	distr$para[1] # parameters are noncentrality (i.e. speed of the mean) and std dev
}

"mu.speed" <- function(tr, time, time.scale=NA) {
	.speed_statistic(tr, time, time.scale, 0, .rice_mean_speed)[,"0",]
}

".rice_mean_speed" <- function(value, distr) {
	# parameters are noncentrality (i.e. speed of the mean) and std dev
	nu <- distr$para[1]
	sigma <- distr$para[2]
	
	# Return the mean of the Rice distribution with parameters nu and sigma
	sigma * sqrt(pi / 2) * lmomco::LaguerreHalf(-0.5 * (nu/sigma)^2)
}

"dspeed" <- function(v, tr, time, time.scale=NA) {
	.speed_statistic(tr, time, time.scale, v, lmomco::pdfrice)
}

"pspeed" <- function(v, tr, time, time.scale=NA) {
	.speed_statistic(tr, time, time.scale, v, lmomco::cdfrice)
}

"qspeed" <- function(p, tr, time, time.scale=NA) {
	.speed_statistic(tr, time, time.scale, p, lmomco::quarice)
}

".speed_statistic" <- function(tr, time, time.scale, value, fn) {
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

	ids <- unique(id(tr))
	
	tr <- na.omit(tr) # This function does not like missing values
		
	# Construct an array indexed by time stamp and an ID.
	# Each element holds the requested statistic for this group at that time
	res <- sapply(time, function(t) {
		params <- velocity(tr, t, time.scale)[,,1]
		# If there is only one ID, its dimension was just dropped; restore it
		dim(params) <- c(length(ids), 3)
		dimnames(params) <- list(ids, c('x','y','var'))
		
		res <- array(rep(NA, length(value)*length(ids)), dim=c(length(ids),length(value)))
		for (i in 1:length(ids)) {
			mu <- sqrt(sum((params[i,c('x','y')])^2))
			s2 <- params[i, 'var']
			
			if (all(!is.na(c(mu, s2)))) {
				distr <- list(type="rice", para=c(mu, sqrt(s2)), source="bbtraj", IFAIL=0, IFAILTEXT="Successful parameter estimation")
				res[i,] <- do.call(fn, list(value, distr))
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
