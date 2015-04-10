"ddistance" <- function(d, tr, time) {
	.distance_statistic(tr, time, d, lmomco::pdfrice, ifelse(distance==0, Inf, 0))
}

"pdistance" <- function(d, tr, time) {
	.distance_statistic(tr, time, d, lmomco::cdfrice, 1)
}

"qdistance" <- function(p, tr, time) {
	.distance_statistic(tr, time, p, lmomco::quarice, 0)
}

".distance_statistic" <- function(tr, time, value, fn, diagonal_value=NA) {
	ids <- unique(id(tr))
	
	tr <- bbFilterNA(tr) # This function does not like missing values
		
	# Construct an array indexed by time stamp and two IDs.
	# Each element holds the requested statistic for these two groups at that time
	res <- sapply(time, function(t) {
		params <- position(tr, t)[,,1]
		
		res <- array(rep(diagonal_value, length(value)*length(ids)^2), dim=c(length(ids),length(ids),length(value)))
		for (i in 1:(length(ids)-1)) {
			for (j in (i+1):length(ids)) {
				mu <- sqrt(sum((params[i,c('x','y')] - params[j,c('x','y')])^2))
				s2 <- params[i, 'var'] + params[j, 'var']
				
				if (all(!is.na(c(mu, s2)))) {
					distr <- list(type="rice", para=c(mu, sqrt(s2)), source="bbtraj", IFAIL=0, IFAILTEXT="Successful parameter estimation")
					res[i,j,] <- res[j,i,] <- do.call(fn, list(value, distr))
				} else {
					res[i,j,] <- res[j,i,] <- NA
				}
			}
		}
		return(res)
	})
	
	# Res is now a list; it really is a 4 dimensional array
	dim(res) <- c(length(ids), length(ids), length(value), length(time))
	dimnames(res) <- list(ids, ids, value, time)
	
	return(res)
}
