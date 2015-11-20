"ddistance" <- function(d, tr, time, groupBy=NULL) {
	.distance_statistic(tr, time, d, lmomco::pdfrice, ifelse(d==0, Inf, 0), groupBy)
}

"pdistance" <- function(d, tr, time, groupBy=NULL) {
	.distance_statistic(tr, time, d, lmomco::cdfrice, 1, groupBy)
}

"qdistance" <- function(p, tr, time, groupBy=NULL) {
	.distance_statistic(tr, time, p, lmomco::quarice, 0, groupBy)
}

setGeneric(".distance_statistic", function(tr, time, value, fn, diagonal_value=NA, groupBy=NULL) standardGeneric(".distance_statistic"))
setMethod(f = ".distance_statistic",
	signature = c(tr="MoveBBStack", time="POSIXct"),
	definition = function (tr, time, value, fn, diagonal_value=NA, groupBy=NULL) {
		d <- .distance_statistic(tr, as.double(time), value, fn, diagonal_value, groupBy)
		dimnames(d)[[4]] <- format(time, usetz=TRUE)
		d
})

setMethod(f = ".distance_statistic",
	signature = c(tr="MoveBBStack", time="numeric"),
	definition = function (tr, time, value, fn, diagonal_value=NA, groupBy=NULL) {
		ids <- unique(.IDs(tr, groupBy))
	
		# Construct an array indexed by time stamp and two IDs.
		# Each element holds the requested statistic for these two groups at that time
		positions <- position(tr, time, groupBy)
		res <- mclapply(seq_len(length(time)), function(t) {
			params <- positions[,,t]
		
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
		res <- simplify2array(res)
		dimnames(res) <- list(ids, ids, value, time)
	
		return(res)
})
