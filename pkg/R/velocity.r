setGeneric("velocity", function(object, time, time.scale=NA, groupBy=NULL) standardGeneric("velocity"))
setMethod(f = "velocity",
	signature = c(time="POSIXct", time.scale="numeric"),
	definition = function(object, time, time.scale, groupBy) {
		v <- velocity(object, as.double(time), time.scale, groupBy)
		dimnames(v)[[3]] <- format(time, usetz=TRUE)
		v
})

setMethod(f = "velocity",
	signature = c(object="MoveBB", time="numeric", time.scale="numeric"),
	definition = function(object, time, time.scale) {
		res <- sapply(time, function(t) {
			# Check if the required time interval is covered by object
			res <- c(x=NA, y=NA, var=NA)
			t.s <- t - time.scale/2
			t.e <- t + time.scale/2
			if (as.double(object@timestamps[1]) <= t.s &&
					as.double(object@timestamps[nrow(object)]) >= t.e) {
				p.s <- position(object, t.s)
				p.e <- position(object, t.e)

				# The mean is always the velocity of the mean position
				res[1:2] <- (p.e[1:2] - p.s[1:2]) / time.scale

				# find in which links t.s and t.e are
				j.s <- max(which(object@timestamps >= t.s)[1], 2)
				j.e <- max(which(object@timestamps >= t.e)[1], 2)

				# Depending on the number of measurements between t.s and t.e,
				# the variance is computed differently
				od <- object@diffusion[[1]]
				if (j.s == j.e) {
					# t.s and t.e are in the same bridge
					# See: Sijben, Stef. "Computational Movement Analysis Using Brownian Bridges." 
					# 	Master's thesis, Eindhoven University of Technology (2013).
					Diff <- list("2f"=integrate(od, t.e, object@timestamps[j.s]),
							"12"=integrate(od, t.s, t.e),
							"sf"=integrate(od, object@timestamps[j.s-1], object@timestamps[j.s]))
					if (any(sapply(Diff, function(d) { d$message }) != "OK")) {
						stop("Problem integrating diffusion coefficient")
					}
					Diff <- sapply(Diff, function(d) { d$value })
					
					alpha <- 1-((Diff['12']+Diff['2f'])/Diff['sf'])
					beta  <- Diff['12']/(Diff['12']+Diff['2f'])
					res['var'] <- (Diff['12'] / (Diff['sf']*time.scale))^2 * sum(object@variance[j.s+c(-1,0)]) + 
							beta/time.scale^2*(alpha*Diff['12'] + Diff['2f'])
				} else if (j.s == j.e-1) {
					# t.s and t.e are in consecutive bridges of one burst
					# See Sijben, Stef. "Computational Movement Analysis Using Brownian Bridges." 
					# 	Master's thesis, Eindhoven University of Technology (2013).
					Diff <- list("s1"=integrate(od, object@timestamps[j.s-1], t.s),
							"si"=integrate(od, object@timestamps[j.s-1], object@timestamps[j.s]),
							"i2"=integrate(od, object@timestamps[j.e-1], t.e),
							"if"=integrate(od, object@timestamps[j.e-1], object@timestamps[j.e]))
					if (any(sapply(Diff, function(d) { d$message }) != "OK")) {
						stop("Problem integrating diffusion coefficient")
					}
					Diff <- sapply(Diff, function(d) { d$value })
					
					alpha <- Diff['s1']/Diff['si']
					beta  <- Diff['i2']/Diff['if']
					res['var'] <- (p.s[1,'var'] + p.e[1,'var'] - 2*alpha*(1-beta)*object@variance[j.s])/time.scale^2
				} else {
					# t.s and t.e are further apart in time, treat as independent
					res['var'] <- (p.e[1,'var'] + p.s[1,'var']) / (t.e-t.s)^2
				} # /if j.s == j.e
			} # /if t.s and t.e inside time range
			res
		}) # /sapply (time, ...)
		
		dim(res) <- c(3, length(time))
		dimnames(res) <- list(c("x","y","var"), time)
	
		t(res)
})

setMethod(f = "velocity",
	signature = c(object="MoveBBStack", time="numeric", time.scale="numeric"),
	definition= function(object, time, time.scale, groupBy=NULL) {
		.mBBStack_statistic(object, time, groupBy, name="velocity", time.scale=time.scale)
})

