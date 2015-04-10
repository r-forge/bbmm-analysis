"direction" <- function(tr, time, time.scale=NA, groupBy=NULL) {
	d <- .direction_statistic(tr, time, time.scale, 0, .fn_mean_direction, groupBy)
	# Drop the second dimension, which indexes the '0' value
	# This construction prevents dropping other dimensions with unit extent
	dn <- dimnames(d)
	dim(d) <- dim(d)[-2]
	dimnames(d) <- dn[-2]
	d
}

".fn_mean_direction" <- function(alpha, nu, theta) {
	theta # parameters are noncentrality (i.e. speed of the mean) and theta (angle of the mean)
}

"ddirection" <- function(d, tr, time, time.scale=NA, groupBy=NULL) {
	.direction_statistic(tr, time, time.scale, d, .direction_pdf, groupBy)
}

".direction_pdf" <- function(alpha, nu, theta) {
	eta <- theta - alpha
	
	exp(-0.5*nu^2) / (2*pi) +
	nu * cos(eta) / sqrt(2 * pi) * exp(-0.5*(nu * sin(eta))^2) * pnorm(nu * cos(eta))
}

"pdirection" <- function(d, tr, time, time.scale=NA, groupBy=NULL, lower=0) {
	.direction_statistic(tr, time, time.scale, d,
			function(alpha, nu, theta) { .direction_cdf(alpha, nu, theta, lower) }, groupBy)
}

".direction_cdf" <- function(alpha, nu, theta, lower=0) {
	res <- sapply(alpha, function(a) {
		integrate(function(x) { .direction_pdf(x, nu, theta) }, lower, a)
	})
	
	unlist(res["value",])
}

#"qdirection" <- function(p, tr, time, time.scale=NA, groupBy=NULL) {
#	.direction_statistic(tr, time, time.scale, p, lmomco::quarice, groupBy)
#}

setGeneric(".direction_statistic", function(object, time, time.scale, value, fn, groupBy=NULL) standardGeneric(".direction_statistic"))
setMethod(f = ".direction_statistic",
	signature = c(time="POSIXct", time.scale="numeric"),
	definition = function(object, time, time.scale, value, fn, groupBy=NULL) {
		s <- .direction_statistic(object, as.double(time), time.scale, value, fn, groupBy)
		dimnames(s)[[3]] <- format(time, usetz=TRUE)
		s
})

setMethod(f=".direction_statistic",
		signature = c(object="MoveBB", time="numeric", time.scale="numeric", value="numeric"),
		definition = function(object, time, time.scale, value, fn) {
	vl <- velocity(object, time, time.scale)
	
	r <- t(apply(vl, 1, function(v) {
		nu <- sqrt(sum((v[c('x','y')])^2)) / sqrt(v['var'])
		theta <- atan2(v['y'], v['x'])
		
		if (!any(is.na(c(nu, theta)))) {
			#print(value)
			#print(do.call(fn, list(value, nu, theta)))
			res <- do.call(fn, list(value, nu, theta))
		} else {
			res <- rep(NA, length(value))
		}
		names(res) <- as.character(value)
		res
	}))
	matrix(r, nrow(vl), length(value), dimnames=list(time, value))
})

setMethod(f=".direction_statistic",
		signature = c(object="MoveBBStack", time="numeric", time.scale="numeric", value="numeric"),
		definition = function(object, time, time.scale, value, fn, groupBy=NULL) {
	.mBBStack_statistic(object, time, groupBy, name=".direction_statistic", time.scale=time.scale, value=value, fn=fn)
})

