"speed" <- function(tr, time, time.scale=NA, groupBy=NULL) {
	s <- .speed_statistic(tr, time, time.scale, 0, .fn_mean_speed, groupBy)
	# Drop the second dimension, which indexes the '0' value
	# This construction prevents dropping other dimensions with unit extent
	dn <- dimnames(s)
	dim(s) <- dim(s)[-2]
	dimnames(s) <- dn[-2]
	s
}

".fn_mean_speed" <- function(value, distr) {
	distr$para[1] # parameters are noncentrality (i.e. speed of the mean) and std dev
}

"mu.speed" <- function(tr, time, time.scale=NA, groupBy=NULL) {
	s <- .speed_statistic(tr, time, time.scale, 0, .rice_mean_speed, groupBy)
	# Drop the second dimension, which indexes the '0' value
	# This construction prevents dropping other dimensions with unit extent
	dn <- dimnames(s)
	dim(s) <- dim(s)[-2]
	dimnames(s) <- dn[-2]
	s
}

".rice_mean_speed" <- function(value, distr) {
	# parameters are noncentrality (i.e. speed of the mean) and std dev
	nu <- distr$para[1]
	sigma <- distr$para[2]
	
	# Return the mean of the Rice distribution with parameters nu and sigma
	sigma * sqrt(pi / 2) * lmomco::LaguerreHalf(-0.5 * (nu/sigma)^2)
}

"dspeed" <- function(v, tr, time, time.scale=NA, groupBy=NULL) {
	.speed_statistic(tr, time, time.scale, v, lmomco::pdfrice, groupBy)
}

"pspeed" <- function(v, tr, time, time.scale=NA, groupBy=NULL) {
	.speed_statistic(tr, time, time.scale, v, lmomco::cdfrice, groupBy)
}

"qspeed" <- function(p, tr, time, time.scale=NA, groupBy=NULL) {
	.speed_statistic(tr, time, time.scale, p, lmomco::quarice, groupBy)
}

setGeneric(".speed_statistic", function(object, time, time.scale, value, fn, groupBy=NULL) standardGeneric(".speed_statistic"))
setMethod(f = ".speed_statistic",
	signature = c(time="POSIXct", time.scale="numeric"),
	definition = function(object, time, time.scale, value, fn, groupBy=NULL) {
		s <- .speed_statistic(object, as.double(time), time.scale, value, fn, groupBy)
		dimnames(s)[[if(length(dim(s)) == 3) { 3 } else { 1 }]] <- format(time, usetz=TRUE)
		s
})

setMethod(f=".speed_statistic",
		signature = c(object="MoveBB", time="numeric", time.scale="numeric", value="numeric"),
		definition = function(object, time, time.scale, value, fn) {
	vl <- velocity(object, time, time.scale)
	
	r <- t(apply(vl, 1, function(v) {
		mu <- sqrt(sum((v[c('x','y')])^2))
		sd <- sqrt(v['var'])
		
		if (!any(is.na(c(mu, sd)))) {
			distr <- list(type="rice", para=c(mu, sd), source="bbtraj", IFAIL=0, IFAILTEXT="Successful parameter estimation")
			res <- do.call(fn, list(value, distr))
		} else {
			res <- rep(NA, length(value))
		}
		names(res) <- as.character(value)
		res
	}))
	matrix(r, nrow(vl), length(value), dimnames=list(time, value))
})

setMethod(f=".speed_statistic",
		signature = c(object="MoveBBStack", time="numeric", time.scale="numeric", value="numeric"),
		definition = function(object, time, time.scale, value, fn, groupBy=NULL) {
	.mBBStack_statistic(object, time, groupBy, name=".speed_statistic", time.scale=time.scale, value=value, fn=fn)
})

