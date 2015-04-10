"distance" <- function(trajectory, time, ids=NA, simplify=TRUE) {
	.distance_statistic(trajectory, time, NA, .fn_mean_distance, 0, ids=ids, simplify=simplify)
}

".fn_mean_distance" <- function(value, distr) {
	distr$para[1] # parameters are noncentrality (i.e. mean distance) and variance
}

"encounter" <- function(trajectory, time, threshold, ids=NA, simplify=TRUE) {
	return(distance_linear(trajectory,time, ids, simplify) <= threshold)
}

"encounterIntervals" <- function(b1, b2, threshold) {
	useFields <- c("x","y","s12","s22","t")

	# serialize x,y, s1^2, s2^2 and time from the bursts in row-major order
	b1$t <- as.double(b1$date)
	b1 <- b1[!is.na(b1$x), useFields]
	b1 <- b1[!is.na(b1$y),]
	b1 <- b1[!is.na(b1$s22),]
	data1 <- c(t(b1)) # flatten into row-major vector
	
	b2$t <- as.double(b2$date)
	b2 <- b2[!is.na(b2$x), useFields]
	b2 <- b2[!is.na(b2$y),]
	b2 <- b2[!is.na(b2$s22),]
	data2 <- c(t(b2)) # flatten into row-major vector


	cResult <- .C("encounterLinearIntervals", double(1), as.double(threshold),
			as.integer(nrow(b1)), as.double(data1),
			as.integer(nrow(b2)), as.double(data2),
			as.integer(0), double(2*(nrow(b1) + nrow(b2))), PACKAGE="movementAnalysis")
	
	intervals <- matrix(cResult[[8]], nrow=2)
	intervals <- intervals[, !is.nan(intervals[1,]) & !is.nan(intervals[2,]), drop=F]
	
	if (ncol(intervals) > 1) {
		connect <- intervals[2,-ncol(intervals)] >= intervals[1, -1]
		start <- c(TRUE, !connect)
		end   <- c(!connect, TRUE)
		intervals <- matrix(c(intervals[1,start], intervals[2,end]), nrow=2, byrow=T)
	}
	
	rownames(intervals) <- c("start", "end")
	return(intervals)
}
