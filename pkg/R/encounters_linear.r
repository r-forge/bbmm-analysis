"distance" <- function(tr, time) {
	.distance_statistic(tr, time, 0, .fn_mean_distance, 0)[,,"0",]
}

".fn_mean_distance" <- function(value, distr) {
	distr$para[1] # parameters are noncentrality (i.e. mean distance) and std dev
}

"encounter" <- function(d, tr, time) {
	distances <- distance(tr, time)
	
	result <- array(NA, dim=c(dim(distances)[1:2], length(d), dim(distances)[3]),
		dimnames=c(dimnames(distances)[1:2], list(d), dimnames(distances)[3]))
	
	for (dist in d) {
		result[,,as.character(dist),] <- (distances <= dist)
	}
	result
}

"encounterIntervals" <- function(d, tr) {
	tr.data <- lapply(tr, function(burst) {
		burst$t <- as.double(burst$date)
		b <- burst[!is.na(burst$x), c("x","y","diff.coeff","loc.var","t")]
		b <- b[!is.na(b$y),]
		data <- c(t(b)) # flatten into row-major vector
		
		attributes(data) <- list(id=attr(burst, 'id'), nloc=nrow(burst))
		data
	})
	names(tr.data) <- adehabitatLT::burst(tr)
	
	ids <- unique(adehabitatLT::id(tr))
	emptyMatrix <- matrix(NA, nrow=2, ncol=0)
	rownames(emptyMatrix) <- c("start", "end")
	
	result <- replicate(length(ids)^2, emptyMatrix, simplify=FALSE)
	dim(result) <- c(length(ids), length(ids))
	rownames(result) <- colnames(result) <- ids
	rm(emptyMatrix)
	diag(result) <- NA

	bursts <- .overlappingBursts(tr)
	for (i in 1:nrow(bursts)) {
		data1 <- tr.data[[bursts[i, 'burst1']]]
		n1 <- attr(data1, 'nloc')
		data2 <- tr.data[[bursts[i, 'burst2']]]
		n2 <- attr(data2, 'nloc')
		
		cResult <- .C("encounterLinearIntervals", double(1), as.double(d),
				as.integer(n1), as.double(data1),
				as.integer(n2), as.double(data2),
				as.integer(0), double(2*(n1 + n2)), PACKAGE="movementAnalysis")
	
		intervals <- matrix(cResult[[8]], nrow=2)
		intervals <- intervals[, !is.nan(intervals[1,]) & !is.nan(intervals[2,]), drop=F]
	
		# If one interval starts when the previous one stops, merge them into one interval
		if (ncol(intervals) > 1) {
			connect <- intervals[2,-ncol(intervals)] >= intervals[1, -1]
			start <- c(TRUE, !connect)
			end   <- c(!connect, TRUE)
			intervals <- matrix(c(intervals[1,start], intervals[2,end]), nrow=2, byrow=T)
		}
		
		result[[attr(data1, 'id'), attr(data2, 'id')]] <- cbind(
				result[[attr(data1, 'id'), attr(data2, 'id')]], intervals)
	}

	intervals <- lapply(result, function(t) { as.POSIXct(t, origin="1970-01-1")} )
	attributes(intervals) <- attributes(result)

	# Make the result matrix symmetric
	intervals[lower.tri(intervals)] <- t(intervals)[lower.tri(intervals)]
	intervals
}
