# Provides location information for each id in the trajectory at the given time(s)
setGeneric("position", function(object, time, groupBy=NULL) standardGeneric("position"))
setMethod(f = "position",
		signature = c(time="POSIXct"),
		definition = function(object, time, groupBy=NULL) {
	p <- position(object, as.double(time), groupBy)
	if (inherits(object, "MoveBBStack")) {
		dimnames(p)[[3]] <- format(time, usetz=TRUE)
	} else {
		rownames(p) <- format(time, usetz=TRUE)
	}
	p
})

setMethod(f = "position",
		signature = c(object="MoveBB", time="numeric"),
		definition = function(object, time) {
	res <- sapply(time, function(t) {
		# Extract the bursts that contain the selected time	
		res <- c(x=NA, y=NA, var=NA)
		if (as.double(object@timestamps[1]) <= t &&
				as.double(object@timestamps[nrow(object)]) >= t) {
					ts <- as.double(object@timestamps)
					j <- max(which(ts >= t)[1], 2)

					diff0T <- integrate(object@diffusion[[1]], ts[j-1], ts[j])
					diff0t <- integrate(object@diffusion[[1]], ts[j-1], t)
					if (diff0T$message != "OK" || diff0t$message != "OK") {
						stop("Problem integrating diffusion coefficient")
					}
					diff0T <- diff0T$value
					diff0t <- diff0t$value
					alpha <- diff0t / diff0T
					res[c("x","y")] <- (1-alpha)*object@coords[j-1,] + alpha*object@coords[j,]
					res["var"] <- (
							  diff0t*(1-alpha) # 1-alpha = diff_tT/diff_0T
							+ (1-alpha)^2 * object@variance[j-1]
							+ alpha^2 * object@variance[j]
					)
		}
		res
	})
	
	dim(res) <- c(3, length(time))
	dimnames(res) <- list(c("x","y","var"), time)
	
	t(res)
})

setMethod(f = "position",
	signature = c(object="MoveBBStack", time="numeric"),
	definition = function (object, time, groupBy=NULL) {
		.mBBStack_statistic(object, time, groupBy, name="position")
})

## Generic function that splits the object of type moveBBStack,
## applies the function identified by "name" to each element and merges the results.
.mBBStack_statistic <- function (object, time, groupBy=NULL, name=c("position","velocity",".speed_statistic",".direction_statistic"), ...) {
	name <- match.arg(name)
	fun <- match.fun(name)

	if (is.null(groupBy)) {
		trs <- split(object)
	} else {
		trs <- split(object, groupBy)
	}

	## Sort time argument, undo after computation
	to <- order(time)
	time <- time[to]

	res <- lapply(trs, function(x, time, ...) {
		## Split into all bursts
		if (class(x) == "MoveBB") {
			xs <- list(x)
		} else {
			xs <- split(x)
		}
		
		## Sort burst by starting time
		xs <- xs[order(sapply(xs, function(x) { x@timestamps[1] }))]
	
		fieldnames <- colnames(fun(xs[[1]], seq_len(0)))
		res <- matrix(NA, ncol=length(fieldnames), nrow=length(time),
				dimnames=list(time, fieldnames))
	
		ti <- 1
		## go through bursts one by one
		for (x in xs) {
			## Find where burst starts and ends in the list of timestamps
			while (ti < length(time) && time[ti] < x@timestamps[1]) { 
				ti <- ti+1 
			}
			tstart <- ti
			while (ti <= length(time) && time[ti] <= x@timestamps[nrow(x)]) {
				ti <- ti+1
			}
			tfinish <- ti ## Points to element past last used
			ts <- seq_len(tfinish-tstart) + tstart - 1
		
			## Call the target function on a subset of timestamps
			if (length(ts) > 0) {
				res[ts,] <- fun(x, time[ts], ...)
			}
		}
		res
	}, time, ...)

	fieldNames <- colnames(res[[1]])
	# Convert to 3d array, set correct order of dimensions
	res <- simplify2array(res)
	dimnames(res) <- list(as.character(time[to]), fieldNames, names(trs))
	res[to,,] <- res # Undo sorting by time
	aperm(res, c(3,2,1))
}
