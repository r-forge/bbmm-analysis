# Provides location information for each id in the trajectory at the given time(s)
setGeneric("position", function(object, time, groupBy=NULL) standardGeneric("position"))
setMethod(f = "position",
	signature = c(time="POSIXct"),
	definition = function(object, time) {
	p <- position(object, as.double(time))
	dimnames(p)[[3]] <- format(time, usetz=TRUE)
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
	signature = c(object="MoveBBStack"),
	definition = function (object, time, groupBy=NULL) {
		.mBBStack_statistic(object, time, groupBy, name="position")
})

## Generic function that splits the object of type moveBBStack,
## applies the function identified by "name" to each element and merges the results.
.mBBStack_statistic <- function (object, time, groupBy=NULL, name=c("position","velocity",".speed_statistic",".direction_statistic"), ...) {
	name <- match.arg(name)

	tr <- split(object)
	
	tr.id <- .IDs(object, groupBy)
	ids <- unique(tr.id)
	
	# Collect positions for each burst and group them by their ID
	tr.pos <- lapply(tr, name, time, ...)
	id.pos <- split(tr.pos, tr.id)
	
	# Merge the positions of the bursts for each ID
	res <- lapply(id.pos, function(p) {
		p.set <- sapply(p, function(p) {
			apply(p, 1, function(row) { all(!is.na(row)) })
		})
		dim(p.set) <- c(length(time), length(p)) # Make sure it is always a matrix
		dimnames(p.set) <- list(time, names(p))

		# Check that no two bursts contain a position for the same time
		if (any(apply(p.set, 1, sum) > 1)) {
			stop("Multiple bursts for same ID overlap")
		}
		
		# Select all the set values from the different bursts
		p.select <- apply(p.set, 1, which.max)
		pos <- sapply(1:length(p.select), function(i) {
			p[[p.select[i]]][i,]
		})
		matrix(as.double(pos), ncol=length(time), dimnames=list(colnames(p[[1]]) ,time))
	})
	
	fieldNames <- rownames(res[[1]])
	# Convert to 3d array, set correct order of dimensions
	res <- simplify2array(res)
	dim(res) <- c(length(fieldNames), length(time), length(ids))
	dimnames(res) <- list(fieldNames, as.character(time), ids)
	aperm(res, c(3,1,2))
}
