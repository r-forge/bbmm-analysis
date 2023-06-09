setGeneric("utilizationDistribution", function(object, raster=NULL, timestepSize=60,
		grid.dim=100, grid.pad=0.2, groupBy=NULL, cutoff.level=1) {
	standardGeneric("utilizationDistribution")
})

setMethod(f = "utilizationDistribution", 
		signature = c(object = "MoveBB", raster = "RasterLayer", timestepSize = "numeric"), 
		function(object, raster, timestepSize=60, cutoff.level=1) {
	r <- .utilizationDistribution(object, raster, timestepSize, cutoff.level=1)
	new("UD", r/sum(values(r)))
})

setMethod(f = "utilizationDistribution", 
          signature = c(object = "MoveBBStack", raster = "RasterLayer", timestepSize = "numeric"), 
          function(object, raster, timestepSize=60, groupBy=NULL, cutoff.level=1) {
	os <- split(split(object), IDs(object, groupBy))

	# Group the results if necessary
	values(raster) <- 0
	UDs <- lapply(os, function(tr) {
		u <- raster
		for (burst in tr) {
			values(u) <- values(u) + values(.utilizationDistribution(burst, raster, timestepSize, cutoff.level=1))
		}
		u
	})

	# Normalize all UDs, then create stack and return
	new("UDStack", stack(lapply(UDs, function(u) { u/sum(values(u)) })))
})

setMethod(f = "utilizationDistribution", 
		signature = c(object = ".MoveTrack", raster="missing"),
		function(object, timestepSize=60, grid.dim=100, grid.pad=0.2, groupBy=NULL, cutoff.level=1) {
	if (!is.numeric(grid.dim) || !is.numeric(grid.pad)) {
		stop("grid.dim and grid.pad must be numeric")
	}
	raster <- .defaultRaster(extent(object), grid.dim, grid.pad, object@proj4string)
	utilizationDistribution(object=object, raster=raster,
			timestepSize=timestepSize, groupBy=groupBy, cutoff.level=cutoff.level)
})

".utilizationDistribution" <- function(object, raster, timestepSize, cutoff.level) {
	timesteps <- .UDtimesteps(object, timestepSize)

	# extract the grid lines from the provided grid
	ext <- extent(raster)
	xc <- seq(ext@xmin, ext@xmax, length.out=ncol(raster))
	yc <- seq(ext@ymin, ext@ymax, length.out=nrow(raster))
	
	# For each timestep, compute the mean location, variance, and a weighing factor
	gxc <- xc
	gyc <- yc
#		if (cutoff.level < 1) {
#			# Find a box that is guaranteed to contain the requested contour.
#			# The bounding box containing the cutoff.level contour for every time step
#			# is guaranteed to contain this contour for the whole trajectory.
#			vDiff <- sqrt(ts[,'var'] * -2 * log(1-cutoff.level))
#			gRows <- which(xc >= min(ts[,'x'] - vDiff) & xc <= max(ts[,'x'] + vDiff))
#			gxc <- xc[gRows]
#			gCols <- which(yc >= min(ts[,'y'] - vDiff) & yc <= max(ts[,'y'] + vDiff))
#			gyc <- yc[gCols]
#		}
		
	pad <- rep(0, 15)
	
	cResult <- .OpenCL("utilizationDistribution", length(gxc)*length(gyc),
			as.integer(nrow(timesteps)),
				as.double(c(timesteps[,1],pad)),
				as.double(c(timesteps[,2],pad)),
				as.double(c(timesteps[,3],pad+1)),
				as.double(c(timesteps[,4],pad)), 
			as.double(gxc), as.double(gyc),
			as.integer(length(gyc)), as.integer(length(gxc)))

	r <- raster
	values(r) <- c(matrix(cResult, ncol(raster), byrow=T)[,nrow(raster):1])
		
#		if (cutoff.level < 1) {
#			# Remember the original dimensions of the requested grid and where the cut-out box fits in
#			attr(res, 'cutoff') <- list(
#					rows=range(gRows), 
#					cols=range(gCols), 
#					dim=c(length(xc),length(yc)), 
#					dimnames=(if (is.null(grid)) list(xc, yc) else NULL)
#			)
#		}
	
	r@crs <- raster@crs
	r
}

#"encounterDistribution" <- function(tr, threshold, raster=NULL, timestepSize=60,
#		grid.dim=100, grid.pad=0.2) {
#	tr <- bbFilterNA(tr) # Filter out missing measurements, these break the algorithm
#	
#	if (is.null(raster)) {
#		raster <- .defaultRaster(tr, grid.dim, grid.pad)
#	}
#		
#	if (inherits(raster,"RasterLayer")) {
#		# extract the grid lines from the provided grid
#		ext <- extent(raster)
#		xc <- seq(ext@xmin, ext@xmax, length.out=ncol(raster))
#		yc <- seq(ext@ymin, ext@ymax, length.out=nrow(raster))
#	} else {
#		stop("raster must be an instance of 'RasterLayer' if set")
#	}	
#	
#	useFields <- c("x","y","diff.coeff","loc.var","t")
#	
#	timeSteps <- .UDtimesteps(tr, timestepSize)
#	ids <- names(timeSteps)
#	
#	# UDs is two-dimensional list where each entry represents the distribution of encounters for a pair of groups
#	UDs <- as.list(rep(NA, length(ids)^2))
#	dim(UDs) <- c(length(ids), length(ids))
#	rownames(UDs) <- colnames(UDs) <- ids
#	
#	for (id1 in ids) {
#		# It makes no sense to compute encounters between a group and itself, so don't do that
#		for (id2 in ids[-which(ids==id1)]) {
#			# encounters can only be determined if we have a position estimate for both IDs
#			encounterTimes <- intersect(rownames(timeSteps[[id1]]), rownames(timeSteps[[id2]]))
#			
#			weight <- pmin(timeSteps[[id1]][encounterTimes,'weight'], timeSteps[[id2]][encounterTimes,'weight'])
#			
#			ts  <- timeSteps[[id1]][encounterTimes,]
#			ts2 <- timeSteps[[id2]][encounterTimes,]
#			
#			pad <- rep(0, 15)
#			
#			cResult <- .OpenCL("encounterUD", length(xc)*length(yc),
#			as.double(threshold), as.integer(length(encounterTimes)),
#			as.double(c( ts[,1],pad)), as.double(c( ts[,2], pad)), as.double(c( ts[,3],pad+1)), as.double(c(ts[,4], pad)),
#			as.double(c(ts2[,1],pad)), as.double(c(ts2[,2], pad)), as.double(c(sqrt(ts2[,3]),pad+1)),
#			as.double(xc), as.double(yc), as.integer(length(yc)), as.integer(length(xc)))
#			
#			r <- raster
#			values(r) <- c(matrix(cResult, ncol(raster), byrow=T)[,nrow(raster):1])
#			
#			UDs[[id1,id2]] <- r
#		}
#	}
#	diag(UDs) <- unstack(utilizationDistribution(tr, timestepSize=timestepSize, raster=raster))
#	return(UDs)
#}

".UDtimesteps" <- function(tr, timestepSize=60) {
	ts <- as.double(tr@timestamps)
	data <- t(cbind(tr@coords, tr@variance, ts))
				
	# Align the time of each step to a multiple of timestepSize
	timeLimits <- timestepSize*round(ts[c(1,nrow(tr))]/timestepSize)
	nsteps <- as.integer((timeLimits[2] - timeLimits[1]) / timestepSize) + 1
	timeStamps <- seq(timeLimits[1], timeLimits[2], length.out=nsteps)

	timeSteps <- t(.Call("UDTimesteps", data, tr@diffusion, timeStamps))

	dimnames(timeSteps) <- list(timeStamps, c('mu.x', 'mu.y', 'var', 'weight'))
	
	return(timeSteps)
}

".defaultRaster" <- function(ext, grid.dim, padding, proj4string=NA) {
	# TODO: at a later stage, we could possibly remove this in favour of the
	# implementation in the Move package
	padding <- rep(padding, length.out=4)

	xpad <- padding[c(2,4)] * (ext@xmax-ext@xmin)
	xr <- c(ext@xmin-xpad[2], ext@xmax+xpad[1])
	
	ypad <- padding[c(1,3)] * (ext@ymax-ext@ymin)
	yr <- c(ext@ymin-ypad[2], ext@ymax+ypad[1])

	yRange <- diff(yr)
	xRange <- diff(xr)
	raster <- min(xRange, yRange) / grid.dim
	
	# calculation of the coordinates to fit squared raster cells
	ymin <- ext@ymin - (ceiling(yRange/raster) * raster - yRange)/2
	ymax <- ext@ymax + (ceiling(yRange/raster) * raster - yRange)/2
	xmin <- ext@xmin - (ceiling(xRange/raster) * raster - xRange)/2
	xmax <- ext@xmax + (ceiling(xRange/raster) * raster - xRange)/2
	
	nrow <- ((ymax - ymin)/raster)
	ncol <- ((xmax - xmin)/raster)
	raster(extent(xmin, xmax, ymin, ymax), nrow, ncol, crs=proj4string)
}
