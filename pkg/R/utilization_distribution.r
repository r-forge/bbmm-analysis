"utilizationDistribution" <- function(tr, grid=NULL, timestepSize=60, xc=NULL, yc=NULL,
		grid.dim=100, grid.pad=0.2) {
	tr <- bbFilterNA(tr) # Filter out missing measurements, these break the algorithm
		
	if (inherits(grid,"asc")) {
		# extract the grid lines from the provided grid
		xc <- seq(from=attr(grid, "xll"), by=attr(grid, "cellsize"), length.out=attr(grid, "dim")[1])
		yc <- seq(from=attr(grid, "yll"), by=attr(grid, "cellsize"), length.out=attr(grid, "dim")[2])
	} else if (!is.null(grid)) {
		stop("grid must be an instance of 'asc' if set")
	}
	
	if (is.null(xc) || is.null(yc)) {
		g <- .defaultGrid(tr, grid.dim, grid.pad)
		xc <- g$x
		yc <- g$y
	}
	
	timeSteps <- .UDtimesteps(tr, timestepSize)
	
	UDs <- list()
	# For each timestep, compute the mean location, variance, and a weighing factor
	for (id in names(timeSteps)) {
		ts <- timeSteps[[id]]
	
		pad <- rep(0, 15)
			
		cResult <- .OpenCL("utilizationDistribution", length(xc)*length(yc),
				as.integer(nrow(ts)),
				as.double(c(ts[,1],pad)), as.double(c(ts[,2],pad)), as.double(c(ts[,3],pad+1)), as.double(c(ts[,4],pad)), 
				as.double(xc), as.double(yc), as.integer(length(yc)), as.integer(length(xc)))
		
#		cResult <- .OpenCL("getUD", length(xc)*length(yc),
#			as.integer(length(timesteps)/4), as.double(timesteps),
#			as.double(xc), as.double(yc), as.integer(length(yc)), as.integer(length(xc)))

		res <- matrix(cResult, nrow=length(xc), byrow=TRUE)			
		
		if (!is.null(grid)) {
			attributes(res) <- attributes(grid)
		} else {
			rownames(res) <- xc
			colnames(res) <- yc
		}
		
		UDs[[id]] <- res
	}
	
	class(UDs) <- c("utilizationDistribution", "list")
	return(UDs)
}

"encounterDistribution" <- function(tr, threshold, grid=NULL, timestepSize=60, xc=NULL, yc=NULL,
		grid.dim=100, grid.pad=0.2) {
	tr <- bbFilterNA(tr) # Filter out missing measurements, these break the algorithm
	
	if (inherits(grid,"asc")) {
		# extract the grid lines from the provided grid
		xc <- seq(from=attr(grid, "xll"), by=attr(grid, "cellsize"), length.out=attr(grid, "dim")[1])
		yc <- seq(from=attr(grid, "yll"), by=attr(grid, "cellsize"), length.out=attr(grid, "dim")[2])
	} else if (!is.null(grid)) {
		stop("grid must be an instance of 'asc' if set")
	}
	
	if (is.null(xc) || is.null(yc)) {
		g <- .defaultGrid(tr, grid.dim, grid.pad)
		xc <- g$x
		yc <- g$y
	}	
	
	useFields <- c("x","y","diff.coeff","loc.var","t")
	
	timeSteps <- .UDtimesteps(tr, timestepSize)
	ids <- names(timeSteps)
	
	# UDs is two-dimensional list where each entry represents the distribution of encounters for a pair of groups
	UDs <- as.list(rep(NA, length(ids)^2))
	dim(UDs) <- c(length(ids), length(ids))
	rownames(UDs) <- colnames(UDs) <- ids
	
	for (id1 in ids) {
		# It makes no sense to compute encounters between a group and itself, so don't do that
		for (id2 in ids[-which(ids==id1)]) {
			# encounters can only be determined if we have a position estimate for both IDs
			encounterTimes <- intersect(rownames(timeSteps[[id1]]), rownames(timeSteps[[id2]]))
			
			weight <- pmin(timeSteps[[id1]][encounterTimes,'weight'], timeSteps[[id2]][encounterTimes,'weight'])
			
			ts  <- timeSteps[[id1]][encounterTimes,]
			ts2 <- timeSteps[[id2]][encounterTimes,]
			
			pad <- rep(0, 15)
			
			cResult <- .OpenCL("encounterUD", length(xc)*length(yc),
			as.double(threshold), as.integer(length(encounterTimes)),
			as.double(c( ts[,1],pad)), as.double(c( ts[,2], pad)), as.double(c( ts[,3],pad+1)), as.double(c(ts[,4], pad)),
			as.double(c(ts2[,1],pad)), as.double(c(ts2[,2], pad)), as.double(c(sqrt(ts2[,3]),pad+1)),
			as.double(xc), as.double(yc), as.integer(length(yc)), as.integer(length(xc)))
			
			res <- matrix(cResult, nrow=length(xc), byrow=TRUE)			
		
			if (!is.null(grid)) {
				attributes(res) <- attributes(grid)
			} else {
				rownames(res) <- xc
				colnames(res) <- yc
			}
			
			UDs[[id1,id2]] <- res
		}
	}
	diag(UDs) <- utilizationDistribution(tr, timestepSize=timestepSize, xc=xc, yc=yc)
	class(UDs) <- c("utilizationDistribution", "list")
	return(UDs)
}

".UDtimesteps" <- function(tr, timestepSize=60) {
	ids <- unique(id(tr))
	
	result <- list()
	
	for (id in ids) {
		bursts <- tr[sapply(tr, nrow) >= 2][id=id]
		result[[id]] <- do.call("rbind", lapply(bursts, function(burst) {
				burst$t <- as.double(burst$date)
				burst <- burst[!(is.na(burst$x) || is.na(burst$y) || is.na(burst$loc.var)),
				               c("x","y","diff.coeff","loc.var","t")]
				data <- c(t(burst)) # flatten into vector in row-major order
				
				timeLimits <- timestepSize*round(burst$t[c(1,nrow(burst))]/timestepSize)
				nsteps <- as.integer((timeLimits[2] - timeLimits[1]) / timestepSize) + 1
				
				cResult <- .C("UDTimesteps", double(nsteps*4),
						as.integer(nrow(burst)), as.double(data),
						as.integer(nsteps), timeLimits)
				timeSteps <- matrix(cResult[[1]], ncol=4, byrow=T)
				colnames(timeSteps) <- c('x', 'y', 'var', 'weight')
				rownames(timeSteps) <- seq(timeLimits[1], timeLimits[2], length=nrow(timeSteps))
				return(timeSteps)
		}))
	}
	return(result)
}

".defaultGrid" <- function(tr, grid.dim, padding) {
	xr <- range(unlist(sapply(tr, function(b) { b$x })))
	xpad <- padding * (xr[2]-xr[1])
	xr <- c(xr[1]-xpad, xr[2]+xpad)
	
	yr <- range(unlist(sapply(tr, function(b) { b$y })))
	ypad <- padding * (yr[2]-yr[1])
	yr <- c(yr[1]-ypad, yr[2]+ypad)

	cellSize <- min(xr[2] - xr[1], yr[2]-yr[1]) / grid.dim
	
	print(paste("Using default grid: bounding box for trajectory extended by", padding, "on each side."))
	return(list(
		x=seq(xr[1], xr[2], by=cellSize),
		y=seq(yr[1], yr[2], by=cellSize)
	))
}
