"levyUD" <- function(tr, raster=NULL, timestepSize=60,
		grid.dim=100, grid.pad=0.2, byburst=FALSE) {
	tr <- na.omit(tr) # Filter out missing measurements, these break the algorithm
	
	if (is.null(raster)) {
		raster <- .defaultRaster(tr, grid.dim, grid.pad)
	}
		
	if (inherits(raster,"RasterLayer")) {
		# extract the grid lines from the provided grid
		ext <- extent(raster)
		xc <- seq(ext@xmin, ext@xmax, length.out=ncol(raster))
		yc <- seq(ext@ymin, ext@ymax, length.out=nrow(raster))
	} else {
		stop("raster must be an instance of 'RasterLayer' if set")
	}
		
	timeSteps <- .UDlevySteps(tr, timestepSize, byburst)
	
	UDs <- list()
	# For each timestep, compute the mean location, variance, and a weighing factor
	for (id in names(timeSteps)) {
		ts <- timeSteps[[id]]
		# Some rows represent a normally distributed variable, treat these separately
		ts <- lapply(split(as.data.frame(ts), is.infinite(ts[,'scaleFactor'])),
			as.matrix)
		nts <- ts$`TRUE`
		
		pad <- rep(0, 15)
		nResult <- .OpenCL("utilizationDistribution", length(xc)*length(yc),
				as.integer(nrow(nts)),
				as.double(c(nts[,'c3x'],pad)), as.double(c(nts[,'c3y'],pad)), as.double(c(nts[,'c2'],pad+1)), as.double(c(nts[,'weight'],pad)), 
				as.double(xc), as.double(yc), as.integer(length(yc)), as.integer(length(xc)))
		
		# Compute the distribution for the cauchy-ratio-distributed time-steps
		cResult <- .OpenCL("UDLevy", length(xc)*length(yc),
				as.integer(nrow(ts$`FALSE`)),
				as.double(ts$`FALSE`), 
				as.double(xc), as.double(yc), as.integer(length(yc)), as.integer(length(xc)))

		r <- raster
		values(r) <- c(matrix(cResult, ncol(raster), byrow=T)[,nrow(raster):1])
		
		UDs[[id]] <- r
	}
	
	UDs <- stack(UDs)
	if (byburst) {
		attr(UDs, "id") <- attr(timeSteps, "id")
	}
	return(UDs)
}

".UDlevySteps" <- function(tr, timestepSize=60, byburst=FALSE) {
	tr <- na.omit(tr)
	tr.split <- split(tr, sapply(tr, attr, if(byburst) { 'burst' } else { 'id' }))

	result <- lapply(tr.split, function (tr.name) {
		bursts <- tr.name[sapply(tr.name, nrow) >= 2]
		
		prevRows <- 0 # Count the number of relocations in previous bursts
		do.call("rbind", lapply(bursts, function(burst) {
				burst$t <- as.double(burst$date)
				burst <- burst[,c("x","y","diff.coeff","loc.var","t")]
				data <- c(t(burst)) # flatten into vector in row-major order
				
				# Align the time of each step to a multiple of timestepSize
				timeLimits <- timestepSize*round(burst$t[c(1,nrow(burst))]/timestepSize)
				nsteps <- as.integer((timeLimits[2] - timeLimits[1]) / timestepSize) + 1
				
				ts <- seq(timeLimits[1], timeLimits[2], length=nsteps)
				is <- prevRows + sapply(ts, function (t) {
					max(c(1, which(burst$t <= t)))
				})
				
				prevRows <<- prevRows + nrow(burst)
				return(data.frame(time=ts, index=is-1)) # C++ has 0-based indices
		}))
	})
	
	if (byburst) {
		attr(result, 'id') <- as.factor(sapply(tr.split, function (t) {
			attr(t[[1]], 'id')
		}))
	}
	
	return(result)
}
'LTS' <- .UDlevySteps
