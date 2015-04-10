"speedDistribution" <- function(tr, raster=NULL, timestepSize=60, time.scale=timestepSize,
		grid.dim=100, grid.pad=0.2) {
	tr <- na.omit(tr) # Filter out missing measurements, these break the algorithm
	
	# Convert date representations to numeric values
	for (i in 1:length(tr)) {
		burst <- tr[[i]]
		burst$date <- as.double(burst$date)
		tr[[i]] <- burst
	}
		
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
	
	if (is.null(xc) || is.null(yc)) {
		g <- .defaultGrid(tr, grid.dim, grid.pad)
		xc <- g$x
		yc <- g$y
	}
	
	if (time.scale <= 0) {
		stop("time.scale must be positive")
	}
	
	result <- lapply(split(tr), function(tr.id) {
		r <- matrix(0, length(xc), length(yc))
		cResult <- list(r, r)
		for (burst in tr.id) {
			burst <- burst[,c("x","y","diff.coeff","loc.var","date")]
			data <- c(t(burst)) # flatten into vector in row-major order

			# First, compute average speeds over intervals starting at the fixed points
			cResult <- .C("speedDistribution", 
					as.double(cResult[[1]]), as.double(cResult[[2]]),
					as.double(data), as.integer(nrow(burst)),
					as.double(xc), as.double(yc), as.integer(length(xc)), as.integer(length(yc)),
					as.double(timestepSize), as.double(0), as.double(time.scale))
			# Then, compute average speeds over intervals ending at the fixed points
			cResult <- .C("speedDistribution", 
					as.double(cResult[[1]]), as.double(cResult[[2]]),
					as.double(data), as.integer(nrow(burst)),
					as.double(xc), as.double(yc), as.integer(length(xc)), as.integer(length(yc)),
					as.double(timestepSize), as.double(-time.scale), as.double(0))
		}
		values <- matrix(cResult[[1]] / cResult[[2]], ncol(raster))
		r <- raster
		values(r) <- c(values[,nrow(raster):1])
		return(r)
	})
	stack(result)
}

