"speedDistribution" <- function(tr, grid=NULL, timestepSize=60, time.scale=timestepSize,
		xc=NULL, yc=NULL, grid.dim=100, grid.pad=0.2) {
	tr <- na.omit(tr) # Filter out missing measurements, these break the algorithm
	
	# Convert date representations to numeric values
	for (i in 1:length(tr)) {
		burst <- tr[[i]]
		burst$date <- as.double(burst$date)
		tr[[i]] <- burst
	}
		
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
	
	if (time.scale <= 0) {
		stop("time.scale must be positive")
	}
	
	lapply(split(tr), function(tr.id) {
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
		matrix(cResult[[1]] / cResult[[2]], length(xc), length(yc))
	})
}

