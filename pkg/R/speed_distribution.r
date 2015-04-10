"speedDistribution" <- function(tr, raster=NULL, timestepSize=60, time.scale=timestepSize,
		grid.dim=100, grid.pad=0.2, groupBy=NULL) {
			
	if (is.null(raster)) {
		raster <- .defaultRaster(extent(tr), grid.dim, grid.pad)
	}
		
	if (inherits(raster,"RasterLayer")) {
		# extract the grid lines from the provided grid
		ext <- extent(raster)
		xc <- seq(ext@xmin, ext@xmax, length.out=ncol(raster))
		yc <- seq(ext@ymin, ext@ymax, length.out=nrow(raster))
	} else {
		stop("raster must be an instance of 'RasterLayer' if set")
	}
	
	if (time.scale <= 0) {
		stop("time.scale must be positive")
	}
	
	# Split the trajectory into bursts and find out how they should be grouped
	trs <- split(tr)
	resultNames <- names(trs)
	burstNames <- resultNames
	if (!is.null(groupBy)) {
		resultNames <- unique(tr@idData[,groupBy])
		burstNames <- as.character(sapply(names(trs), function (n) { tr@idData[n, groupBy] }))
	}

	emptyResult <- matrix(0, length(xc), length(yc))
	result <- lapply(resultNames, function(id) {
		res <- list(emptyResult, emptyResult)
		for (burst in trs[burstNames==id]) {
			# First, compute average speeds over intervals starting at the fixed points
			cResult <- .Call("speedDistribution",
					list(coords=t(burst@coords),
						ts=as.double(burst@timestamps),
						var=burst@variance,
						diff=burst@diffusion
					),
					list(X=xc, Y=yc),
					as.double(timestepSize), c(0.0, time.scale)
			)
			res <- mapply("+", res, cResult, SIMPLIFY=FALSE)
			
			# Then, compute average speeds over intervals ending at the fixed points
			cResult <- .Call("speedDistribution",
					list(coords=t(burst@coords),
						ts=as.double(burst@timestamps),
						var=burst@variance,
						diff=burst@diffusion
					),
					list(X=xc, Y=yc),
					as.double(timestepSize), c(0.0, time.scale)
			)
			res <- mapply("+", res, cResult, SIMPLIFY=FALSE)
		}
		values <- matrix(res[[1]] / res[[2]], ncol(raster))
		r <- raster
		values(r) <- c(values[,nrow(raster):1])
		return(r)
	})
	names(result) <- resultNames
	stack(result)
}

