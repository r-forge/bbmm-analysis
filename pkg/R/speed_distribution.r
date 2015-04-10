setGeneric("speedDistribution", function(object, raster=NULL, timestepSize=60,
		time.scale=timestepSize, grid.dim=100, grid.pad=0.2, groupBy=NULL) {
	standardGeneric("speedDistribution")
})

setMethod(f = "speedDistribution", 
		signature = c(object = ".MoveTrack", raster="missing"),
		function(object, timestepSize=60, time.scale=timestepSize, 
			grid.dim=100, grid.pad=0.2, groupBy=NULL) {
	if (!is.numeric(grid.dim) || !is.numeric(grid.pad)) {
		stop("grid.dim and grid.pad must be numeric")
	}
	raster <- .defaultRaster(extent(object), grid.dim, grid.pad, object@proj4string)
	speedDistribution(object=object, raster=raster,
			timestepSize=timestepSize, time.scale=time.scale,
			groupBy=groupBy)
})

setMethod(f = "speedDistribution", 
		signature = c(object = ".MoveTrack", raster="RasterLayer", timestepSize="numeric", time.scale="missing"),
		function(object, raster, timestepSize=60, groupBy=NULL) {
	speedDistribution(object, raster, timestepSize, time.scale, groupBy)
})


setMethod(f = "speedDistribution", 
		signature = c(object = "MoveBB", raster = "RasterLayer", timestepSize = "numeric", time.scale="numeric"), 
		function(object, raster, timestepSize=60, time.scale) {
	r <- .speedDistribution(object, raster, timestepSize, time.scale)
	
	# r contains the sum of expected speeds and total weight for each cell;
	# Divide to obtain expected speed
	# Raster has flipped coordinates?
	values(raster) <- c((r[[1]]/r[[2]])[,ncol(r[[1]]):1])
	raster
})

setMethod(f = "speedDistribution", 
		signature = c(object = "MoveBBStack", raster = "RasterLayer", timestepSize = "numeric", time.scale="numeric"), 
		function(object, raster, timestepSize=60, time.scale, groupBy=NULL) {
	os <- split(split(object), IDs(object, groupBy))

	# Group the results if necessary
	SDs <- lapply(os, function(tr) {
		# Create a list of two zeroed matrices to add up the results in
		u <- rep(list(matrix(0, ncol(raster), nrow(raster))), 2)
		for (burst in tr) {
			spdist <- .speedDistribution(burst, raster, timestepSize, time.scale)
			u <- mapply("+", u, spdist, SIMPLIFY=FALSE)
		}
		u
	})

	# Normalize all UDs, then create stack and return
	stack(lapply(SDs, function(r) {
		# r contains the sum of expected speeds and total weight for each cell;
		# Divide to obtain expected speed.
		# Raster has flipped coordinates?
		values(raster) <- c((r[[1]]/r[[2]])[,ncol(r[[1]]):1])
		raster
	}))
})

".speedDistribution" <- function(tr, raster, timestepSize, time.scale) {
	# extract the grid lines from the provided grid
	ext <- extent(raster)
	xc <- seq(ext@xmin, ext@xmax, length.out=ncol(raster))
	yc <- seq(ext@ymin, ext@ymax, length.out=nrow(raster))
	
	# First, compute average speeds over intervals starting at the fixed points
	res <- .Call("speedDistribution",
			list(coords=t(tr@coords),
				ts=as.double(tr@timestamps),
				var=tr@variance,
				diff=tr@diffusion
			),
			list(X=xc, Y=yc),
			as.double(timestepSize), as.double(c(0.0, time.scale))
	)
	
	# Then, compute average speeds over intervals ending at the fixed points
	cResult <- .Call("speedDistribution",
			list(coords=t(tr@coords),
				ts=as.double(tr@timestamps),
				var=tr@variance,
				diff=tr@diffusion
			),
			list(X=xc, Y=yc),
			as.double(timestepSize), as.double(c(0.0, time.scale))
	)
	# Sum the results (sum of speeds and total weights) and return
	mapply("+", res, cResult, SIMPLIFY=FALSE)
}

