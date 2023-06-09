# TODO: We should probably use S4 generics
# TODO: can we somehow connect this to move::getVolumeUD?
contour.ud <- function (x, levels=0.99, col=hcl(1:nlayers(x) * 360/nlayers(x), l=35),
		labels=levels, add=FALSE, ...) {
	
	#if (!inherits(x, "utilizationDistribution"))
    #    stop("x should be an object of class utilizationDistribution")
	
	for (i in 1:nlayers(x)) {
		# Find the requested levels
		alphaLevel <- .UD.alphaLevels(x[[i]])
		
		contour(alphaLevel, levels=levels,
				labels=labels, add=(add | i > 1), col=col[i], ...)
	}
}

setGeneric("contourPolygons", function(ud, levels=0.99, ...) standardGeneric("contourPolygons"))
setMethod(f="contourPolygons",
	signature=c(levels="missing"), 
	definition = function(ud, levels, ...) {

	contourPolygons(ud, levels=0.99, ...)
})

setMethod(f="contourPolygons",
	signature=c(ud="RasterLayer", levels="numeric"), 
	definition = function(ud, levels, ...) {

	coords <- list(
			x=seq(extent(ud)@xmin, extent(ud)@xmax, length.out=ncol(ud)),
			y=seq(extent(ud)@ymin, extent(ud)@ymax, length.out=nrow(ud))
	)
	
	## Surround the working matrix by a border of values outside any contour
	## This avoids weird effects when the countour touches the raster border
	v <- matrix(2, ncol(ud)+2, nrow(ud)+2)
	v[2:(ncol(ud)+1),(nrow(ud)+1):2] <- matrix(values(.UD.alphaLevels(ud)), ncol(ud))
	
	## Make up coordinates for the border rows/columns
	coords <- lapply(coords, function (cs) {
		delta <- cs[2]-cs[1]
		c(cs[1]-delta, cs, cs[length(cs)]+delta)
	})
	
	lines <- contourLines(coords$x, coords$y, v, levels=levels)
	# Get a separate list of line segments for each level
	lines <- split(lines, as.factor(sapply(lines, "[[", "level")))
		
	# The result is a SpatialPolygons object with one set of Polygons for each level
	SpatialPolygons(lapply(lines, function(ls) {
		ps <- Polygons(lapply(ls, function(l) {
			# Repeat the first coordinate to make sure the ring is closed
			Polygon(rbind(cbind(l$x, l$y), c(l$x[1], l$y[1])))
		}), ls[[1]]$level)
		maptools::checkPolygonsHoles(ps)
	}), proj4string=ud@crs)
})

setMethod(f="contourPolygons",
	signature=c(ud="RasterStack", levels="numeric"), 
	definition = function(ud, levels, ...) {

	cp <- lapply(unstack(ud), "contourPolygons", levels, ...)
	names(cp) <- names(ud)
	cp
})

# Given a RasterLayer representing a UD, compute for each cell the minimum alpha level for which it would be inside the contour
".UD.alphaLevels" <- function(ud) {
	# Find the requested levels
	# Flatten the matrix into a vector and label each element with its current index
	v <- values(ud)
	names(v) <- 1:length(v)
	
	# Determine the minimum alpha level at which it is inside the contour for each cell
	v <- sort(v, decreasing=TRUE)
	alphaLevel <- cumsum(v)/sum(v)
	alphaLevel <- alphaLevel[order(as.integer(names(v)))]
	
	# Reassemble the matrix of alpha values
	values(ud) <- alphaLevel
	ud
}

