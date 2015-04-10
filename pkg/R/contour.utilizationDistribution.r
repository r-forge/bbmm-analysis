contour.utilizationDistribution <- function (x, levels=0.99, col=hcl(1:length(x) * 360/length(x), l=35),
		xlim=NULL, ylim=NULL, labels=levels, add=FALSE, ...) {
	coords <- .UD.coords(x)
	
	if (is.null(xlim)) { xlim <- range(coords$x) }
	if (is.null(ylim)) { ylim <- range(coords$y) }
		
	if (!inherits(x, "utilizationDistribution"))
        stop("x should be an object of class utilizationDistribution")
	
	if (add != TRUE) {
		plot(1, 1, type = "n", xlim = xlim, ylim = ylim, ...)
	}
	
	for (i in 1:length(x)) {
		# Find the requested levels
		alphaLevel <- .UD.alphaLevels(x[[i]])

		contour(coords$x, coords$y, alphaLevel, levels=levels,
				labels=labels, add=TRUE, col=col[i], ...)
	}
}

contourPolygons <- function(ud, levels=0.99) {
	if (!inherits(ud, "utilizationDistribution")) {
		stop("ud must be of class utilizationDistribution")
	}
	
	coords <- .UD.coords(ud)
	lapply(ud, function(r) {
		lines <- contourLines(coords$x, coords$y, .UD.alphaLevels(r), levels=levels)
		# Get a separate list of line segments for each level
		lines <- split(lines, as.factor(sapply(lines, "[[", "level")))
		
		# The result is a SpatialPolygons object with one set of Polygons for each level
		SpatialPolygons(lapply(lines, function(ls) {
			ps <- Polygons(lapply(ls, function(c) {
				Polygon(cbind(c$x, c$y))
			}), ls[[1]]$level)
			maptools::checkPolygonsHoles(ps)
		}))
	})
}

# Given a matrix representing a UD, compute for each cell the minimum alpha level for which it would be inside the contour
".UD.alphaLevels" <- function(ud) {
	# Find the requested levels
	# Flatten the matrix into a vector and label each element with its current index
	v <- c(ud)
	names(v) <- 1:length(v)
	
	# Determine the minimum alpha level at which it is inside the contour for each cell
	v <- sort(v, decreasing=TRUE)
	alphaLevel <- cumsum(v)/sum(v)
	alphaLevel <- alphaLevel[order(as.integer(names(v)))]
	
	# Reassemble the matrix of alpha values
	matrix(alphaLevel, nrow(ud), ncol(ud))
}

".UD.coords" <- function(ud) {
	if (inherits(ud[[1]], "asc")) {
		return(adehabitat::getXYcoords(ud[[1]]))
	} else {
		xc <- as.double(rownames(ud[[1]]))
		yc <- as.double(colnames(ud[[1]]))
		
		if (!is.null(xc) && !is.null(yc)) {
			return(list(x=xc, y=yc))
		}
	}
	
	list(x=seq(0,1,length.out=nrow(ud[[1]])), y=seq(0,1,length.out=ncol(ud[[1]])))
}
