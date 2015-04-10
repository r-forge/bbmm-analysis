contour.utilizationDistribution <- function (x, levels=0.99, col=hcl(1:length(x) * 360/length(x), l=35),
		xlim=NULL, ylim=NULL, labels=levels, add=FALSE, ...) {
	if (inherits(x[[1]], "asc")) {
		coords <- sapply(x, adehabitat::getXYcoords)
		xc <- coords['x',1][[1]]
		yc <- coords['y',1][[1]]
	} else {
		xc <- as.double(rownames(x[[1]]))
		yc <- as.double(colnames(x[[1]]))
	}
	
	if (is.null(xlim)) { xlim <- range(xc) }
	if (is.null(ylim)) { ylim <- range(yc) }
		
	if (!inherits(x, "utilizationDistribution"))
        stop("x should be an object of class utilizationDistribution")
	
	if (add != TRUE) {
		plot(1, 1, type = "n", xlim = xlim, ylim = ylim, ...)
	}
	
	for (i in 1:length(x)) {
		# Find the requested levels
		alphaLevel <- .UD.alphaLevels(x[[i]])

		contour(xc, yc, alphaLevel, levels=levels,
				labels=labels, add=TRUE, col=col[i], ...)
	}
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
