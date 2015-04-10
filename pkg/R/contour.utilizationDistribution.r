contour.utilizationDistribution <- function (ud, levels=0.99, col=hcl(1:length(ud) * 360/length(ud), l=35),
		xlim=NULL, ylim=NULL, labels=levels, add=FALSE, ...) {
	if (inherits(ud[[1]], "asc")) {
		coords <- sapply(ud, getXYcoords)
		xc <- coords['x',1][[1]]
		yc <- coords['y',1][[1]]
	} else {
		xc <- as.double(rownames(ud[[1]]))
		yc <- as.double(colnames(ud[[1]]))
	}
		
	if (is.null(xlim)) { xlim <- range(xc) }
	if (is.null(ylim)) { ylim <- range(yc) }
		
	if (!inherits(ud, "utilizationDistribution"))
        stop("ud should be an object of class utilizationDistribution")
	
	if (add != TRUE) {
		plot(1, 1, type = "n", xlim = xlim, ylim = ylim, ...)
	}
	
	for (i in 1:length(ud)) {
		# Find the requested levels
		# Flatten the matrix into a vector and label each element with its current index
		v <- c(ud[[i]])
		names(v) <- 1:length(v)
	
		# Determine the minimum alpha level at which it is inside the contour for each cell
		v <- sort(v, decreasing=TRUE)
		alphaLevel <- cumsum(v)/sum(v)
		alphaLevel <- alphaLevel[order(as.integer(names(v)))]
	
		# Reassemble the matrix of alpha values
		alphaLevel <- matrix(alphaLevel, nrow(ud[[i]]), ncol(ud[[i]]))

		#image(x=as.double(colnames(x[[i]])), y=as.double(rownames(x[[i]])),t(alphaLevel), col=rgb(1,0,0,255:0/255))
		contour(xc, yc, alphaLevel, levels=levels,
				labels=labels, add=TRUE, col=col[i], ...)
	}
}
