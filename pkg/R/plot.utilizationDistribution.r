plot.utilizationDistribution <- function (ud, col=hcl(1:length(ud) * 360/length(ud), 50, 70),
		xlim=NULL, ylim=NULL, add=FALSE, ...) {
	if (length(ud) == 0) { return(NULL) }
		
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

	if (!is.list(col)) {
		col <- lapply(as.data.frame(col2rgb(col)), function(col) {
				rgb(col[1], col[2], col[3], sqrt(0:(255^2)), maxColorValue=255)
		})
	}
	# recycle the provided colours if necessary
	col <- rep(col, length.out=length(ud))
	
	if (add != TRUE) {
		plot(1, 1, type = "n", xlim=xlim, ylim = ylim, ...)
	}
	for (i in 1:length(ud)) {
		u <- ud[[i]]
		if (max(u) > 0) {
			image(xc, yc, unclass(u), col=col[[i]], add=T)
		}
	}
}
