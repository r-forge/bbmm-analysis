plot.utilizationDistribution <- function (x, col=hcl(1:length(x) * 360/length(x), 50, 70),
		xlim=NULL, ylim=NULL, add=FALSE, ...) {
	if (length(x) == 0) { return(NULL) }
		
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

	if (!is.list(col)) {
		col <- lapply(as.data.frame(col2rgb(col)), function(col) {
				rgb(col[1], col[2], col[3], sqrt(0:(255^2)), maxColorValue=255)
		})
	}
	# recycle the provided colours if necessary
	col <- rep(col, length.out=length(x))
	
	if (add != TRUE) {
		plot(1, 1, type = "n", xlim=xlim, ylim = ylim, ...)
	}
	for (i in 1:length(x)) {
		u <- x[[i]]
		if (max(u, na.rm=TRUE) > 0) {
			image(xc, yc, unclass(u), col=col[[i]], add=T)
		}
	}
}
