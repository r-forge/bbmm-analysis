plot.utilizationDistribution <- function (x, col=hcl(1:length(x) * 360/length(x), 50, 70),
		xlim=NULL, ylim=NULL, add=FALSE, plot.range=sapply(x, range), ...) {
	if (length(x) == 0) { return(NULL) }
		
	if (inherits(x[[1]], "asc")) {
		coords <- sapply(x, adehabitat::getXYcoords)
		xc <- coords['x',1][[1]]
		yc <- coords['y',1][[1]]
	} else {
		xc <- as.double(rownames(x[[1]]))
		yc <- as.double(colnames(x[[1]]))
	}
	if (length(xc) == 0) { xc <- 1:nrow(x[[1]]) }
	if (length(yc) == 0) { yc <- 1:ncol(x[[1]]) }
		
	if (is.null(xlim)) { xlim <- range(xc) }
	if (is.null(ylim)) { ylim <- range(yc) }

	if (!inherits(x, "utilizationDistribution"))
        stop("x should be an object of class utilizationDistribution")

	if (!is.list(col)) {
		col <- as.list(col)
	}
	# recycle the provided colours if necessary
	col <- rep(col, length.out=length(x))
	
	if (add != TRUE) {
		# Define the dimensions of the canvas
		plot(1, 1, type = "n", xlim=xlim, ylim = ylim, ...)
	}
	for (i in 1:length(x)) {
		u <- x[[i]]
		
		if (is.matrix(plot.range)) {
			pr <- plot.range[,i]
		} else {
			pr <- plot.range
		}
		
		if (max(u, na.rm=TRUE) > 0) {
			# Normalize the values before plotting
			u <- (u-pr[1]) / (pr[2] - pr[1])
			# rasterImage seems to assume another matrix orientation than image
			u <- t(u)[ncol(u):1,]
			if (typeof(col[[i]])=="builtin") {
				# A function was provided that maps the colours
				stop("Sorry, functions as colour maps are not yet supported")
			} else if (length(col[[i]]) > 1) {
				# A list of colours was provided, map to them
				img <- ifelse(is.na(u),
							"#FFFFFF00", # Transparent white for missing values
							col[[i]][ceiling(u*length(col[[i]]))])
			} else {
				d <- dim(u)
				img <- array(c(
						rep(col2rgb(col[[i]])/255, each=prod(d)), u), dim=c(d,4))
			}

			xyc <- if(inherits(x[[i]], 'asc')) {
					getXYcoords(x[[i]])
				} else {
					lapply(dimnames(x[[i]]), as.double)
				}
			names(xyc) <- c('x','y')
			xyr <- sapply(xyc, range)
			
			rasterImage(img, xyr[1,'x'], xyr[1,'y'], xyr[2,'x'], xyr[2,'y'])
		}
	}
}
