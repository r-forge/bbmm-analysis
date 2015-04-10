barplot.encounterDuration <- function(height, col=NULL, units="auto", ...) {
	if (is.null(col)) {
		nColours <- ifelse(is.data.frame(height), ncol(height)-4, dim(height)[3])
		col <- hcl(1:nColours * 360/nColours)
	}
		
	if (is.data.frame(height)) {
		# duration by burst
		
		# The first 2 columns do not contain values for movement models
		models <- colnames(height)[3:length(colnames(height))]
		
		if (!any(is.na(suppressWarnings(as.numeric(rownames(height)))))) {
			# All row names are numbers; override the row names using the burst names
			rownames(height) <- paste(height$burst1, height$burst2, sep=" - ")
		}

		plotMatrix <- t(sapply(height[,models], function(d) { as.numeric(d, units=units) }))
		colnames(plotMatrix) <- rownames(height)
	} else {
		# duration by id, we need to convert to a matrix suitable for plotting
		d <- dim(height)
		IDs <- dimnames(height)[[1]]
		
		plotMatrix <- matrix(nrow=d[3], ncol=d[1]*(d[1]-1)/2)
		rownames(plotMatrix) <- dimnames(height)[[3]]

		c <- rep("Unknown", ncol(plotMatrix))

		k <- 0
		for (i in 1:(d[1]-1)) {
			id1 <- IDs[i]
			for (j in (i+1):d[1]) {
				id2 <- IDs[j]
				
				k <- k+1
				c[k] <- paste(id1, id2, sep="-");
				plotMatrix[,k] <- as.numeric(height[id1, id2,], units=units)
			}
		}
		colnames(plotMatrix) <- c
	}
	
	barplot(plotMatrix, beside=TRUE, legend=rownames(plotMatrix), col=col, ...)
}
