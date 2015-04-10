barplot.encounterDuration <- function(ed,
		col=rainbow(ifelse(is.data.frame(ed), ncol(ed)-4, dim(ed)[3])),
		...) {
	if (is.data.frame(ed)) {
		# duration by burst
		
		# The first 4 columns do not contain values for models
		models <- colnames(ed)[5:length(colnames(ed))]
		
		if (!any(is.na(as.numeric(rownames(ed))))) {
			# All row names are numbers; override the row names using the burst names
			rownames(ed) <- paste(ed$burst1, ed$burst2, sep=" - ")
		}
		
		plotMatrix <- t(ed[,models])
	} else {
		# duration by id, we need to convert to a matrix suitable for plotting
		d <- dim(ed)
		IDs <- dimnames(ed)[[1]]
		
		plotMatrix <- matrix(nrow=d[3], ncol=d[1]*(d[1]-1)/2)
		rownames(plotMatrix) <- dimnames(ed)[[3]]

		c <- rep("Unknown", ncol(plotMatrix))

		k <- 0
		for (i in 1:(d[1]-1)) {
			id1 <- IDs[i]
			for (j in (i+1):d[1]) {
				id2 <- IDs[j]
				
				k <- k+1
				c[k] <- paste(id1, id2, sep="-");
				plotMatrix[,k] <- ed[id1, id2,]
			}
		}
		colnames(plotMatrix) <- c
	}
	
	barplot(plotMatrix, beside=TRUE, legend=rownames(plotMatrix), col=col, ...)
}
