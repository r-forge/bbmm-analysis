#visibility <- function(DEM, maxDist=Inf) {
#	xc <- seq(attr(DEM, 'xll'), by=attr(DEM, 'cellsize'), length.out=dim(DEM)[1])
#	yc <- seq(attr(DEM, 'yll'), by=attr(DEM, 'cellsize'), length.out=dim(DEM)[2])
#	
#	result <- as.list(rep(NA, prod(dim(DEM))))
#	dim(result) <- dim(DEM)
#	
#	for (i in 1:length(xc)) {
#		iNear <- which(abs(xc[i]-xc) <= maxDist)
#		for (j in 1:length(yc)) {
#			jNear <- which(abs(yc[j]-yc) <= maxDist)
#			xy <- c(xc[i], yc[j])
#			
#			visible <- data.frame(x=NULL, y=NULL)
#			for (ni in iNear) {
#				for (nj in jNear) {
#					xyN <- c(xc[ni], yc[nj])
#				
#					# xNear and yNear define a bounding box for the relevant locations
#					# now restrict ourselves to the circle that is actually in range
#					if (sqrt(sum((xy - xyN)^2)) <= maxDist && visible(DEM, c(i,j), c(ni, nj), xc, yc)) {
#						visible <- rbind(visible, list(x=ni, y=nj))
#					}
#				}
#			}
#			
#			# visible is now using indices, convert to coordinates
#			visible$x <- xc[visible$x]
#			visible$y <- yc[visible$y]
#			result[[i,j]] <- visible
#		}
#	}
#	
#	dimnames(result) <- list(xc, yc)
#	result
#}

#visible <- function(DEM, ij1, ij2,
#		xc=seq(attr(DEM, 'xll'), by=attr(DEM, 'cellsize'), length.out=dim(DEM)[1]),
#		yc=seq(attr(DEM, 'yll'), by=attr(DEM, 'cellsize'), length.out=dim(DEM)[2])) {
#	cellSize <- attr(DEM, 'cellsize')
#	
#	# Walk along the x axis and investigate crossings with x=c edges
#	x1 <- xc[ij1[1]]
#	y1 <- yc[ij1[2]]
#	z1 <- DEM[ij1[1], ij1[2]]
#	x2 <- xc[ij2[1]]
#	y2 <- yc[ij2[2]]
#	z2 <- DEM[ij2[1], ij2[2]]
#	
#	dx <- x2 - x1
#	dy <- y2 - y1
#	dz <- z2 - z1

#	# Walk along the x axis and investigate crossings with vertical edges	
#	if (abs(ij1[1] - ij2[1]) > 1) { # If there is some x extent to walk along
#		yi <- ij1[2]
#		direction <- sign(ij2[1] - ij1[1])
#		for (xi in (ij1[1]+direction):(ij2[1]-direction)) {
#			x <- xc[xi]
#			y <- y1 + (x-x1)/dx * dy
#		
#			while (yc[yi] <= y && yi < length(yc)) { yi <- yi+1 }
#		
#			zLoS <- z1 + (x-x1)/dx * dz
#			zDEM <- DEM[xi, yi-1] + (y-yc[yi-1])/cellSize * (DEM[xi, yi] - DEM[xi, yi-1])
#		
#			if (zLoS < zDEM) { # If the LoS exactly touches the terrain, allow it through
#				return(FALSE)
#			}
#		}
#	}
#	
#	# Walk along the y axis and investigate crossings with horizontal edges
#	if (abs(ij1[2] - ij2[2]) > 1) { # If there is some y extent to walk along
#		xi <- ij1[1]
#		direction <- sign(ij2[2] - ij1[2])
#		for (yi in (ij1[2]+direction):(ij2[2]-direction)) {
#			y <- yc[yi]
#			x <- x1 + (y-y1)/dy * dx
#			
#			while (xc[xi] <= x && xi < length(xc)) { xi <- xi+1 }
#		
#			zLoS <- z1 + (y-y1)/dy * dz
#			zDEM <- DEM[xi-1, yi] + (x-xc[xi-1])/cellSize * (DEM[xi, yi] - DEM[xi-1, yi])
#			
#			if (zLoS < zDEM) { # If the LoS exactly touches the terrain, allow it through
#				return(FALSE)
#			}
#		}
#	}
#	
#	return(TRUE)
#}

#visibility_dPos <- function(maxCells) {
#	dPos <- matrix(NA, nrow=2, ncol=0)
#	for (i in (-maxCells:maxCells)) {
#		for (j in (-maxCells:maxCells)) {
#			if (i^2+j^2 <= maxCells^2) {
#				dPos <- cbind(dPos, c(i, j))
#			}
#		}
#	}
#	rownames(dPos) <- c('dx', 'dy')
#	
#	dPos
#}

#visibility_OCL <- function(DEM, maxDist) {
#	# Compose a list of all cell offsets within a radius of maxDist from the origin
#	dPos <- visibility_dPos(floor(maxDist/attr(DEM, 'cellsize')))

#	result <- .OpenCL('DEMVisibility', as.integer(length(DEM)*ncol(dPos)),
#			dim=as.integer(dim(DEM)), as.double(DEM),
#			as.integer(dPos), as.integer(ncol(dPos)))
#	
#	dim(result) <- c(ncol(dPos), dim(DEM))
#	
#	result
#}

#visibility_OCL_convert_results <- function(clResult, DEM, maxDist) {
#	dPos <- visibility_dPos(floor(maxDist/attr(DEM, 'cellsize')))

#	xc <- seq(attr(DEM, 'xll'), by=attr(DEM, 'cellsize'), length.out=dim(DEM)[1])
#	yc <- seq(attr(DEM, 'yll'), by=attr(DEM, 'cellsize'), length.out=dim(DEM)[2])
#	
#	result <- as.list(rep(NA, length(DEM)))
#	dim(result) <- dim(DEM)
#	dimnames(result) <- list(xc, yc)
#	
#	for (i in 1:nrow(DEM)) {
#		for (j in 1:ncol(DEM)) {
#			vis <- data.frame(t(dPos), visibility=clResult[,i,j])
#			vis <- vis[!is.na(vis$visibility) & vis$visibility > 0,]
#			vis$x <- xc[i+vis$dx]
#			vis$y <- yc[j+vis$dy]
#			
#			result[[i,j]] <- vis[,c('visibility','x','y')]
#		}
#	}
#	
#	result
#}

#visibility_fast <- function(DEMs, xc=1:(dim(DEMs)[2]), yc=1:(dim(DEMs)[3])) {
#	cResult <- .Call('visibility', DEMs)
#	result <- lapply(cResult, function(d) {
#		res <- as.data.frame(t(matrix(d, nrow=3, dimnames=list(c('x','y','visibility')))))
#		res$x <- xc[res$x+1]
#		res$y <- yc[res$y+1]
#		res$visibility <- res$visibility / dim(DEMs)[1]
#		res
#	})
#	
#	dim(result) <- dim(DEMs)[2:3]
#	dimnames(result) <- list(xc,yc)
#	
#	result
#}
