setGeneric("summary")
setMethod("summary", 
          signature=".MoveTrackSingle", 
          definition=function(object){
            if (!require(circular)) 
              stop("You need to install the circular package to proceed") #for angle
            if (!isLonLat(object)) {
              object <- spTransform(object, CRSobj="+proj=longlat")
              warning("The projeciton of the object was changed to longlat inside this function")
            }
            data.frame(distanceSummary(object), timeSummary(object), speedSummary(object), angleSummary(object))
          })
          
setMethod("summary", 
          signature=".MoveTrackStack", 
          definition=function(object){
            do.call(rbind, lapply(split(object), summary))
          })
          
setGeneric("dataSummary", function(object, ...){standardGeneric("dataSummary")})
setMethod("dataSummary",
		signature=".MoveTrackStack",
		definition=function(object, ..., units.direction="radian") {
			so <- split(object)
			djl <- djl(so)
			dist <- displacement.distance(so)
			dir <- displacement.direction(so, units.direction)
			
			data.frame(object@idData,
					nobs=sapply(so, nrow),
					DJL=djl,
					displ.dist=dist,
					displ.dir=dir)
})

"djl" <- function(tr) {
	res <- sapply(tr, function(burst) {
		if (nrow(burst) > 0) {
			# Compute the distance between each pair of coordinates, then sum these
			sum(apply(burst@coords[-1,] - burst@coords[-nrow(burst),], 1,
				function(d) { sqrt(sum(d^2)) }))
		} else {
			NA
		}
	})
	names(res) <- names(tr)
	res
}

"displacement.distance" <- function(tr) {
	res <- sapply(na.omit(tr, drop.bursts=FALSE), function(burst) {
		if (nrow(burst) > 0) {
			sqrt(sum((burst@coords[nrow(burst),] - burst@coords[1,])^2))
		} else {
			NA
		}
	})
	names(res) <- names(tr)
	res
}

"displacement.direction" <- function(tr, units=c("radian","degree","compass")) {
	units <- match.arg(units)
	res <- sapply(na.omit(tr, drop.bursts=FALSE), function(burst) {
		if (nrow(burst) > 0) {
			atan2(burst@coords[nrow(burst),2] - burst@coords[1,2],
				burst@coords[nrow(burst),1] - burst@coords[1,1])
		} else {
			NA
		}
	})
	if (units == "degree" || units == "compass") {
		res <- res * 180 / pi
	}
	if (units == "compass") {
		# A compass has reversed (clockwise) direction and is rotated by 90 degrees
		# Rotate by 450 degrees to make sure every value is positive before the mod operation
		res <- (450 - res) %% 360
	}
	names(res) <- names(tr)
	res
}
