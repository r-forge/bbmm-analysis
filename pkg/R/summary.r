setGeneric("summary")
setMethod("summary", 
          signature=".MoveTrackSingle", 
          definition=function(object){
            
	data <- data.frame(distanceSummary(object),
				timeSummary(object),
				speedSummary(object),
				angleSummary(object),
				displ.dist=displacement.distance(object),
				displ.dir=displacement.direction(object),
				nobs=nrow(object),
				idData(object),
			row.names="summary")
	line <- Line(object@coords)
	SpatialLinesDataFrame(
			SpatialLines(
				list(Lines(list(line), ID="summary")), 
				object@proj4string),
			data=data)
})
          
setMethod("summary", 
          signature=".MoveTrackStack", 
          definition=function(object){
          	os <- split(object)
            do.call(rbind, mapply(function(object, name) {
            	s <- summary(object)
            	spChFIDs(s, name)
            }, os, names(os)))
          })

"displacement.distance" <- function(object) {
	if (nrow(object) > 0) {
		sqrt(sum((object@coords[nrow(object),] - object@coords[1,])^2))
	} else {
		NA
	}
}

"displacement.direction" <- function(object, units=c("radian","degree","compass")) {
	units <- match.arg(units)
	if (nrow(object) > 0) {
		res <- atan2(object@coords[nrow(object),2] - object@coords[1,2],
				object@coords[nrow(object),1] - object@coords[1,1])
	} else {
		res <- NA
	}
	
	if (units == "degree" || units == "compass") {
		res <- res * 180 / pi
	}
	if (units == "compass") {
		# A compass has reversed (clockwise) direction and is rotated by 90 degrees
		# Rotate by 450 degrees to make sure every value is positive before the mod operation
		res <- (450 - res) %% 360
	}
	res
}
