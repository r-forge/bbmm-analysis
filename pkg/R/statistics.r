setGeneric("statistics", function(object, add=TRUE) {
	standardGeneric("statistics")
})

setMethod("statistics", 
          signature=c(object=".MoveTrackSingle"), 
          definition=function(object, add=TRUE) {
    
    dp <- rbind(object@coords[-1,] - object@coords[-nrow(object),], c(NA,NA))
    angle <- mapply(atan2, dp[,2], dp[,1]) #signature is atan2(y,x)

	data <- data.frame(
				dt=as.difftime(c(as.numeric(diff(object@timestamps), units="secs"), NA), units="secs"),
				dx=dp[,1],
				dy=dp[,2],
				dist=apply(dp^2, 1, function(d) { sqrt(sum(d)) }),
				displ=apply(object@coords, 1, function(d) {
					sqrt(sum((d-object@coords[1,])^2)) }),
				abs.angle=angle,
				turn.angle=c(NA, angle[-1] - angle[-length(angle)]),
			row.names=rownames(object))
	if (add) {
		object@data <- cbind(object@data, data)
		object
	} else {
		cbind(object@coords, data) # Add the coordinates if we return a separate data.frame
	}
})
          
setMethod("statistics", 
		signature=".MoveTrackStack", 
		definition=function(object, add=TRUE){
	data <- do.call(rbind, lapply(split(object), statistics, add=FALSE))
	if (add) {
		object@data <- cbind(object@data, data[,-c(1,2)])
		object
	} else {
		data
	}
})
