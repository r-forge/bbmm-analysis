setGeneric("moveBB", function(x, y, var, time, data, proj=NA, ...) standardGeneric("moveBB"))
setMethod(f="moveBB",
	signature=c(var="missing"), 
	definition = function(x,y,var,time,data,proj,...){
		moveBB(x=x, y=y, var=0, time=time, data=data, proj=proj,...)
	}
)
setMethod(f="moveBB",
	signature=c(x=".MoveTrack", y="missing", var="numeric", time="missing", data="missing", proj="missing"),
	definition = function (x, y=NULL, var, time=NULL, data=NULL, proj=NULL, ...) {
		.moveBB(x, var)
})
setMethod(f="moveBB",
	signature=c(var="numeric", data="missing"), 
	definition = function(x,y,var,time,data,proj,...){
		data <- data.frame(x,y,var,time)
		moveBB(x=x, y=y, var=var, time=time, data=data, proj=proj,...)
	}
)
setMethod(f="moveBB",
	signature=c(x="numeric", y="numeric", var="numeric", time="POSIXct", data="data.frame", proj="ANY"), 
	definition = function(x,y,var,time,data,proj,sensor='unknown',animal='unnamed', ...){
		# Create a move object, then compose it with the extra BB properties
		if (length(var) == length(x)) {
			# Drop variance for missing observations to align correctly
			var <- var[!is.na(x)]
		}
		.moveBB(move(x, y, time, data, proj, sensor, animal, ...), var)
	}
)

setAs("MoveBBStack", "MoveBB", function(from) {
	m <- from
	class(m) <- "MoveBB"
	## MoveBB is sorted by time instead of ID,time as in MoveBBStack
	m <- m[order(m@timestamps)]
	for (n in names(m@idData)) {
		if (length(m@idData[[n]]) > 1) {
			## Not idData for the coerced object
			m@idData[[n]] <- NULL
			#TODO: insert as regular data?
		} else if (is.factor(m@idData[[n]])) {
			m@idData[[n]] <- droplevels(m@idData[[n]])	
		}
	}
	m
})

.moveBB <- function(moveObj, var) {
	if (length(var) == 1) {
		var <- rep(var, nrow(moveObj@coords))
	}

	emptyDiff <- lapply(rownames(moveObj@idData), function(b) { NULL })
	names(emptyDiff) <- rownames(moveObj@idData)

	if (is(moveObj, "Move")) {
		res <- new("MoveBB",
				variance=var,
				diffusion=emptyDiff,
				moveObj
		)
	} else {
		res <- new ("MoveBBStack",
				variance=var,
				diffusion=emptyDiff,
				moveObj
		)
	}
	
	# Set the default diffusion coefficient if the trajectory is projected
	if (!isLonLat(res)) {
		diffusion(res) <- diffusionCoefficient(res)
	}
	res
}

IDs <- function(moveObj, groupBy=NULL) {
	if (is.null(groupBy)) {
		ids <- rownames(moveObj@idData)
	} else {
		ids <- moveObj@idData[[groupBy]]
		if (is.null(ids)) {
			stop(paste("Invalid field name '", groupBy, "' for groupBy", sep=""))
		}
	}
	ids <- as.factor(ids)
	names(ids) <- rownames(moveObj@idData)
	ids
}
.IDs <- IDs # Deprecated hidden version of the IDs function
