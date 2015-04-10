setGeneric("moveBB", function(x, y, var, time, data, proj=NA, ...) standardGeneric("moveBB"))
setMethod(f="moveBB",
	signature=c(var="missing"), 
	definition = function(x,y,var,time,data,proj,...){
		moveBB(x=x, y=y, var=0, time=time, data=data, proj=proj,...)
	}
)
setMethod(f="moveBB",
	signature=c(var="numeric",data="missing"), 
	definition = function(x,y,var,time,data,proj,...){
		data <- data.frame(x,y,var,time)
		moveBB(x=x, y=y, var=var, time=time, data=data, proj=proj,...)
	}
)
setMethod(f="moveBB",
	signature=c(x="numeric", y="numeric", var="numeric", time="POSIXct", data="data.frame", proj="ANY"), 
	definition = function(x,y,var,time,data,proj,sensor='unknown',animal='unnamed', ...){
		# Create a move object, then compose it with the extra BB properties
		.moveBB(move(x, y, time, data, proj, sensor, animal, ...), var[!is.na(x)])
	}
)

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
	
	# Set the diffusion coefficient if the trajectory is projected
	if (!isLonLat(res)) {
		diffusion(res) <- diffusionCoefficient(res)
	}
	res
}

.IDs <- function(moveObj, groupBy=NULL) {
	if (is.null(groupBy)) {
		ids <- rownames(moveObj@idData)
	} else {
		ids <- moveObj@idData[[groupBy]]
	}
	ids <- as.character(ids)
	names(ids) <- rownames(moveObj@idData)
	ids
}