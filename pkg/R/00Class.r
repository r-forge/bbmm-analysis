setClass(Class = ".BBInfo",
	representation = representation(
			variance = "numeric",
			diffusion = "list"
	),
	prototype = prototype(
			variance = numeric(0),
			diffusion = list(function(t) { rep(0, length(t)) })
	),
	validity = function(object) {
		if (any(is.na(object@variance) | object@variance < 0)) {
			stop("Variance must be set and nonnegative")
		}
		if (any(sapply(object@diffusion, class) != "function"
				& !sapply(object@diffusion, is.null))) {
			stop("Diffusion coefficient mappings must be functions")
		}
	}
)

setClass(Class = "MoveBB",contains=c("Move", ".BBInfo"), 
	representation=representation(),
	prototype=prototype(),
	validity = function(object){
		if (length(object@variance) != nrow(object@coords)) {
			stop(paste("Number of location variances does not match number of coordinates: ", length(object@variance), nrow(object@coords)))
		}
		if (length(object@diffusion) != 1) {
			stop("There can only be one diffusion coefficient mapping")
		}
	}
)

setClass(Class = "MoveBBStack", contains=c("MoveStack", ".BBInfo"), 
	representation=representation(),
	prototype=prototype(),
	validity = function(object){
		if (length(object@variance) != nrow(object@coords)) {
			stop("Number of location variances does not match number of coordinates")
		}
		if (length(object@diffusion) != nrow(object@idData)) {
			stop("There must be exactly one diffusion coefficient mapping for each track")
		}
	}
)
