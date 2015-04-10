setGeneric("moveBBStack", function(x) standardGeneric("moveBBStack"))
setMethod(f = "moveBBStack", 
	signature = c(x="list"),
	definition = function(x){
		if (any((as.character(lapply(x, class)))!="MoveBB")) 
			stop("One or more objects in the list are not from class MoveBB")
		
		# Let the parent implementation deal with most work, then add the custom fields
		ms <- moveStack(lapply(x, as, "Move"))
		var <- unlist(sapply(x, slot, "variance"))
		diff <- unlist(lapply(lapply(x, slot, "diffusion"), unname))
		
		new("MoveBBStack",
				ms,
				variance=var,
				diffusion=diff
		)
})

