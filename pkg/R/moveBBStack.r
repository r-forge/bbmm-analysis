setGeneric("moveBBStack", function(x) standardGeneric("moveBBStack"))
setMethod(f = "moveBBStack", 
	signature = c(x="list"),
	definition = function(x){
		if (any((as.character(lapply(x, class)))!="MoveBB")) 
			stop("One or more objects in the list are not from class MoveBB")
		
		ms <- moveStack(lapply(x, as, "Move"))

		var <- sapply(x, slot, "variance")
print(var)
stop("Stop")
})
