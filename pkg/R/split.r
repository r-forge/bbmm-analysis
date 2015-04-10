###create a list of Move objects from a Move Stack (hand over additional arguments!)
setGeneric("split") 
setMethod(f = "split",
	signature = c(x="MoveBBStack", f="missing"),
	definition = function(x, f, ...){
		vs <- split(x@variance, x@trackId)
		s<-split(as(x,'.MoveTrackStack'))
		
		moveBBList <- list()
		for (ID in names(s)) {
			moveBBList[[ID]] <- new('MoveBB',
				s[[ID]],
				as(x, '.MoveGeneral'),
				variance=vs[[ID]],
				diffusion=x@diffusion[ID]
			)
		}
		moveBBList
})

setGeneric("unstack")
setMethod(f = "unstack",
		signature = c(x="UDStack"),
		definition = function(x, ...) {
	us <- callNextMethod()
	lapply(us, function(u) { new("UD", u) })
})
