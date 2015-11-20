###create a list of MoveBB objects from a MoveBBStack (hand over additional arguments!)
setGeneric("split") 
setMethod(f = "split",
	signature = c(x="MoveBBStack", f="missing"),
	definition = function(x, f, ...) {
		xs <- split(x, x@trackId)
		lapply(xs, as, "MoveBB")
})
#		vs <- split(x@variance, x@trackId)
#		s<-split(as(x,'.MoveTrackStack'))
#		
#		moveBBList <- list()
#		for (ID in names(s)) {
#			moveBBList[[ID]] <- new('MoveBB',
#				s[[ID]],
#				as(x, '.MoveGeneral'),
#				variance=vs[[ID]],
#				diffusion=x@diffusion[ID]
#			)
#		}
#		moveBBList
#})

setMethod(f = "split",
	signature = c(x="MoveBBStack", f="character"),
	definition = function(x, f, ...) {
		## Split by idData field named by f
		## Find the associated f value for each obs by trackId
		split(x, x@idData[[f]][x@trackId])
})

setMethod(f = "split",
	signature = c(x="MoveBBStack", f="factor"),
	definition = function(x, f, ...) {
		# First generate lists of observations for each subset
		obs.id <- split(seq_len(nrow(x)), f)
		# Then extract the requested observations using subsetting
		mclapply(obs.id, function(i) { x[i] })
})


setGeneric("unstack")
setMethod(f = "unstack",
		signature = c(x="UDStack"),
		definition = function(x, ...) {
	us <- callNextMethod()
	lapply(us, function(u) { new("UD", u) })
})
