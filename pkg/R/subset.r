setMethod("[", 
		signature(x="MoveBB"),
		definition=function(x,i,j,...) {
	if(!missing(i)) {
		x@variance <- x@variance[i]
	} else { i<-T }
	if (missing(j)) { j<-T }
	
	callNextMethod(x=x,i=i,j=j,...)
})

setMethod("[", 
		signature(x="MoveBBStack"),
		definition=function(x,i,j,...) {
	if(!missing(i)) {
		x@variance <- x@variance[i]
	} else { i<-T }
	if (missing(j)) { j<-T }
	
	x<-callNextMethod(x=x,i=i,j=j,...)
	
	x@diffusion <- x@diffusion[rownames(x@idData)]
	print(names(x@diffusion))
	x
})
