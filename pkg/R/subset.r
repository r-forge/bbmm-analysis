setGeneric("[")
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
		x@diffusion <- x@diffusion[unique(x@trackId[i])]
	} else { i<-T }
	if (missing(j)) { j<-T }
	callNextMethod(x=x,i=i,j=j,...)
})

setMethod("[[",
		signature(x="UDStack"),
		definition=function(x, i) {
	u <- callNextMethod()
	new("UD", u)
})

