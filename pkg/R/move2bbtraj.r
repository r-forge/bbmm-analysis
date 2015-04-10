## define bbtrajs when needed
#if(!isClass("bbtraj"))
#    setClass("bbtraj")


#setAs("Move", "bbtraj", function(from){
#      if(isLonLat(from))
#	      warning('Converting a long lat projected object while the ltraj does not deal with long lat projected data')
#      adehabitatLT::as.ltraj(as.data.frame(coordinates(from)),date=timestamps(from), id=rownames(from@idData), typeII=T, infolocs=data.frame(sensor=from@sensor,from@data))
#})
#setAs("MoveStack", "bbtraj", function(from){
#      if(isLonLat(from))
#	      warning('Converting a long lat projected object while the ltraj does not deal with long lat projected data')
#      adehabitatLT::as.ltraj(as.data.frame(coordinates(from)),date=timestamps(from), id=from@trackId, typeII=T, infolocs=data.frame(sensor=from@sensor,from@data))
#})

#setAs("bbtraj", "Move", function(from) {
#    if (!inherits(from, "bbtraj"))
#        stop("from should be of class \"bbtraj\"")
#    if(length(from)!=1)
#	    stop("Can only convert one individual to a move object")
#    if(!attr(from,"typeII"))
#	    stop('Can only work on typeII objects')
#	from <- na.omit(from)
#	missing <- if (is.null(attr(from[[1]], 'na.action'))) {
#			data.frame(date=NULL)
#		} else {
#			as.data.frame(attr(from[[1]], 'na.action'))
#		}
#    spdf<-adehabitatLT::ltraj2spdf(from)
#    new("Move", spdf, sensor=rep(factor("unknown"), nrow(spdf)), timestamps=spdf$date, 
#    	timestampsUnUsedRecords=missing$date, sensorUnUsedRecords=rep(factor("unknown"),nrow(missing)), dataUnUsedRecords=missing,
#    	idData=data.frame(row.names=paste(attr(from[[1]], 'id'),'_', attr(from[[1]],'burst'), sep=""),
#    	burst=attr(from[[1]],'burst'), id=attr(from[[1]],'id')))
#})

#setAs("bbtraj", "MoveStack", function(from) {
#    if (!inherits(from, "bbtraj"))
#        stop("from should be of class \"bbtraj\"")
#    res<-list()
#    for(i in 1:length(from))
#    {
#	    res[[i]]<-as(from[i], 'Move')
#    }
#    moveStack(res)


#})

