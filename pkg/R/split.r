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


setMethod(f = "split",
	  signature = c(x="move:::.MoveTrackStack", f="missing"),
	  definition = function(x, f, ...){
		  moveList <- list()
		  unUsed<-as(x,".unUsedRecordsStack")
		  sls <- split(x, x@trackId)
		  tss <- split(x@timestamps, x@trackId)
		  sns <- split(x@sensor, x@trackId)
		  uurs <- split(unUsed, unUsed@trackIdUnUsedRecords)
		  for (ID in names(sls)) {
			  s<-sls[[ID]]
			  spdf <- SpatialPointsDataFrame(coords = matrix(s@coords, ncol=2),
							 data=s@data,
							 proj4string=x@proj4string)
			  mt <- new(Class=".MoveTrack",
				    spdf,
				    timestamps=tss[[ID]],
				    sensor=sns[[ID]])
			  unUsedSub<-as(uurs[[ID]],'.unUsedRecords')
			  moveObj <- new(Class=".MoveTrackSingle", 
					 mt,
					 idData=x@idData[row.names(x@idData)==ID, ,drop=F],
	#				 dateCreation=x@dateCreation,
	#				 study=x@study,
	#				 citation=x@citation,
	#				 license=x@license,
					 unUsedSub)
			  moveList[[ID]]  <- moveObj
		  }
		  return(moveList)
	  }) 
