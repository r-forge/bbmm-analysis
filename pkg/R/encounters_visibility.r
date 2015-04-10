#'encounter.visible' <- function(vis, tr, time) {
#	ids <- unique(adehabitatLT::id(tr))
#	
#	tr <- bbFilterNA(tr) # This function does not like missing values
#		
#	# Construct an array indexed by time stamp and two IDs.
#	# Each element holds the requested statistic for these two groups at that time
#	res <- sapply(time, function(t) {
#		params <- position(tr, t)[,,1]
#		
#		res <- array(rep(1, length(ids)^2), dim=c(length(ids),length(ids)))
#		
#		v <- lapply(vis, function(v) { c(t(as.matrix(v))) })
#		dim(v) <- dim(vis)
#		dimnames(v) <- dimnames(vis)
#		
#		for (i in 1:(length(ids)-1)) {
#			for (j in (i+1):length(ids)) {
#				res[i,j] <- res[j,i] <- .Call('encounterProbability',
#					v, params[i,], params[j,], as.double(rownames(v)), as.double(colnames(v)))
#			}
#		}
#		
#		res
#	})

#	dim(res) <- c(length(ids), length(ids), length(time))
#	dimnames(res) <- list(ids, ids, time)
#	res
#}
