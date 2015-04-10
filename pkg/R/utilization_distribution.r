"utilizationDistribution" <- function(tr, xc, yc, timestepSize=60)
{
	
	timeSteps <- .UDtimesteps(tr, timestepSize)
	
	UDs <- list()
	# For each timestep, compute the mean location, variance, and a weighing factor
	for (id in names(timeSteps)) {
		ts <- timeSteps[[id]]
	
		pad <- rep(0, 15)
			
		cResult <- .OpenCL("utilizationDistribution", length(xc)*length(yc),
				as.integer(nrow(ts)),
				as.double(c(ts[,1],pad)), as.double(c(ts[,2],pad)), as.double(c(ts[,3],pad+1)), as.double(c(ts[,4],pad)), 
				as.double(xc), as.double(yc), as.integer(length(yc)), as.integer(length(xc)))
		
#		cResult <- .OpenCL("getUD", length(xc)*length(yc),
#			as.integer(length(timesteps)/4), as.double(timesteps),
#			as.double(xc), as.double(yc), as.integer(length(yc)), as.integer(length(xc)))

		res <- matrix(cResult, nrow=length(yc))
		rownames(res) <- yc
		colnames(res) <- xc
		
		UDs[[id]] <- res
	}
	
	return(UDs)
}

"encounterDistribution" <- function(tr, threshold, xc, yc, timestepSize=60) {
	useFields <- c("x","y","diff.coeff","loc.var","t")
	
	timeSteps <- .UDtimesteps(tr, timestepSize)
	ids <- names(timeSteps)
	
	# UDs is two-dimensional list where each entry represents the distribution of encounters for a pair of groups
	UDs <- as.list(rep(NA, length(ids)^2))
	dim(UDs) <- c(length(ids), length(ids))
	rownames(UDs) <- colnames(UDs) <- ids
	
	for (id1 in ids) {
		# It makes no sense to compute encounters between a group and itself, so don't do that
		for (id2 in ids[-which(ids==id1)]) {
			# encounters can only be determined if we have an estimate for both IDs
			encounterTimes <- intersect(rownames(timeSteps[[id1]]), rownames(timeSteps[[id2]]))
			
			weight <- pmin(timeSteps[[id1]][encounterTimes,'weight'], timeSteps[[id2]][encounterTimes,'weight'])
			
			ts  <- timeSteps[[id1]][encounterTimes,]
			ts2 <- timeSteps[[id2]][encounterTimes,]
			
			pad <- rep(0, 15)
			
			cResult <- .OpenCL("encounterUD", length(xc)*length(yc),
			as.double(threshold), as.integer(length(encounterTimes)),
			as.double(c( ts[,1],pad)), as.double(c( ts[,2], pad)), as.double(c( ts[,3],pad+1)), as.double(c(ts[,4], pad)),
			as.double(c(ts2[,1],pad)), as.double(c(ts2[,2], pad)), as.double(c(sqrt(ts2[,3]),pad+1)),
			as.double(xc), as.double(yc), as.integer(length(yc)), as.integer(length(xc)))
			
			res <- matrix(cResult, nrow=length(yc))
			rownames(res) <- yc
			colnames(res) <- xc
			UDs[[id1,id2]] <- res
		}
	}
	diag(UDs) <- utilizationDistribution(tr, xc, yc, timestepSize)
	return(UDs)
}

".UDtimesteps" <- function(tr, timestepSize=60) {
	ids <- unique(id(tr))
	
	result <- list()
	
	for (id in ids) {
		bursts <- tr[id=id]
		result[[id]] <- do.call("rbind", lapply(bursts, function(burst) {
				burst$t <- as.double(burst$date)
				burst <- burst[!(is.na(burst$x) || is.na(burst$y) || is.na(burst$loc.var)),
				               c("x","y","diff.coeff","loc.var","t")]
				data <- c(t(burst)) # flatten into vector in row-major order
				
				timeLimits <- timestepSize*round(burst$t[c(1,nrow(burst))]/timestepSize)
				nsteps <- as.integer((timeLimits[2] - timeLimits[1]) / timestepSize) + 1
				
				cResult <- .C("UDTimesteps", double(nsteps*4),
						as.integer(nrow(burst)), as.double(data),
						as.integer(nsteps), timeLimits)
				timeSteps <- matrix(cResult[[1]], ncol=4, byrow=T)
				colnames(timeSteps) <- c('x', 'y', 'var', 'weight')
				rownames(timeSteps) <- seq(timeLimits[1], timeLimits[2], length=nrow(timeSteps))
				return(timeSteps)
		}))
	}
	return(result)
}

