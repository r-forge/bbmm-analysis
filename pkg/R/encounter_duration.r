"encounterDurationById" <- function(encounterDuration, groupBy=NULL) {
	# Compute the total duration of encounters for each pair of IDs
	

	# Make sure all IDs in the encounterList will be present in the matrix
	burstIDs <- .IDs(attr(encounterDuration, "trajectory"), groupBy=groupBy)
	useIDs <- unique(burstIDs)
	
	models <- names(encounterDuration)[3:length(encounterDuration)]

	encounterMatrix <- array(rep(0, length(models) * length(useIDs)^2), 
		dim=c(length(useIDs),length(useIDs),length(models)),
		dimnames=list(useIDs, useIDs, models))
	# No encounters between ID and itself, set diagonal entries
	for (id in useIDs) {
		encounterMatrix[id, id,] <- NA
	}
	
	duration <- sapply(encounterDuration[,models], as.numeric)
	if (nrow(encounterDuration) > 0) {
		for (i in 1:nrow(encounterDuration)) {
			id1 <- burstIDs[encounterDuration$burst1[i]]
			id2 <- burstIDs[encounterDuration$burst2[i]]
			
			encounterMatrix[id1, id2, models] <- encounterMatrix[id1, id2, models] + c(duration[i, models])
		}
	}
	
	# We only have the part above the diagonal now, make the matrix symmetric
	encounterMatrix <- encounterMatrix + aperm(encounterMatrix, c(2,1,3))
	encounterMatrix <- as.difftime(encounterMatrix, units=attr(encounterDuration[[models[1]]], "units"))
	class(encounterMatrix) <- c("encounterDuration", class(encounterMatrix))
	encounterMatrix
}

"encounterDuration" <- function(tr, threshold, model=c("BBMM", "Linear"), groupBy=NULL, timestepSize=60)
{
# If groupBy=NULL, produces encounter duration for each pair of distinct elements in the MoveBBStack that overlap in time
# Else, produces encounter duration for each pair of distinct values of the field indicated by groupBy
	model <- match.arg(model, several.ok=TRUE)

	# Allocate maximum-size data frame to collect the results in
	trs <- split(tr)
	n <- length(trs)
	emptyCol <- list(rep(NA, n*(n-1)/2)) # size n choose 2, TODO: this can probably be optimized
	result <- as.data.frame(rep(emptyCol, 2+length(model)))
	names(result) <- c("burst1", "burst2", model)
	rm(emptyCol)
	
	resultSize <- 0 # Insert the next result into row resultSize+1

	# Sort the bursts in tr by their starting time
	trs <- trs[sort.list(sapply(trs, function(b) { min(b@timestamps) }))]
	activeBursts <- list()

	# For each pair of bursts, compute overlap and add to proper result
	for (burst in names(trs)) {
		# Filter out the bursts that end before this one starts
		activeBursts <- activeBursts[sapply(activeBursts, function(b) {
			max(b@timestamps)
			}) >= min(trs[[burst]]@timestamps)]
		
		for (b in names(activeBursts)) {
			# b overlaps with burst in time, so compute encounter durations
			resultSize <- resultSize+1
			# Make sure lexicographically smaller burst name appears first
			result[resultSize,1:2] <- range(b, burst)
			
			result[resultSize,model] <- sapply(model, function(m) {
				fn <- paste(".burstEncounterDuration", m, sep="")
				do.call(fn, list(trs[[burst]], activeBursts[[b]], threshold, timestepSize))
			})
		}
		
		activeBursts <- c(activeBursts, trs[burst])
	}

	# Store the results in difftime format	
	result[,model] <- lapply(result[,model], function(d) { as.difftime(d, units="secs") })
	
	attr(result, "trajectory") <- tr
	if (is.null(groupBy)) {
		result <- result[1:resultSize,] # return only the filled in rows
		class(result) <- c("encounterDuration", "data.frame")
		result
	} else {
		encounterDurationById(result[1:resultSize,], groupBy=groupBy)
	}
}

".burstEncounterDurationLinear" <- function (b1, b2, dist, timestepSize) {
	# serialize the fields we need from the bursts in row-major order
	# Substitute 0 for diffusion coefficient, linear model doesn't care
	d1 <- cbind(b1@coords, 0, b1@variance, as.double(b1@timestamps))
	d1 <- d1[apply(d1, 1, function(r) { !any(is.na(r)) }),]

	data1 <- c(t(d1)) # flatten into row-major vector
	
	d2 <- cbind(b2@coords, 0, b2@variance, as.double(b2@timestamps))
	d2 <- d2[apply(d2, 1, function(r) { !any(is.na(r)) }),]
	
	data2 <- c(t(d2)) # flatten into row-major vector

	cResult <- .C("encounterLinear",
			double(1), as.double(dist),
			as.integer(nrow(d1)), as.double(data1),
			as.integer(nrow(d2)), as.double(data2),
			as.double(timestepSize), PACKAGE="movementAnalysis")
	cResult[[1]] # The result is collected in the first C argument
}

".burstEncounterDurationBBMM" <- function (b1, b2, dist, timestepSize) {
	timeLimits <- c(range(as.double(b1@timestamps)), range(as.double(b2@timestamps)))
	timeLimits <- c(max(timeLimits[c(1,3)]), min(timeLimits[c(2,4)]))
	# Align the time of each step to a multiple of timestepSize
	timeStamps <- unique(c(seq(timeLimits[1], timeLimits[2], by=timestepSize), timeLimits[2]))

	p1 <- position(b1, timeStamps)
	p2 <- position(b2, timeStamps)
	
	# TODO: fix the stepsize of the first and two last rows
	data <- cbind(p1[,c('x','y')] - p2[,c('x','y')], p1[,'var']+p2[,'var'], timestepSize)
	data[1, "timestepSize"] <- timestepSize / 2
	data[nrow(data),   "timestepSize"] <-
			(timeStamps[length(timeStamps)] - timeStamps[length(timeStamps)-1]) / 2
	data[nrow(data)-1, "timestepSize"] <-
			data[1, "timestepSize"] + data[nrow(data), "timestepSize"]

	cResult <- .C("encounterBBMM",
			double(1), as.double(dist),
			as.integer(nrow(data)), as.double(t(data)),
			PACKAGE="movementAnalysis")
	cResult[[1]] # The result is collected in the first C argument
}
