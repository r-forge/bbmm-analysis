"encounterDurationById" <- function(encounterDurationByBurst, IDs=NULL) {
	# Compute the total duration of encounters for each pair of IDs

	# Make sure all IDs in the encounterList will be present in the matrix
	useIDs <- sort(unique(c(IDs, encounterDurationByBurst$id1, encounterDurationByBurst$id2)))

	models <- names(encounterDurationByBurst)[5:length(encounterDurationByBurst)]

	encounterMatrix <- array(rep(0, length(models) * length(useIDs)^2), 
		dim=c(length(useIDs),length(useIDs),length(models)),
		dimnames=list(useIDs, useIDs, models))
	for (id in useIDs) {
		encounterMatrix[id, id,] <- NA
	}
	
	duration <- sapply(encounterDurationByBurst[,models], as.numeric)
	if (nrow(encounterDurationByBurst) > 0) {
		for (i in 1:nrow(encounterDurationByBurst)) {
			id1 <- encounterDurationByBurst$id1[i]
			id2 <- encounterDurationByBurst$id2[i]
			
			encounterMatrix[id1, id2, models] <- encounterMatrix[id1, id2, models] + c(duration[i, models])
		}
	}
	
	# We only have the part above the diagonal now, make the matrix symmetric
	encounterMatrix <- encounterMatrix + aperm(encounterMatrix, c(2,1,3))
	encounterMatrix <- as.difftime(encounterMatrix, units=attr(encounterDurationByBurst[[models[1]]], "units"))
	class(encounterMatrix) <- c("encounterDuration", class(encounterMatrix))
	encounterMatrix
}

"encounterDuration" <- function(tr, threshold, model=c("BBMM", "linear"), byburst=FALSE, timestepSize=60)
{
# If byburst, produces encounter duration for each pair of bursts that overlap in time
# Else, produces encounter duration for each pair of distinct IDs
	tr <- bbFilterNA(tr) # This function does not like missing values

	# Allocate maximum-size data frame to collect the results in
	n <- length(tr)
	emptyCol <- list(rep(NA, n*(n-1)/2)) # size n choose 2, TODO: this can probably be optimized
	result <- as.data.frame(rep(emptyCol, 4+length(model)))
	names(result) <- c("id1", "id2", "burst1", "burst2", model)
	rm(emptyCol)
	
	resultSize <- 0 # Insert the next result into row resultSize+1

	# Sort the bursts in tr by their starting time
	tr <- tr[sort.list(sapply(tr, function(b) {min(b$date)}))]
	activeBursts <- list()

	# For each pair of bursts, compute overlap and add to proper result
	for (burst in tr) {
		# Filter out the bursts that end before this one starts
		activeBursts <- activeBursts[sapply(activeBursts, function(b) {
			max(b$date)
			}) >= min(burst$date)]
		
		for (b in activeBursts) {
			# b overlaps with burst in time, so compute encounter durations
			resultSize <- resultSize+1
			result[resultSize,1:4] <- c(attr(b, "id"), attr(burst, "id"), attr(b, "burst"), attr(burst, "burst"))
			result[resultSize,model] <- sapply(model, function(m) {
				.burstEncounterDuration(burst, b, threshold, timestepSize, m) 
			})
		}
		
		activeBursts <- c(activeBursts, list(burst))
	}

	# Store the results in difftime format	
	result[,model] <- lapply(result[,model], function(d) { as.difftime(d, units="secs") })
	
	if (byburst) {
		result <- result[1:resultSize,] # return only the filled in rows
		class(result) <- c("encounterDuration", "data.frame")
		result
	} else {
		ids <- unique(sapply(tr, function(x) { attr(x, "id") }))
		encounterDurationById(result[1:resultSize,], ids)
	}
}

".burstEncounterDuration" <- function (b1, b2, dist, timestepSize, model) {
	implementations <- c(linear="encounterLinear", BBMM="encounterBBMM")
	# serialize the fields we need from the bursts in row-major order
	b1$t <- as.double(b1$date)
	b1 <- b1[!is.na(b1$x), c("x","y","diff.coeff","loc.var","t")]
	b1 <- b1[!is.na(b1$y),]
	b1 <- b1[!is.na(b1$loc.var),]
	data1 <- c(t(b1)) # flatten into row-major vector
	
	b2$t <- as.double(b2$date)
	b2 <- b2[!is.na(b2$x), c("x","y","diff.coeff","loc.var","t")]
	b2 <- b2[!is.na(b2$y),]
	b2 <- b2[!is.na(b2$loc.var),]
	data2 <- c(t(b2)) # flatten into row-major vector

	cResult <- .C("encounterLinearIntervals",
			double(1), as.double(dist),
			as.integer(nrow(b1)), as.double(data1),
			as.integer(nrow(b2)), as.double(data2),
			as.double(timestepSize), double(2*(nrow(b1) + nrow(b2))), PACKAGE="movementAnalysis")
	cResult[[1]] # The result is collected in the first C argument
}
