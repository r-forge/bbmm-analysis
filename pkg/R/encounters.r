# Computes all pairs of bursts that overlap in time
".overlappingBursts" <- function(tr) {
	tr <- bbFilterNA(tr) # This function does not like missing values

	# Allocate maximum-size data frame to collect the results in
	n <- length(tr)
	emptyCol <- rep(NA, n*(n-1)/2) # size n choose 2, TODO: this can probably be optimized
	result <- data.frame(burst1=emptyCol, burst2=emptyCol)
	names(result) <- c("burst1", "burst2")
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
			
			# The call to range() sorts the burst names
			result[resultSize,] <- range(c(attr(burst,'burst'), attr(b, 'burst')))
		}
		
		activeBursts <- c(activeBursts, list(burst))
	}
	
	result[1:resultSize,]
}
