velocity <- function(tr, time, time.scale=NA) {
	tr <- na.omit(tr) # This function does not like missing values

	if (inherits(time, "POSIXct")) {
		time <- as.double(time)
	}
	if (is.na(time.scale)) {
		if (length(time) > 1) {
			# By default, use the average time between two elements of time argument
			time.scale <- (time[length(time)]-time[1])/(length(time)-1)
		} else {
			stop('Cannot determine appropriate time.scale')
		}
	}
	
	ids <- unique(adehabitatLT::id(tr))

	res <- sapply(time, function(t) {
		t.s <- t - time.scale/2
		t.e <- t + time.scale/2
		
		# Extract the bursts that contain the selected times
		tr.s <- tr[which(sapply(tr, function(b) {
			b$date[1] <= t.s && b$date[nrow(b)] >= t.s
		}))]
		tr.e <- tr[which(sapply(tr, function(b) {
			b$date[1] <= t.e && b$date[nrow(b)] >= t.e
		}))]
		
		p.s <- position(tr, t.s)[,,as.character(t.s)]
		p.e <- position(tr, t.e)[,,as.character(t.e)]
		# If there is only one ID, its dimension was just dropped; restore it
		dim(p.s) <- dim(p.e) <- c(length(ids), 3)
		dimnames(p.s) <- dimnames(p.e) <- list(ids, c('x','y','var'))
		
		res <- matrix(NA, nrow=length(ids), ncol=3)
		for (i in 1:length(ids)) {
			id <- ids[i]
		
			if (length(tr.s[id=id]) == 0 || length(tr.e[id=id]) == 0) { next }
			b.s <- tr.s[id=id][[1]]
			b.e <- tr.e[id=id][[1]]

			b.s$date <- as.double(b.s$date)
			b.e$date <- as.double(b.e$date)
			
			j.s <- max(which(b.s$date >= t.s)[1], 2)
			j.e <- max(which(b.e$date >= t.e)[1], 2)
			
			# The mean is always the velocity of the mean position
			res[i,1] <- (p.e[id,'x'] - p.s[id,'x']) / (t.e-t.s)
			res[i,2] <- (p.e[id,'y'] - p.s[id,'y']) / (t.e-t.s)
			
			# Depending on the number of measurements between t.s and t.e,
			# the variance is computed differently
			if (attr(b.s, 'burst') == attr(b.e, 'burst') && j.s == j.e) {
				# t.s and t.e are in the same bridge
				# See .... TODO: Add reference when available
				dt <- b.s$dt[j.s-1]
				res[i,3] <- sum(b.s$loc.var[c(j.s-1,j.s)]) / dt^2 +
					(1 / (t.e-t.s) - 1 / dt) * b.s$diff.coeff[j.s-1]
			} else if (attr(b.s, 'burst') == attr(b.e, 'burst') && j.s == j.e-1) {
				# t.s and t.e are in consecutive bridges of one burst
				# See .... TODO: Add reference when available
				res[i,3] <- (p.e[id,'var'] + p.s[id,'var'] - (2
							* (t.s - b.s$date[j.s-1]) / (b.s$date[j.s] - b.s$date[j.s-1])
							* (b.e$date[j.e]   - t.e) / (b.e$date[j.e] - b.e$date[j.e-1])
							* b.s$loc.var[j.s])
						) / (t.e-t.s)^2
			} else {
				# t.s and t.e are further apart in time, treat as independent
				res[i,3] <- (p.e[id,'var'] + p.s[id,'var']) / (t.e-t.s)^2
			}
		}
		return(res)
	})
	
	dim(res) <- c(length(ids),3,length(time))
	dimnames(res) <- list(ids, c("x","y","var"), time)
	return(res)
}
