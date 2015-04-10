# Provides location information for each id in the trajectory at the given time(s)
position <- function(tr, time) {
	if (inherits(time, "POSIXct")) {
		time <- as.double(time)
	}
	
	tr <- bbFilterNA(tr)
	ids <- unique(id(tr))

	res <- sapply(time, function(t) {
	
		# Extract the bursts that contain the selected time	
		tr <- tr[which(sapply(tr, function(b) {
			b$date[1] <= t && b$date[nrow(b)] >= t
		}))]

		if (length(unique(id(tr))) != length(tr))
			stop("Multiple bursts contain requested time for some ID")

		res <- matrix(nrow=length(ids), ncol=3)
		if (length(tr) > 0) {
			for (i in 1:length(tr)) {
				b <- tr[[i]]

				b$date <- as.double(b$date)			
				j <- max(which(b$date >= t)[1], 2)

		
				start <- b[j-1,]
				end <- b[j,]
			
				T <- end$date-start$date
				alpha <- (t-start$date)/T
				res[i, 1:2] <- (1-alpha)*c(start$x, start$y) + alpha*c(end$x, end$y)
				res[i, 3] <- T*alpha*(1-alpha) * end$diff.coeff + (1-alpha)^2 * start$loc.var + alpha^2 * end$loc.var
			}
		}
		return(res)
	})
	
	dim(res) <- c(length(ids),3,length(time))
	dimnames(res) <- list(ids, c("x","y","var"), time)
	
	return(res)
}
