# Provides location information for each id in the trajectory at the given time(s)
setGeneric("position", function(object, time, groupBy=NULL) standardGeneric("position"))
setMethod(f = "position",
	signature = c(object="MoveBBStack", time="POSIXct"),
	definition = function (object, time, groupBy=NULL) {
		p <- position(object, as.double(time), groupBy)
		dimnames(p)[[3]] <- format(time, usetz=TRUE)
		p
})

setMethod(f = "position",
	signature = c(object="MoveBBStack", time="numeric"),
	definition = function (object, time, groupBy=NULL) {
		tr <- split(object)
		
		ids <- unique(.IDs(object, groupBy))
		
		res <- sapply(time, function(t) {
			# Extract the bursts that contain the selected time	
			tr <- tr[sapply(tr, function(b) {
				as.double(b@timestamps[1]) <= t &&
				as.double(b@timestamps[nrow(b)]) >= t
			})]

			if (length(tr) > length(ids))
				stop("Multiple bursts contain requested time for some ID")

			res <- matrix(nrow=length(ids), ncol=3, dimnames=list(ids, c("x","y","var")))
			if (length(tr) > 0) {
				for (i in names(tr)) {
					b <- tr[[i]]
					i <- .IDs(b, groupBy) # get the proper ID for the burst
					
					ts <- as.double(b@timestamps)
					j <- max(which(ts >= t)[1], 2)

					diff0T <- integrate(b@diffusion[[1]], ts[j-1], ts[j])
					diff0t <- integrate(b@diffusion[[1]], ts[j-1], t)
					if (diff0T$message != "OK" || diff0t$message != "OK") {
						stop("Problem integrating diffusion coefficient")
					}
					diff0T <- diff0T$value
					diff0t <- diff0t$value
					alpha <- diff0t / diff0T
					res[i, c("x","y")] <- (1-alpha)*b@coords[j-1,] + alpha*b@coords[j,]
					res[i, "var"] <- (
							  diff0t*(1-alpha)
							+ (1-alpha)^2 * b@variance[j-1]
							+ alpha^2 * b@variance[j]
					)
				}
			}
			return(res)
		})
			
		dim(res) <- c(length(ids),3,length(time))
		dimnames(res) <- list(ids, c("x","y","var"), time)
	
		return(res)
})
