as.bbtraj <- function(xys, date, id, burst=id, typeII = TRUE,
                     slsp =  c("remove", "missing"),
                     infolocs = data.frame(pkey = paste(id, date, sep=".")),
                     make.regular=FALSE)
{
	## Various verifications
	if (typeII) {
		if (!inherits(date,"POSIXct"))
			stop("For objects of type II,\n date should be of class \"POSIXct\"")
	} else {
		date <- 1:nrow(xys)
	}
	if (length(date) != nrow(xys))
		stop("date should be of the same length as xys")
	## Length of infolocs, if provided
	if (!is.null(infolocs)) {
		if (nrow(infolocs)!=nrow(xys))
			stop("infolocs should have the same number of rows as xys")
	}

	slsp <- match.arg(slsp)

	## length of id
	if (length(id)==1)
		id <- rep(as.character(id), nrow(xys))
	if (length(id)!=nrow(xys))
		stop("id should be of the same length as xys, or of length 1")
	id <- as.character(id)

	## length of burst
	if (length(burst)==1)
		burst <- rep(as.character(burst), nrow(xys))
	if (length(burst)!=nrow(xys))
		stop("burst should be of the same length as xys, or of length 1")
	burst <- as.character(burst)

	## Verification that there is only one burst per id
	id1 <- factor(id)
	burst1 <- factor(burst)
	if (!all(apply(table(id1,burst1)>0,2,sum)==1))
		stop("one burst level should belong to only one id level")

	x <- xys[,1]
	y <- xys[,2]
	s2 <- xys[,3]
	res <- split(data.frame(x=x,y=y,diff.coeff=rep(-1,nrow(xys)),loc.var=s2, date=date), burst)
	liid <- split(id, burst)
	
	if (!is.null(infolocs))
		linfol <- split(infolocs, burst)

	## sort the dates
	if (!is.null(infolocs))
		linfol <- lapply(1:length(linfol),
						 function(j) linfol[[j]][order(res[[j]]$date),,drop=FALSE])

	## sort the dates
	res <- lapply(res, function(y) y[order(y$date),, drop=FALSE])

	## Unique dates?
	rr <- any(unlist(lapply(res,
	                        function(x) { length(unique(x$date))!=length(x$date) })))
	if (rr)
		stop("non unique dates for a given burst")



	## Descriptive parameters
	foo <- function(x) {
		x1 <- x[-1, ]
		x2 <- x[-nrow(x), ]
		dist <- c(sqrt((x1$x - x2$x)^2 + (x1$y - x2$y)^2),NA)
		R2n <- (x$x - x$x[1])^2 + (x$y - x$y[1])^2
		dt <- c(unclass(x1$date) - unclass(x2$date), NA)
		dx <- c(x1$x - x2$x, NA)
		dy <- c(x1$y - x2$y, NA)
		abs.angle <- ifelse(dist<1e-07,NA,atan2(dy,dx))
		## absolute angle = NA if dx==dy==0
		so <- cbind.data.frame(dx=dx, dy=dy, dist=dist,
		                       dt=dt, R2n=R2n, abs.angle=abs.angle)
		return(so)
	}

	speed <- lapply(res, foo)
	res <- lapply(1:length(res), function(i) cbind(res[[i]],speed[[i]]))

	## The relative angle
	ang.rel <- function(df,slspi=slsp) {
		ang1 <- df$abs.angle[-nrow(df)] # angle i-1
		ang2 <- df$abs.angle[-1] # angle i

		if(slspi=="remove"){
			dist <- c(sqrt((df[-nrow(df),"x"] - df[-1,"x"])^2 +
			               (df[-nrow(df),"y"] - df[-1,"y"])^2),NA)
			wh.na <- which(dist<1e-7)
			if(length(wh.na)>0){
				no.na <- (1:length(ang1))[!(1:length(ang1)) %in% wh.na]
				for (i in wh.na){
					indx <- no.na[no.na<i]
					ang1[i] <- ifelse(length(indx)==0,NA,ang1[max(indx)])
				}
			}
		}
		res <- ang2-ang1
		res <- ifelse(res <= (-pi), 2*pi+res,res)
		res <- ifelse(res > pi, res -2*pi,res)
		return(c(NA,res))
	}

	rel.angle <- lapply(res, ang.rel)
	res <- lapply(1:length(res),
	              function(i) data.frame(res[[i]], rel.angle=rel.angle[[i]]))
	res <- lapply(1:length(res), function(i) {
		x <- res[[i]]
		attr(x, "id") <- as.character(liid[[i]][1])
		attr(x,"burst") <- levels(factor(burst))[i]
		return(x)
	})

	## The diffusion coefficient
	class(res) <- c("bbtraj","ltraj","list")
	attr(res,"typeII") <- typeII
	attr(res,"regular") <- is.regular(res)
	
	
	
	s2 <- diffusionCoefficient(res)
	res <- lapply(res, function(b) {
		b$diff.coeff <- s2[attr(b, "id")]
		return(b)
	})
	
	class(res) <- c("bbtraj","ltraj","list")
	attr(res,"typeII") <- typeII
	attr(res,"regular") <- is.regular(res)
	
	if (is.list(make.regular)) {
		# call setNA and sett0 to make the trajectory regular
		make.regular$ltraj <- res
		make.regular$ltraj <- do.call(adehabitatLT::setNA, make.regular)
		ltraj <- do.call(adehabitatLT::sett0, make.regular)
		
		# ltraj has some dropped columns, just copy the date fields from there to res
		res <- lapply(res, function(b) {
			ltb <- ltraj[burst=attr(b, "burst")][[1]]
			b$date <- ltb$date
			b$dt <- ltb$dt
			b
		})
		attributes(res) <- attributes(ltraj)
	}

	## And possibly, the data.frame infolocs
	if (!is.null(infolocs)) {
		res <- lapply(1:length(res), function(i) {
			x <- res[[i]]
			y <- linfol[[i]]
			row.names(y) <- row.names(x)
			attr(x, "infolocs") <- y
			return(x)
		})
	}

	## Output
	class(res) <- c("bbtraj","ltraj","list")
	attr(res,"typeII") <- typeII
	attr(res,"regular") <- is.regular(res)

	return(res)
}

"[.bbtraj" <- function(x, i, id, burst)
{
	if (!inherits(x, "bbtraj"))
		stop("x should be of class \"bbtraj\"")
	
	y <- "[.ltraj"(x, i, id, burst) # use the selection operator from ltraj to select the right bursts

	class(y) <- class(x)
	return(y)
}

"[<-.bbtraj" <- function(x, i, id, burst, value)
{
	if (!inherits(x, "ltraj"))
		stop("x should be of class \"ltraj\"")
	
	y <- "[<-.ltraj"(x, i, id, burst, value)
	
	class(y) <- class(x)
	return(y)
}

"diffusionCoefficient" <- function (tr, byburst = FALSE, nsteps = 1000) {
	tr <- bbFilterNA(tr)

	resultNames <- unique(id(tr))
	if (byburst) {
		resultNames <- burst(tr)
	}
	
	inputData <- list()
	for (n in resultNames) { inputData[[n]] <- matrix(0, nrow=3, ncol=0) }
	maxDiffCoeff <- rep(0, length(resultNames)) # For each element of the result, find a proper range to search
	names(maxDiffCoeff) <- resultNames

	# For each even numbered measurement in a burst (except the last if burst has even length):
	#  - compute the relevant parameters for the estimation of the diffusion coefficient:
	#      - the distance between the observation and the mean derived from the neighbouring observations
	#      - parameters used to compute the location variance from the diffusion coefficient
	#  - compute the ML value for the diffusion coefficient looking only at one bridge.
	#      When the diffusion coefficient gets larger than the maximum of these,
	#      the likelihood of the observations becomes a decreasing function of the diffusion coefficient.
	for (b in 1:length(tr)) {
		burst <- tr[[b]]

		if (nrow(burst) >= 3) { # We can't estimate the likelihood for shorter bursts
			burst$date <- as.double(burst$date) - min(as.double(burst$date))
		
			bdata <- matrix(NA, nrow=3, ncol=floor((nrow(burst)-1)/2))
		
			fac <- attr(burst, ifelse(byburst, "burst", "id"))
			# i runs over all even numbers that have a measurement before and after them
			for (i in seq(2, nrow(burst)-1, by=2)) {
				alpha <- (burst$date[i] - burst$date[i-1]) / (burst$date[i+1] - burst$date[i-1])
			
				# Squared distance from measured loc to estimated mean
				bdata[1, i/2] <- 
						(burst$x[i] - (1-alpha) * burst$x[i-1] - alpha * burst$x[i+1])^2 +
						(burst$y[i] - (1-alpha) * burst$y[i-1] - alpha * burst$y[i+1])^2
				# The coefficients for a linear function mapping diffusion coefficient to variance
				bdata[2, i/2] <- (burst$date[i+1]-burst$date[i-1]) * alpha * (1-alpha)
				bdata[3, i/2] <- (1-alpha)^2*burst$loc.var[i-1] + alpha^2*burst$loc.var[i+1]
			
				# Find the variance at which this bridge reaches its maximum likelihood
				maxVar <- 0.5 * bdata[1, i/2]
				maxDiffCoeff[fac] <- max(maxDiffCoeff[fac], (maxVar-bdata[3, i/2])/bdata[2, i/2])
			}
			
			inputData[[fac]] <- cbind(inputData[[fac]], bdata)
		}
	}
	
	result <- rep(NA, length(resultNames))
	names(result) <- resultNames
	for(fac in names(result)) {
		if (maxDiffCoeff[fac] <= 0) {
			# We already know that the ML is achieved for diff.coeff 0
			result[fac] <- 0
		} else {
			candidates <- seq(0,maxDiffCoeff[fac], length.out=nsteps)
		
			cResult <- .C("diffusion_static", double(1),
					c(inputData[[fac]]), ncol(inputData[[fac]]),
					candidates, length(candidates), PACKAGE="movementAnalysis")
			result[fac] <- cResult[[1]]
		}
	}
	
	return(result)
}

"bbFilterNA" <- function(tr) {
	# Filter out all rows with NA in all bursts
	trNew <- lapply(tr, function(burst) {
		filter <- (!is.na(burst$x)
			& !is.na(burst$y)
			& !is.na(burst$diff.coeff)
			& !is.na(burst$loc.var))
	
		newBurst <- burst[filter,]
		atts <- attributes(burst)
		attributes(newBurst) <- atts[which(names(atts) != "row.names")]
		names <- atts[["row.names"]]
		rownames(newBurst) <- names[filter]
		
		return(newBurst)
	})
	trNew <- trNew[unlist(sapply(trNew, function(b) { nrow(b) > 0}))]
	attributes(trNew) <- attributes(tr)
	class(trNew) <- class(tr)
	return(trNew)
}
