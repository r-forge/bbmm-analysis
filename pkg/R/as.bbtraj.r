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

	## Verification that there is only one id per burst
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
	# Add extra information about the relocations
	res <- lapply(res, .bbtraj.extra.info, slsp)
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
	names(res) <- sapply(res, function(burst) { attr(burst, "burst") })

	return(res)
}

# Compute the derived quantities for the relocations in a bridge.
# Derived quantities include: distance, daily displacement, duration, direction
".bbtraj.extra.info" <- function(burst, slspi) {
	if (nrow(burst) == 0) { return(burst) }

	## Descriptive parameters
	x1 <- burst[-1, ]
	x2 <- burst[-nrow(burst), ]
	dist <- c(sqrt((x1$x - x2$x)^2 + (x1$y - x2$y)^2),NA)        # Length of bridge
	
	first.row <- head(rbind(
			na.omit(burst[,c("x","y","loc.var")]),
			data.frame(x=NA, y=NA, loc.var=NA)), 1)  # First row with an actual measurement, if any exists
	R2n <- (burst$x - first.row$x)^2 + (burst$y - first.row$y)^2 # Squared displacement from burst start
	dt <- c(unclass(x1$date) - unclass(x2$date), NA)             # Duration of last bridge
	dx <- c(x1$x - x2$x, NA)                                     # Displacement on x-axis of bridge
	dy <- c(x1$y - x2$y, NA)                                     # Displacement on y-axis of bridge
	abs.angle <- ifelse(dist<1e-07,NA,atan2(dy,dx))              # Bridge direction: Radians from positive x-axis
	## absolute angle = NA if dx==dy==0
	## The relative angle, i.e. change in direction in the range (-pi, pi]
	ang1 <- abs.angle[-nrow(burst)] # angle i-1
	ang2 <- abs.angle[-1] # angle i

	if(slspi=="remove"){
		wh.na <- which(dist<1e-7)
		if(length(wh.na)>0){
			no.na <- (1:length(ang1))[!(1:length(ang1)) %in% wh.na]
			for (i in wh.na){
				indx <- no.na[no.na<i]
				ang1[i] <- ifelse(length(indx)==0,NA,ang1[max(indx)])
			}
		}
	}
	ra <- ang2-ang1
	ra <- ifelse(ra <= (-pi), 2*pi+ra,ra)
	ra <- ifelse(ra > pi, ra -2*pi,ra)
	ra <- c(NA, ra)

	so <- cbind.data.frame(dx=dx, dy=dy, dist=dist,
                       dt=dt, R2n=R2n, abs.angle=abs.angle, rel.angle=ra)
	b <- cbind(burst[,!colnames(burst) %in% colnames(so)], so)
	attr(b, "id") <- attr(burst,"id")
	attr(b, "burst") <- attr(burst,"burst")
	
	return(b)
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

#"[[.bbtraj" <- function(x, i, exact=TRUE) {
#	# If i is a string, regard it as a burst name
#	if (is.character(i)) {
#		burst.name <- i
#		i <- which(sapply(x, function(b) { attr(b, 'burst') }) == i)
#		if (length(i) == 0) { stop(paste("Unknown burst name", burst.name)) }
#	}

#	NextMethod("[[")
#}

"na.omit.bbtraj" <- function(object, slsp=c("remove", "missing"), ...) {
	slsp <- match.arg(slsp)

	# Filter out all rows with NA in all bursts
	trNew <- lapply(object, function(burst) {
		filter <- (!is.na(burst$x)
			& !is.na(burst$y)
			& !is.na(burst$diff.coeff)
			& !is.na(burst$loc.var))
	
		# Don't bother recomputing everything if we don't delete any rows
		if (all(filter)) { return(burst) }
		
		newBurst <- .bbtraj.extra.info(burst[filter,], slsp)
		rn <- split(rownames(burst), filter)
		rownames(newBurst) <- rn$`TRUE`
		
		# Store the information about the removed rows
		rm.rows <- rn$`FALSE`
		names(rm.rows) <- rn$`FALSE`
		attr(rm.rows, "class") <- "omit"
		attr(newBurst, "na.action") <- rm.rows
		return(newBurst)
	})
	keepBursts <- unlist(sapply(trNew, function(b) { nrow(b) > 0}))
	trNew <- trNew[keepBursts]
	attributes(trNew) <- attributes(object)[names(attributes(object)) != "names"]
	names(trNew) <- names(object)[keepBursts]
	class(trNew) <- class(object)
	return(trNew)
}
"bbFilterNA" <- na.omit.bbtraj # Keep this one for backwards compatibility

summary.bbtraj <- function(object, ..., units.direction="radian") {
	if (!inherits(object, "bbtraj"))
      stop("object should be of class \"bbtraj\"")
    res <- summary.ltraj(object,...)
    if (attr(object,"typeII")) {
    	object <- na.omit(object) # Remove the NAs already, then the call inside the functions has hardly any overhead
    	djl <- djl(object)
    	dist <- displacement.distance(object)
    	dir <- displacement.direction(object, units.direction)
    	res <- cbind(res, DJL=djl, displ.dist=dist, displ.dir=dir)
    }
    res
}

"djl" <- function(tr) {
	res <- sapply(na.omit(tr), function(burst) {
		sum(burst$dist, na.rm=T)
	})
	names(res) <- names(tr)
	res
}

"displacement.distance" <- function(tr) {
	res <- sapply(na.omit(tr), function(burst) {
		sqrt(burst$R2n[nrow(burst)])
	})
	names(res) <- names(tr)
	res
}

"displacement.direction" <- function(tr, units=c("radian","degree","compass")) {
	units <- match.arg(units)
	res <- sapply(na.omit(tr), function(burst) {
		atan2(burst$y[nrow(burst)] - burst$y[1], burst$x[nrow(burst)] - burst$x[1])
	})
	if (units == "degree" || units == "compass") {
		res <- res * 180 / pi
	}
	if (units == "compass") {
		# A compass has reversed direction and is rotated by 90 degrees
		# Rotate by 450 degrees to make sure every value is positive before the mod operation
		res <- (450 - res) %% 360
	}
	names(res) <- names(tr)
	res
}
