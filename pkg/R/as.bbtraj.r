as.bbtraj <- function(xys, date=NULL, id, burst=id, typeII = TRUE,
                     slsp =  c("remove", "missing"))
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
    res <- split(data.frame(x=x,y=y,diff.coeff=rep(NA,nrow(xys)),loc.var=s2, date=date), burst)
    liid <- split(id, burst)

    ## sort the dates
    res <- lapply(res, function(y) y[order(y$date),])

    ## Unique dates?
    rr <- any(unlist(lapply(res,
                            function(x) (length(unique(x$date))!=length(x$date)))))
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

    ## Output
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
    
    
    s2 <- .diffusion_coefficient(res, c(0,10)); # TODO: make sure the range for the coefficient works out
	res <- lapply(res, function(tr) {
		tr$diff.coeff <- s2[[attr(tr, "id")]]
		return(tr)
	})

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

".diffusion_coefficient" <- function (tr, rangesig1, le=1000, byburst = FALSE)
{
    x <- adehabitatLT:::.ltraj2traj(tr)
    if (!inherits(x, "traj"))
        stop("tr should be of class \"ltraj\"")

    sorties <- list()
    gr <- grid
    x <- x[!is.na(x$x), ]
    x <- x[!is.na(x$y), ]

    fac <- x$burst
    if (!byburst)
        fac <- x$id
    fac <- factor(fac)
    lixy <- split(x, fac)
    so <- list()
    for (i in 1:length(lixy)) {
        dft <- lixy[[i]]
        df <- dft[, c("x", "y", "loc.var")]
        vsig <- seq(rangesig1[1], rangesig1[2],
                    length=le)
        date <- as.double(dft$date) - min(as.double(dft$date))
        huhu <- .C("diffusion", as.double(t(as.matrix(df))),
                   as.double(date), as.integer(nrow(df)),
                   double(length(vsig)), as.double(vsig^2),
                   as.integer(length(vsig)),
                   PACKAGE="movementAnalysis")
                   
       # C implementation signature:
       # void diffusion_static(double *xyr, double *Tr, 
       #	 int *nloc, double *Lr, double *sigma, int *nsig)
       
       # Maximisation of the likelihood for the Brownian bridge
       
       # Meaning of the params:
       # xyr    The location data
       # Tr     Relative time of the measurements
       # nloc   Number of location measurements
       # Lr     Log-likelihood for each considered value of the diffusion coefficient
       # sigma  The values of the diffusion coefficient to consider
       # nsig   The number of different diffusion coefficients to consider
       # sigma2 The location uncertainty variance

       so[[i]] <- vsig[which.max(huhu[[4]])] # select the sigma1 that maximizes the likelihood
    }
    names(so) <- names(lixy)
    #class(so) <- "liker"
    return(so)
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
	trNew <- trNew[sapply(trNew, function(b) { nrow(b) > 0})]
	attributes(trNew) <- attributes(tr)
	class(trNew) <- class(tr)
	return(trNew)
}
