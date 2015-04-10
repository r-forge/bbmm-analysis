#plot.enc.ltraj <- function (x, id = unique(unlist(lapply(x, attr, which="id"))),
#                        burst = unlist(lapply(x, attr, which="burst")),
#                        asc = NULL, area = NULL, xlim = NULL,
#                        ylim = NULL, colasc = gray((240:1)/256),
#                        colpol = "green", addpoints = TRUE,
#                        addlines = TRUE, perani = TRUE, final = TRUE,
#                        threshold = 100, ...)
#{
#    polygon <- area
#    if (!is.null(area)) {
#        if (!inherits(area, "area"))
#            stop("area should be an object of class area")
#    }
#    if (!inherits(x, "ltraj"))
#        stop("x should be an object of class ltraj")

#    # Select id and burst
#    x <- x[id=id]
#    x <- x[burst=burst]
#    typeII <- attr(x, "typeII")

#    # Remove NA
#    x <- lapply(x, function(i) {
#        jj <- i[!is.na(i$x),]
#        attr(jj, "id") <- attr(i,"id")
#        attr(jj, "burst") <- attr(i,"burst")
#        return(jj)
#    })
#    class(x) <- c("ltraj","list")
#    attr(x, "typeII") <- typeII
#    attr(x,"regular") <- is.regular(x)

#    # keeps only the coordinates
#    uu <- lapply(x, function(i) {
#        i[,c("x","y")]
#    })

#    # Which kind of plot?
#    idc <- "id"

#    id <- unique(unlist(lapply(x, function(i) {
#        attr(i, idc)
#    })))

#    if (length(id)>1)
#        opar <- par(mar = c(0.1, 0.1, 2, 0.1), mfrow = n2mfrow(length(id)))

#    # xlim is the limit of the graph
#    if (is.null(xlim)) {
#        idtt <- unique(adehabitatLT::id(x))
#        oo <- lapply(idtt,
#                     function(i) unlist(lapply(x[id=i], function(j) j$x)))
#        maxxl <- max(unlist(lapply(oo, function(kk) diff(range(kk)))))
#        xlim <- lapply(oo, function(ki) c(min(ki), min(ki)+maxxl))
#    }  else {
#        xlim <- split(rep(xlim, length(id)), gl(length(id),2))
#    }

#    if (is.null(ylim)) {
#        idtt <- unique(adehabitatLT::id(x))
#        oo <- lapply(idtt,
#                     function(i) unlist(lapply(x[id=i], function(j) j$y)))
#        maxyl <- max(unlist(lapply(oo, function(kk) diff(range(kk)))))
#        ylim <- lapply(oo, function(ki) c(min(ki), min(ki)+maxyl))
#    } else {
#        ylim <- split(rep(ylim, length(id)), gl(length(id),2))
#    }
#    names(xlim) <- id
#    names(ylim) <- id

#    for (i in id) {
#        if (!is.null(asc)) {
#            if (length(id)==1) {
#                image(asc, col = colasc, xlim = xlim[i][[1]],
#                      ylim = ylim[i][[1]], ...)
#            } else {
#                image(asc, col = colasc, xlim = xlim[i][[1]],
#                      ylim = ylim[i][[1]], main = i,
#                      axes = (length(id)==1), ...)
#            }
#        } else {
#            if (length(id)==1) {
#                plot(1, 1, type = "n", asp = 1,
#                     xlim = xlim[i][[1]],
#                     ylim = ylim[i][[1]], ...)
#            } else {
#                plot(1, 1, type = "n", asp = 1,
#                     xlim = xlim[i][[1]],
#                     ylim = ylim[i][[1]], axes = (length(id)==1),
#                     main = i, ...)
#            }
#        }
#        if (length(id)>1)
#            box()
#        if (!is.null(polygon)) {
#            pol <- split(polygon[, 2:3], factor(polygon[, 1]))
#            for (j in 1:length(pol)) polygon(pol[[j]], col = colpol)
#        }
#        if (addlines) {
#            xtmp <- x[id=i]
#             TODO: This is big time cheating. Only two bursts will have their lines plotted, but they will have encounters
#            coordinates <- .encounterSegments(xtmp[[1]], xtmp[[2]], threshold)
#            for (j in 1:2) {
#            	xc <- coordinates[[j]]$x
#            	yc <- coordinates[[j]]$y
#            	enc <- coordinates[[j]]$encounter
#            	for (k in 1:(nrow(coordinates[[j]])-1)) {
#            		lines(xc[k:(k+1)], yc[k:(k+1)], 
#            				col=ifelse(enc[k], c("green","red")[j], "black"),
#            				lwd=ifelse(enc[k], 3, 1))
#            	}
#            }
#        }
#        if (addpoints) {
#            xtmp <- x[id=i]
#            for (j in 1:length(xtmp)) {
#                points(xtmp[[j]]$x, xtmp[[j]]$y,
#                       pch = 21, col = "black", bg = "white")
#            }
#        }
#        if (final) {
#            xtmp <- x[id=i]
#            for (j in 1:length(xtmp)) {
#                points(xtmp[[j]]$x[c(1, length(xtmp[[j]]$x))],
#                       xtmp[[j]]$y[c(1, length(xtmp[[j]]$x))],
#                       pch = c(2,14), col = c("blue","red"), cex=2, lwd=2)
#            }
#        }
#    }
#    if (length(id)>1)
#        par(opar)
#}


#".encounterSegments" <- function(b1, b2, threshold) {
#	intervals <- as.vector(t(encounterIntervals(b1, b2, threshold)))
#	
#	return(list(.burst_encounter_segments(b1, intervals),
#			.burst_encounter_segments(b2, intervals)))
#}

#".burst_encounter_segments" <- function(burst, intervals) {
#	burst$t <- as.double(burst$date)
#	nrows <- nrow(burst) + length(intervals)
#	res <- data.frame(x=rep(NA, nrows), y=rep(NA, nrows), encounter=rep(FALSE, nrows))
#	encounter <- FALSE
#	j <- 0
#	for (i in 1:nrow(burst)) {
#		while (j < length(intervals) && burst$t[i] > intervals[j+1]) {
#			alpha <- (intervals[j+1] - burst$t[i-1]) / (burst$t[i] - burst$t[i-1])
#			x <- (1-alpha) * burst$x[i-1] + alpha * burst$x[i]
#			y <- (1-alpha) * burst$y[i-1] + alpha * burst$y[i]
#			
#			encounter <- !encounter			
#			res[i+j,] <- list(x, y, encounter)
#			j <- j+1
#		}
#	
#		res[i+j,] <- list(burst$x[i], burst$y[i], encounter)
#	}
#	
#	return(res[1:(nrow(burst)+j),])
#}

