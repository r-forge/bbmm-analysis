## Model-based segmentation (discrete), by default based on the BBMM

# TODO: Implement fixed number of segments using exponential-binary search over p?

#setGeneric("segment.dscr", function(object, LL, p, n.dc) standardGeneric("segment.dscr"))

#setMethod(f = "segment.dscr",
#          signature = c(object="MoveBB", LL="missing"),
#          definition = function (object, LL, p) {
#            LL <- function(tr) {
#              d <- diffusionCoefficient(tr)
#              result <- attr(d, "LL")
#              attr(result, "dc") <- d
#              result
#            }
#            attr(LL, "step.size") <- 2
#            attr(LL, "monotone")  <- TRUE
#            segment.dscr(object, LL, p)
#          })

#setMethod(f = "segment.dscr",
#          signature = c(object="MoveBB", LL="function", p="missing", n.dc="numeric"),
#          definition = function (object, LL, p) {
#            .segment.dscr(object, LL, log(nrow(object))) ## Use BIC by default
#          })

#setMethod(f = "segment.dscr",
#          signature = c(object="MoveBB", LL="matrix", p="numeric"),
#          definition = function (object, LL, p) {
#            s <- .segment.dscr(object, LL, p)
#            ## TODO: map s to subtrajectories,
#            ## TODO: deal with step size (i.e. ncol(LL) vs nrow(object))
#          })

segment.dscr <- function (tr, p=NA, n.dc=100) {
	## TODO: signature: tr, LL (matrix listing LL values), p
	## TODO: move preprocessing out to Methods
	## TODO: deal with dc=0 and position error = 0
	
	step.size = 2 #TODO: variable step size
 	if (nrow(tr) %% step.size != 1) {
		tr <- tr[-nrow(tr)] ## Make sure tr has odd length
	}

	if (is.na(p)) {
		p <- log(nrow(tr))
	}
	
	dc.alpha <- 3 ## Exponent for power-law distribution of values
	dc.range <- range(sapply(1:(nrow(tr)-2), function(i) {
		diffusionCoefficient(tr[i:(i+2)])
	}))^(1/dc.alpha)
	dc.values <- seq(dc.range[1], dc.range[2], length.out=n.dc)^dc.alpha
		
	opt.IC             <- rep(Inf, nrow(tr))
	opt.pred <- opt.dc <- rep(NA,  nrow(tr))
	opt.IC[1]          <- 0
	
	cur.IC   <- rep(Inf, n.dc)
	cur.pred <- rep(NA,  n.dc)
	for (i in seq(3,nrow(tr),by=step.size)) {
		prev.IC   <- cur.IC
		prev.pred <- cur.pred
		
		LL <- diffusion.LL(list(tr[(i-step.size):i]), dc.values)
		for (v in 1:n.dc) {
			## Test two options: append and extend
			append <- opt.IC[i-step.size] + p
			extend <- prev.IC[v]
			
			
			if (append < extend) {
				cur.IC[v]   <- append - 2*LL[v] #LL[i,v]
				cur.pred[v] <- i - step.size
			} else {
				cur.IC[v]   <- extend - 2*LL[v] #LL[i,v]
				cur.pred[v] <- prev.pred[v]
			}
			
			if (!is.nan(cur.IC[v]) && cur.IC[v] < opt.IC[i]) {
				opt.IC[i]   <- cur.IC[v]
				opt.pred[i] <- cur.pred[v]
				opt.dc[i]   <- dc.values[v]
			}
		}
	}
	
	## Walk the predecessor pointers to construct the segmentation
	i <- nrow(tr) ## TODO: take step.size into account
	res <- list()
	res.names <- character(0)
	while (i > 1) {
		prev.i <- opt.pred[i]
		
		cur.segment <- tr[prev.i:i]
		diffusion(cur.segment) <- opt.dc[i]
		res <- c(res, cur.segment)
		res.names <- c(res.names, paste(prev.i, i, sep="-"))
#		res <- c(res, i)
				
		i <- prev.i
	}
	names(res) <- res.names
	#moveBBStack(rev(res))
	rev(res)
}

