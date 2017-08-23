## Model-based classification, by default based on the BBMM

setGeneric("classify", function(object, LL, p) standardGeneric("classify"))

setMethod(f = "classify",
		signature = c(object="MoveBBStack", LL="missing"),
		definition = function (object, LL, p) {
	classify(object, 1000, p)
})

setMethod(f = "classify",
		signature = c(object="MoveBBStack", LL="numeric"),
		definition = function (object, LL, p) {
	dc.alpha <- 3 ## Exponent for power-law distribution of values
	dc.range <- range(diffusionCoefficient(object))^(1/dc.alpha)
	dc.values <- seq(dc.range[1], dc.range[2], length.out=LL)^dc.alpha
	tr.ll <- diffusion.LL(split(object), dc.values)
	classify(object, tr.ll, p)
})

setMethod(f = "classify",
		signature = c(object="MoveBBStack", LL="matrix", p="missing"),
		definition = function (object, LL, p) {
	.classify(LL, log(ncol(LL))) ## Use BIC by default
})

setMethod(f = "classify",
		signature = c(object="MoveBBStack", LL="matrix", p="numeric"),
		definition = function (object, LL, p) {
	.classify(LL, p)
})



".classify" <- function(LL, p) {
	
	## First, ensure the LL functions are sorted by increasing arg max
	M <- apply(LL, 2, which.max)
	order.inp <- order(M)
	LL[,seq_len(ncol(LL))] <- LL[,order.inp]
	M <- M[order(M)] ## In which iteration a trajectory gets included
	fn.until <- findInterval(0:(nrow(LL)-1), M) ## Opt_i classifies functions 1 to fn.until[i]
	
	## Transform the LL to the form used in the IC
	nLL <- -2 * LL
	
	## The classification with everything at x_1 to initialize with
	C <- rep(NA, ncol(LL))
	attr(C, "IC") <- Inf
		
	Opt <- matrix(NA, nrow=nrow(LL)+1, ncol=ncol(LL))
	for (i in 1:nrow(LL)) {
		## Everything at x_i as extension of Opt_0
		Opt[i,seq_len(fn.until[i])] <- i
		ICi.Opt <- sum(nLL[i,seq_len(fn.until[i])]) + p
		
		for (j in seq_len(i-1)) {
			## Extend Opt_j to a candidate for Opt_i
			O <- Opt[j,]
			if (fn.until[j] < fn.until[i]) {
				## Assign the added functions to either x_j or x_i depending on LL
				l.upd <- (fn.until[j]+1):fn.until[i]
				O[l.upd] <- c(j,i)[apply(LL[c(j,i),l.upd,drop=F], 2, which.max)]
			}
			
			## Compute IC conditioned on using x_i
			nLL.O <- apply(cbind(1:ncol(LL), O), 1, function(idx) { nLL[idx[2],idx[1]] })
			ICi.O <- sum(nLL.O[!is.na(nLL.O)]) + p*length(unique(c(O[1:fn.until[i]],i)))
			
			if (ICi.O < ICi.Opt) {
				Opt[i,] <- O
				ICi.Opt <- ICi.O
			}
		}
		## Construct a complete classification from Opt_i
		Ci <- Opt[i,]
		if (fn.until[i] < ncol(LL)) {
			Ci[(fn.until[i]+1):ncol(LL)] <- i
		}
		
		## If Ci is better than C, replace it
		nLL.Ci <- apply(cbind(1:ncol(LL), Ci), 1, function(idx) { nLL[idx[2],idx[1]] })
		IC.Ci <- sum(nLL.Ci) + p*length(unique(Ci))
		if (IC.Ci < attr(C, "IC")) {
			C <- Ci
			attr(C, "IC") <- IC.Ci
		}
	}
	
	## Undo the sorting operation at the start of the function
	C[order.inp] <- C
	if (!is.null(attr(LL, "dc.values"))) {
		attr(C, "dc.values") <- attr(LL, "dc.values")[C]
	}
	C
}

