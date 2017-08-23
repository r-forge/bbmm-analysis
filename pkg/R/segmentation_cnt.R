## Model-based segmentation (continuous), by default based on the BBMM

setGeneric("segment.cnt", function(object, LL, p) standardGeneric("segment.cnt"))

setMethod(f = "segment.cnt",
          signature = c(object="MoveBB", LL="missing"),
          definition = function (object, LL, p) {
            LL <- function(tr) {
              d <- diffusionCoefficient(tr)
              result <- attr(d, "LL")
              attr(result, "dc") <- d
              result
            }
            attr(LL, "step.size") <- 2
            attr(LL, "monotone")  <- TRUE
            segment.cnt(object, LL, p)
          })

setMethod(f = "segment.cnt",
          signature = c(object="MoveBB", LL="function", p="missing"),
          definition = function (object, LL, p) {
            .segment.cnt(object, LL, log(nrow(object))) ## Use BIC by default
          })

setMethod(f = "segment.cnt",
          signature = c(object="MoveBB", LL="function", p="numeric"),
          definition = function (object, LL, p) {
            s <- .segment.cnt(object, LL, p)
            ## TODO: map s to subtrajectories,
            ## TODO: deal with step size (i.e. ncol(LL) vs nrow(object))
          })

.segment.cnt <- function (tr, LL, p) {
  opt.pred  <- opt.DC <- rep(NA, nrow(tr))
  opt.IC              <- rep(Inf, nrow(tr))
  opt.IC[1]           <- 0
  
  step.size <- attr(LL, "step.size")
  if (is.null(step.size)) { step.size <- 1 }
  
  ## Start computing possible breakpoints after first observation
  for (i in seq(1+step.size, nrow(tr), by=step.size)) {
    for (j in seq(i-step.size, 1, by=-step.size)) {
      sll <- p - (2 * LL(tr[j:i]))
      ## End the search when there can be no better solution
      ## TODO: can we be more aggressive? Maybe based on previous iteration of outer loop?
      if (attr(LL, "monotone") && sll >= opt.IC[i]) {
        print(paste("break", i, j))
        break
      }
      
      sll <- opt.IC[j] + sll
      if (sll < opt.IC[i]) {
        opt.IC[i] <- sll
        opt.pred[i] <- j
        opt.DC[i] <- attr(sll, "dc")
      }
    }
  }
  return(list(opt.pred, opt.IC, opt.DC))
}

