\name{encounterDistribution}
\alias{encounterDistribution}
\title{
Spatial distribution of encounters
}
\description{
Computes the expected duration of encounters at each location for every pair of IDs.
}
\usage{
encounterDistribution(tr, threshold, xc, yc, timestepSize = 60)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tr}{The trajectory for which to compute the UD}
  \item{threshold}{The maximum distance at which an encounter occurs}
  \item{xc}{The x coordinates of the vertical grid lines}
  \item{yc}{The y coordinates of the horizontal grid lines}
  \item{timestepSize}{The difference between consecutive time steps, in seconds}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The return value is a list, indexed by two IDs. Each element of the list is
a matrix representing the expected duration of encounters at each cell of the
specified grid.
The diagonal entries of the result list contain the utilization distribution of
each ID, since an entity is always at a distance 0 from itself.
}
\note{
%%  ~~further notes~~
The \code{\link{image}} function has ugly colours, use the \code{col} attribute to define a better colour map.

Also note that this function may take a rather long time to complete, so please be patient,
specify a sufficiently small grid or use a larger time step.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Warning}{
There seems to be some problem with the result being transposed, this needs further investigation. Until then, you can plot the transpose of the result using \code{image(t(ud[["BD","NH"]]))}.
}
\examples{
data("example_data", package="movementAnalysis")
tr <- as.bbtraj(data.frame(x=example_data$X, y=example_data$Y, var=example_data$StdDev^2),
	date=example_data$DateTime, id=example_data$GroupID, burst=example_data$GroupDayNo)
tr <- bbFilterNA(tr) # Some operations in the following do not like NAs

# Define grid lines: equally spaced between the min and max coordinate in tr
xmin <- min(unlist(sapply(tr, function(b) { b$x })))
xmax <- max(unlist(sapply(tr, function(b) { b$x })))
xc <- seq(xmin, xmax, length.out=30)

ymin <- min(unlist(sapply(tr, function(b) { b$y })))
ymax <- max(unlist(sapply(tr, function(b) { b$y })))
yc <- seq(ymin, ymax, length.out=30)

# Compute the UD and plot the result for one ID
ud <- encounterDistribution(tr, 100, xc, yc)
image(ud[["BD","NH"]])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
