\name{utilizationDistribution}
\alias{utilizationDistribution}
\title{
Compute utilization distribution for a trajectory
}
\description{
This function computes the utilization distribution (UD) for each of the IDs
present in the trajectory. The user specifies the grid on which the UD is
evaluated and may also specify the size of the time step in the numerical
integration.
}
\usage{
utilizationDistribution(tr, xc, yc, timestepSize = 60)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tr}{The trajectory for which to compute the UD}
  \item{xc}{The x coordinates of the vertical grid lines}
  \item{yc}{The y coordinates of the horizontal grid lines}
  \item{timestepSize}{The difference between consecutive time steps, in seconds}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list, indexed by the IDs in \code{tr}.
Each element of the list is a matrix, indexed by the coordinates specified
in \code{xc} and \code{yc}.
}
\note{
%%  ~~further notes~~
The \code{\link{image}} function has ugly colours, use the \code{col} attribute to define a better colour map.
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
xc <- seq(xmin, xmax, length.out=100)

ymin <- min(unlist(sapply(tr, function(b) { b$y })))
ymax <- max(unlist(sapply(tr, function(b) { b$y })))
yc <- seq(ymin, ymax, length.out=100)

# Compute the UD and plot the result for one ID
ud <- utilizationDistribution(tr, xc, yc)
image(ud[["BD"]])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
