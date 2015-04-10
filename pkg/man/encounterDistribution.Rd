\name{encounterDistribution}
\alias{encounterDistribution}
\title{
Spatial distribution of encounters
}
\description{
Computes the expected duration of encounters at each location for every pair of IDs.
}
\usage{
encounterDistribution(tr, threshold, grid=NULL, timestepSize = 60, xc=NULL, yc=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tr}{The trajectory for which to compute the UD}
  \item{threshold}{The maximum distance at which an encounter occurs}
  \item{grid}{a grid of class \code{\link{asc}}}
  \item{timestepSize}{The difference between consecutive time steps, in seconds}
  \item{xc}{The x coordinates of the vertical grid lines}
  \item{yc}{The y coordinates of the horizontal grid lines}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
One can specify the grid in two ways. If \code{grid} is set, the coordinates of
the grid lines are derived from there. Otherwise, these coordinates have to be
specified in the parameters \code{xc} and \code{yc}.
}
\value{
The return value is a list, indexed by two IDs.
If \code{grid} is given, each element of the result list is an object of class
\code{asc}, representing the same grid. Otherwise, each element of the list is
a matrix, indexed by the coordinates specified in \code{xc} and \code{yc}.

Element \code{result[["id1", "id2"]]} of the result represents the distribution
of the position of \code{id1} while it had encounters with \code{id2}. This is
not the same as \code{result[["id2","id1"]]}, since that is the distribution of
\code{id2}'s location during its encounters with \code{id1}.

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

# Define grid lines: equally spaced between the min and max coordinate in tr
xlim <- range(unlist(sapply(tr, function(b) { b$x })), na.rm=TRUE)
xc <- seq(xlim[1], xlim[2], length.out=20)

ylim <- range(unlist(sapply(tr, function(b) { b$y })), na.rm=TRUE)
yc <- seq(ylim[1], ylim[2], length.out=20)

# Compute the UD and plot the result for one ID
ud <- encounterDistribution(tr, 100, xc=xc, yc=yc)
image(ud[["BD","NH"]])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
