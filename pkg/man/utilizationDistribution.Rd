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
utilizationDistribution(tr, grid=NULL, timestepSize = 60,
    xc=NULL, yc=NULL, grid.dim = 100, grid.pad = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tr}{The trajectory for which to compute the UD}
  \item{grid}{a grid of class \code{\link{asc}}}
  \item{timestepSize}{The difference between consecutive time steps, in seconds}
  \item{xc}{The x coordinates of the vertical grid lines}
  \item{yc}{The y coordinates of the horizontal grid lines}
  \item{grid.dim}{If all of \code{grid}, \code{xc} and \code{yc} are omitted,
  this specifies the dimension of the shortest axis of the grid used.}
  \item{grid.pad}{If the grid is automatically generated, its range is the range
  of the relocations extended \code{grid.pad} times the length of the range on each side.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
One can specify the grid in two ways. If \code{grid} is set, the coordinates of
the grid lines are derived from there. Otherwise, these coordinates have to be
specified in the parameters \code{xc} and \code{yc}.
}
\value{
Returns a list, indexed by the IDs in \code{tr}.
If \code{grid} is given, each element of the result list is an object of class
\code{asc}, representing the same grid. Otherwise, each element of the list is
a matrix, indexed by the coordinates specified in \code{xc} and \code{yc}.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data("example_data", package="movementAnalysis")
tr <- as.bbtraj(data.frame(x=example_data$X, y=example_data$Y, var=example_data$StdDev^2),
	date=example_data$DateTime, id=example_data$GroupID, burst=example_data$GroupDayNo)

# Define grid lines: equally spaced between the min and max coordinate in tr
xlim <- range(unlist(sapply(tr, function(b) { b$x })), na.rm=TRUE)
xc <- seq(xlim[1], xlim[2], length.out=100)

ylim <- range(unlist(sapply(tr, function(b) { b$y })), na.rm=TRUE)
yc <- seq(ylim[1], ylim[2], length.out=100)

# Compute the UD and plot the result for one ID
ud <- utilizationDistribution(tr, xc=xc, yc=yc)
plot(ud)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
