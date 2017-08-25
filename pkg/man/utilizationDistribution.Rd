\name{utilizationDistribution}
\alias{utilizationDistribution}
\alias{[[.utilizationDistribution}
\alias{utilizationDistribution.byID}
\title{
Compute utilization distribution for a trajectory
}
\description{
This function computes the utilization distribution (UD) for each of the IDs
present in the trajectory. The user may specify the grid on which the UD is
evaluated and may also specify the size of the time step in the numerical
integration.
}
\usage{
utilizationDistribution(tr, grid=NULL, timestepSize = 60,
    xc=NULL, yc=NULL, grid.dim = 100, grid.pad = 0.2,
    byburst = FALSE, cutoff.level = 1)
utilizationDistribution.byID(ud)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tr}{The trajectory for which to compute the UD}
  \item{grid}{a grid of class \code{\link[adehabitat]{asc}}}
  \item{timestepSize}{The difference between consecutive time steps, in seconds}
  \item{xc}{The x coordinates of the vertical grid lines}
  \item{yc}{The y coordinates of the horizontal grid lines}
  \item{grid.dim}{If all of \code{grid}, \code{xc} and \code{yc} are omitted,
  this specifies the dimension of the shortest axis of the grid used.}
  \item{grid.pad}{If the grid is automatically generated, its range is the range
  of the relocations extended \code{grid.pad} times the length of the range on each side.}
  \item{byburst}{Whether the result should be a UD for each burst.}
  \item{cutoff.level}{What fraction of the UD volume must be actually stored.}
  \item{ud}{The result of \code{utilizationDistribution} with \code{byburst=TRUE}}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
One can specify the grid in three ways:
\itemize{
	\item If \code{grid} is set, the coordinates of
	the grid lines are derived from there. 
	\item These coordinates can also be specified in the parameters \code{xc} and \code{yc}.
	\item Otherwise, the grid is determined automatically. The grid then ranges over
	a bounding box of all measurements, extended on each side by a specified fraction of the range.
	The number of grid cells can be controlled via \code{grid.dim}. The amount by which
	the grid is extended is controlled via \code{grid.pad}, which may be a vector specifying
	the extension on top,right,bottom,left respectively. It is recycled as usual.
}

If \code{cutoff.level} is set to a value less than 1, the result will not store
a utilization distribuiton for the full grid, but instead an axis-aligned box
that contains the full contour at the requested level.

When a UD is extracted from a list containing such partial UDs, the cells that
were left out are restored and filled with zeros.

\code{utilizationDistribution.byID} can be used to convert the result of \code{utilizationDistribution}
with \code{byburst=TRUE} to a result that is aggregated by ID.
}
\value{
Returns a list, indexed by the IDs in \code{tr}.
If \code{grid} is given, each element of the result list is an object of class
\code{asc}, representing the same grid. Otherwise, each element of the list is
a matrix, indexed by the coordinates specified in \code{xc} and \code{yc}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data("vervet_monkeys", package="moveBB")

# Compute the UD and plot the result
#ud <- utilizationDistribution(monkey.tr, byburst=TRUE)
#plot(ud)

# Convert to one UD for each group
#ud <- utilizationDistribution.byID(ud)
#plot(ud)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
