\name{encounterDistribution}
\alias{encounterDistribution}
\title{
Spatial distribution of encounters
}
\description{
Computes the expected duration of encounters at each location for every pair of IDs.
}
\usage{
encounterDistribution(tr, threshold, grid=NULL, timestepSize = 60,
    xc=NULL, yc=NULL, grid.dim = 100, grid.pad = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tr}{The trajectory for which to compute the UD}
  \item{threshold}{The maximum distance at which an encounter occurs}
  \item{grid}{a grid of class \code{\link[adehabitat]{asc}}}
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
data("vervet_monkeys", package="moveBB")

## Define grid lines: equally spaced between the min and max coordinate in monkey.tr
#TODO: define from monkey.tr@extent
xlim <- range(monkey.data$X, na.rm=TRUE)
xc <- seq(xlim[1], xlim[2], length.out=20)

ylim <- range(monkey.data$Y, na.rm=TRUE)
yc <- seq(ylim[1], ylim[2], length.out=20)

## Compute the UD and plot the result for one ID
#ud <- encounterDistribution(monkey.tr, 100, xc=xc, yc=yc)
#image(ud[["BD","NH"]])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
