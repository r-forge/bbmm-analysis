\name{speedDistribution}
\alias{speedDistribution}
\title{
Compute average speed at each location
}
\description{
For each grid cell, this function computes the average speed over time conditioned on the trajectory being in that cell.
}
\usage{
speedDistribution(tr, grid = NULL, timestepSize = 60, time.scale = timestepSize, xc = NULL, yc = NULL, grid.dim = 100, grid.pad = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tr}{The trajectory for which to compute the UD}
  \item{grid}{a grid of class \code{\link[adehabitat]{asc}}}
  \item{timestepSize}{The difference between consecutive time steps, in seconds.}
  \item{time.scale}{The length of the intervals over which to compute the speed.}
  \item{xc}{The x coordinates of the vertical grid lines}
  \item{yc}{The y coordinates of the horizontal grid lines}
  \item{grid.dim}{If all of \code{grid}, \code{xc} and \code{yc} are omitted,
  this specifies the dimension of the shortest axis of the grid used.}
  \item{grid.pad}{If the grid is automatically generated, its range is the range
  of the relocations extended \code{grid.pad} times the length of the range on each side.}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Returns a list, indexed by the IDs in \code{tr}.
If \code{grid} is given, each element of the result list is an object of class
\code{asc}, representing the same grid. Otherwise, each element of the list is
a matrix, indexed by the coordinates specified in \code{xc} and \code{yc}, or by
the coordinates computed from \code{tr}.
}

\seealso{
\code{\link{velocity}}, \code{\link{speed}}, \code{\link{bbtraj}}
}
\examples{
data("vervet_monkeys", package="moveBB")

## Compute the speed distribution and plot the result
#sd <- speedDistribution(monkey.tr)
#image(sd[[1]]) # Plot only the distribution of the first entity
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
