\name{contourPolygons}
\alias{contourPolygons}
\title{
Compute home range polygons from utilization distribution
}
\description{
Given a utilization distribution, compute the polygons that contain a certain fraction of the UD.
}
\usage{
contourPolygons(ud, levels = 0.99)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ud}{The utilization distribution.}
  \item{levels}{The level(s) at which to compute contours.}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list of \code{\linkS4class{SpatialPolygons}} as defined in package \code{sp}. Each element of the list represents the result for one utilization distribution from \code{ud}. If multiple levels are requested, these are separate \code{\linkS4class{Polygons}} in one element of the result list.
}

\seealso{
\code{\link{utilizationDistribution}}
}
\examples{
data("vervet_monkeys", package="moveBB")

ud <- utilizationDistribution(monkey.tr)
cp <- contourPolygons(ud, levels=c(0.95, 0.99))
str(cp[[1]])

# Draw the home range for the first animal
plot(cp[[1]], col=1:2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
