\name{burst-statistics}
\alias{displacement.distance}
\alias{displacement.direction}
\alias{djl}
\title{
Collect various statistics about individual bursts
}
\description{
These functions compute a certain statistic of every burst in a given trajectory. The statistics are net displacement distance, direction of net displacement and day journey length (DJL), which is the total length of all relocations.
}
\usage{
displacement.distance(tr)
displacement.direction(tr, units = c("radian", "degree", "compass"))
djl(tr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tr}{The trajectory object}
  \item{units}{The units in which to return the result}
}
\details{
All functions use the linear motion model to compute the statistics.

\code{displacement.direction} can return its result in different units. The default is in radians from the positive \eqn{x} axis, but it can also report degrees from this axis or compass heading (i.e. degrees from north, measured clockwise).
}
\value{
A vector containing the requested statistic for each burst, indexed by the burst name.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{bbtraj}}
}
\examples{
data("vervet_monkeys", package="movementAnalysis")

djl(monkey.tr)

displacement.direction(monkey.tr, units="compass")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
