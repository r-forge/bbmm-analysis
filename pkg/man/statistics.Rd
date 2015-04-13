\name{statistics}
\title{
Collect statistics about relocations
}
\description{
%These functions compute a certain statistic of every burst in a given trajectory. The statistics are net displacement distance, direction of net displacement and day journey length (DJL), which is the total length of all relocations.

}
\usage{
statistics(object, add=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{The trajectory object}
  \item{add}{If \const{TRUE}, the statistics will be added to the object's data attribute.}
}
\details{
%All functions use the linear motion model to compute the statistics.
%
%\code{displacement.direction} can return its result in different units. The default is in radians from the positive \eqn{x} axis, but it can also report degrees from this axis or compass heading (i.e. degrees from north, measured clockwise).
The statistics computed are: displacement in the x and y dimensions and total displacement since last observation, total displacement since start of burst, direction of relocation and turning angle at this observation.
}
\value{
A \code{data.frame} containing the statistics for each relocation, or a copy of \code{object} with the statistics added.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data("vervet_monkeys", package="movementAnalysis")

statistics(monkey.tr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
