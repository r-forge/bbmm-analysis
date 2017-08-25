\name{direction}
\alias{direction}
\alias{ddirection}
\alias{pdirection}
\title{
Compute the (distribution of) movement direction.
}
\description{
Direction (in the linear motion model) and density of the direction distribution (in the BBMM) over specific time intervals.
}
\usage{
direction(tr, time, time.scale = NA)
ddirection(d, tr, time, time.scale = NA)
pdirection(d, tr, time, time.scale = NA, lower = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tr}{The trajectory object}
  \item{time}{A vector of times at which to report the direction}
  \item{time.scale}{The length of the intervals over which to report the direction}
  \item{d}{The direction for which to compute the density}
  \item{lower}{The lower end of the interval for the CDF}
}
\details{
The direction at a time \eqn{t} is defined as the direction of the net displacement over the time interval \eqn{[t-time.scale/2, t+time.scale/2]}. Directions are in radians from the positive \eqn{x} axis.

Since direction has a circular distribution, the distribution function has to be computed with respect to some base value. Using the \code{lower} argument, one can change this base value. Thus, when calling \code{ddirection(x, ..., lower=l)}, the density function will be integrated from \code{l} to \code{x}. The latter may be a vector.
}
\value{
\code{direction} returns a matrix indexed by the ID in the trajectory and time.
\code{ddirection} returns an array indexed by ID, the requested direction and time.
}
\references{
Buchin, K., Sijben, S., Arseneau, T. J. M., AND Willems, E. P.
Detecting movement patterns using brownian bridges. \emph{In Proceedings of
the 20th International Conference on Advances in Geographic Information
Systems}, 119--128.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{velocity}}, \code{\link{bbtraj}}
}
\examples{
data("vervet_monkeys", package="moveBB")

## Compute the (linear model) direction over two half-hour intervals
#d <- direction(monkey.tr, as.POSIXct(c("2011-01-18 15:15:15", "2011-01-19 16:30:00")), 1800)
#d

## Compute the density of four directions over the same time intervals
#ddirection(d, monkey.tr, as.POSIXct(c("2011-01-18 15:15:15", "2011-01-19 16:30:00")), 1800)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
