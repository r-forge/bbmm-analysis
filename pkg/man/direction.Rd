\name{direction}
\alias{direction}
\alias{ddirection}
\title{
Compute the (distribution of) movement direction.
}
\description{
Direction (in the linear motion model) and density of the direction distribution (in the BBMM) over specific time intervals.
}
\usage{
direction(tr, time, time.scale = NA)
ddirection(d, tr, time, time.scale = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tr}{The trajectory object}
  \item{time}{A vector of times at which to report the direction}
  \item{time.scale}{The length of the intervals over which to report the direction}
  \item{d}{The direction for which to compute the density}
}
\details{
The direction at a time \eqn{t} is defined as the direction of the net displacement over the time interval \eqn{[t-time.scale/2, t+time.scale/2]}. Directions are in radians from the positive \eqn{x} axis.
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

\section{Warning}{
    \bold{The definition of the velocity variance may still contain errors, leading to incorrect results. Use at your own risk.}
}

\seealso{
\code{\link{velocity}}, \code{\link{bbtraj}}
}
\examples{
data("example_data", package="movementAnalysis")
tr <- as.bbtraj(data.frame(x=example_data$X, y=example_data$Y, var=example_data$StdDev^2),
	date=example_data$DateTime, id=example_data$GroupID, burst=example_data$GroupDayNo)

# Compute the (linear model) direction over two half-hour intervals
d <- direction(tr, as.POSIXct(c("2011-01-18 15:15:15", "2011-01-19 16:30:00")), 1800)
d

# Compute the density of four directions over the same time intervals
ddirection(d, tr, as.POSIXct(c("2011-01-18 15:15:15", "2011-01-19 16:30:00")), 1800)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
