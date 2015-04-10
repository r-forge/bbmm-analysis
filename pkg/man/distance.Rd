\name{distance}
\alias{distance}
\alias{ddistance}
\alias{pdistance}
\alias{qdistance}
\alias{encounter}
\alias{encounterIntervals}
\title{
Compute distance statistics
}
\description{
Density, distribution function, quantile function for the distances between
entities at the given time(s), using the Brownian bridge movement model.
Distance between entities and whether encounters occur in the linear movement model.
}
\usage{
distance(tr, time)

ddistance(d, tr, time)
pdistance(d, tr, time)
qdistance(p, tr, time)

encounter(d, tr, time)
encounterIntervals(d, tr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{Vector of distances}
  \item{p}{Vector of probabilities}
  \item{tr}{The trajectory object}
  \item{time}{Vector of times}
}
\details{
\code{encounterIntervals} accepts only a single distance.
}
\value{
\code{distance} returns a three dimensional array indexed by two IDs and the time,
each value representing the distance in the linear movement model between these two IDs at that time.

\code{encounter} returns whether an encounter occurs for each value of \code{d}
and for each of the requested times. The value is computed between each pair of IDs in the trajectory.
This means that the result is a 4 dimensional array, indexed by the IDs involved,
the value of \code{d} and the time.

\code{encounterIntervals} computes the time intervals in which there is an encounter
between each pair of IDs. In the linear model this can be solved analytically,
which is often faster than requesting whether an encounter occurs for many different
times.

\code{ddistance} computes the density, \code{pdistance} the distribution function and
\code{qdistance} evaluates the quantile function for the requested parameters. \code{pdistance} is equivalent to the probability that an encounter occurs, given a threshold distance.

The functions return the requested values for each value of the first parameter
(\code{d} or \code{p}) and for each of the requested times. The value is computed
between each pair of IDs in the trajectory.

This means that the result is a 4 dimensional array, indexed by the IDs involved,
the value of \code{d} or \code{p} and the time.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{bbtraj}}
}
\examples{
data("example_data", package="movementAnalysis")
tr <- as.bbtraj(data.frame(x=example_data$X, y=example_data$Y, var=example_data$StdDev^2),
	date=example_data$DateTime, id=example_data$GroupID, burst=example_data$GroupDayNo)

# Compute the distance at two distinct times
distance(tr, as.POSIXct(c("2011-01-18 15:15:15", "2011-01-19 16:30:00")))

# See whether an encounter occurs at the 50m or 100m thesholds
encounter(c(10, 100), tr, as.POSIXct(c("2011-01-18 15:15:15", "2011-01-19 16:30:00")))

# Compute when there are encounters
ei <- encounterIntervals(100, tr)
ei[["BD","NH"]]

# Compute the 5th and 95th percentile of the distance at two distinct times
qdistance(c(0.05, 0.95), tr, as.POSIXct(c("2011-01-18 15:15:15", "2011-01-19 16:30:00")))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
