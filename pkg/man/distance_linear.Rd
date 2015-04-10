\name{distance (Linear movement model)}
\alias{distance}
\alias{encounter}
\alias{encounterIntervals}
\title{
Compute distances using the Linear movement model
}
\description{
Compute the distance between trajectories at the requested times. Report when
encounters occur.
}
\usage{
distance(tr, time)
encounter(d, tr, time)
encounterIntervals(d, tr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{Vector of distances}
  \item{tr}{The trajectory object}
  \item{time}{Vector of times}
}
\details{
\code{encounterIntervals} accepts only a single distance.
}
\value{
\code{distance} returns a three dimensional array indexed by two IDs and the time,
each value representing the distance between these two IDs at that time.

\code{encounter} returns whether an encounter occurs for each value of \code{d}
and for each of the requested times. The value is computed between each pair of IDs in the trajectory.
This means that the result is a 4 dimensional array, indexed by the IDs involved,
the value of \code{d} and the time.

\code{encounterIntervals} computes the time intervals in which there is an encounter
between each pair of IDs. In the linear model this can be solved analytically,
which is often faster than requesting whether an encounter occurs for many different
times.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{ddistance}}, \code{\link{pdistance}} and \code{\link{qdistance}}
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
ei
ei[["BD", "NH"]]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
