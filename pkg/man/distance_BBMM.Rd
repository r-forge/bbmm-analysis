\name{distance (BBMM)}
\alias{ddistance}
\alias{pdistance}
\alias{qdistance}
\title{
Compute distance statistics using the BBMM
}
\description{
Density, distribution function, quantile function for the distances between
entities at the given time(s), using the Brownian bridge movement model.
}
\usage{
ddistance(d, tr, time)
pdistance(d, tr, time)
qdistance(p, tr, time)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{Vector of distances}
  \item{p}{Vector of probabilities}
  \item{tr}{The trajectory object}
  \item{time}{Vector of times}
}
\details{
}
\value{
\code{ddistance} computes the density, \code{pdistance} the distribution function and
\code{qdistance} evaluates the quantile function for the requested parameters.

The functions return the requested values for each value of the first parameter
(\code{d} or \code{p}) and for each of the requested times. The value is computed
between each pair of IDs in the trajectory.

This means that the result is a 4 dimensional array, indexed by the IDs involved,
the value of \code{d} or \code{p} and the time.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data("example_data", package="movementAnalysis")
tr <- as.bbtraj(data.frame(x=example_data$X, y=example_data$Y, var=example_data$StdDev^2),
	date=example_data$DateTime, id=example_data$GroupID, burst=example_data$GroupDayNo)

# Compute the 5th and 95th percentile of the distance at two distinct times
qdistance(c(0.05, 0.95), tr, as.POSIXct(c("2011-01-18 15:15:15", "2011-01-19 16:30:00")))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
