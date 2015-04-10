\name{velocity}
\alias{velocity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute distribution parameters for velocity
}
\description{
The velocity over any time interval is a bivariate normal distribution. This function computes the parameters for this distribution. The velocity in the linear model is identical to the mean velocity, so ignore the variance in that case.
}

\usage{
velocity(tr, time, time.scale = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tr}{The trajectory object}
  \item{time}{Vector of times}
  \item{time.scale}{The length of the intervals over which to report the direction}
}
\details{
The velocity at a time \eqn{t} is defined as the net rate of displacement over the time interval \eqn{[t-time.scale/2, t+time.scale/2]}, i.e. \deqn{v(t, \Delta t) = \frac{x(t + \Delta t / 2) - x(t - \Delta t / 2)}{\Delta t}.}{v(t, \Delta t) = (x(t + \Delta t / 2) - x(t - \Delta t / 2)) / \Delta t.}
}
\value{
Returns a three dimensional array indexed by ID in the trajectory, variable name and time. The variables are named \code{x}, \code{y} and \code{var} for the \eqn{x} and \eqn{y} components of the mean and the variance respectively.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{bbtraj}}
}
\examples{
data("vervet_monkeys", package="movementAnalysis")

# Compute the velocity distribution parameters over two half-hour intervals
velocity(monkey.tr, as.POSIXct(c("2011-01-18 15:15:15", "2011-01-19 16:30:00")), 1800)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
