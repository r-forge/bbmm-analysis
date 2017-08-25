\name{speed}
\alias{speed}
\alias{mu.speed}
\alias{dspeed}
\alias{pspeed}
\alias{qspeed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get statistics about speed
}
\description{
Get speed over any time intervals in the linear model. Get density, distribution
function, quantiles and expected value in BBMM.
}
\usage{
speed(tr, time, time.scale = NA)
mu.speed(tr, time, time.scale=NA)

pspeed(v, tr, time, time.scale=NA)
dspeed(v, tr, time, time.scale=NA)
qspeed(p, tr, time, time.scale=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{Vector of speeds}
  \item{p}{Vector of probabilities}
  \item{tr}{The trajectory object}
  \item{time}{Vector of times}
  \item{time.scale}{The length of the intervals over which to report the speed or direction}
}
\details{
The speed at a time \eqn{t} is defined as the absolute value of the \code{\link{velocity}} over the time interval \eqn{[t-time.scale/2, t+time.scale/2]}.
}
\value{
\code{speed} and \code{mu.speed} return a matrix indexed by the ID in the trajectory and time.
\code{dspeed}, \code{pspeed} and \code{qspeed} return an array indexed by ID, the requested speed or probability and time.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{velocity}}, \code{\link{bbtraj}}
}
\examples{
data("vervet_monkeys", package="moveBB")

# Compute the (linear model) speed over two half-hour intervals
v <- speed(monkey.tr, as.POSIXct(c("2011-01-18 15:15:15", "2011-01-19 16:30:00")), 1800)
v

# Compute the density of four speeds over the same time intervals
dspeed(v, monkey.tr, as.POSIXct(c("2011-01-18 15:15:15", "2011-01-19 16:30:00")), 1800)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
