\name{position}
\alias{position}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Evaluate position parameters at the given times
}
\description{
Given a trajectory of class \code{\link{bbtraj}} and a list of date/times
(as timestamps or objects of class \code{\link{POSIXct}}), this function
computes the distribution parameters of the position for each group in the
trajectory at each of the requested times. The distribution parameters are
the mean location and the variance of the location.
}
\usage{
position(tr, time)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tr}{An object of class \code{\link{bbtraj}}}
  \item{time}{A list of date/times at which to evaluate the position parameters.}
}
\details{
If a requested time does not have a relocation, the values are interpolated
in the usual way for the Brownian bridge movement model.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
An array indexed by ID, parameter name, time.
}
\references{
%% ~put references to the literature/web site here ~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{bbtraj}}
}
\examples{
data("example_data", package="movementAnalysis")
tr <- as.bbtraj(data.frame(x=example_data$X, y=example_data$Y, var=example_data$StdDev^2),
	date=example_data$DateTime, id=example_data$GroupID, burst=example_data$GroupDayNo)
	
position(tr, as.POSIXct(c("2011-01-19 08:30:00")))

# , , 1295422200
#
#           x       y      var
# BD 322487.6 6900484 1133.653
# NH 323012.9 6900720 1192.534
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line