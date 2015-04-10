\name{barplot.encounterDuration}
\alias{barplot.encounterDuration}
\title{
Draw a barplot representing the duration of encounters.
}
\description{
Given the result from a call to \code{\link{encounterDuration}}, this function draws a barplot. The plot contains one group of bars for each pair of IDs or bursts in the data (depending on the type of data). Each group contains a bar for each model for which the encounter duration was computed.
}
\usage{
barplot.encounterDuration(ed,
    col = rainbow(ifelse(is.data.frame(ed), ncol(ed) - 4, dim(ed)[3])),
    ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ed}{The object containing the duration of encounters.}
  \item{col}{The colours to use in the plot.}
  \item{\dots}{Other \link{graphical parameters}, e.g. from \code{\link{barplot.default}}.}
}
\details{
If \code{ed} is an array, i.e. it was the result from a call to \code{\link{encounterDuration}} with \code{byburst=FALSE}, the name displayed under each pair of bars are the names of the IDs involved, separated by a hyphen.
Otherwise, you may specify the labels yourself in the row names of \code{ed} (which is a data.frame). If all labels are numeric, this function replaces them by the names of the bursts involved, separated by a hyphen.
}
\value{
	The return value is the same as what \code{\link{barplot.default}} returns, i.e. a vector or matrix indicating the coordinates of bar midpoints.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data("example_data", package="movementAnalysis")
tr <- as.bbtraj(data.frame(x=example_data$X, y=example_data$Y, var=example_data$StdDev^2),
	date=example_data$DateTime, id=example_data$GroupID, burst=example_data$GroupDayNo)
	
ed <- encounterDuration(tr, 100, byburst=TRUE)
barplot(ed)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
