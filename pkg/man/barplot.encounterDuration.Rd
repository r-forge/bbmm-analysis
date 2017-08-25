\name{barplot.encounterDuration}
\alias{barplot.encounterDuration}
\title{
Draw a barplot representing the duration of encounters.
}
\description{
Given the result from a call to \code{\link{encounterDuration}}, this function draws a barplot. The plot contains one group of bars for each pair of IDs or bursts in the data (depending on the type of data). Each group contains a bar for each model for which the encounter duration was computed.
}
\usage{
\method{barplot}{encounterDuration}(height, col = NULL, units="auto", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{height}{The object containing the duration of encounters.}
  \item{col}{The colours to use in the plot. By default, takes colours with equally spaced hues in the \code{\link{hcl}} space.}
  \item{units}{The unit used on the y axis of the plot. Must be a unit accepted by \code{\link{difftime}}.}
  \item{\dots}{Other \link{graphical parameters}, e.g. from \code{\link{barplot.default}}.}
}
\details{
If \code{height} is an array, i.e. it was the result from a call to \code{\link{encounterDuration}} with \code{groupBy=NULL}, the name displayed under each pair of bars are the names of the IDs involved, separated by a hyphen.
Otherwise, you may specify the labels yourself in the row names of \code{height} (which is a \code{data.frame}). If all labels are numeric, this function replaces them by the names of the bursts involved, separated by a hyphen.
}
\value{
	The return value is the same as what \code{\link{barplot.default}} returns, i.e. a vector or matrix indicating the coordinates of bar midpoints.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data("vervet_monkeys", package="moveBB")
	
ed <- encounterDuration(monkey.tr, 100, groupBy="GroupID")
barplot(ed)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
