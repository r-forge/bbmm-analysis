\name{plot.utilizationDistribution}
\alias{plot.utilizationDistribution}
\alias{contour.utilizationDistribution}
\title{
Plot or draw contours for utilization distributions.
}
\description{
Plots the given utilization distribution, either by applying a colourmap to the values or by drawing contours where the value passes certain levels.
}
\usage{
\method{plot}{utilizationDistribution}(x, col = hcl(1:length(x) * 360/length(x), 50, 70),
    xlim = NULL, ylim = NULL, add = FALSE, ...)
\method{contour}{utilizationDistribution}(x, levels = 0.99, col = hcl(1:length(x) * 360/length(x), l=35), 
    xlim = NULL, ylim = NULL, labels = levels, add = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The utilization distribution to plot.}
  \item{col}{The colours to use in the plot. May be a vector of colours or a list containing the colour scheme for each ID.}
  \item{xlim}{The minimum and maximum x coordinate to plot.}
  \item{ylim}{The minimum and maximum y coordinate to plot.}
  \item{add}{Logical specifying if plot should be added to an already existing plot; defaults to FALSE.}
  \item{levels}{The level(s) at which to draw contours.}
  \item{labels}{The labels to print with the contours. Use \code{drawlabels=FALSE} to hide them altogether.}
  \item{\dots}{Other \link{graphical parameters}, e.g. from \code{\link{plot.default}} or \code{\link{contour.default}}.}
}
\details{
If \code{col} is a vector, it is assumed to contain a colour code for each grid to be plotted. The colour code for each value of the UD is derived by varying the opacity of the provided colour. If the length of \code{col} does not match the number of UDs to plot, it is recycled as usual.

If \code{col} is a list, it is assumed to contain one element for each UD to plot, specifying exactly which colours to use in that plot.
}
%\value{
%
%}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{utilizationDistribution}}, \code{\link{SpatialPolygons}}
}
\examples{
data("vervet_monkeys", package="movementAnalysis")

# Compute the UD and plot the result for one ID
ud <- utilizationDistribution(monkey.tr)
plot(ud)

# Specify nicer colours
plot(ud, col=c("red","green"))
# Add a contour to the existing plot in the same colours
contour(ud, col=c("red","green"), add=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
