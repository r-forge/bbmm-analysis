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
\method{plot}{utilizationDistribution}(ud, col = hcl(1:length(ud) * 360/length(ud), 50, 70),
    xlim = NULL, ylim = NULL, add = FALSE, ...)
\method{contour}{utilizationDistribution}(ud, levels = 0.99, col = col=hcl(1:length(ud) * 360/length(ud), l=35), 
    xlim = NULL, ylim = NULL, labels = levels, add = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ud}{The utilization distribution to plot.}
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
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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

# Define grid lines: equally spaced between the min and max coordinate in tr
xlim <- range(unlist(sapply(tr, function(b) { b$x })), na.rm=TRUE)
xc <- seq(xlim[1], xlim[2], length.out=100)

ylim <- range(unlist(sapply(tr, function(b) { b$y })), na.rm=TRUE)
yc <- seq(ylim[1], ylim[2], length.out=100)

# Compute the UD and plot the result for one ID
ud <- utilizationDistribution(tr, xc=xc, yc=yc)
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
