\name{MoveBB-class}
\Rdversion{1.1}
\docType{class}
\alias{MoveBB-class}
\alias{[,MoveBB,ANY,ANY-method}
\alias{diffusionCoefficient,MoveBB-method}
\alias{position,MoveBB,numeric-method}
\alias{speedDistribution,MoveBB,RasterLayer,numeric,numeric-method}
\alias{utilizationDistribution,MoveBB,RasterLayer,numeric-method}
\alias{velocity,MoveBB,numeric,numeric-method}

\title{Class \code{"MoveBB"}}
\description{
The \code{MoveBB} object contains at least time and coordinate information of an animal, along with location variance for each relocation and information about the diffusion coefficient. It can contain further data as for class \code{\linkS4class{Move}}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MoveBB", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
  	\item{idData}{Object of class \code{"data.frame"}: additional (one row) data;  }
    \item{dateCreation}{Object of class \code{"numeric"}: time stamp when the file was downloaded;  }
    \item{study}{Object of class \code{"character"}: name of the study;  }
    \item{citation}{Object of class \code{"character"}: how to cite the study;  }
    \item{license}{Object of class \code{"character"}: the license under which the data were published;  }
    \item{timesMissedFixes}{Object of class \code{"POSIXct"}: stores the timestamps of lines of the data set that were removed because they included NA locations }
    \item{bbox}{belongs to the SpatialPointsDataFrame}
    \item{coords}{coordinates of the track, belongs to the SpatialPointsDataFrame}
    \item{coords.nrs}{belongs to the SpatialPointsDataFrame}
    \item{data}{additional data of that object that is stored in the SpatialPointsDataFrame}
    \item{proj4string}{projection of the coordinates}
    \item{timestamps}{timestamps according to the coordinates}
    %
    \item{\code{variance}:}{location variance for each relocation}
    \item{\code{diffusion}:}{function providing diffusion coefficient over time}
  }
}
\section{Extends}{
Class \code{"\linkS4class{Move}"}
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "MoveBB", i = "ANY", j = "ANY")}: Select specific relocations }
    \item{diffusionCoefficient}{\code{signature(tr = "MoveBB")}: Optimize diffusion coefficient for trajectory }
    \item{position}{\code{signature(object = "MoveBB", time = "numeric")}: Compute position distribution at given times }
    \item{velocity}{\code{signature(object = "MoveBB", time = "numeric", time.scale = "numeric")}: Compute velocity distribution at given times }
    \item{utilizationDistribution}{\code{signature(object = "MoveBB", raster = "RasterLayer", timestepSize = "numeric")}: Compute utilization distribution for trajectory }
    \item{speedDistribution}{\code{signature(object = "MoveBB", raster = "RasterLayer", timestepSize = "numeric", time.scale = "numeric")}: Compute spatial distribution of expected speed }
	 }
}

\details{
	The main way in which \code{MoveBB} extends \code{\linkS4class{Move}} is that it provides explicit support for location uncertainty, e.g. using the Brownian bridge movement model (BBMM).
	In particular, every relocation has not only a position and a timestamp, but also stores the variance of the location, e.g. GPS error.
	Further, an object of type \code{MoveBB} stores the \emph{diffusion coefficient} of the Brownian motion, or a similar parameter for other movement models.
	The diffusion coefficient is allowed to vary with time and is therefore modeled as a function that maps time to the value of the diffusion coefficient. The default estimation function just computes a constant.
}

\references{
Horne, J., Garton, E., Krone, S. and Lewis, J. Analyzing animal movements
using Brownian bridges. \emph{Ecology 88}, 9 (2007), 2354--2363.
}

\author{
Stef Sijben
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\linkS4class{MoveBBStack}}, \code{\linkS4class{move:Move}}
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("MoveBB")

data(vervet_monkeys)
split(monkey.tr)
}
\keyword{classes}
