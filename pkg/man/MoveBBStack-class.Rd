\name{MoveBBStack-class}
\Rdversion{1.1}
\docType{class}
\alias{MoveBBStack-class}
\alias{[,MoveBBStack,ANY,ANY-method}
\alias{diffusionCoefficient,MoveBBStack-method}
\alias{position,MoveBBStack,ANY-method}
\alias{speedDistribution,MoveBBStack,RasterLayer,numeric,numeric-method}
\alias{split,MoveBBStack,missing-method}
\alias{utilizationDistribution,MoveBBStack,RasterLayer,numeric-method}
\alias{velocity,MoveBBStack,numeric,numeric-method}

\title{Class \code{"MoveBBStack"}}
\description{
A collection of trajectories, which can separately be stored using the \code{\linkS4class{Move}} class.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MoveBBStack", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{trackId}:}{Object of class \code{"factor"} ~~ }
    \item{\code{timestamps}:}{Object of class \code{"POSIXct"} ~~ }
    \item{\code{idData}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{sensor}:}{Object of class \code{"factor"} ~~ }
    \item{\code{data}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{coords.nrs}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{coords}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{bbox}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{proj4string}:}{Object of class \code{"CRS"} ~~ }
    \item{\code{trackIdUnUsedRecords}:}{Object of class \code{"factor"} ~~ }
    \item{\code{timestampsUnUsedRecords}:}{Object of class \code{".OptionalPOSIXct"} ~~ }
    \item{\code{sensorUnUsedRecords}:}{Object of class \code{"factor"} ~~ }
    \item{\code{dataUnUsedRecords}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{dateCreation}:}{Object of class \code{"POSIXct"} ~~ }
    \item{\code{study}:}{Object of class \code{"character"} ~~ }
    \item{\code{citation}:}{Object of class \code{"character"} ~~ }
    \item{\code{license}:}{Object of class \code{"character"} ~~ }
    \item{\code{variance}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{diffusion}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{MoveStack}"}, directly.
Class \code{"\linkS4class{.BBInfo}"}, directly.
Class \code{"\linkS4class{.MoveTrackStack}"}, by class "MoveStack", distance 2.
Class \code{"\linkS4class{.MoveGeneral}"}, by class "MoveStack", distance 2.
Class \code{"\linkS4class{.MoveTrack}"}, by class "MoveStack", distance 3.
Class \code{"\linkS4class{.unUsedRecordsStack}"}, by class "MoveStack", distance 3.
Class \code{"\linkS4class{SpatialPointsDataFrame}"}, by class "MoveStack", distance 4.
Class \code{"\linkS4class{.unUsedRecords}"}, by class "MoveStack", distance 4.
Class \code{"\linkS4class{SpatialPoints}"}, by class "MoveStack", distance 5.
Class \code{"\linkS4class{Spatial}"}, by class "MoveStack", distance 6.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "MoveBBStack", i = "ANY", j = "ANY")}: ... }
    \item{diffusionCoefficient}{\code{signature(tr = "MoveBBStack")}: ... }
    \item{position}{\code{signature(object = "MoveBBStack", time = "ANY")}: ... }
    \item{speedDistribution}{\code{signature(object = "MoveBBStack", raster = "RasterLayer", timestepSize = "numeric", time.scale = "numeric")}: ... }
    \item{split}{\code{signature(x = "MoveBBStack", f = "missing")}: ... }
    \item{utilizationDistribution}{\code{signature(object = "MoveBBStack", raster = "RasterLayer", timestepSize = "numeric")}: ... }
    \item{velocity}{\code{signature(object = "MoveBBStack", time = "numeric", time.scale = "numeric")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("MoveBBStack")

data(vervet_monkeys)
monkey.tr
}
\keyword{classes}
