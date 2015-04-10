\name{moveBB}
\docType{methods}
\alias{moveBB-methods}
\alias{moveBB,ANY,ANY,missing,ANY,ANY-method}
\alias{moveBB,ANY,ANY,numeric,ANY,missing-method}
\alias{moveBB,numeric,numeric,numeric,POSIXct,data.frame-method}

\title{Create a \code{\linkS4class{MoveBB}} object}
\description{
This method creates an object of type \code{\linkS4class{MoveBB}} or \code{\linkS4class{MoveBBStack}}.
}

\usage{
\S4method{moveBB}{ANY,ANY,missing,ANY,ANY,ANY}(x,y,var,time,data,proj,...)
\S4method{moveBB}{ANY,ANY,numeric,ANY,missing,ANY}(x,y,var,time,data,proj,...)
\S4method{moveBB}{numeric,numeric,numeric,POSIXct,data.frame,ANY}(x,y,var,time,data,proj,sensor='unknown',animal='unnamed',...)
}

\arguments{
\item{x}{vector of x coordinates(e.g. \code{data$x})}
\item{y}{vector of y coordinates}
\item{var}{vector of location variances}
\item{time}{column indicator for non-Movebank data for the time stamps, with POSIXct conversion, i.e. \code{as.POSIXct(data$timestamp, format="\%Y-\%m-\%d \%H:\%M:\%S", tz="UTC")}}
\item{data}{Optional extra data associated with the relocations, if empty it is filled with the coordinates and timestamps}
\item{proj}{projection method; requires a valid CRS (see \code{\link{CRS-class}}) object, like CRS("+proj=longlat +ellps=WGS84"); default is NA}
\item{sensor}{sensor name, either single character or a vector with length of the number of coordinates}
\item{animal}{animal ID or name, either single character or a vector with length of the number of coordinates}
\item{...}{Additional arguments}
}

\author{Stef Sijben}

\seealso{
	\code{\link{moveBBStack}}, \code{\linkS4class{MoveBB}}, \code{\linkS4class{MoveBBStack}}
}
\examples{
	data(vervet_monkeys)
	
	head(monkey.data)
	
	## Generates an object of type MoveBBStack, where the individual trajectories
	## are identified by the value of the "GroupDayNo" field.
	tr <- moveBB(monkey.data$X, monkey.data$Y, monkey.data$StdDev^2,
	        monkey.data$DateTime, animal=monkey.data$GroupDayNo, 
	        data=monkey.data, proj="+proj=utm +zone=36 +south ellps=WGS84")
	
	split(tr) ## Get a list of MoveBB objects
}

\keyword{methods}
%\keyword{ ~~ other possible keyword(s) ~~ }
