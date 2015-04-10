\name{monkey.data}
\alias{monkey.data}
\docType{data}
\title{
Trajectory data of two monkeys over two days
}
\description{
This data set contains GPS data obtained from two vervet monkeys over two days.
The data is sampled at one sample per hour from 5:00 to 19:00, although there
are missing values.
}
\usage{data(vervet_monkeys)}
\format{
	A \code{data.frame} containing the raw data. The fields are:
	\describe{
		\item{GroupID}{The identifier of the monkey.}
		\item{GroupDayNo}{The identifier of the monkey and the day the measurement was obtained.}
		\item{DateTime}{The time of the measurement.}
		\item{X}{X coordinate of the location.}
		\item{Y}{Y coordinate of the location.}
		\item{StdDev}{The standard deviation of the location error.}
	}
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\examples{
data(vervet_monkeys)
head(monkey.data)
}
\keyword{datasets}
