\name{monkey.tr}
\alias{monkey.tr}
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
	An object of type \code{\link{bbtraj}} created from the data in \code{\link{monkey.data}},
	containing four bursts, one for each \code{(monkey, day)} pair.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\examples{
data(vervet_monkeys)

monkey.tr
plot(monkey.tr)
}
\keyword{datasets}
