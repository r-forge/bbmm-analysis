\name{as.bbtraj}
\alias{as.bbtraj}
\alias{bbtraj}
\alias{bbFilterNA}
%\alias{print.bbtraj}
%\alias{summary.bbtraj}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Working with Trajectories in 2D Space: the Class bbtraj
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
The class \code{bbtraj} stores trajectories of animals, similar to the class \code{\link{ltraj}} in package \code{\link{adehabitatLT}}. The difference is that a bbtraj object also stores the variances used when dealing with Brownian bridges. When converting a data set to an object of this class (using the function \code{as.bbtraj}) the most likely value for the diffusion coefficient is calculated and stored in the resulting object.
\code{bbFilterNA} filters the missing measurements from a trajectory, removing bursts that are empty after this filtering.
}
\usage{
as.bbtraj(xys, date = NULL, id, burst = id, typeII = TRUE,
          slsp = c("remove", "missing"), make.regular=FALSE)
bbFilterNA(tr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xys}{a data.frame containing the x and y coordinates and variance
    of the relocations}
  \item{date}{for trajectories of type II, a vector of class \code{POSIXct}
    giving the date for each relocation. For trajectories of type I, this
    argument is not taken into account.}
  \item{id}{either a character string indicating the identity of the
    animal or a factor with length equal to \code{nrow(xy)} }
  \item{burst}{either a character string indicating the identity of the
    burst of relocations or a factor with length equal to
    \code{nrow(xy)}}
  \item{typeII}{logical. \code{TRUE} indicates a trajectory of type II
    (time recorded, e.g. radio-tracking), whereas \code{FALSE} indicates
    a trajectory of type I (time not recorded, e.g. sampling of tracks in
    the snow)}
  \item{slsp}{a character string used for the computation of the
    turning angles (for details, see \code{\link{as.ltraj}})}
  \item{make.regular}{Either \code{FALSE}, or a list specifying the arguments to \code{\link{setNA}} and \code{\link{sett0}}.}
  \item{tr}{A trajectory of class \code{ltraj} or \code{bbtraj}}
}
\details{
If \code{make.regular} is not \code{FALSE}, it must be a list containing parameters
accepted by \code{setNA} and \code{sett0}, by the same names or in the proper order.
as these functions are used to make the resulting trajectory regular.
}
\value{
  An object of class \code{bbtraj} is a list with one component per burst of
  relocations.  Each component is a data frame with two attributes:
  the attribute \code{"id"} indicates the identity of the animal, and
  the attribute \code{"burst"} indicates the identity of the
  burst. Each data frame stores the following columns:
\item{x}{the x coordinate for each relocation}
\item{y}{the y coordinate for each relocation}
\item{diff.coeff}{The diffusion coefficient for the move. This has no
  meaning if the previous relocation is missing, but usually it is set
  to the same value for many relocations.}
\item{loc.var}{The variance of the measured position}
\item{date}{the date for each relocation (type II) or a vector of
  integer giving the order of the relocations in the trajectory.}
\item{dx}{the increase of the move in the x direction.  At least two
  successive relocations are needed to compute \code{dx}.  Missing
  values are returned otherwise.}
\item{dy}{the increase of the move in the y direction.  At least two
  successive relocations are needed to compute \code{dy}.  Missing
  values are returned otherwise.}
\item{dist}{the length of each move.  At least two
  successive relocations are needed to compute \code{dist}.  Missing
  values are returned otherwise.}
\item{dt}{the time interval between successive relocations}
\item{R2n}{the squared net displacement between the current relocation
  and the first relocation of the trajectory}
\item{abs.angle}{the angle between each move and the x axis.  At least two
  successive relocations are needed to compute \code{abs.angle}.  Missing
  values are returned otherwise.}
\item{rel.angle}{the turning angles between successive moves.  At least
  three successive relocations are needed to compute \code{rel.angle}.
  Missing values are returned otherwise.}
}
\references{
%% ~put references to the literature/web site here ~
Horne, J., Garton, E., Krone, S. and Lewis, J. Analyzing animal movements
using Brownian bridges. \emph{Ecology 88}, 9 (2007), 2354--2363.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{ltraj}}, \code{\link{adehabitatLT}}
}
\examples{
data("example_data", package="movementAnalysis")
example_data

tr <- as.bbtraj(data.frame(x=example_data$X, y=example_data$Y, var=example_data$StdDev^2),
    date=example_data$DateTime, id=example_data$GroupID, burst=example_data$GroupDayNo)
tr
is.regular(tr)

# Now, make the trajectory regular while constructing it
tr <- as.bbtraj(data.frame(x=example_data$X, y=example_data$Y, var=example_data$StdDev^2),
    date=example_data$DateTime, id=example_data$GroupID, burst=example_data$GroupDayNo,
    make.regular=list(date.ref=strptime("00:00", "\%H:\%M"), dt=1, units="hour"))
is.regular(tr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
