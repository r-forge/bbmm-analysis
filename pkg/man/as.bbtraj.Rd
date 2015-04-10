\name{as.bbtraj}
\alias{as.bbtraj}
\alias{bbtraj}
\alias{diffusionCoefficient}

\alias{[.bbtraj}
\alias{[<-.bbtraj}
\alias{na.omit.bbtraj}
\alias{bbFilterNA}
%\alias{print.bbtraj}
\alias{summary.bbtraj}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Working with Trajectories in 2D Space: the Class bbtraj
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
The class \code{bbtraj} stores trajectories of animals, similar to the class \code{\link{ltraj}} in package \code{adehabitatLT}. The difference is that a bbtraj object also stores the variances used when dealing with Brownian bridges. When converting a data set to an object of this class (using the function \code{as.bbtraj}) the most likely value for the diffusion coefficient is calculated and stored in the resulting object.
\code{na.omit} filters the missing measurements from a trajectory, removing bursts that are empty after this filtering. All information that is derived from multiple measurements is recomputed. \code{bbFilterNA} is a deprecated alias for \code{na.omit}.
}
\usage{
as.bbtraj(xys, date, id, burst = id, typeII = TRUE,
          slsp =  c("remove", "missing"),
          infolocs = data.frame(pkey = paste(id, date, sep=".")),
          make.regular=FALSE)
\method{na.omit}{bbtraj}(object, slsp=c("remove", "missing"), ...)
diffusionCoefficient(tr, byburst = FALSE, nsteps = 1000)
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
  \item{infolocs}{ If not \code{NULL}, a data.frame containing extra information about the relocations.}
  \item{make.regular}{Either \code{FALSE}, or a list specifying the arguments to \code{\link{setNA}} and \code{\link{sett0}}.}
  \item{object}{The trajectory from which to filter missing values.}
  \item{...}{Further arguments that might be required by specialized methods.}
  \item{tr}{A trajectory of class \code{ltraj} or \code{bbtraj}}
  \item{byburst}{A logical indicating whether a different diffusion coefficient should be estimated for each burst.}
  \item{nsteps}{The number of different values for the diffusion coefficient to test. The range in which the value can be is determined automatically.}
}
\details{
If \code{make.regular} is not \code{FALSE}, it must be a list containing parameters
accepted by \code{setNA} and \code{sett0}, by the same names or in the proper order.
as these functions are used to make the resulting trajectory regular.

\code{diffusionCoefficient} can be used to estimate the diffusion coefficient for an already existing \code{bbtraj} object. When calling \code{as.bbtraj}, it is invoked automatically. However, one might want to re-estimate the diffusion coefficient after modifying some values or constructing a trajectory from scratch.
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
  meaning if the next relocation is missing, but usually it is set
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
  
  \code{diffusionCoefficient} returns a vector containing the optimal values
  for each ID or burst.
}
\references{
%% ~put references to the literature/web site here ~
Horne, J., Garton, E., Krone, S. and Lewis, J. Analyzing animal movements
using Brownian bridges. \emph{Ecology 88}, 9 (2007), 2354--2363.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link[adehabitatLT]{ltraj}}
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

# Extract a single burst from the trajectory
tr[["BD_1417"]]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
