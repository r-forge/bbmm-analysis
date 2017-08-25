\name{encounterDuration}
\alias{encounterDuration}
\alias{encounterDurationById}
\title{
Compute duration of encounters between groups
}
\description{
Computes the duration of encounters between each pair of groups over the whole
measurement period of a trajectory. The user can select what movement model(s)
to apply and whether the result should be on the level of bursts or IDs.
\code{encounterDurationById} returns the duration of encounters for each pair of
IDs, given the duration of encounters for each pair of bursts. This allows to
obtain both types of result without recomputing everything.
}
\usage{
encounterDuration(tr, threshold, model = c("BBMM", "linear"), byburst = FALSE, timestepSize = 60)
encounterDurationById(encounterDurationByBurst)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tr}{The trajectory to analyze}
  \item{threshold}{The maximum distance at which an encounter is detected}
  \item{model}{The movement models for which to compute the encounter duration}
  \item{byburst}{If \code{TRUE}, the result contains durations for pairs of bursts.
  	If \code{FALSE}, the result contains encounter durations between pairs of IDs.}
  \item{timestepSize}{If a movement model requires numerical integration in the time
  	dimension, this is the size of each time step.}
  \item{encounterDurationByBurst}{The result of calling \code{encounterDuration} with \code{byburst=TRUE}.}
}
\details{
Since the duration of encounters is a random variable in the Brownian bridge
movement model, this function cannot give exact results there. Instead, it
reports the expected duration of encounters in the BBMM.

If you already have the encounter duration between bursts and you also want the
encounter duration between IDs, you should use \code{encounterDurationById}
instead of calling \code{encounterDuration} again with different parameters. They
give identical results, but the former is much faster since it does not recompute
all relevant durations.
}
\value{
If \code{byburst == FALSE}, the result is a \code{data.frame} with one row
for each interesting pair of bursts. A pair of bursts is interesting if they
overlap in time, since otherwise the duration is always zero.
The result contains the following fields:
\item{id1,id2}{The IDs of the bursts involved}
\item{burst1,burst2}{The names of the bursts involved}
\item{}{In addition there is one column for each model requested, named after
	the model. These columns are of type \code{\link{difftime}}, with unit \code{secs}.
	They contain the encounter duration according to that model.}

If \code{byburst == TRUE}, the result is a 3 dimensional array, indexed by
the two IDs and the movement model. The entries in the array are again of class
\code{\link{difftime}} with unit \code{secs}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{bbtraj}}
}
\examples{
data("vervet_monkeys", package="moveBB")

ed <- encounterDuration(monkey.tr, 100, byburst=TRUE)
barplot(ed)

barplot(encounterDurationById(ed))

# This one has the same result, but needs to recompute everything
barplot(encounterDuration(monkey.tr, 100))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
