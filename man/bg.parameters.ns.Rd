\name{bg.parameters.ns}
\alias{bg.parameters.ns}
\alias{left.sigma}
\alias{PAV}
\alias{average.for.PAV}
\title{Estimation of non-specific Binding Background Parameters}
\description{
An internal function to be used by \code{\link{gcrma}}
}
\usage{
bg.parameters.ns(x,affinities,order.aff=NULL,Q=.25,nbreaks=40,
                             monotonize.mu=TRUE,
                             monotonize.sigma=FALSE)
}
\arguments{
  \item{x}{PM or MM intensities after optical background correction, before non-specific-binding correction.}
  \item{affinities}{Probe affinities for probes with known sequences.}
  \item{order.aff}{The order of affinities. If \code{NULL}, it will be
    computed in this function.}
  \item{Q}{Between 0 and 1. For given affinity level, the Q-th quantile
    of log intensities is taken as the mean of the background
    distribution for this affinity level.}
  \item{nbreaks}{The number of strata of affinity levels.}
  \item{monotonize.mu} {Logical value, if \code{TRUE}, estimated mean
    background if forced to be monotonic with probe affinity level.}
   \item{monotonize.sigma} {Logical value, if \code{TRUE}, estimated
     background standard deviation if forced to be monotonic with probe
     affinity level.}
}
%\details{}
\value{a vector of same length as x.}
\seealso{\code{\link{gcrma}}}
\author{Rafeal Irizarry, Zhijin (Jean) Wu}
\keyword{manip}


