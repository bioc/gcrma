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
bg.parameters.ns(x,affinities,affinities2=NULL,affinities3=NULL,span=.2)
}
\arguments{
  \item{x}{PM or MM intensities after optical background correction, before non-specific-binding correction.}
  \item{affinities}{Probe affinities for probes with known
    sequences.Used to estimate the function between non-specific binding
  and affinities.}
\item{affinities2}{
   Probe affinities for the probes whoes expected non-specific binding
   intensity is to be predicted.
}  
 \item{affinities3}{
   Probe affinities for another extra group of probes whoes expected non-specific binding
   intensity is to be predicted.
 }
 \item{span}{The span parameter passed to loess function}
}
%\details{}
\value{a vector of same length as x.}
\seealso{\code{\link{gcrma}}}
\author{Rafeal Irizarry, Zhijin (Jean) Wu}
\keyword{manip}


