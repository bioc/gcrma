\name{bg.adjust.affinities}
\alias{bg.adjust.affinities}
\alias{bg.adjust.mm}
\alias{bg.adjust.fullmodel}
\alias{bg.adjust.constant}
\alias{bg.adjust.optical}

\title{Background adjustment with sequence information (internal function)}
\description{
An internal function to be used by \code{\link{gcrma}}.
}
\usage{
bg.adjust.affinities(x,affinities,index=seq(along=x),k,Q=0.25,fast=TRUE)
}

\arguments{
  \item{x}{PM intensities after optical background correction, before non-specific-binding correction.}
  \item{index}{The index of x with known sequences. (For some types of
    arrays the sequences of a small subset of probes are not provided by
  Affymetrix.)}
  \item{affinities}{Probe affinities for probes with known sequences.}
  \item{k}{A tuning parameter.}
  \item{Q}{}
  \item{fast}{Logicalvalue. If \code{TRUE} a faster add-hoc algorithm is
    used.}
}
\details{ Assumes PM=background1+signal,mm=background2, (log(background1),log(background2))'
follow bivariate normal distribution, signal distribution follows power law.
\code{bg.parameters.gcrma} and \code{sg.parameters.gcrma}
provide adhoc estimates of the parameters.}
\value{a vector of same length as x.}
\seealso{\code{\link{gcrma}}}
\author{Rafeal Irizarry, Zhijin(Jean) Wu}
\keyword{manip}% at least one, from doc/KEYWORDS

