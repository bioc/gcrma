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
  \item{k}{A tuning parameter. See details}
  \item{Q}{A number between 0 and 1 that determines what quantile to use
  as an estimate of the mean background noise.}
  \item{fast}{Logical value. If \code{TRUE} a faster add-hoc algorithm is
    used.}
}
\details{ Assumes PM=background1+signal,mm=background2,
  (log(background1),log(background2))' 
  follow bivariate normal distribution, signal distribution follows power
  law. 
  \code{bg.parameters.gcrma} and \code{sg.parameters.gcrma} 
  provide adhoc estimates of the parameters.

  the original gcrma uses an emprical bayes estimate. this requiers a
  complicated numerical integration. An add-hoc method tries to immitate
  the empirical bayes estimate with a PM-B but values of PM-B<\code{k}
  going to \code{k}. This can be thought as a shrunken MVUE. For more
  details see Wu et al. (2003). 
}


\value{a vector of same length as x.}
\seealso{\code{\link{gcrma}}}
\author{Rafeal Irizarry, Zhijin(Jean) Wu}
\keyword{manip}% at least one, from doc/KEYWORDS

