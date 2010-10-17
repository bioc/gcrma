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
bg.adjust.fullmodel(pms,mms,ncs=NULL,apm,amm,anc=NULL,index.affinities,k=6
* fast + 0.25 * (1 - fast),rho=.7,fast=FALSE)
bg.adjust.affinities(pms,ncs,apm,anc,index.affinities,k=6
* fast + 0.25 * (1 - fast),fast=FALSE,nomm=FALSE)
}
\arguments{
  \item{pms}{PM intensities after optical background correction, before
    non-specific-binding correction.}
  \item{mms}{MM intensities after optical background correction, before
    non-specific-binding correction.}
  \item{ncs}{Negative control probe intensities after optical background correction, before
    non-specific-binding correction. If \code{ncs=NULL}, the MM probes
    are considered the negative control probes.}
  \item{index.affinities}{The index of pms with known sequences. (For some types of
    arrays the sequences of a small subset of probes are not provided by
    Affymetrix.)}
  \item{apm}{Probe affinities for PM probes with known sequences.}
  \item{amm}{Probe affinities for MM probes with known sequences.}
  \item{anc}{Probe affinities for Negative control probes with known
    sequences. This is ignored when \code{ncs=NULL}.}
  \item{rho}{correlation coefficient of log background intensity in a pair of pm/mm probes. Default=.7}
  \item{k}{A tuning parameter. See details.}
%  \item{Q}{A number between 0 and 1 that determines what quantile to use
%  as an estimate of the mean background noise.}
  \item{fast}{Logical value. If \code{TRUE} a faster add-hoc algorithm is used.}
  \item{nomm}{Logical value indicating if MM intensities are available
    and will to be used to estimate background.}
}
\details{Assumes PM=background1+signal,mm=background2,
  (log(background1),log(background2))' 
  follow bivariate normal distribution, signal distribution follows power
  law. 
  \code{bg.parameters.gcrma} and \code{sg.parameters.gcrma} 
  provide adhoc estimates of the parameters.

  the original gcrma uses an empirical Bayes estimate. this requires a
  complicated numerical integration. An add-hoc method tries to imitate
  the empirical Bayes estimate with a PM-B but values of PM-B<\code{k}
  going to \code{k}. This can be thought as a shrunken MVUE. For more
  details see Wu et al. (2003). 
}
\value{a vector of same length as x.}
\seealso{\code{\link{gcrma}}}
\author{Rafeal Irizarry, Zhijin(Jean) Wu}
\keyword{manip}% at least one, from doc/KEYWORDS
