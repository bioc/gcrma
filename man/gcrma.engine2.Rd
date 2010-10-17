\name{gcrma.engine2}
\alias{gcrma.engine2}
\title{GCRMA background adjust engine(internal function)}
\description{
  This function adjust for non-specific binding when each array has its
  own probe affinity information. It takes an AffyBatch object of probe
  intensities and an AffyBatch of probe  affinity,
  returns one matrix of non-specific binding corrected PM probe intensities.
}
\usage{
gcrma.engine2(object,pmIndex=NULL,mmIndex=NULL,
              NCprobe=NULL,affinity.info,
              type=c("fullmodel","affinities","mm","constant"),
              k=6*fast+0.5*(1-fast),
              stretch=1.15*fast+1*(1-fast),correction=1,GSB.adjust=TRUE,rho=0.7,
              verbose=TRUE,fast=TRUE)
}
\arguments{
  \item{object}{an \code{\link[affy:AffyBatch-class]{AffyBatch}}. Note: this is an internal
    function. Optical noise should have been corrected for. }
  \item{pmIndex}{Index of PM probes.This will be computed within the
    function if left \code{NULL}}
  \item{mmIndex}{Index of MM probes.This will be computed within the
    function if left \code{NULL}}
  \item{NCprobe}{Index of negative control probes. When set as
    \code{NULL},the MM probes will be used. These probes
    are used to estimate parameters of non-specific binding on each
    array. These will be also used to estimate probe affinity profiles when
    affinity.info is not provided.}
  \item{affinity.info}{\code{NULL} or an \code{AffyBatch} containing the
    affinities in the \code{exprs} slot. This object can be created
    using the function \code{\link{compute.affinities}}.}
  \item{type}{"fullmodel" for sequence and MM model. "affinities" for
    sequence information only. "mm" for using MM without sequence
    information.}
  \item{k}{A tuning factor.}
  \item{stretch}{.}
  \item{correction}{.}
  \item{GSB.adjust}{Logical value. If \code{TRUE}, probe effects in specific binding will
    be adjusted.}
  \item{rho}{correlation coefficient of log background intensity in a pair of pm/mm probes. Default=.7}
  \item{verbose}{Logical value. If \code{TRUE} messages about the progress of
    the function is printed.} 
  \item{fast}{Logicalvalue. If \code{TRUE} a faster add-hoc algorithm is
    used.}
}
\details{
  Note that this expression measure is given to you in log base 2
  scale. This differs from most of the other expression measure
  methods.

  The tunning factor \code{k} will have different meainngs if one uses
  the fast (add-hoc) algorithm or the empirical bayes approach. See Wu
  et al. (2003)
     
}
\value{
 A matrix of PM intensties. 
}
\seealso{gcrma.engine}
\author{Rafeal Irizarry & Zhijin Wu}
\keyword{manip}

