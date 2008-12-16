\name{gcrma}
\alias{gcrma}
\alias{gcrma.bg.transformation.fast}
\alias{gcrma.bg.transformation}
\alias{GSB.adj}
\title{Robust Multi-Array expression measure using sequence information}
\description{
  This function converts an \code{AffyBatch} into an \code{ExpressionSet}
  using the robust multi-array average (RMA) expression measure with help of probe sequence.
}
\usage{
gcrma(object,affinity.info=NULL,
      affinity.source=c("reference","local"),NCprobe=NULL,
      type=c("fullmodel","affinities","mm","constant"),
      k=6*fast+0.5*(1-fast),stretch=1.15*fast+1*(1-fast),correction=1,
      GSB.adjust=TRUE,
      rho=.7,optical.correct=TRUE,verbose=TRUE,fast=TRUE,
      subset=NULL,normalize=TRUE,\dots)
}
\arguments{
  \item{object}{an \code{\link[affy:AffyBatch-class]{AffyBatch}}}
  \item{affinity.info}{\code{NULL} or an \code{AffyBatch} containing the
    affinities in the \code{exprs} slot. This object can be created
    using the function \code{\link{compute.affinities}}.}
  \item{affinity.source}{\code{reference}: use the package internal
    Non-specific binding data or \code{local}: use the experimental
    data in \code{object}. If \code{local} is chosen, either MM probes or a user-defined
    list of probes (see \code{NCprobes}) are used to estimate
    affinities.}
  \item{NCprobe}{Index of negative control probes. When set as
    \code{NULL},the MM probes will be used. These probes
    are used to estimate parameters of non-specific binding on each
    array. These will be also used to estimate probe affinity profiles when
    affinity.info is not provided.}
  \item{type}{"fullmodel" for sequence and MM model. "affinities" for
    sequence information only. "mm" for using MM without sequence
    information.}
  \item{k}{A tuning factor.}
  \item{stretch}{.}
  \item{correction}{.}
  \item{GSB.adjust}{Logical value. If \code{TRUE}, probe effects in specific binding will
    be adjusted.}
  \item{rho}{correlation coefficient of log background intensity in a pair of pm/mm probes. Default=.7}
  \item{optical.correct}{Logical value. If \code{TRUE}, optical
    background correction is performed.}
  \item{verbose}{Logical value. If \code{TRUE} messages about the progress of
    the function is printed.} 
  \item{fast}{Logical value. If \code{TRUE} a faster ad hoc algorithm is
    used.}
  \item{subset}{a character vector with the the names of the probesets to be
          used in expression calculation.}
  \item{normalize}{logical value. If 'TRUE' normalize data using quantile
    normalization.}
   \item{\dots}{further arguments to be passed (not currently implemented -
          stub for future use).}
}
\details{
  Note that this expression measure is given to you in log base 2
  scale. This differs from most of the other expression measure
  methods.

  The tuning factor \code{k} will have different meanings if one uses
  the fast (add-hoc) algorithm or the empirical Bayes approach. See Wu
  et al. (2003)    
}
\value{
 An \code{ExpressionSet}.
}
\author{Rafeal Irizarry}
\examples{
if(require(affydata) & require(hgu95av2probe) & require(hgu95av2cdf)){
     data(Dilution)
     ai <- compute.affinities(cdfName(Dilution))
     Dil.expr<-gcrma(Dilution,affinity.info=ai,type="affinities")
}
}
\keyword{manip}
