\name{bg.adjust.gcrma}
\alias{bg.adjust.gcrma}
\title{GCRMA background adjust (internal function)}
\description{
  This function performs background adjustment (optical noise and
  non-specific binding on an \code{AffyBatch}
  project and returns an \code{AffyBatch} object in which the PM
  intensities are adjusted.
}
\usage{
bg.adjust.gcrma(object,affinity.info=NULL,
      affinity.source=c("reference","local"),
      NCprobe=NULL,
      type=c("fullmodel","affinities","mm","constant"),
      k=6*fast+0.5*(1-fast),stretch=1.15*fast+1*(1-fast),correction=1,
      GSB.adjust=TRUE,
      rho=.7,optical.correct=TRUE,verbose=TRUE,fast=TRUE)
}
\arguments{
  \item{object}{an \code{\link[affy:AffyBatch-class]{AffyBatch}}}
  \item{affinity.info}{\code{NULL} or an \code{AffyBatch} containing the
    affinities in the \code{exprs} slot. This object can be created
    using the function \code{\link{compute.affinities}}.}
  \item{affinity.source}{\code{reference}: use the package internal
    Non-specific binding data or \code{local}: use the experimental
    data in \code{object}. If \code{local} is chosen, either MM probes or a user-defined
    list of probes (see \code{NCprobes}) are used to estimate affinities.}
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
}
\details{
  The returned value is an \code{AffyBatch} object, in which the PM probe intensities
  have been background adjusted. The rest is left the same as the
  starting \code{AffyBatch} object.
  
  The tunning factor \code{k} will have different meainngs if one uses
  the fast (ad hoc) algorithm or the empirical bayes approach. See Wu
  et al. (2003)
    
}
\value{
 An \code{AffyBatch}.
}
\author{Rafeal Irizarry}
\examples{
 if(require(affydata) & require(hgu95av2probe) & require(hgu95av2cdf)){
          data(Dilution)
          ai <- compute.affinities(cdfName(Dilution))
          Dil.adj<-bg.adjust.gcrma(Dilution,affinity.info=ai,type="affinities")
     }
}
\keyword{manip}

