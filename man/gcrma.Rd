\name{gcrma}
\alias{gcrma}
\alias{gcrma.bg.transformation.fast}
\alias{gcrma.bg.transformation}
\alias{gcrma.engine}
\title{Robust Multi-Array expression measure using sequence information}
\description{
  This function converts an \code{AffyBatch} into an \code{exprSet}
  using the robust multi-array average (RMA) expression measure with help of probe sequence.
}
\usage{
gcrma(object,affinity.info=NULL,
                  type=c("fullmodel","affinities","mm","constant"),
                  k=6*fast+0.5*(1-fast),stretch=1.15*fast+1*(1-fast),correction=1,
                  rho=.7,optical.correct=TRUE,verbose=TRUE,fast=TRUE)
}
\arguments{
  \item{object}{an \code{\link[affy:AffyBatch-class]{AffyBatch}}}
  \item{affinity.info}{\code{NULL} or an \code{AffyBatch} containing the
    affinities in the \code{exprs} slot. This object can be created
    using the function \code{\link{compute.affinities}}.}
  \item{type}{"fullmodel" for sequence and MM model. "affinities" for
    sequence information only. "mm" for using MM without sequence
    information.}
  \item{k}{A tuning factor.}
  \item{rho}{correlation coefficient of log background intensity in a pair of pm/mm probes. Default=.7}
  \item{stretch} {.}
  \item{correction}{.}
  \item{optical.correct}{Logical value. If \code{TRUE}, optical
    background correction is performed.}
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
 An \code{exprSet}.
}
\author{Rafeal Irizarry}
\examples{
if(require(affydata) & require(hgu95av2probe)){
     data(Dilution)
     ai <- compute.affinities(cdfName(Dilution))
     Dil.expr<-gcrma(Dilution,affinity.info=ai,type="affinities")
}
}
\keyword{manip}

