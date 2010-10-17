\name{gcrma.engine}
\alias{gcrma.engine}
\title{GCRMA background adjust engine(internal function)}
\description{
  This function adjust for non-specific binding when all arrays in the
  dataset share the same probe affinity information. It takes matrices of PM
  probe intensities, MM probe intensities, other negative control probe
  intensities(optional) and the associated probe affinities, and return
  one matrix of non-specific binding corrected PM probe intensities.
}
\usage{
gcrma.engine(pms,mms,ncs=NULL,
                         pm.affinities=NULL,mm.affinities=NULL,anc=NULL,
                         type=c("fullmodel","affinities","mm","constant"),
                         k=6*fast+0.5*(1-fast),
                         stretch=1.15*fast+1*(1-fast),correction=1,GSB.adjust=TRUE,rho=0.7,
                         verbose=TRUE,fast=FALSE)
}
\arguments{
  \item{pms}{The matrix of PM intensities}
  \item{mms}{The matrix of MM intensities}
  \item{ncs}{The matrix of negative control probe intensities. When left
as\code{NULL}, the MMs are considered the negative control probes.}
  \item{pm.affinities}{The vector of PM probe affinities. Note: This can be
shorter than the number of rows in \code{pms} when some probes do not
have sequence information provided.}
 \item{mm.affinities}{The vector of MM probe affinities.}
 \item{anc}{The vector of Negative Control probe affinities. This is
ignored if MMs are used as negative controls (\code{ncs=NULL})}
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
\seealso{gcrma.engine2}
\author{Rafeal Irizarry & Zhijin Wu}
\keyword{manip}

