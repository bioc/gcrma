\name{justGCRMA}
\alias{justGCRMA}
\alias{just.gcrma}
\title{Compute GCRMA Directly from CEL Files}
\description{
  This function converts CEL files into an \code{exprSet}
  using the robust multi-array average (RMA) expression measure with help of probe sequences.
}
\usage{
just.gcrma <- function(..., filenames=character(0),
                       phenoData=new("phenoData"),
                       description=NULL,
                       notes="", background=FALSE,
                       compress=getOption("BioC")$affy$compress.cel,
                       normalize=TRUE, bgversion=2, affinity.info=NULL,
                       type=c("fullmodel","affinities","mm","constant"),
                       k=6*fast+0.5*(1-fast), stretch=1.15*fast+1*(1-fast),
                       correction=1, rho=0.7, optical.correct=TRUE,
                       verbose=TRUE, fast=TRUE, minimum=1)

justGCRMA <- function(..., filenames=character(0),
                     widget=getOption("BioC")$affy$use.widgets,
                     compress=getOption("BioC")$affy$compress.cel,
                     celfile.path=getwd(),
                     sampleNames=NULL,
                     phenoData=NULL,
                     description=NULL,
                     notes="",
                     normalize=TRUE, background=FALSE, 
                     bgversion=2, affinity.info=NULL,
                     type=c("fullmodel","affinities","mm","constant"),
                     k=6*fast+0.5*(1-fast), stretch=1.15*fast+1*(1-fast),
                     correction=1, rho=0.7, optical.correct=TRUE,
                     verbose=TRUE, fast=TRUE, minimum=1)
}
\arguments{
  \item{\dots}{ file names separated by comma.}
  \item{filenames}{ file names in a character vector.}
  \item{widget}{a logical specifying if widgets should be used.}
  \item{compress}{are the CEL files compressed ?}
  \item{phenoData}{ a \code{\link[Biobase]{phenoData}} object}
  \item{description}{ a \code{\link[Biobase]{MIAME}} object}
  \item{notes}{ notes }
  \item{affinity.info}{\code{NULL} or a list of three components:
   apm,amm and index, for PM probe
    affinities, MM probe affinities, the index of probes with known
    sequence, respectively.}
  \item{type}{"fullmodel" for sequence and MM model. "affinities" for
    sequence information only. "mm" for using MM without sequence
    information.}
  \item{k}{A tuning factor.}
  \item{rho}{correlation coefficient of log background intensity in a pair of pm/mm probes. Default=.7}
  \item{stretch} {.}
  \item{correction}{.}
  \item{normalize}{logical value. If \code{TRUE} normalize data using
    quantile normalization}
  \item{optical.correct}{Logical value. If \code{TRUE}, optical
    background correction is performed.}
  \item{verbose}{Logical value. If \code{TRUE} messages about the progress of
    the function is printed.} 
  \item{fast}{Logicalvalue. If \code{TRUE} a faster add-hoc algorithm is
    used.}
}
\details{
  This method should require much less RAM than the conventional
  method of first creating an \code{AffyBatch} and then running
  \code{\link[gcrma]{gcrma}}.
  
  Note that this expression measure is given to you in log base 2
  scale. This differs from most of the other expression measure
  methods.

  The tuning factor \code{k} will have different meanings if one uses
  the fast (add-hoc) algorithm or the empirical bayes approach. See Wu
  et al. (2003)
     
}
\value{
 An \code{exprSet}.
}
\author{James W. MacDonald}
\examples{

}
\keyword{manip}

