\name{justGCRMA}
\alias{justGCRMA}
\alias{just.gcrma}
\title{Compute GCRMA Directly from CEL Files}
\description{
  This function converts CEL files into an \code{ExpressionSet}
  using the robust multi-array average (RMA) expression measure with help of probe sequences.
}
\usage{
            just.gcrma(\dots, filenames=character(0),
                       phenoData=new("AnnotatedDataFrame"),
                       description=NULL,
                       notes="", compress=getOption("BioC")$affy$compress.cel,
                       normalize=TRUE, bgversion=2, affinity.info=NULL,
                       type=c("fullmodel","affinities","mm","constant"),
                       k=6*fast+0.5*(1-fast), stretch=1.15*fast+1*(1-fast),
                       correction=1, rho=0.7, optical.correct=TRUE,
                       verbose=TRUE, fast=TRUE, minimum=1, optimize.by =
                       c("speed","memory"), 
                       cdfname = NULL, read.verbose = FALSE)

            justGCRMA(\dots, filenames=character(0),
                     widget=getOption("BioC")$affy$use.widgets,
                     compress=getOption("BioC")$affy$compress.cel,
                     celfile.path=getwd(),
                     sampleNames=NULL,
                     phenoData=NULL,
                     description=NULL,
                     notes="",
                     normalize=TRUE, 
                     bgversion=2, affinity.info=NULL,
                     type=c("fullmodel","affinities","mm","constant"),
                     k=6*fast+0.5*(1-fast), stretch=1.15*fast+1*(1-fast),
                     correction=1, rho=0.7, optical.correct=TRUE,
                     verbose=TRUE, fast=TRUE, minimum=1,
                     optimize.by = c("speed","memory"),
                     cdfname = NULL, read.verbose = FALSE)
}
\arguments{
  \item{\dots}{file names separated by comma.}
  \item{filenames}{file names in a character vector.}
  \item{widget}{a logical specifying if widgets should be used.}
  \item{compress}{are the CEL files compressed?}
  \item{phenoData}{a \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}} object.}
  \item{description}{a \code{\link[Biobase:class.MIAME]{MIAME}} object.}
  \item{notes}{notes.}
  \item{affinity.info}{\code{NULL} or a list of three components:
    apm,amm and index, for PM probe affinities, MM probe affinities,
    the index of probes with known sequence, respectively.}
  \item{type}{"fullmodel" for sequence and MM model. "affinities" for
    sequence information only. "mm" for using MM without sequence
    information.}
  \item{k}{A tuning factor.}
  \item{rho}{correlation coefficient of log background intensity in a pair of pm/mm probes.
    Default=.7.}
  \item{stretch}{.}
  \item{correction}{.}
  \item{normalize}{Logical value. If \code{TRUE}, then normalize data using
    quantile normalization.}
  \item{optical.correct}{Logical value. If \code{TRUE}, then optical
    background correction is performed.}
  \item{verbose}{Logical value. If \code{TRUE}, then messages about the progress of
    the function is printed.} 
  \item{fast}{Logical value. If \code{TRUE}, then a faster add-hoc algorithm is used.}
  \item{optimize.by}{"speed" will use a faster algorithm but more RAM, and "memory"
    will be slower, but require less RAM.}
  \item{bgversion}{integer value indicating which RMA background to use 1: use
          background similar to pure R rma background given in affy
          version 1.0 - 1.0.2 2: use background similar to pure R rma
          background given in affy version 1.1 and above.} 
  \item{minimum}{.}
  \item{celfile.path}{a character denoting the path 'ReadAffy' should look for
          cel files.}
  \item{sampleNames}{a character vector of sample names to be used in the
    'AffyBatch'.}
  \item{cdfname}{Used to specify the name of an alternative cdf package. If set to
     \code{NULL}, the usual cdf package based on Affymetrix' mappings
     will be used. Note that the name should not include the 'cdf' on
     the end, and that the corresponding probe package is also required
     to be installed. If either package is missing an error will
     result.}
   \item{read.verbose}{Logical value. If \code{TRUE}, then messages will be
     printed as each celfile is read in.}
}
\details{
  This method should require much less RAM than the conventional
  method of first creating an \code{AffyBatch} and then running
  \code{\link[gcrma]{gcrma}}.

  This is a simpler version than \code{gcrma}, so some of the arguments
  available in \code{gcrma} are not available here. For example, it is
  not possible to use the MM probes to estimate background. Instead, the
  internal NSB estimates are used (which is also the default for \code{gcrma}).
  
  Note that this expression measure is given to you in log base 2
  scale. This differs from most of the other expression measure
  methods.

  The tuning factor \code{k} will have different meanings if one uses
  the fast (add-hoc) algorithm or the empirical Bayes approach. See Wu
  et al. (2003)

  \code{fast.bkg} and \code{mem.bkg} are two internal functions.
}
\value{
 An \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} object.
}
\author{James W. MacDonald}
\keyword{manip}
