\name{gcrma}
\alias{gcrma}
\title{Robust Multi-Array expression measure using sequence information}
\description{
  This function converts an `AffyBatch' into an `exprSet' using the
     robust multi-array average (RMA) expression measure with help of probe sequence.
}
\usage{
gcrma(object,estimate="eb",summary.method = "rlm",summary.subset=geneNames(object)[1:10]),...)
}
\arguments{
   \item{object}{an \code{\link[affy:AffyBatch-class]{AffyBatch}}}
  \item{estimate}{"eb" or "mle". An character indicating background correction method.
"eb" for empirical Bayes, "mle" for maximum likelihood} 
  \item{summary.subset}{An vector of gene names with the the probe sets do to be
    converted to expression. }
 \item{normalize}{logical value. If \code{TRUE} normalize data using
    normalize.method (default as quantile normalization)}
  \item{rho}{correlation coefficient of log background intensity in a pair of pm/mm probes. Default=.8}
  \item{baseline}{For estimate="mle" only. The lowest log(expression) value allowed. Default=.25}
  \item{lower.bound}{For estimate="eb" only. Lower bound of signal intensity allowed a priori. Default=1.}
  \item{step}{a parameter in approximating "eb" estimates.}
}
\details{
  Note that this expression measure is given to you in log base 2
     scale. This differs from most of the other expression measure
     methods.
}
\value{
 An `exprSet'
}
\author{Jean Wu, Rafeal Irizarry}


\examples{
require(affydata)
data(Dilution)
Dil.expr<-gcrma(Dilution,estimate="mle",summary.method = "medianpolish",summary.subset=geneNames(Dilution)[1:10])

}
\keyword{manip}
