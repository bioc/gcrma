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
  \item{summary.method}{"medianpolish" or "rlm". Method to extract the expression value of a probeset from  background adjusted intensities. }
 \item{normalize}{logical value. If \code{TRUE} normalize data using
    normalize.method (default as quantile normalization)}
  \item{rho}{correlation coefficient of log background intensity in a pair of pm/mm probes. Default=.8}
  \item{triple.goal} {logical value. If \code{TRUE} Shen and Louis'
 triple-goal adjustment is used.}
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
\references{Wei Shen and Thomas Louis.  Triple-goal estimates in
  two-stage hierarchical models. J.R.Statist.Soc.B(1998)}
   


\examples{
require(affydata)
data(Dilution)
Dil.expr<-gcrma(Dilution,estimate="mle",summary.method = "medianpolish",summary.subset=geneNames(Dilution)[1:10])

}
\keyword{manip}
