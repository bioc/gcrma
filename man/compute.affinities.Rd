\name{compute.affinities}
\alias{compute.affinities}
\title{Probe Affinity computation}
\description{
An internal function to calculate probe affinities from their
sequences.}
\usage{compute.affinities(cdfname,pmonly=FALSE,verbose=TRUE)}
\arguments{
  \item{cdfname}{ Object of class \code{character} representing the name of
          \code{CDF} file associated with the arrays in the \code{AffyBatch}.}
  \item{pmonly}{Logical value. If \code{TRUE}, only affinities of PM
probes are calculated. If \code{FALSE}, both PM and MM probe affinities
are calculated.}
 \item{verbose}{Logical value. If \code{TRUE} messages about the progress of
    the function is printed.} 

}
\details{
The affinity of a probe is described as the sum of position-dependent base affinities.
Each base at each position contributes to the total affinity of a probe
in an additive fashion. For a given type of base, the positional effect
is modeled as a spline function with 5 degrees of freedom.
}

\value{A list of the following components:

  \code{apm}: affinities for PM probes.

  \code{amm}: affinities for MM probes.

  \code{index}: index of the subset of PM probes with sequence available.
}
\references{
Hekstra, D., Taussig, A. R., Magnasco, M., and Naef, F. (2003) Absolute mRNA concentrations from sequence-specific calibration of oligonucleotide array. Nucleic Acids Research, 31. 1962-1968.
}
\seealso{\code{\link{gcrma},\link{affinity.spline.coefs}}}
\author{Rafeal Irizarry}
\keyword{manip}% at least one, from doc/KEYWORDS

