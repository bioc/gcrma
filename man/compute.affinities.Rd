\name{compute.affinities}
\alias{compute.affinities}
\alias{compute.affinities2}
\alias{check.probes}
\alias{base.profiles}
\alias{base.profiles.mm}
\alias{base.profiles.nc}
\alias{compute.affinity.coef}
\alias{compute.affinities.local}
\alias{plotBaseProfiles}

\title{Probe Affinity computation}
\description{
An internal function to calculate probe affinities from their
sequences.}
\usage{
compute.affinities(cdfname,verbose=TRUE)
compute.affinities2(cdfname,verbose=TRUE)
check.probes(probepackage,cdfname)
}
\arguments{
  \item{cdfname}{ Object of class \code{character} representing the name of
          \code{CDF} file associated with the arrays in the \code{AffyBatch}.}
\item{probepackage}{\code{character} representing the name of the
  package with the probe sequence information.}
 \item{verbose}{Logical value. If \code{TRUE} messages about the progress of
    the function is printed.} 

}
\details{
The affinity of a probe is described as the sum of position-dependent base affinities.
Each base at each position contributes to the total affinity of a probe
in an additive fashion. For a given type of base, the positional effect
is modeled as a spline function with 5 degrees of freedom.

Use \code{compute.affinities2} if there are no MM probes.

\code{check.probes} makes sure things are matching as they should.
}

\value{\code{compute.affinities} returns an \code{AffyBatch} with the
  affinities for PM probes in the pm locations and the affinities for MM
  probes in the mm locations. NA  will be added for probes with no
  sequence information.
}
\references{
Hekstra, D., Taussig, A. R., Magnasco, M., and Naef, F. (2003) Absolute mRNA concentrations from sequence-specific calibration of oligonucleotide array. Nucleic Acids Research, 31. 1962-1968.
}
\seealso{\code{\link{gcrma},\link{affinity.spline.coefs}}}
\author{Rafeal Irizarry}
\keyword{manip}% at least one, from doc/KEYWORDS

