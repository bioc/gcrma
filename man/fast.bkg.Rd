\name{fast.bkg}
\alias{fast.bkg}
\alias{mem.bkg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Internal functions for justGCRMA}
\description{
  These are internal functions for justGCRMA that are called based on
  memory or speed constraints.
}
\usage{
fast.bkg(filenames, pm.affinities, mm.affinities, index.affinities,
type, minimum, optical.correct, verbose, k, rho, correction, stretch,
fast, cdfname, read.verbose)
mem.bkg(filenames, pm.affinities, mm.affinities, index.affinities, type,
minimum, optical.correct, verbose, k, rho, correction, stretch, fast,
cdfname, read.verbose)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filenames}{A list of cel files.}
  \item{pm.affinities}{Values passed from \code{compute.affinities}.}
  \item{mm.affinities}{Values passed from \code{compute.affinities}.}
  \item{index.affinities}{Values passed from \code{compute.affinities}.}
  \item{type}{"fullmodel" for sequence and MM model. "affinities" for
          sequence information only. "mm" for using MM without sequence
          information.}
  \item{minimum}{A minimum value to be used for \code{optical.correct}.}
  \item{optical.correct}{Logical value. If \code{TRUE}, optical
    background correction is performed.}
  \item{verbose}{Logical value. If \code{TRUE}, messages about the
    progress of the function are printed.}
  \item{k}{A tuning factor.}
  \item{rho}{correlation coefficient of log background intensity in a pair
    of pm/mm probes. Default=.7}
  \item{correction}{  }
  \item{stretch}{  }
  \item{fast}{Logical value. If \code{TRUE}, then a faster ad hoc
    algorithm is used.}
  \item{cdfname}{Used to specify the name of an alternative cdf package. If set to
    \code{NULL}, the usual cdf package based on Affymetrix' mappings
    will be used. Note that the name should not include the 'cdf' on
    the end, and that the corresponding probe package is also required
    to be installed. If either package is missing an error will result.}
  \item{read.verbose}{Logical value. If \code{TRUE}, a message is
    returned as each celfile is read in.}
}
\details{
  Note that this expression measure is given to you in log base 2
  scale. This differs from most of the other expression measure
  methods.

  The tuning factor 'k' will have different meanings if one uses
  the fast (add-hoc) algorithm or the empirical Bayes approach. See
  Wu et al. (2003)
}
\value{
  An \code{ExpressionSet}.
}
\author{James W. MacDonald <jmacdon@med.umich.edu>}
\seealso{\code{\link{gcrma}}}
\keyword{internal}% at least one, from doc/KEYWORDS
