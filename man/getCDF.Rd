\name{getCDF}
\alias{getCDF}
\alias{getProbePackage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Functions for Automatic Download of Packages }
\description{
  These are internal functions that are called by \code{justGCRMA} and
  \code{GCRMA} in order to automatically download and install cdf
  environments and probe packages.
}
\usage{
getCDF(cdfname, lib = .libPaths()[1], verbose = TRUE)
getProbePackage(probepackage, lib = .libPaths()[1], verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cdfname}{ Name of the cdfenv to install. }
  \item{probepackage}{ Name of the probe package to install. }
  \item{lib}{ Directory of the R library where packages will be installed. }
  \item{verbose}{ Output informative comments? Defaults to \code{TRUE}}
}

\value{
  Nothing is returned. These functions are called simply to install
  environments. 
}
\author{ James W. MacDonald, based on \code{getCDFinfo}, written by Jeff
Gentry.}

\seealso{ \code{getCDFinfo} 
  %\code{\link[affy]{getCDFinfo}} there is no such function in affy
  %called getCDFinfo 
}

\keyword{ internal }
