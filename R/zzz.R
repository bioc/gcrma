.First.lib <- function(libname, pkgname, where) {
  library.dynam("gcrma", pkgname, libname)  

  where <- match(paste("package:", pkgname, sep=""), search())

  require(affy, quietly=TRUE) ##Biobase uses methods

  cacheMetaData(as.environment(where))

}
