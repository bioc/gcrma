getCDF <- function(cdfname, lib=.libPaths()[1], verbose=TRUE){

  ## 1-23-2006 JWM: removed dependency on reposTools
  
  options(show.error.messages = FALSE)
  attempt <- try(do.call(library, list(cdfname, lib.loc = lib)))
  options(show.error.messages = TRUE)
  if (inherits(attempt, "try-error")){
    
    if(testBioCConnection()){
      if(verbose)
        print("Checking to see if your internet connection works...")
      ## Check for file permissions
      if(file.access(lib, mode = 0) < 0) stop(paste("Directory", lib,"does not",
                            "seem to exist.\n", "Please check your 'lib' parameter",
                            "and try again."))
      
      if(file.access(lib, mode = 2) < 0) stop(paste("You do not have write access to", lib,
                            "\nPlease check your permissions or provide",
                            "a different 'lib' parameter."))
      
      biocContribUrl <- sapply(repositories(), contrib.url)
      biocPkgs <- available.packages(biocContribUrl)
      if(!cdfname %in% biocPkgs[,"Package"]){
        if(verbose)
          print(paste("Environment",cdfname,
                      "was not found in the Bioconductor",
                      "repository."))
        return(list(paste("Bioconductor -",cdfname,"not available")))
      }else{
        install.packages(cdfname, lib=lib,
                         repos=repositories(),
                         dependencies=TRUE)
      }
    }else{
      stop(paste("The current operation could not access",
                 "the Bioconductor repository. Please",
                 "check your internet connection, and",
                 "report further problems to",
                 "bioconductor@stat.math.ethz.ch"))
    }
  }
  ## Now load the library
  do.call(library, list(cdfname, lib.loc=lib))
  ##Check that library is loaded
  if(!cdfname %in% .packages()) stop(paste("The package", cdfname,
                                           "could not be loaded."))
}


getProbePackage <- function(probepackage, lib=.libPaths()[1], verbose=TRUE){

  ## 1-23-2006 JWM: Removed dependency on reposTools
  options(show.error.messages = FALSE)
  attempt <- try(do.call(library, list(probepackage, lib.loc = lib)))
  options(show.error.messages = TRUE)
  if (inherits(attempt, "try-error")){
    
    if(testBioCConnection()){
      if(verbose)
        print("Checking to see if your internet connection works...")
      ## Check for file permissions
      if(file.access(lib, mode = 0) < 0) stop(paste("Directory", lib,"does not",
                            "seem to exist.\n", "Please check your 'lib' parameter",
                            "and try again."))
      
      if(file.access(lib, mode = 2) < 0) stop(paste("You do not have write access to", lib,
                            "\nPlease check your permissions or provide",
                            "a different 'lib' parameter."))
      
      biocContribUrl <- sapply(repositories(), contrib.url)
      biocPkgs <- available.packages(biocContribUrl)
      if(!probepackage %in% biocPkgs[,"Package"]){
        if(verbose)
          print(paste("Environment",probepackage,
                      "was not found in the Bioconductor",
                      "repository."))
        return(list(paste("Bioconductor -",probepackage,"not available")))
      }else{
        install.packages(probepackage, lib=lib,
                         repos=repositories(),
                         dependencies=TRUE)
      }
    }else{
      stop(paste("The current operation could not access",
                 "the Bioconductor repository. Please",
                 "check your internet connection, and",
                 "report further problems to",
                 "bioconductor@stat.math.ethz.ch"))
    }
  }
  
  ## Now load the library
  do.call(library, list(probepackage, lib.loc=lib))
  ##Check that library is loaded
  if(!probepackage %in% .packages()) stop(paste("The package", probepackage,
                                                "could not be loaded."))
}







