getCDF <- function(cdfname, lib=.libPaths()[1], verbose=TRUE){

  options(show.error.messages = FALSE)
  attempt <- try(do.call("library", list(cdfname, lib.loc=lib)))
  options(show.error.messages = TRUE)
  if( inherits(attempt, "try-error")){
    require(reposTools) || stop("Package 'reposTools' is required",
                              " for this operation.")
    if(testBioCConnection()){
      if(verbose)
        print("Checking to see if your internet connection works...")
                                        # Check for file permissions
      if(file.access(lib, mode = 0) < 0) stop(paste("Directory", lib,"does not",
                            "seem to exist.\n", "Please check your 'lib' parameter",
                            "and try again."))
      
      if(file.access(lib, mode = 2) < 0) stop(paste("You do not have write access to", lib,
                            "\nPlease check your permissions or provide",
                            "a different 'lib' parameter."))

      z <- install.packages2(cdfname, lib=lib)
      if(! cdfname %in% updatedPkgs(z)) {
        stop(paste("Environment", cdfname,"was not found in the Bioconductor",
                                                   "repository."))
      }else{
        if(verbose)
          print(paste("Installation of environment", cdfname,
                      "was successful."))
      }
    }else{
      stop(paste("The current operation could not access",
                 "the Bioconductor repository. Please",
                 "check your internet connection, and",
                 "report further problems to",
                 "bioconductor@stat.math.ethz.ch"))
    }
  }
    # Now load the library
    do.call("library", list(cdfname, lib.loc=lib))
    #Check that library is loaded
    if(!cdfname %in% .packages()) stop(paste("The package", cdfname,
                                             "could not be loaded."))
}


getProbePackage <- function(probepackage, lib=.libPaths()[1], verbose=TRUE){

  options(show.error.messages = FALSE)
  attempt <- try(do.call("library", list(probepackage, lib.loc=lib)))
  options(show.error.messages = TRUE)
  if( inherits(attempt, "try-error")){
    require(reposTools) || stop("Package 'reposTools' is required",
                              " for this operation.")
    if(testBioCConnection()){
      if(verbose)
        print("Checking to see if your internet connection works...")
                                        # Check for file permissions
      if(file.access(lib, mode = 0) < 0) stop(paste("Directory", lib,"does not",
                            "seem to exist.\n", "Please check your 'lib' parameter",
                            "and try again."))
      
      if(file.access(lib, mode = 2) < 0) stop(paste("You do not have write access to", lib,
                            "\nPlease check your permissions or provide",
                            "a different 'lib' parameter."))

      z <- install.packages2(probepackage, lib=lib)
      if(! probepackage %in% updatedPkgs(z)){
        stop(paste("Environment", probepackage,"was not found in the Bioconductor",
                                                   "repository."))
      }else{
        if(verbose)
          print(paste("Installation of environment", probepackage,
                      "was successful."))
      }
    }else{
      stop(paste("The current operation could not access",
                 "the Bioconductor repository. Please",
                 "check your internet connection, and",
                 "report further problems to",
                 "bioconductor@stat.math.ethz.ch"))
    }
  }
    # Now load the library
    do.call("library", list(probepackage, lib.loc=lib))
    #Check that library is loaded
    if(!probepackage %in% .packages()) stop(paste("The package", probepackage,
                                             "could not be loaded."))
}


      
    
        

  
