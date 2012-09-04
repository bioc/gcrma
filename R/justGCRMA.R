### A user friendly wrapper for just.gcrma
justGCRMA <- function(..., filenames=character(0),
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
                     verbose=TRUE, fast=TRUE, minimum=1, optimize.by = c("speed","memory"),
                      cdfname = NULL, read.verbose = FALSE){

   l <- AllButCelsForReadAffy(..., filenames=filenames,
                             widget=widget,
                             celfile.path=celfile.path,
                             sampleNames=sampleNames,
                             phenoData=phenoData,
                             description=description)
  

  ##and now we are ready to read cel files
  return(just.gcrma(filenames=l$filenames,
                    phenoData=l$phenoData,
                    description=l$description,
                    notes=notes,
                    compress=compress,
                    verbose=verbose,
                    normalize=normalize,
                    bgversion=bgversion,
                    affinity.info=affinity.info,
                    type=type, k=k, stretch=stretch,
                    correction=correction, rho=rho,
                    optical.correct=optical.correct,
                    fast=fast, minimum=minimum,
                    optimize.by=optimize.by,
                    cdfname = cdfname,
                    read.verbose = read.verbose))

}


just.gcrma <- function(..., filenames=character(0),
                       phenoData=new("AnnotatedDataFrame"),
                       description=NULL,
                       notes="", compress=getOption("BioC")$affy$compress.cel,
                       normalize=TRUE, bgversion=2, affinity.info=NULL,
                       type=c("fullmodel","affinities","mm","constant"),
                       k=6*fast+0.5*(1-fast), stretch=1.15*fast+1*(1-fast),
                       correction=1, rho=0.7, optical.correct=TRUE,
                       verbose=TRUE, fast=TRUE, minimum=1, optimize.by = c("speed","memory"),
                       cdfname = NULL, read.verbose = FALSE) {

  auxnames <- unlist(list(...))
  filenames <- c(filenames, auxnames)

  checkValidFilenames(filenames)
  
  n <- length(filenames)
  
  ## error if no file name !
  if (n == 0)
    stop("No file name given !")
  
  pdata <- pData(phenoData)
  ##try to read sample names from phenoData. if not there use CEL filenames
  if(dim(pdata)[1]!=n){#if empty pdata filename are samplenames
    #warning("Incompatible phenoData object. Created a new one.\n")
    
    samplenames <- gsub("^/?([^/]*/)*", "", unlist(filenames))
    pdata <- data.frame(sample=1:n,row.names=samplenames)
    phenoData <- new("AnnotatedDataFrame",
                     data = pdata,
                     varMetadata = data.frame(
                     labelDescription = "arbitrary numbering",
                     row.names = "sample"))
  }
  else samplenames <- rownames(pdata)
  

  if (is.null(description))
  {
    description <- new("MIAME")
    description@preprocessing$filenames <- filenames
    description@preprocessing$affyversion <- library(help=affy)$info[[2]][[2]][2]
  }

  ## get information from cdf environment

  headdetails <- read.celfile.header(filenames[[1]])
  dim.intensity <- headdetails[[2]]
  if(is.null(cdfname))
    cdfName <- headdetails[[1]]
  else
    cdfName <- cdfname
  
  type <- match.arg(type)

  pmonly <- (type=="affinities"|type=="constant")
  needaff <- (type=="fullmodel"|type=="affinities")

  if( needaff ){
    if(is.null (affinity.info)){
      affinity.info <- compute.affinities(cdfName,verbose=verbose)
    }
      
    pm.affinities <- pm(affinity.info)
    mm.affinities <- mm(affinity.info)
    
    index.affinities <- which(!is.na(pm.affinities))
    
    pm.affinities <- pm.affinities[index.affinities]
    mm.affinities <- mm.affinities[index.affinities]
    
    ##Recover memory
    rm(affinity.info)
    gc()
    
  }

  speed <- match.arg(optimize.by)
  if(speed == "speed"){
    pms <- fast.bkg(filenames = filenames, pm.affinities = pm.affinities,
                    mm.affinities = mm.affinities, index.affinities = index.affinities,
                    type = type, minimum = minimum, optical.correct = optical.correct,
                    verbose = verbose, k = k, rho = rho, correction = correction,
                    stretch = stretch, fast = fast, cdfname = cdfname,
                    read.verbose = read.verbose)
  }
  if(speed == "memory"){
    pms <- mem.bkg(filenames = filenames, pm.affinities = pm.affinities,
                    mm.affinities = mm.affinities, index.affinities = index.affinities,
                    type = type, minimum = minimum, optical.correct = optical.correct,
                    verbose = verbose, k = k, rho = rho, correction = correction,
                   stretch = stretch, fast = fast, cdfname = cdfname,
                   read.verbose = read.verbose)
  }

 
  tmp <- new("AffyBatch",
             cdfName=cdfName,
             annotation=cleancdfname(cdfName, addcdf=FALSE))
  
  ngenes <- length(featureNames(tmp))

  pNList <- probeNames(tmp)
  pNList <- split(0:(length(pNList) -1), pNList)


  
  exprs <- .Call("rma_c_complete",pms, pNList, ngenes, normalize, FALSE, bgversion, verbose, PACKAGE="affy")

  colnames(exprs) <- samplenames
  se.exprs <- array(NA, dim(exprs))
  
  annotation <- annotation(tmp)
  
  new("ExpressionSet",
      phenoData = phenoData,
      annotation = annotation,
      experimentData = description,
      exprs = exprs, se.exprs = se.exprs)
}

## A function to do background correction fast, but taking more RAM

fast.bkg <- function(filenames, pm.affinities, mm.affinities,
                     index.affinities, type, minimum, optical.correct,
                     verbose, k, rho, correction, stretch, fast, cdfname, read.verbose){
  
  pms <- read.probematrix(filenames=filenames, which="pm", cdfname = cdfname,
                          verbose = read.verbose)$pm
  mms <- read.probematrix(filenames=filenames, which="mm", cdfname = cdfname)$mm

  if(optical.correct){
     if(verbose) cat("Adjusting for optical effect.")
     for (i in 1:ncol(pms)){
       if(verbose) cat(".")
       tmp <- min(c(pms[,i], mms[,i]), na.rm=TRUE)
       pms[,i] <- pms[,i]- tmp + minimum
       mms[,i] <- mms[,i]- tmp + minimum
    }
     if(verbose) cat("Done.\n")
  }
  if(type=="fullmodel" | type=="affinities"){
    set.seed(1)
    Subset <- sample(1:length(pms[index.affinities,]),25000)
    y <- log2(pms)[index.affinities,][Subset]
    Subset <- (Subset-1)%%nrow(pms[index.affinities,])+1
    x <- pm.affinities[Subset]
    fit1 <- lm(y~x)
  }

  if(verbose) cat("Adjusting for non-specific binding")
  for(i in 1:ncol(pms)){
    if(verbose) cat(".")

          
    if(type=="fullmodel"){
      pms[,i] <- bg.adjust.fullmodel(pms[,i],mms[,i],ncs=NULL,
                                     pm.affinities,mm.affinities,anc=NULL,
                                     index.affinities,k=k,
                                     rho=rho,fast=fast)
      pms[,i] <- GSB.adj(pms[,i], index.affinities, pm.affinities, fit1$coef, k)
   #   pms[index.affinities,i] <- 2^(log2(pms[index.affinities,i])-
    #                                fit1$coef[2]*pm.affinities+mean(fit1$coef[2]*pm.affinities))
    }
    if(type=="affinities"){
      pms[,i] <- bg.adjust.affinities(pms[,i],ncs=mms[,i],pm.affinities,mm.affinities,
                                      index.affinities, k=k,
                                      fast=fast)
      pms[,i] <- GSB.adj(pms[,i], index.affinities, pm.affinities, fit1$coef, k)
      
      ## pms[index.affinities,i] <- 2^(log2(pms[index.affinities,i])-
##                                     fit1$coef[2]*pm.affinities + 
##                                    mean(fit1$coef[2]*pm.affinities))
    }
    if(type=="mm") pms[,i] <- bg.adjust.mm(pms[,i],correction*mms[,i],k=k,fast=fast)
    if(type=="constant"){
      pms[,i] <- bg.adjust.constant(pms[,i],k=k,Q=correction*mean(pms[,i]<mms[,i]),fast=fast)
    }
    if(stretch!=1){
      mu <- mean(log(pms[,i]))
      pms[,i] <- exp(mu + stretch*(log(pms[,i])-mu))
    }
  }
  if(verbose) cat("Done.\n")
  return(pms)
}

## A function to do background correction using less RAM

mem.bkg <- function(filenames, pm.affinities, mm.affinities,
                     index.affinities, type, minimum, optical.correct,
                     verbose, k, rho, correction, stretch, fast, cdfname, read.verbose){
  
  pms <- read.probematrix(filenames=filenames, which="pm", cdfname = cdfname,
                          verbose = read.verbose)$pm

  ## tmps used to carry optical correct value to bkg correction loop
   if(optical.correct){
     if(verbose) cat("Adjusting for optical effect.")
     tmps <- NULL
     for (i in 1:ncol(pms)){
       if(verbose) cat(".")
       mm <- read.probematrix(filenames=filenames[i], which="mm",
                              cdfname = cdfname)$mm[,1]
       tmp <-  min(c(pms[,i], mm), na.rm=TRUE)
       pms[,i] <- pms[,i]- tmp + minimum
       tmps <- c(tmps, tmp)
    }
   }
  if(verbose) cat("Done.\n")

  
  if(type=="fullmodel" | type=="affinities"){
    set.seed(1)
    Subset <- sample(1:length(pms[index.affinities,]),25000)
    y <- log2(pms)[index.affinities,][Subset]
    Subset <- (Subset-1)%%nrow(pms[index.affinities,])+1
    x <- pm.affinities[Subset]
    fit1 <- lm(y~x)
  }

  if(verbose) cat("Adjusting for non-specific binding")
  for(i in 1:ncol(pms)){
    if(verbose) cat(".")

    mm <- read.probematrix(filenames=filenames[i], which="mm", cdfname = cdfname)$mm[,1]
    
    if(optical.correct)
      mm <- mm - tmps[i] + minimum
          
    if(type=="fullmodel"){
      pms[,i] <- bg.adjust.fullmodel(pms[,i],mm,ncs=NULL,
                                     pm.affinities,mm.affinities,anc=NULL,
                                     index.affinities,k=k,
                                     rho=rho,fast=fast)
      pms[index.affinities,i] <- 2^(log2(pms[index.affinities,i])-
                                    fit1$coef[2]*pm.affinities+mean(fit1$coef[2]*pm.affinities))
    }
    if(type=="affinities"){
      pms[,i] <- bg.adjust.affinities(pms[,i],ncs=mm,pm.affinities,mm.affinities,
                                      index.affinities, k=k,
                                      fast=fast)
      pms[index.affinities,i] <- 2^(log2(pms[index.affinities,i])-
                                    fit1$coef[2]*pm.affinities + 
                                    mean(fit1$coef[2]*pm.affinities))
    }
    if(type=="mm") pms[,i] <- bg.adjust.mm(pms[,i],correction*mm,k=k,fast=fast)
    if(type=="constant"){
      pms[,i] <- bg.adjust.constant(pms[,i],k=k,Q=correction*mean(pms[,i]<mm),fast=fast)
    }
    if(stretch!=1){
      mu <- mean(log(pms[,i]))
      pms[,i] <- exp(mu + stretch*(log(pms[,i])-mu))
    }
  }
  if(verbose) cat("Done.\n")
  return(pms)
}


  




