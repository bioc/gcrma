### A user friendly wrapper for just.gcrma
justGCRMA <- function(..., filenames=character(0),
                     widget=getOption("BioC")$affy$use.widgets,
                     compress=getOption("BioC")$affy$compress.cel,
                     celfile.path=getwd(),
                     sampleNames=NULL,
                     phenoData=NULL,
                     description=NULL,
                     notes="", normalize=TRUE, 
                     bgversion=2, affinity.info=NULL,
                     type=c("fullmodel","affinities","mm","constant"),
                     k=6*fast+0.5*(1-fast), stretch=1.15*fast+1*(1-fast),
                     correction=1, rho=0.7, optical.correct=TRUE,
                     verbose=TRUE, fast=TRUE, minimum=1){
  ##first figure out filenames
  auxnames <- unlist(as.list(substitute(list(...)))[-1])
  
  if (widget){
    require(tkWidgets)
    widgetfiles <- fileBrowser(textToShow="Choose CEL files",
                               testFun=hasSuffix("[cC][eE][lL]"))
  }
  else
    widgetfiles <- character(0)
  
  filenames <- .Primitive("c")(filenames, auxnames, widgetfiles)
  
  if(length(filenames)==0) filenames <- list.celfiles(celfile.path,full.names=TRUE)
  
  if(length(filenames)==0) stop("No cel filenames specified and no cel files in specified directory:",celfile.path,"\n")
  
  
  ##now assign sampleNames if phenoData not given
  if(is.null(phenoData)){
    if(is.null(sampleNames)){
      if(widget){
        require(tkWidgets)
        tksn <- tkSampleNames(filenames=filenames)
        sampleNames <- tksn[,1]
        ##notice that a description of the files is ingored for now
        ##soon to go into MIAME
      }
      else{
        sampleNames <- sub("^/?([^/]*/)*", "", filenames, extended=TRUE)
      }
    }
    else{
      if(length(sampleNames)!=length(filenames)){
        warning("sampleNames not same length as filenames. Using filenames as sampleNames instead\n")
        sampleNames <- sub("^/?([^/]*/)*", "", filenames, extended=TRUE)
      }
    }
  }
  
  ##now get phenoData
  if(is.character(phenoData)) ##if character read file
    phenoData <- read.phenoData(filename=phenoData)
  else{
    if(class(phenoData)!="phenoData"){
      if(widget){
        require(tkWidgets)
        phenoData <- read.phenoData(sampleNames=sampleNames,widget=TRUE)
      }
      else
        phenoData <- read.phenoData(sampleNames=sampleNames,widget=FALSE)
    }
  }
  
  ##get MIAME information
  if(is.character(description)){
    description <- read.MIAME(filename=description,widget=FALSE)
  }
  else{
    if(class(description)!="MIAME"){
      if(widget){
        require(tkWidgets)
        description <- read.MIAME(widget=TRUE)
      }
      else
        description <- new("MIAME")
    }
  }
  
  ##MIAME stuff
  description@preprocessing$filenames <- filenames
  if(exists("tksn")) description@samples$description <- tksn[,2]
  description@preprocessing$affyversion <- library(help=affy)$info[[2]][[2]][2]

  ##and now we are ready to read cel files
  return(just.gcrma(filenames=filenames,
                    phenoData=phenoData,
                    description=description,
                    notes=notes,
                    compress=compress,
                    verbose=verbose,
                    normalize=normalize,
                    bgversion=bgversion,
                    affinity.info=affinity.info,
                    type=type, k=k, stretch=stretch,
                    correction=correction, rho=rho,
                    optical.correct=optical.correct,
                    fast=fast, minimum=minimum))
}


just.gcrma <- function(..., filenames=character(0),
                       phenoData=new("phenoData"),
                       description=NULL,
                       notes="", background=FALSE,
                       compress=getOption("BioC")$affy$compress.cel,
                       normalize=TRUE, bgversion=2, affinity.info=NULL,
                       type=c("fullmodel","affinities","mm","constant"),
                       k=6*fast+0.5*(1-fast), stretch=1.15*fast+1*(1-fast),
                       correction=1, rho=0.7, optical.correct=TRUE,
                       verbose=TRUE, fast=TRUE, minimum=1) {

  require(affy, quietly=TRUE)

  auxnames <- as.list(substitute(list(...)))[-1]
  filenames <- .Primitive("c")(filenames, auxnames)
  
  n <- length(filenames)
  
  ## error if no file name !
  if (n == 0)
    stop("No file name given !")
  
  pdata <- pData(phenoData)
  ##try to read sample names from phenoData. if not there use CEL filenames
  if(dim(pdata)[1]!=n){#if empty pdata filename are samplenames
    warning("Incompatible phenoData object. Created a new one.\n")
    
    samplenames <- gsub("^/?([^/]*/)*", "", unlist(filenames), extended=TRUE	)
    pdata <- data.frame(sample=1:n,row.names=samplenames)
    phenoData <- new("phenoData",pData=pdata,varLabels=list(sample="arbitrary numbering"))
  }
  else samplenames <- rownames(pdata)
  

  if (is.null(description))
  {
    description <- new("MIAME")
    description@preprocessing$filenames <- filenames
    description@preprocessing$affyversion <- library(help=affy)$info[[2]][[2]][2]
  }

  ## get information from cdf environment

  headdetails <- .Call("ReadHeader", filenames[[1]], compress)
  dim.intensity <- headdetails[[2]]
  cdfName <- headdetails[[1]]
 
  type <- match.arg(type)

  pmonly <- (type=="affinities"|type=="constant")
  needaff <- (type=="fullmodel"|type=="affinities")

  if( needaff & is.null(affinity.info)){
    if(verbose) cat("Computing affinities.")
    affinity.info <- compute.affinities(cdfName,verbose=verbose)
    if(verbose) cat("Done.\n")
    
    pm.affinities <- pm(affinity.info)
    mm.affinities <- mm(affinity.info)

    index.affinities <- which(!is.na(pm.affinities))

    ##Recover memory
    rm(affinity.info)
    gc()
    
  }

  pms <- read.probematrix(filenames=filenames, which="pm")$pm
  mms <- read.probematrix(filenames=filenames, which="mm")$mm

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
      pms[,i] <- bg.adjust.fullmodel(pms[,i],mms[,i],
                                     pm.affinities,mm.affinities,
                                     index.affinities,k=k,
                                     Q=correction*mean(pms[,i]<mms[,i]),
                                     Qmm=correction*0.5,rho=rho,fast=fast)
      pms[index.affinities,i] <- 2^(log2(pms[index.affinities,i])-
                                    fit1$coef[2]*pm.affinities+mean(fit1$coef[2]*pm.affinities))
    }
    if(type=="affinities"){
      pms[,i] <- bg.adjust.affinities(pms[,i],pm.affinities,
                                      index.affinities, k=k,
                                      Q=correction*mean(pms[,i]<mms[,i]),
                                      fast=fast)
      pms[index.affinities,i] <- 2^(log2(pms[index.affinities,i])-
                                    fit1$coef[2]*pm.affinities + 
                                    mean(fit1$coef[2]*pm.affinities))
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

 
  tmp <- new("AffyBatch",
             cdfName=cdfName,
             annotation=cleancdfname(cdfName, addcdf=FALSE))
  pmIndex <- pmindex(tmp)
  probenames <- rep(names(pmIndex), unlist(lapply(pmIndex,length)))
  pmIndex <- unlist(pmIndex)

  ngenes <- length(geneNames(tmp))
  
  ##background correction - not used, but need to pass to .Call
  
  bg.dens <- function(x){density(x,kernel="epanechnikov",n=2^14)}
  
  exprs <- .Call("rma_c_complete",pms,pms,probenames,ngenes,body(bg.dens),new.env(),normalize,background=FALSE,bgversion)

  colnames(exprs) <- filenames
  se.exprs <- array(NA, dim(exprs))
  
  annotation <- annotation(tmp)
  
  new("exprSet", exprs = exprs, se.exprs = se.exprs, phenoData = phenoData, 
      annotation = annotation, description = description, notes = notes)
}


  



