compute.affinities <- function(cdfname,verbose=TRUE){
  if(verbose) cat("Computing affinities")   
  data(affinity.spline.coefs, package = "gcrma") ###needs to change to data(something)
  affinity.basis.matrix <- ns(1:25,df=length(affinity.spline.coefs)/3)
  
  cleancdf <- cleancdfname(cdfname,addcdf=FALSE)
  cdfpackagename <- paste(cleancdf,"cdf",sep="")
  probepackagename <- paste(cleancdf,"probe",sep="")
  
  getCDF(cdfpackagename)
  getProbePackage(probepackagename)
  p <- get(probepackagename)
  
  p <- check.probes(p, cdfname)
  
  prlen <- unique(nchar(p$sequence))
  stopifnot(length(prlen)==1)
  
  A13 <- sum(affinity.basis.matrix[13,]*affinity.spline.coefs[1:5])
  T13 <- 0
  C13 <- sum(affinity.basis.matrix[13,]*affinity.spline.coefs[6:10])
  G13 <- sum(affinity.basis.matrix[13,]*affinity.spline.coefs[11:15])
  
  if(verbose) cat(".")
  
  apm <- vector("numeric",length(p$sequence))
  amm <- vector("numeric",length(p$sequence))
  
  for(i in seq(along=apm)) {
    charMtrx <- .Call("gcrma_getSeq", p$sequence[i],
                      PACKAGE="gcrma")
    A <- cbind(charMtrx[1,] %*% affinity.basis.matrix,
               charMtrx[2,] %*% affinity.basis.matrix,
               charMtrx[3,] %*% affinity.basis.matrix)
    
    apm[i] <- A %*% affinity.spline.coefs
    
    if (charMtrx[1,13] == 1) {
      amm[i] <- apm[i] + T13 - A13
    }
    else {
      if (charMtrx[4,13] == 1) {
        amm[i] <- apm[i] + A13 - T13
      }
      else{
        if (charMtrx[3,13]) {
          amm[i] <- apm[i] + C13 - G13
        }
        else {
          amm[i] <- apm[i] + G13 - C13
        }
      }
    }
  }
  
  ##put it in an affybatch
  #tmp <- get("xy2i",paste("package:",cdfpackagename,sep=""))
  affinity.info <- new("AffyBatch",cdfName=cdfname)
  pmIndex <-  unlist(indexProbes(affinity.info,"pm"))
  mmIndex <-  unlist(indexProbes(affinity.info,"mm"))
  subIndex <- match(xy2indices(p$x,p$y, cdf=cdfpackagename),pmIndex)
  tmp.exprs=matrix(NA,nrow=max(cbind(pmIndex,mmIndex)),ncol=1)
  tmp.exprs[pmIndex[subIndex]]=apm
  if(!is.null(amm)){ tmp.exprs[mmIndex[subIndex]]=amm }
  exprs(affinity.info)=tmp.exprs
  if(verbose) cat("Done.\n")
  return(affinity.info)
}

check.probes <- function(probepackage, cdfname){
  cdfnames <- names(pmindex(new("AffyBatch", cdfName=cdfname)))
  ppnames <- as.character(probepackage$Probe.Set.Name)
  
  if (sum(!(ppnames %in% cdfnames)) != 0){
    Index <- ppnames %in% cdfnames
    probepackage <- probepackage[Index,]
  }
  return(probepackage)
}


