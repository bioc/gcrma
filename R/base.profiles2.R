complementSeq <- function(seq,start=1,stop=0){
  if (is.numeric(stop) && length(stop) == 1L && !is.na(stop) && stop == 0)
    stop <- -1
  seq <- DNAStringSet(seq)
  subseq(seq, start=start, end=stop) <-
    complement(subseq(seq, start=start, end=stop))
  as.character(seq)
}

# compute affinity.spline.coefs from MMs
base.profiles.mm <- function(object,verbose=TRUE){
  
  cleancdf <- cleancdfname(cdfName(object),addcdf=FALSE)
  ##  object <- bg.adjust.optical(object) #should already been done
  cdfpackagename <- paste(cleancdf,"cdf",sep="")
  #myxy2i <- get("xy2i",paste("package:",cdfpackagename,sep=""))
  probepackagename <- paste(cleancdf,"probe",sep="")
  getCDF(cdfpackagename)
  getProbePackage(probepackagename)
  p <- get(probepackagename)
  seqs=p$seq#the PM sequences
  seqs=complementSeq(seqs, start=13, stop=13)#the MM sequences

  pmIndex <-  unlist(indexProbes(object,"pm"))
  mmIndex <-  unlist(indexProbes(object,"mm"))
  subIndex <- match(xy2indices(p$x,p$y, cdf=cdfpackagename),pmIndex)
  bgy <- as.matrix(intensity(object)[mmIndex[subIndex],])
  affinity.spline.coefs <- base.profiles(bgy,seqs)
  affinity.spline.coefs
}
# compute affinity.spline.coefs from a vector of negative control probes

base.profiles.nc <- function(object, NCprobe,verbose=TRUE){
  cleancdf  <- cleancdfname(cdfName(object),addcdf=FALSE)
  ##object <- bg.adjust.optical(object)
  cdfpackagename <- paste(cleancdf,"cdf",sep="")
  #myxy2i <- get("xy2i",paste("package:",cdfpackagename,sep=""))
  probepackagename <- paste(cleancdf,"probe",sep="")
  getCDF(cdfpackagename)
  getProbePackage(probepackagename)
  p <- get(probepackagename)
  PMseq=p$seq#the PM sequences
  MMseq=complementSeq(PMseq, start=13, stop=13)#the MM sequences

  pmIndex <-  unlist(indexProbes(object,"pm"))
  mmIndex <-  unlist(indexProbes(object,"mm"))
  subIndex1 <- match(NCprobe,c(xy2indices(p$x,p$y, cdf=cdfpackagename),
                               xy2indices(p$x,p$y+1, cdf=cdfpackagename)))
  seqs <-c(PMseq,MMseq)[subIndex1[!is.na(subIndex1)]]

  bgy <- as.matrix(intensity(object)[NCprobe,])
  if(length(seqs)<length(NCprobe)){
    cat("\nNote: some of your negative control probes do not have sequence information\n")
    subIndex2 <- match(c(xy2indices(p$x,p$y, cdf=cdfpackagename),
                         xy2indices(p$x,p$y+1, cdf=cdfpackagename)),NCprobe)
    subIndex2 <- subIndex2[!is.na(subIndex2)]
    bgy <- bgy[!is.na(subIndex1),]}

  affinity.spline.coefs <- base.profiles(bgy,seqs)
  affinity.spline.coefs
}

##seqs: probe sequences of those used to estimate base profiles
##intensities: observed intensity on these probes, a vector or a matrix
base.profiles <- function(intensities,seqs){
  prlen <- unique(nchar(seqs))
  stopifnot(length(prlen)==1)
  mS <- .Call("gcrma_getSeq2",paste(seqs,collapse=""),length(seqs),PACKAGE="gcrma")

  ##Create a basis matrix
  B <- ns(1:25,df=5)
  ##Create the prediction matrix in the model
  X <- cbind(mS[,1:25]%*%B,mS[,26:50]%*%B,mS[,51:75]%*%B)
  
  Y=as.matrix(intensities)
  affinity.spline.coefs <- apply(Y,2,function(y)
                                 lm(log2(y) ~ X)$coef[-1])
  rownames(affinity.spline.coefs) <- c(rep("A",5),rep("C",5),rep("G",5))
  affinity.spline.coefs
}

plotBaseProfiles <- function(affinity.spline.coefs,prlen=25){
  if(!is.vector(affinity.spline.coefs)) stop("affinity.spline.coefs must be a vector")
  P <- length(affinity.spline.coefs)/3 ##number of parameters per letter
  B <- ns(1:25,df=5)
  coeff <- matrix(0, nrow=prlen, ncol=4)
  colnames(coeff) <- c("A", "C", "G", "T")
  for (i in 0:2)
    coeff[, i+1] <- B%*%affinity.spline.coefs[i*P+(1:P)]
  coeff <- coeff-rowSums(coeff[,1:3])/4
  matplot(coeff,  type="b", col=1:4, lwd=2, pch=colnames(coeff))
}

##this is almost the same as compute.affinities
##except that affinity.spline.coefs is not loaded with data()
##but computed using the data provided by the user
compute.affinities.local <- function(object,Array=1,NCprobe=NULL,verbose=TRUE,
                                     affinity.spline.coefs=NULL,optical.correct=TRUE){
  if(optical.correct) object <- bg.adjust.optical(object)
  if(is.null(affinity.spline.coefs))
    affinity.spline.coefs <- compute.affinity.coef(object,Array,NCprobe,verbose=TRUE)
  if(verbose) cat("Computing affinities")
  affinity.basis.matrix <- ns(1:25,df=nrow(as.matrix(affinity.spline.coefs))/3)

  cdfname <- cdfName(object)
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
  
  APM=AMM=matrix(NA,length(p$sequence),length(object))
  for(K in 1:ncol(APM)){
    if(verbose) cat(".")
    apm <- vector("numeric",length(p$sequence))
    amm <- vector("numeric",length(p$sequence))
    
    for(i in seq(along=apm)) {
      charMtrx <- .Call("gcrma_getSeq", p$sequence[i],
                        PACKAGE="gcrma")
      A <- cbind(charMtrx[1,] %*% affinity.basis.matrix,
                 charMtrx[2,] %*% affinity.basis.matrix,
                 charMtrx[3,] %*% affinity.basis.matrix)
      
      apm[i] <- A %*% affinity.spline.coefs[,K]
      
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
    APM[,K]=apm;AMM[,K]=amm
  }
  ##put it in an affybatch
  #tmp <- get("xy2i",paste("package:",cdfpackagename,sep=""))
  affinity.info <- new("AffyBatch",cdfName=cdfname)
  pmIndex <-  unlist(indexProbes(affinity.info,"pm"))
  mmIndex <-  unlist(indexProbes(affinity.info,"mm"))
  subIndex <- match(xy2indices(p$x,p$y, cdf=cdfpackagename),pmIndex)
  tmp.exprs=matrix(NA,nrow=max(cbind(pmIndex,mmIndex)),ncol=ncol(APM))
  tmp.exprs[pmIndex[subIndex],]=APM
  if(!is.null(amm)){ tmp.exprs[mmIndex[subIndex],]=AMM }
  exprs(affinity.info)=tmp.exprs
  if(verbose) cat("Done.\n")
  return(affinity.info)
}

compute.affinity.coef <- function(object,Array,NCprobe,verbose=TRUE){
  if(verbose) cat("Computing base-position profiles for probe affinities")
  if(is.null(NCprobe)){ #use MMs
    if(is.null(Array))  #"use all arrays"
      affinity.spline.coefs <- base.profiles.mm(object,verbose)
    else if(is.numeric(Array) & length(Array)==1) #use one Array
      affinity.spline.coefs <- base.profiles.mm(object[,Array],verbose)
    else stop("\nParameter Array should be either NULL(use all arrays)
 or a number indexing the array to be used\n")
  }

  else if(is.numeric(NCprobe)){ #"use NCs from each array"
    if(is.null(Array))  #"use all arrays"
      affinity.spline.coefs <- base.profiles.nc(object,NCprobe,verbose)
    else if(is.numeric(Array) & length(Array)==1) #use one Array
      affinity.spline.coefs <- base.profiles.nc(object[,Array],NCprobe,verbose)
    else stop("\nParameter Array should be either NULL(use all arrays)
 or a number indexing the array to be used\n")
  }

  else
    stop("\nNCprobe should be either NULL (use MMs)
or a vector of probe indexes for negative control probes\n")
  if(verbose) cat("Done.\n")
  affinity.spline.coefs
}

