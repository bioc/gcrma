compute.affinities <- function(cdfname,pmonly=FALSE,
                               verbose=TRUE){
  require(splines,quietly = TRUE)
  require(matchprobes,quietly = TRUE)
  data(affinity.spline.coefs) ###needs to change to data(something)
  affinity.basis.matrix <- ns(1:25,df=length(affinity.spline.coefs)/3)

  cleancdf <- cleancdfname(cdfname,addcdf=FALSE)
  cdfpackagename <- paste(cleancdf,"cdf",sep="")
  probepackagename <- paste(cleancdf,"probe",sep="")
  
  library(cdfpackagename,character.only=TRUE)
  library(probepackagename,character.only=TRUE)
  p <- get(probepackagename)
  
  p <- check.probes(p, cdfname) #missing line
  
  prlen <- unique(nchar(p$sequence))
  stopifnot(length(prlen)==1)
  
  if(verbose) cat(".")
  mat <- matrix(unlist(strsplit(p$sequence, "")), byrow=TRUE, ncol=prlen)
  stopifnot(all(mat %in% c("A", "C", "G", "T")))
  S <- cbind(mat=="A", mat=="C", mat=="G")

  if(verbose) cat(".")
  S <- array(as.numeric(S), dim=dim(S))
  colnames(S) <- paste ("S", 1:(3*prlen), sep="")
  Xpm <- cbind(S[,1:25]%*%affinity.basis.matrix,S[,26:50]%*%affinity.basis.matrix,S[,51:75]%*%affinity.basis.matrix)
  apm <- Xpm%*%affinity.spline.coefs
  
  if(!pmonly){
    cat(".")
    for(i in 1:nrow(S))
      if(S[i,38]==1 | S[i,63]==1) S[i,c(38,63)]=1-S[i,c(38,63)] else S[i,13]=1-S[i,13]

    Xmm <- cbind(S[,1:25]%*%affinity.basis.matrix,S[,26:50]%*%affinity.basis.matrix,S[,51:75]%*%affinity.basis.matrix)
    amm <- Xmm%*%affinity.spline.coefs
  }
  else amm <- NULL
  
  tmp <- get("xy2i",paste("package:",cdfpackagename,sep=""))
  pmIndex <-  unlist(indexProbes(new("AffyBatch",cdfName=cdfname),"pm"))
  Index <- match(tmp(p$x,p$y),pmIndex)
  
  return(list(pm=apm,mm=amm,index=Index))
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

  
