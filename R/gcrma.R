gcrma <- function(object,affinity.info=NULL,
                  affinity.source=c("reference","local"),
                  NCprobe=NULL,
                  type=c("fullmodel","affinities","mm","constant"),
                  k=6*fast+0.5*(1-fast),stretch=1.15*fast+1*(1-fast),correction=1,
                  GSB.adjust=TRUE,rho=0.7,
                  optical.correct=TRUE,verbose=TRUE,fast=TRUE,
                  subset=NULL,normalize=TRUE,...){
  object <- bg.adjust.gcrma(object,
                            affinity.info=affinity.info,
                            affinity.source=affinity.source, NCprobe=NCprobe,
                            type=type,k=k,stretch=stretch,correction=correction,
                            GSB.adjust=GSB.adjust,rho=rho,
                            optical.correct=optical.correct,verbose=verbose,fast=fast)
  return(rma(object,subset=subset,background=FALSE,normalize=normalize,verbose=verbose))
}



bg.adjust.gcrma <- function(object,affinity.info=NULL,
                            affinity.source=c("reference","local"),
                            NCprobe=NULL,
                            type=c("fullmodel","affinities","mm","constant"),
                            k=6*fast+0.5*(1-fast),stretch=1.15*fast+1*(1-fast),correction=1,
                            GSB.adjust=TRUE,rho=0.7,
                            optical.correct=TRUE,verbose=TRUE,fast=TRUE){
  ##OPTICAL BG
  if(optical.correct)
    object <- bg.adjust.optical(object,verbose=verbose)
  type <- match.arg(type)
  pmonly <- (type=="affinities"|type=="constant")
  needaff <- (type=="fullmodel"|type=="affinities")
  ##OBTAIN AFFINITY.INFO

  Source=match.arg(affinity.source)
  if( needaff & is.null(affinity.info) & Source=="reference" )  ##use the GCRMA internal NSB data 
    affinity.info <- compute.affinities(cdfName(object),verbose=verbose)
  if( needaff & is.null(affinity.info) & Source=="local" )  ##use experimental data
    affinity.info <- compute.affinities.local(object,Array=NULL,NCprobe=NCprobe,
                                              verbose=verbose,optical.correct=FALSE)

  ##clean NCprobe
  if(!is.null(NCprobe)){
    if(length(affinity.info)==1){
      anc=as.matrix(intensity(affinity.info)[NCprobe,])
      index2=which(!is.na(anc))
      anc=anc[index2,]
      ncs=as.matrix(intensity(object)[NCprobe[index2],])}
    else {
      anc=intensity(affinity.info)[NCprobe,]
      index2=which(!is.na(anc[,1]))
      anc=anc[index2,]
      ncs=as.matrix(intensity(object)[NCprobe[index2],])
    }
  }
  else anc=ncs=NULL
  
  pmIndex=unlist(indexProbes(object,"pm"))
  mmIndex=unlist(indexProbes(object,"mm"))
  
  if(!is.null(affinity.info)){
    if(length(affinity.info)==1){ #one array affinity.info, NCprobe=MM
      pm(object) <- gcrma.engine(pms=pm(object),
                                 mms=mm(object),ncs,
                                 pm.affinities=pm(affinity.info),
                                 mm.affinities=mm(affinity.info),anc,
                                 type=type,k=k,
                                 stretch=stretch,
                                 correction=correction,GSB.adjust=GSB.adjust,rho=rho,
                                 verbose=verbose,fast=fast)
    }
    else if (length(affinity.info)==length(object)){#multi-array affinity.info
      pm(object) <- gcrma.engine2(object,pmIndex,mmIndex,
                                  NCprobe=NCprobe,
                                  affinity.info,
                                  type=type,k=k,
                                  stretch=stretch,
                                  correction=correction,GSB.adjust=GSB.adjust,rho=rho,
                                  verbose=verbose,fast=fast)}
  }
  else{#local affinity.info, not provided.
    
  }
  return(object)
}
