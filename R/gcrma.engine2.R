##multi-array affinity.info
gcrma.engine2 <- function(object,pmIndex=NULL,mmIndex=NULL,
                          NCprobe=NULL,
                          affinity.info,
                          type=c("fullmodel","affinities","mm","constant"),
                          k=6*fast+0.5*(1-fast),
                          stretch=1.15*fast+1*(1-fast),correction=1,GSB.adjust=TRUE,rho=0.7,
                          verbose=TRUE,fast=TRUE){

  if(is.null(pmIndex)) pmIndex=unlist(indexProbes(object, "pm"))
  if(is.null(mmIndex)) mmIndex=unlist(indexProbes(object, "mm"))
  pm.affinities=intensity(affinity.info)[pmIndex,]
  mm.affinities=intensity(affinity.info)[mmIndex,]
  pms=intensity(object)[pmIndex,]
  mms=intensity(object)[mmIndex,]
  if(!is.null(pm.affinities)){
    index.affinities <- which(!is.na(pm.affinities[,1]))
    pm.affinities <- pm.affinities[index.affinities,]
    if(!is.null(mm.affinities)){
      mm.affinities <- mm.affinities[index.affinities,]
    }
  }
  
###get GSB.adjust parameters #multiple arrays
  if(GSB.adjust &(type=="fullmodel" | type=="affinities")){
    set.seed(1)
    Subset <- sample(1:length(pms[index.affinities,]),25000)
    y <- log2(pms)[index.affinities,][Subset]
    x <- pm.affinities[Subset]
    fit1 <- lm(y~x)$coef
  }

  if(verbose) cat("Adjusting for non-specific binding")
  for(i in 1:ncol(pms)){
    if(verbose) cat(".")
    ##fullmodel
    if(type=="fullmodel") {
      if(is.null(NCprobe)) 
        pms[,i] <- bg.adjust.fullmodel(pms[,i],mms[,i],ncs=NULL,
                                       apm=pm.affinities[,i],amm=mm.affinities[,i],
                                       anc=NULL,
                                       index.affinities,k=k,
                                       rho=rho,fast=fast)
    
      else
        pms[,i] <- bg.adjust.fullmodel(pms[,i],mms[,i],ncs=intensity(object)[NCprobe,i],
                                       apm=pm.affinities[,i],amm=mm.affinities[,i],
                                       anc=intensity(affinity.info)[NCprobe,i],
                                       index.affinities,k=k,
                                       rho=rho,fast=fast)
    
      if(GSB.adjust)
        pms[,i] <- GSB.adj(Yin=pms[,i],subset=index.affinities,aff=pm.affinities[,i],fit1=fit1,k=k)

    }  
    ##affinities
    if(type=="affinities"){
      if(is.null(NCprobe))
        pms[,i] <-  bg.adjust.affinities(pms[,i],ncs=mms[,i],
                                         apm=pm.affinities[,i],anc=mm.affinities[,i],
                                         index.affinities,k=k,fast=fast)
      else
        pms[,i] <- bg.adjust.affinities(pms[,i],ncs=intensity(object)[NCprobe,i],
                                        apm=pm.affinities[,i],anc=intensity(affinity.info)[,i],
                                        index.affinities,k=k,fast=fast,nomm=TRUE)
      if(GSB.adjust)
        pms[,i] <- GSB.adj(Yin=pms[,i],subset=index.affinities,aff=pm.affinities,fit1=fit1,k=k)
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
