##default, one array affinity.info
gcrma.engine <- function(pms,mms,ncs=NULL,
                         pm.affinities=NULL,mm.affinities=NULL,anc=NULL,
                         type=c("fullmodel","affinities","mm","constant"),
                         k=6*fast+0.5*(1-fast),
                         stretch=1.15*fast+1*(1-fast),correction=1,GSB.adjust=TRUE,rho=0.7,
                         verbose=TRUE,fast=FALSE){
  
  if(!is.null(pm.affinities)){
    index.affinities <- which(!is.na(pm.affinities))
    pm.affinities <- pm.affinities[index.affinities]
    if(!is.null(mm.affinities)){
      mm.affinities <- mm.affinities[index.affinities]
    }
  }
  
###get GSB.adjust parameters
  if(GSB.adjust &(type=="fullmodel" | type=="affinities")){
    set.seed(1)
      Subset <- sample(1:length(as.matrix(pms)[index.affinities,]),25000)
    y <- log2(pms)[index.affinities,][Subset]
    Subset <- (Subset-1)%%nrow(as.matrix(pms)[index.affinities,])+1
    x <- pm.affinities[Subset]
    fit1 <- lm(y~x)$coef
  }
  
  
  if(verbose) cat("Adjusting for non-specific binding")
  for(i in 1:ncol(pms)){
    if(verbose) cat(".")
    if(type=="fullmodel") {
      pms[,i] <- bg.adjust.fullmodel(pms[,i],mms[,i],ncs=ncs[,i],
                                     apm=pm.affinities,amm=mm.affinities,anc=anc,
                                     index.affinities,k=k,
                                     rho=rho,fast=fast)
      if(GSB.adjust) 
        pms[,i] <- GSB.adj(Yin=pms[,i],subset=index.affinities,aff=pm.affinities,fit1=fit1,k=k)
    }
    if(type=="affinities"){
      if(is.null(ncs))
        pms[,i] <-  bg.adjust.affinities(pms[,i],ncs=mms[,i],apm=pm.affinities,anc=mm.affinities,
                                         index.affinities,k=k,fast=fast)
      else
        pms[,i] <- bg.adjust.affinities(pms[,i],ncs=ncs[,i],apm=pm.affinities,anc=anc,
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


GSB.adj <- function(Yin,subset,aff,fit1,k=k){ #subset are the index of Yin with the affinities aff
    y0=Yin[subset]
    y0.adj=k+ 2^(-fit1[2]*(aff-mean(aff)))*(y0-k)
    Yin[subset]=y0.adj
    Yin
     }
