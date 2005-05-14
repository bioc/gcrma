gcrma <- function(object,affinity.info=NULL,
                  type=c("fullmodel","affinities","mm","constant"),
                  k=6*fast+0.5*(1-fast),stretch=1.15*fast+1*(1-fast),
                  correction=1,rho=0.7,
                  optical.correct=TRUE,verbose=TRUE,fast=TRUE){
  
  type <- match.arg(type)
  
  
  pmonly <- (type=="affinities"|type=="constant")
  needaff <- (type=="fullmodel"|type=="affinities")
  if( needaff & is.null(affinity.info)){
    if(verbose) cat("Computing affinities")
    affinity.info <- compute.affinities(cdfName(object),verbose=verbose)
    if(verbose) cat("Done.\n")
  }
  
  if(optical.correct)
    object <- bg.adjust.optical(object,verbose=verbose)
  
  pm(object) <- gcrma.engine(pms=pm(object),
                             mms=mm(object),
                             pm.affinities=pm(affinity.info),
                             mm.affinities=mm(affinity.info),
                             type=type,k=k,
                             stretch=stretch,
                             correction=correction,rho=rho,
                             verbose=verbose,fast=fast)
  
  return(rma(object,background=FALSE,verbose=verbose))
}

##for now we need the mms for everything
gcrma.engine <- function(pms,mms,pm.affinities=NULL,mm.affinities=NULL,
                         type=c("fullmodel","affinities","mm","constant"),
                         k=6*fast+0.25*(1-fast),
                         stretch=1.15*fast+1*(1-fast),correction=1,rho=0.7,
                         verbose=TRUE,fast=TRUE){
  
  type <- match.arg(type)
  
  if(!is.null(pm.affinities)){
    index.affinities <- which(!is.na(pm.affinities))
    pm.affinities <- pm.affinities[index.affinities]
    if(!is.null(mm.affinities)){
      mm.affinities <- mm.affinities[index.affinities]
    }
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
                                     rho=rho,fast=fast)
      pms[index.affinities,i] <- 2^(log2(pms[index.affinities,i])-
                                    fit1$coef[2]*pm.affinities+mean(fit1$coef[2]*pm.affinities))
    }
    if(type=="affinities"){
      pms[,i] <- bg.adjust.affinities(pms[,i],mms[,i],
                                      pm.affinities,mm.affinities,
                                      index.affinities, k=k,
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

  return(pms)
}


