
bg.adjust.optical <- function(abatch,minimum=1,verbose=TRUE){
  Index <- unlist(indexProbes(abatch,"both"))
  Index <- Index[!is.na(Index)]
  if(verbose) cat("Adjusting for optical effect")
  for(i in 1:length(abatch)){
    if(verbose) cat(".")
    exprs(abatch)[Index,i] <- exprs(abatch)[Index,i] -
      min(exprs(abatch)[Index,i],na.rm=TRUE) + minimum
  }
  if(verbose) cat("Done.\n")
  
  abatch
}


##########################################################################################
bg.adjust.fullmodel<- function(pms,mms,ncs=NULL,apm,amm,anc=NULL,index.affinities,k=6*fast+0.25*(1-fast),rho=.7,fast=FALSE){
  if(is.null(ncs)){
    parameters <- bg.parameters.ns(mms[index.affinities],amm,apm)
    mu.pm <- vector("numeric",length(pms))
    mu.mm <-  vector("numeric",length(pms))
    mu.pm[index.affinities] <- parameters$bg.mu2
    mu.mm[index.affinities] <- parameters$bg.mu
    sigma<- parameters$bg.sigma
  }
  
  else{
    parameters <- bg.parameters.ns(ncs,anc,apm,amm)
    mu.pm <- vector("numeric",length(pms))
    mu.mm <-  vector("numeric",length(pms))
    mu.pm[index.affinities] <- parameters$bg.mu2
    mu.mm[index.affinities] <- parameters$bg.mu3
    sigma<- parameters$bg.sigma
  }
  
  ##fill in the pms for which we dont have affinities
  if(length(index.affinities)<length(pms)){
    mu.pm[-index.affinities] <- median(mu.pm[index.affinities])
    mu.mm[-index.affinities] <- median(mu.mm[index.affinities])
  } 
  
  if(fast){
    ##this is the unbiased  
    bhat <- exp(mu.pm + rho*(log(mms)-mu.mm) + 1/2*(1 - rho^2)*sigma^2)
    var.y=exp(2*mu.pm+sigma^2)*(exp(sigma^2)-exp(sigma^2*rho^2))
    return(gcrma.bg.transformation.fast(pms,bhat,var.y,k=k))
  }
  else return(gcrma.bg.transformation(pms,mu.pm + rho*(log(mms)-mu.mm),sqrt(1 - rho^2)*sigma,k=k))
}

##########################################################################################
bg.adjust.affinities<- function(pms,ncs,apm,anc,index.affinities,k=6*fast+0.25*(1-fast),fast=FALSE,nomm=FALSE){
  if(!nomm) parameters <- bg.parameters.ns(ncs[index.affinities],anc,apm)
  else parameters <- bg.parameters.ns(ncs,anc,apm)
  mu.pm <- vector("numeric",length(pms))
  mu.pm[index.affinities] <- parameters$bg.mu2
  sigma<- parameters$bg.sigma
  ##fill in the pms for which we dont have affinities
  if(length(index.affinities)<length(pms)){
    mu.pm[-index.affinities] <- median(mu.pm[index.affinities])
  } 
  
  if(fast){
    ##this is the unbiased  
    bhat <- exp(mu.pm +  1/2*sigma^2)
    var.y=exp(2*mu.pm+sigma^2)*(exp(sigma^2)-1)
    return(gcrma.bg.transformation.fast(pms,bhat,var.y,k=k))
  }
  else return(gcrma.bg.transformation(pms,mu.pm,sigma,k=k))
} 
##########################################################################################
bg.adjust.mm <- function(pms,mms,k=6*fast+0.25*(1-fast),fast=FALSE){
  mu <- log(mms)
  Index <- which(pms<mms)
  sigma <- sqrt(mean((log(pms)[Index]-mu[Index])^2))

  bhat <- exp(mu+0.5*sigma^2)
  var.y <- exp(2*mu+sigma^2)*(exp(sigma^2)-1)

  if(fast) return(gcrma.bg.transformation.fast(pms,bhat,var.y,k=k))
  else return(gcrma.bg.transformation(pms,bhat,var.y,k=k))
}

##########################################################################################
bg.adjust.constant <- function(x,k=6*fast+0.25*(1-fast),Q=0.25,fast=FALSE){
  mu <- log(quantile(x,Q))
  sigma <- left.sigma(log(x),mu)
 
  if(fast){
     bhat <- exp(mu+1/2*sigma^2)
     var.y <- rep(exp(2*mu+sigma^2)*(exp(sigma^2)-1),length(x))
     return(gcrma.bg.transformation.fast(x,bhat,var.y,k=k))
   }
  else return(gcrma.bg.transformation(x,rep(mu,length(x)),sigma,k=k))
}
left.sigma <- function(x,mu) sqrt(mean((x[x<mu]-mu)^2))
