bg.adjust.optical <- function(abatch,minimum=1,verbose=TRUE){
  Index <- unlist(indexProbes(abatch,"both"))

  if(verbose) cat("Adjusting for optical effect")
  for(i in 1:length(abatch)){
    if(verbose) cat(".")
    exprs(abatch)[Index,i] <- exprs(abatch)[Index,i] -
      min(exprs(abatch)[Index,i],na.rm=TRUE) + minimum
  }
  if(verbose) cat("Done.\n")
  
  abatch
}

bg.adjust.mm <- function(pms,mms,k=6*fast+0.25*(1-fast),fast=TRUE){
  mu <- log(mms)
  Index <- which(pms<mms)
  sigma <- sqrt(mean((log(pms)[Index]-mu[Index])^2))

  bhat <- exp(mu+0.5*sigma^2)
  var.y <- exp(2*mu+sigma^2)*(exp(sigma^2)-1)

  if(fast) return(gcrma.bg.transformation.fast(pms,bhat,var.y,k=k))
  else return(gcrma.bg.transformation(pms,bhat,var.y,k=k))
}

bg.adjust.constant <- function(x,k=6*fast+0.25*(1-fast),Q=0.25,fast=TRUE){
  mu <- log(quantile(x,Q))
  sigma <- left.sigma(log(x),mu)
 
  if(fast){
     bhat <- exp(mu+1/2*sigma^2)
     var.y <- rep(exp(2*mu+sigma^2)*(exp(sigma^2)-1),length(x))
     return(gcrma.bg.transformation.fast(x,bhat,var.y,k=k))
   }
  else return(gcrma.bg.transformation(x,mu,sigma,k=k))
}

bg.adjust.affinities <- function(x,affinities,index=seq(along=x),
                                 k=6*fast+0.25*(1-fast),Q=0.25,fast=TRUE){
  parameters <- bg.parameters.ns(x[index],affinities,Q=Q)
  mu <- vector("numeric",length(x))
  sigma <- vector("numeric",length(x))
  mu[index] <- parameters$bg.mu
  sigma[index] <- parameters$bg.sigma
  
  ##fill in the pms for which we dont have affinities
  if(length(index)<length(x)){
    mu[-index] <- median(mu[index])
    sigma[-index] <- median(sigma[index])
  }
  
  if(fast){
      bhat <- exp(mu + 1/2*sigma^2)
      var.y <- exp(2*mu+sigma^2)*(exp(sigma^2)-1)
      return(gcrma.bg.transformation.fast(x,bhat,var.y,k=k))
    }
  else return(gcrma.bg.transformation(x,mu,sigma,k=k))
}

bg.adjust.fullmodel <- function(pms,mms,pm.affinities,mm.affinities,
                                index.affinities,k=6*fast+0.25*(1-fast),
                                Q=0.25,Qmm=0.5,rho=0.7,fast=TRUE){
  
  parameters <- bg.parameters.ns(pms[index.affinities],pm.affinities,Q=Q)
  mu.pm <- vector("numeric",length(pms))
  sigma <- vector("numeric",length(pms))
  mu.pm[index.affinities] <- parameters$bg.mu
  sigma[index.affinities] <- parameters$bg.sigma

  parameters <- bg.parameters.ns(mms[index.affinities],mm.affinities,Q=Qmm)
  mu.mm <-  vector("numeric",length(pms))
  mu.mm[index.affinities] <- parameters$bg.mu
  
  ##fill in the pms for which we dont have affinities
  if(length(index.affinities)<length(pms)){
    mu.pm[-index.affinities] <- median(mu.pm[index.affinities])
    mu.mm[-index.affinities] <- median(mu.mm[index.affinities])
    sigma[-index.affinities] <- median(sigma[index.affinities])
  } 
  
  ##mu.mm <- mu.mm-mean(mu.mm)+mean(mu.pm) ##correction in case not the same
  ##this is the unbiased
  
  if(fast){
    bhat <- exp(mu.pm + rho*(log(mms)-mu.mm) + 1/2*(1 - rho^2)*sigma^2)
    var.y=exp(2*mu.pm+sigma^2)*(exp(sigma^2)-exp(sigma^2*rho^2))
    return(gcrma.bg.transformation.fast(pms,bhat,var.y,k=k))
  }
  else return(gcrma.bg.transformation(pms,mu.pm + rho*(log(mms)-mu.mm),sqrt(1 - rho^2)*sigma,k=k))
}
