require(affy)
require(MASS)
bg.correct.gcrma <- function(object,...)
{
  mms=mm(object)
  pms=pm(object)
  pm(object) <- apply(array(c(pms,mms),dim=c(dim(pms),2)),2,bg.adjust.gcrma,...)
  return(object)
  }

bg.adjust.gcrma <- function(Data,gcgroup,estimate=c("eb","mle"),rho=0.8,step=60,lower.bound=1,baseline=.25,...){
  estimate=match.arg(estimate)
  pm=Data[,1];mm=Data[,2]
  pars.bg <-bg.parameters.gcrma(log(mm),gcgroup=gcgroup)
  pars.sg <- sg.parameters.gcrma(pm,pars.bg,gcgroup=gcgroup)
  base=exp(log(2^16)/step)#defult~=1.2
  mylog=function(x) log(x,base)
    ##cat("background correction method:",estimate)
  if (estimate=="eb") {
    K0=max(0,round(mylog(lower.bound))+1)
    for ( k in seq(along=gcgroup)) {
      a=pars.sg[k]
      tau=pars.bg[2,k]*sqrt(1-rho^2)
      mus=(1-rho)*pars.bg[1,k]+rho*log(mm[gcgroup[[k]]])
      Ks=floor(mylog(pm[gcgroup[[k]]]))
      posty <- function(xxx){
        p=xxx[1];mu=xxx[2];K=xxx[3]	
        pnorms=pnorm(log(p-base^(K0:K)),mu,tau)
        diff.pnorms=pnorms[-length(pnorms)]-pnorms[-1]
        g=(base^a+1)/base^(a*(K0:(K-1)+1))/2*diff.pnorms
        g0=(1/lower.bound^a+1/base^(a*K0))/2*(pnorm(log(p-lower.bound),mu,tau)-pnorm(log(p-base^K0),mu,tau))
        gn=(base^a+1)/base^(a*(K+1))*pnorm(log(p-base^K),mu,tau) 
        f=g*log(base^(K0:(K-1))*(base+1)/2)
        f0=g0*log(lower.bound/2+base^K0/2)
        fn=gn*log(base^K/2+p/2)
        (sum(f)+f0+fn)/(sum(g)+g0+gn)
      }
      pm[gcgroup[[k]]]=apply(cbind(pm[gcgroup[[k]]],mus,Ks),1,posty)
    }
    pm=exp(pm)
  }

  if (estimate=="mle") { 
    for ( k in seq(along=gcgroup)) {
      pm[gcgroup[[k]]]=pm[gcgroup[[k]]]-exp(log(mm[gcgroup[[k]]])*rho+pars.bg[1,k]*(1-rho)-(1-rho^2)*pars.bg[2,k]^2)
    }
    pm[pm<baseline]=baseline
  }
    pm[is.na(pm)]=lower.bound #with extreme mm>>>pm it's possible to get NaN
   pm
}




############################################################
bg.parameters.gcrma <- function(mm,gcgroup,n.pts=2^10){
  max.density <- function(x, n.pts) {
        aux <- density(x, kernel = "epanechnikov", n = n.pts, 
            na.rm = TRUE)
        aux$x[order(-aux$y)[1]]
    }
  pars.bg <- sapply(1:25,function(k){
    x=mm[gcgroup[[k]]]
    bg.mu <- max.density(x,n.pts)
    bg.data <- x[x < bg.mu]
    bg.sd <- sqrt(sum((bg.data-bg.mu)^2)/(length(bg.data) - 1)) * sqrt(2)
        c(bg.mu,bg.sd)})
  ##last group is probe-without-seq-info
  cbind(pars.bg,apply(pars.bg,1,mean))
}

sg.parameters.gcrma <- function(pm,pars.bg,gcgroup){
   pars.sg <- sapply(1:25,function(k){
     x=log(round(pm[gcgroup[[k]]],-1))
     y=log(table(x))
     x=as.numeric(names(y))
     if (sum(x>pars.bg[1,k]+3*pars.bg[2,k]&y>log(2))>40) 
        return(-lm(y~x,subset=x>pars.bg[1,k]+3*pars.bg[2,k]&y>log(2))$coef[2])
      else( return(3.5))
   })
   pars.sg <-c(pars.sg,mean(pars.sg))
   pars.sg[pars.sg<1]=1
   pars.sg
 }
        

##########################################################################################
##########################################################################################
generateExprVal.method.rlm <- function(probes,...){
  gcrma.rlm(probes,...)
}

gcrma.rlm<-function(x,maxit=50,k=.5)#x is #of probes by # of arrays
  { X.model=data.frame(X1=gl(ncol(x),nrow(x)),X2=gl(nrow(x),1,length(x)))
    X.model=model.matrix(~X1+X2,X.model,contrasts = list(X2="contr.sum"))[,-1]  
    rlm0=summary(rlm(as.vector(log2(x))~X.model,maxit=maxit,k=k))
    exprs=c(rlm0$coef[1,1],rlm0$coef[2:ncol(x),1]+rlm0$coef[1,1])
    se.exprs=sqrt(sum(rlm0$cov.unscaled[1:2,1:2]))*rlm0$stddev
    list(exprs=exprs,se.exprs=rep(se.exprs,ncol(x)))
  }

#bgcorrect.methods <- c(bgcorrect.methods,"gcrma")
#express.summary.stat.methods <- c(express.summary.stat.methods,"rlm")

getGroupInfo <- function(object){
  allprobes<-probeNames(object)
  matchaffyName=paste(cleancdfname(object@cdfName,addcdf=FALSE),"probe",sep="")
  require(matchaffyName,character.only = TRUE)
  seqData <- get(matchaffyName)
  mm.seq<- complementSeq(seqData$seq,start=13,stop=13)
  ATCGmm=basecontent(mm.seq)
  CGadj=ATCGmm[,3:4]
  CGadj[CGadj<=4]=4
  CGadj[CGadj>=8]=8
  Index2=which(allprobes%in%seqData$Probe.Set.Name) #probes with seq 
  group2d=as.list(NULL)
  for (C in 4:8) {
    for(G in 4:8) {group2d=c(group2d,list(Index2[which(CGadj[,1]==C & CGadj[,2]==G)]))}}
  c(group2d,list(which(!allprobes%in%seqData$Probe.Set.Name)))
}

gcrma <- function (object,estimate="eb",summary.method = "rlm",normalize = TRUE,normalize.method = "quantiles", rho=.8,step=60,lower.bound=1,baseline=.25,...) {
  old.bgcorrect.methods <- bgcorrect.methods
  on.exit(assign("bgcorrect.methods",old.bgcorrect.methods,env=.GlobalEnv))
  assign("bgcorrect.methods",c(bgcorrect.methods,"gcrma"), env=.GlobalEnv)
  old.express.summary.stat.methods <- express.summary.stat.methods
  on.exit(assign("express.summary.stat.methods",old.express.summary.stat.methods,env=.GlobalEnv))
  assign("express.summary.stat.methods",c(express.summary.stat.methods,"rlm"), env=.GlobalEnv)
  
  gcgroup <- getGroupInfo(object)
  res <- expresso(object,bgcorrect.method = "gcrma",bgcorrect.param =list(gcgroup,estimate,rho,step,lower.bound,baseline),pmcorrect.method = "pmonly",normalize=normalize,normalize.method =normalize.method ,summary.method = summary.method,...)
  return(res)
}
