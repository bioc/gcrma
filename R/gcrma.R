require(affy)
require(MASS)
bg.correct.gcrma <- function(object,...)
{
  mms=mm(object)
  pms=pm(object)
  pm(object) <- apply(array(c(pms,mms),dim=c(dim(pms),2)),2,bg.adjust.gcrma,...)
  return(object)
  }
bg.adjust.gcrma <- function(Data,gcgroup,estimate=c("eb","mle"),rho=0.8,step=60,lower.bound=1,baseline=.25,triple.goal,...){
  estimate=match.arg(estimate)
  pm=Data[,1];mm=Data[,2]
  pars.bg <-bg.parameters.gcrma(log(mm),gcgroup=gcgroup)
  pars.sg <- sg.parameters.gcrma(pm,pars.bg,gcgroup=gcgroup)
  base=exp(log(2^16)/step)#defult~=1.2
  mylog=function(x) log(x,base)
    ##cat("background correction method:",estimate)
  if (estimate=="eb") {
    K0=max(0,round(mylog(lower.bound))+1)
    for ( k in  which(c(lapply(gcgroup, length),recursive=TRUE)>0)) {
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
        f=g*log(base^(K0:(K-1))*(base+1)/2);f2=g*log(base^(K0:(K-1))*(base+1)/2)^2
        f0=g0*log(lower.bound/2+base^K0/2);f02=g0*log(lower.bound/2+base^K0/2)^2
        fn=gn*log(base^K/2+p/2);fn2=gn*log(base^K/2+p/2)^2
        c((sum(f)+f0+fn)/(sum(g)+g0+gn),(sum(f2)+f02+fn2)/(sum(g)+g0+gn))
      }
       tmp=apply(cbind(pm[gcgroup[[k]]],mus,Ks),1,posty)
      pm[gcgroup[[k]]]=tmp[1,];mm[gcgroup[[k]]]=tmp[2,]-tmp[1,]^2
    }
    if (triple.goal)

pm=mean(pm,na.rm=TRUE)+sqrt(1+mean(mm,na.rm=TRUE)/var(pm,na.rm=TRUE))*(pm-mean(pm,na.rm=TRUE))
    pm=exp(pm)
  }

  if (estimate=="mle") { 
    for ( k in which(c(lapply(gcgroup, length),recursive=TRUE)>0)) {
      pm[gcgroup[[k]]]=pm[gcgroup[[k]]]-exp(log(mm[gcgroup[[k]]])*rho+pars.bg[1,k]*(1-rho)-(1-rho^2)*pars.bg[2,k]^2)
    }
    pm[pm<baseline]=baseline
  }
    pm[is.na(pm)]=lower.bound #with extreme mm>>>pm it's possible to get NaN
   pm
}


############################################################
bg.parameters.gcrma<-function(mm,gcgroup,n.pts=2^10,adjust=1){
   max.density <- function(x,n.pts,bw.c) {
   bw=bw.nrd0(x)*adjust
    aux <- density(x, kernel = "epanechnikov", n = n.pts,bw=bw, 
            na.rm = TRUE)
        median(aux$x[aux$y>max(aux$y)*.95])
    }
  pars.bg <- sapply(1:25,function(k){
    x=mm[gcgroup[[k]]]
    bg.mu <- max.density(x,n.pts,bw.c)
    bg.data <- x[x < bg.mu]
    bg.sd <- sqrt(sum((bg.data-bg.mu)^2)/(length(bg.data) - 1)) * sqrt(2)
        c(bg.mu,bg.sd)})
  ##last group is probe-without-seq-info
  cbind(pars.bg,apply(pars.bg,1,mean))
}
############################################################
        
sg.parameters.gcrma <- function(pm,pars.bg,gcgroup){
   pars.sg <- sapply(1:25,function(k){
     x=log(round(pm[gcgroup[[k]]],-1))
     y=log(table(x))
     x=as.numeric(names(y))
     subset=x>quantile(x[x>pars.bg[1,k]],.05)&y>log(2)
     if (sum(subset)>25)
        return(-lm(y~x,subset=subset)$coef[2])
      else( return(3))
   })
   pars.sg=apply(matrix(pars.sg,5,5),2,function(x) PAV(x)$y)
   pars.sg <-c(pars.sg,mean(pars.sg))
   pars.sg[pars.sg<1]=1;pars.sg[pars.sg>3]=3
   pars.sg
 }
  average<-function(y, wt = rep(1, length(y)))
{
# compute a weighted average of a vector, y
        if(any(is.na(wt))) stop("NA's not allowed for wt")
        if(any(wt < 0))
                stop("wt must be a vector of NON-NEGATIVE weights")
        if(length(wt) != length(y)) stop(
                "y and wt must be vectors of the same length")
# if any observations have Infinite weight, return the simple
# (unweighted) average of only those observations (giving no
# weight to observations with finite weight)
        if(any(wt == Inf)) {
                wt[wt < Inf] <- 0
                wt[wt == Inf] <- 1
}
# if all weights are zero, return the simple (unweighted)
# average of y
        if(sum(wt) == 0)
                wt <- rep(1, length(wt))
        return(sum((y * wt)/sum(wt)))
}

PAV <- function(y, wt = rep(1,length(y)))
{
# This is a modification of Derick's PAV program
#
# (Weighted) Pool-Adjacent-Violators (PAV) algorithm
# for non-parametric monotonic (decreasing) regression of y on x
        n <- length(y)
        if(n != length(wt))
                stop("y, and wt must be vectors of equal length")
yhat <- y       # initialize while loop
        j <- count <- 1
        k <- 2
        support <- vector("numeric", n)
        support[count] <- j
        while(k <= n) {
                while(yhat[j] < yhat[k]) {
                        yhat[j:k] <- average(y[j:k], wt[j:k])
                        if(yhat[support[count]] < yhat[k]) {
                                j <- support[count]
                                if(count > 1)
                                        count <- count - 1
                        }
                        else {
                                k <- ifelse(k == n, k, k + 1)
                        }
                }

 count <- count + 1
                support[count] <- j
                j <- k
                k <- k + 1
        }
        return(y = yhat, wt)
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
getGroupInfo<- function(object){
  allprobes<-probeNames(object)
  pms=pm(object)
  CDF=getCdfInfo(object)
  loc.pm=rep(0,length(allprobes));i=1
  for( x in unique(allprobes)) {
    tmp=get(x,CDF)[,1]
    loc.pm[i:(i+length(tmp)-1)]=tmp
    i=i+length(tmp)} #faster than multiget

  matchaffyName=paste(cleancdfname(object@cdfName,addcdf=FALSE),"probe",sep="")
 if (identical(.find.package(matchaffyName, quiet=TRUE),character(0)))
     stop(paste("Probe sequence information not available. Please download and install the",matchaffyName,"library from:\n http://www.bioconductor.org/data/metaData.html\n See the vignette of matchprobes and gcrma for more details"))
  require(matchaffyName,character.only = TRUE)
  seqData <- get(matchaffyName)
  seq.pm=rep(NA,length(seqData$x))
  loc.seq=xy2i(seqData$x,seqData$y)
##check seqData
  Index2=which(loc.pm%in%loc.seq) #probes with seq
  seq.pm[order(loc.pm[Index2])]=seqData$seq[order(loc.seq)]
  mm.seq<- complementSeq(seq.pm,start=13,stop=13)
  ATCGmm=basecontent(mm.seq)
  CGadj=ATCGmm[,3:4]
  CGadj[CGadj<=4]=4
  CGadj[CGadj>=8]=8
    group2d=as.list(NULL)
  for (C in 4:8) {
    for(G in 4:8) {group2d=c(group2d,list(Index2[which(CGadj[,1]==C & CGadj[,2]==G)]))}}
  c(group2d,list((1:length(allprobes))[-Index2]))
}


gcrma <- function (object,estimate="eb",summary.method = "medianpolish",normalize = TRUE,normalize.method = "quantiles", triple.goal=TRUE,rho=.8,step=60,lower.bound=1,baseline=.25,...) {
  old.bgcorrect.methods <- bgcorrect.methods
  on.exit(assign("bgcorrect.methods",old.bgcorrect.methods,env=.GlobalEnv))
  assign("bgcorrect.methods",c(bgcorrect.methods,"gcrma"), env=.GlobalEnv)
  old.express.summary.stat.methods <- express.summary.stat.methods
  on.exit(assign("express.summary.stat.methods",old.express.summary.stat.methods,env=.GlobalEnv))
  assign("express.summary.stat.methods",c(express.summary.stat.methods,"rlm"), env=.GlobalEnv)
  
  gcgroup <- getGroupInfo(object)
  res <- expresso(object,bgcorrect.method = "gcrma",bgcorrect.param =list(gcgroup,estimate,rho,step,lower.bound,baseline,triple.goal),pmcorrect.method = "pmonly",normalize=normalize,normalize.method =normalize.method ,summary.method = summary.method,...)
  return(res)
}
