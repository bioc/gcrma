bg.parameters.ns=function (x, affinities,affinities2=NULL,affinities3=NULL)
{   
       set.seed(1989);sample1=sample(length(x),5000)
       suppressWarnings(lo1<-loess(log(x)~affinities,
                        subset=sample1,degree=1,family="symmetric",span=.2))
       bg.mu=predict(lo1,affinities)
       bg.mu[affinities>max(affinities[sample1])]=max(lo1$fitted)
       bg.mu[affinities<min(affinities[sample1])]=min(lo1$fitted)
       if (is.null(affinities2)) bg.mu2=NULL
       else {bg.mu2=predict(lo1,affinities2)
             bg.mu2[affinities2>max(affinities[sample1])]=max(lo1$fitted)
             bg.mu2[affinities2<min(affinities[sample1])]=min(lo1$fitted)
           }
       res=lo1$res[lo1$res<0];res=c(res,-res)
       bg.sigma=mad(res)

       if (is.null(affinities3))        return(list(bg.mu=bg.mu,bg.mu2=bg.mu2,bg.sigma=bg.sigma))
       
       else {bg.mu3=predict(lo1,affinities3)
             bg.mu3[affinities3>max(affinities[sample1])]=max(lo1$fitted)
             bg.mu3[affinities3<min(affinities[sample1])]=min(lo1$fitted)
             return(list(bg.mu=bg.mu,bg.mu2=bg.mu2,bg.mu3=bg.mu3,bg.sigma=bg.sigma))
           }
     }


# sg.parameters=function(pms,apm,Source){
#   apm=as.matrix(apm)
#   index.affinities <- which(!is.na(apm[,1]))
#   apm=apm[index.affinities,]
#   pms=pms[index.affinities,]
#   set.seed(1)
#   if(Source=="reference"){ #one array
#     Subset <- sample(1:length(pms),25000)
#     y <- log2(pms)[Subset]
#     Subset <- (Subset-1)%%nrow(pms)+1
#     x <- apm[Subset]}
#   else {#multiple arrays
#     set.seed(1)
#     Subset <- sample(1:length(pms),25000)
#     y <- log2(pms)[Subset]
#     x <- apm[Subset]}
#   fit1 <- lm(y~x)
#   return(fit1$coef)
# }
  

  
  
