bg.parameters.ns=function (x, affinities,affinities2=NULL)
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

       return(list(bg.mu=bg.mu,bg.mu2=bg.mu2,bg.sigma=bg.sigma))

}


