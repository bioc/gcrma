gcrma.bg.transformation <- function(pms,mu,tau,k,lower.bound=.5,a=1,...){
  step=60
  base=exp(log(2^16)/step)#defult~=1.2
  mylog=function(x) log(x,base)
  K0=max(0,round(mylog(lower.bound))+1)
  Ks=floor(mylog(pms))
  posty <- function(xxx){
    p=xxx[1];mu=xxx[2];tau=xxx[3];K=xxx[4]	
    pnorms=pnorm(log(p-base^(K0:K)),mu,tau)
    diff.pnorms=pnorms[-length(pnorms)]-pnorms[-1]
    g=(base^a+1)/base^(a*(K0:(K-1)+1))/2*diff.pnorms
    g0=(1/lower.bound^a+1/base^(a*K0))/2*(pnorm(log(p-lower.bound),mu,tau)-pnorm(log(p-base^K0),mu,tau))
    gn=(base^a+1)/base^(a*(K+1))*pnorm(log(p-base^K),mu,tau) 
    f=g*log(base^(K0:(K-1))*(base+1)/2);f2=g*log(base^(K0:(K-1))*(base+1)/2)^2
    f0=g0*log(lower.bound/2+base^K0/2);f02=g0*log(lower.bound/2+base^K0/2)^2
    fn=gn*log(base^K/2+p/2);fn2=gn*log(base^K/2+p/2)^2
    (sum(f)+f0+fn)/(sum(g)+g0+gn)
  }
  pm=apply(cbind(pms,mu,tau,Ks),1,posty)
  pm=exp(pm)
  pm[is.na(pm)]=lower.bound #with extreme mm>>>pm it's possible to get NaN
  pm
}

