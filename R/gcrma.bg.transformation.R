#####################################################################################
##########  LOGNORMAL + LOG UNIFORM /as old code but in C                  ##########
#####################################################################################
gcrma.bg.transformation <- function (pms, mu, tau, k, a = 1, ...) 
{
#  dyn.load("~/Year5/GCrelease1.6/up2/callpnorm.so")
  sub=which(pms>exp(mu-4*tau))
  y1=.C("Rposty", as.double(pms[sub]),as.double(mu[sub]),as.double(tau),as.integer(length(sub)),as.double(k),
    ans=as.double(rep(1,length(sub))),PACKAGE="gcrma")$ans
  pms[sub]=exp(y1)
  pms[-sub]=k
  pms[is.na(pms)]=k
  pms
  }
    
