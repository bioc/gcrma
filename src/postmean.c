#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <stdlib.h>
#include <stdio.h>
/*#include "nmath.h"*/
/*#include "dpq.h"*/
#define min(x1,x2) ((x1) > (x2))? (x2):(x1)
#define max(x1,x2) ((x1) > (x2))? (x1):(x2)

void posty(double p, double mu,double tau,double lower_bound, double *ans)
{ 
  /*  double   lower_bound=.5;*/
  
  int a=1; 
  int  step = 60,K0,K,i;
  double base;
  double G1=0,F1=0;
  double tmp;
  double *pnorms, diff_pnorms;
  /* double *g; */

  base = exp(log(pow(2,16))/step);
  /*  mylog = function(x) log(x, base);*/
  K0 = max(0, floor(log(lower_bound)/log(base)) + 1.0);
  K = floor(log(p)/log(base));
  
  pnorms = (double *) R_alloc(K+1, sizeof(double));
  /* g = (double *) R_alloc(K+2, sizeof(double)); */
  pnorms[0]=pnorm(log(p - pow(base,K0)), mu, tau,1,0);
  G1=(1/pow(lower_bound,a) + 1/pow(base,(a * K0)))/2*(pnorm(log(p -lower_bound), mu, tau,1,0) - pnorm(log(p - pow(base,K0)), mu, tau,1,0));
  F1= G1 * log(lower_bound/2 + pow(base,K0)/2);

  for (i=1;i<=K-K0;i++){
    pnorms[i] = pnorm(log(p - pow(base,(K0+i))), mu, tau,1,0);
    diff_pnorms = pnorms[i-1] - pnorms[i];
    tmp = (pow(base,a) + 1)/pow(base,(a * (K0+i)))/2 * diff_pnorms;
    G1+=tmp;
    F1+=tmp * log(pow(base,(K0+i-1)) * (base + 1)/2);
  }

  tmp= (pow(base,a) + 1)/pow(base,(a * (K + 1))) * pnorm(log(p - pow(base,K)), mu, tau,1,0);
  G1+=tmp;
  F1+= tmp * log(pow(base,K)/2 + p/2);
  *ans=F1/G1;
  
}

void Rposty1(double *p, double *mu, double *tau, double *k,double *ans){
  posty(*p,*mu,*tau,*k,ans);
}


void Rposty(double *p, double *mu, double *tau,int *G, double *k, double *ans){
  int i;
  for(i=0;i<*G;i++)    posty(p[i],mu[i],tau[0],*k,ans+i);
}
