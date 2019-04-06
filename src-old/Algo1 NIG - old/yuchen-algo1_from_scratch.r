#############################################

# Pricing European Plain-Vanilla options via Levy Processes

# created by Joseph Loss
# Copyright 2019

#############################################
library(BH)
library(BSDA)
library(crayon)
library(DistributionUtils)
library(dplyr)
library(Ecfun)
library(fBasics)
library(fda)
library(flexsurv)
library(jrvFinance)
library(fAssets)
library(FNN)
library(Formula)
library(FinancialMath)
library(GeneralizedHyperbolic)
library(GGally)
library(highr)
library(Hmisc)
library(hms)
library(latticeExtra)
library(lintr)
library(MASS)
library(mclust)
library(moments)
library(muhaz)
library(mvtnorm)
library(numDeriv)
library(pillar)
library(plotly)
library(PortfolioAnalytics)
library(prettyunits)
library(progress)
library(quantmod)
library(Rcpp)
library(RcppArmadillo)
library(readr)
library(sde)
library(stabledist)
library(stats)
library(stats4)
library(styler)
library(survival)
library(tbl2xts)
library(tibble)
library(tidyr)
library(fGarch)
library(boot)
library(e1071)
library(bootstrap)
library(urca)
library(quadprog)
library(tseries)
library(xts)
library(Metrics)
library(AER)
library(caret)
library(fExoticOptions)
library(fOptions)
library(TTR)
library(kernlab)
library(qrmtools)
library(quantreg)
library(derivmkts)
library(fOptions)
library(QRM)
library(ragtop)



##############################################
## set up inputs
alpha=15;
beta=-5;
delta=0.5;
r=0.05;
q=0.02;
S0=100;
K=S0;
T0=0.5;
N=100;
## compute mu
mu=r-q+delta*sqrt((alpha*alpha)-((beta+1)*(beta+1)))-sqrt((alpha*alpha)-(beta*beta));
u=0.3
#BSiN to get stander normal 
G22=function(u){
  a0=2.5662823884;
  a1=-18.61500062529;
  a2=41.39119773534;
  a3=-2544106049637;
  b0=-8.47351093090;
  b1=23.08336743743;
  b2=-21.06224101826;
  b3=3.13082909833;
  c0=0.3374754822726147;
  c1=0.9761690190917186;
  c2=0.1607979714918209;
  c3=0.0276438810333863;
  c4=0.0038405729373609;
  c5=0.0003951896511919;
  c6=0.0000321767881768;
  c7=0.0000002888167364;
  c8=0.0000003960315187;
  
  y=u-0.5;
  if(abs(y)<0.42){
    r=y*y;
    x=y*(((a3*r+a2)*r+a1)*r+a0)/((((b3*r+b2)*r+b1)*r+b0)*r+1.0);
  }
  else{
    r=u;
    if(y>0){
      r=1-u;
    }
      r=log(-log(r));
      x=c0+r*(c1+r*(c2+r*(c3+r*c4+r*(c5+r*(c6+r*(c7+r*c8))))));
      if(y<0){
        x=-x;
      }
  }
  return (x);
}
## NIG function
NIG=function(alpha,beta,delta,mu, T, N) {
  ## compute Gamma
  gamma=sqrt(alpha*alpha-beta*beta)
  ## three uniform
  U2=runif(N,0,1);
  U3=runif(N,0,1);
  U=runif(N,0,1);
  
  ## init zeta, zt, X (not real)
  zeta=runif(N,0,1);
  zt=runif(N,0,1);
  X=runif(N,0,1);
  for (i in 1:N) {
    ## two stander normal
    G1=G22(U2[i]);
    G2=G22(U3[i]);
    Z=(G1*G1)/gamma;
    zeta[i]=(1.0/gamma)*(delta*i+0.5*Z-sqrt(delta*i*Z+Z*Z/4));
    #I[i] = G2(a*h, b)
    #X[i+1] =  X[i]+mu/N+beta*delta*delta*I[i] + 
     # delta*sqrt(I[i])*rnorm(1)
    if(U[i]<(delta *i/(delta*i+zeta[i]))){
      zt[i]=zeta[i];
    }else{
      zt[i]=delta*delta*i*i/(gamma*gamma*zeta[i])
    }
    if(zt[i]<0){
      zt[i]=-zt[i]
    }
    ##compute X
  X[i]=mu*i+beta*zt[i]+sqrt(zt[i])*G2;
  }
  return(X)}


# run NIG using input ouput Xt
X=NIG(alpha,beta,delta,mu,T0,N)
X
# compute St
S=S0*exp(sum(X))
S=S0*exp((X))
X
S
#check stander normal
G22(runif(1,0,1))
