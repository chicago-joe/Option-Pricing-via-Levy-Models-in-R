

K=22


Fx = seq(0,1,by=(1)/K)
x = seq(-0.477,0,by=(0.477+0)/K)

N=(Fx[K+1]-Fx[1])/K

Fx
inp=-0.1

bf<-function(inp){
  out=0
  for (variable in x) {
    if(variable>inp){
      return(out)
    }
    out=out+1
  }
}

bff<-function(inp){
  out=0
  for (variable in Fx) {
    if(variable>inp){
      return(out)
    }
    out=out+1
  }
}

Fff <- function(inp) {
  out=0
  if(inp<x[1]){
    out=0
    return(out)
  }else if(inp>=x[K+1]){
    out=1
    return(out)
  }else{
    #xk=bs(inp)
    xk_1=bf(inp)
    out=Fx[xk_1]+(Fx[xk_1+1]-Fx[xk_1])/N*(inp-x[xk_1])
  }
  return(out)
}
Fff(-0.47)
U=(runif(1,0,1))
Fiff<-function(U){
  xk_1=bff(U)
  if(xk_1<(K-1)){
    return(x[xk_1+1]+(x[xk_1+2]-x[xk_1+1])/(Fx[xk_1+2]-Fx[xk_1+1])*(U-Fx[xk_1+1]))
  }else{
    return(0)
  }
}
#=============================
Xt=Fiff(runif(1,0,1))
stock=100*exp(0.05*1)*Fiff(0.41)
put=max(0,100-stock)
put
call=max(0,stock-100)
call
#==============================

bff(0.5)
s0=100
s0*exp(Xt)
T=0.5
simulations= 64*10^3; 
r= 0.05;
strike=100;
out <- rep(0,simulations)
P <- rep(0,simulations)
C <- rep(0,simulations)

for (j in 1:simulations) {
  Xt=Fiff(runif(1,0,1))
  
  out[j]=s0*exp(Xt)
  P[j]=exp(-r*T)*max(0,strike-out[j])
  C[j]=exp(-r*T)*max(0,out[j]-strike)
}
sum(out)/simulations
sum(P)/simulations
sum(C)/simulations


V <- rep(0,simulations)


for (j in 1:simulations) {
  V[j]=100*exp(-r*T)*max(0,100/100-exp(Fiff(runif(1,0,1))))
  
}
sum(V)/simulations



