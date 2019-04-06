

K=22


Fx = seq(0,1,by=(1)/K)
x = seq(-0.477,0,by=(0+0.477)/K)

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

Fiff<-function(U){
  xk_1=bff(U)
  if(xk_1>K-1){
  return(x[xk_1+1]+(x[xk_1+2]-x[xk_1+1])/(Fx[xk_1+2]-Fx[xk_1+1])*(U-Fx[xk_1+1]))
  }else{
    return(0)
  }
}
#=============================
Xt=Fiff(runif(1,0,1))
stock=100*exp(0.05*1)*Fff(-0.11)
put=max(0,100-stock)
put
call=max(0,stock-100)
call
#==============================

bff(0.5)
s0=100
s0*exp(Xt)

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

# -------------------------------------------------------------------------


# random functions --------------------------------------------------------
# while (_K > 0) {
#   N = (xk - x0) / _K
#   
#   }
# 
# while ((_K-1) >= k >= 0) {
#   xk = x0 + k*N
# 
#   }  
# simulate U on (0,1)
# U<-runif(100,0,1)

# -------------------------------------------------------------------------


# Monte Carlo Riemann Sums ------------------------------------------------
# m = 100000
# a = 0; b = 3/2; w = (b - a)/m;
# x = seq(a + w/2, b-w/2, length = m)
# h = x^2; rect.areas = w*h
# 
# sum(rect.areas)      # Riemann
# u = runif(m, a, b)
# h = u^2; y = (b - a)*h
# mean(y)              # Monte Carlo
# 2*sd(y)/sqrt(m)      # MC margin of error
# sum(rect.areas)      # Riemann
# mean(y)              # Monte Carlo
# 
# 2*sd(y)/sqrt(m)      # MC margin of error


# -------------------------------------------------------------------------


# Binary Search Using library(data.tables) --------------------------------
set.seed(2L)
N = 2e7L

DT = data.table(x = sample(letters, N, TRUE), y = sample(1000L, N, TRUE), 
                val = runif(N), key = c("x", "y"))
print(object.size(DT), units = "Mb")
key(DT)

## (1) Usual way of subsetting - vector scan approach
t1 <- system.time(ans1 <- DT[x == "g" & y == 877L])
t1
head(ans1)
dim(ans1)

## (2) Subsetting using keys
t2 < -system.time(ans2<-DT[.("g", 877L)])
t2
head(ans2)
dim(ans2)

## Test if == to each other
identical(ans1$val, ans2$val)


# -------------------------------------------------------------------------


