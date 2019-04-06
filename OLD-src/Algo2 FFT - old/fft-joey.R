# Inverse Fourier Transform to compute CDF
# using a table of probabilities 
#
# Created by Joseph Loss on 01/15/2019
# 
#

library(data.table)
library(stats)
library(stats4)

alpha = 15
beta = -5
d.min = beta - alpha
d.max = beta + alpha
x0 = -0.477
#x[22] = 0
K = 22
bias = .01*10^-2
N = 64*10^3
s0 = 100
k = 100
r = 0.05
T = 0.5
delta = 0.5
d = 0.02

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
    out=Fx[xk_1]+(Fx[xk_1+1]-Fx[xk_1])/N(inp-x[xk_1])
  }
  return(out)
}
Fff(-0.47)

Fiff<-function(U){
  xk_1=bff(U)
  return(x[xk_1+1]+(x[xk_1+2]-x[xk_1+1])/(Fx[xk_1+2]-Fx[xk_1+1])(U-Fx[xk_1+1]))
}
Xt=Fiff(runif(1,0,1))
bff(0.5)
s0=100
s0exp(Xt)

simulations= 6410^3; 
r= 0.05;
strike=100;
out <- rep(0,simulations)
P <- rep(0,simulations)
C <- rep(0,simulations)

for (j in 1:simulations) {
  Xt=Fiff(runif(1,0,1))
  
  out[j]=s0exp(Xt)
  P[j]=exp(-rT)max(0,strike-out[j])
  C[j]=exp(-rT)*max(0,out[j]-strike)
}
sum(out)/simulations
sum(P)/simulations
sum(C)/simulations



DT <- data.table(x = rep(0,K),
                 Fx = rep(0,K),
                 key=c("x","Fx"))

DT$x[K]=0
DT$x[1]=-0.477
x[22] = 0

K=22


x = seq(-0.477,0,by=(0+0.477)/K)
N=(x[K]-x[1])/K

x
inp=0.3

bf<-function(inp){
  out=1
  for (variable in x) {
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
  }else if(inp>=x[K]){
    out=1
    return(out)
  }else{
    #xk=bs(inp)
    xk_1=bf(inp)
    out=x[xk_1]+(x[xk_1+1]-x[xk_1])/N*(inp-xk_1)
  }
  return(out)
}
Fff(-0.1)
print(bf(-0.1))


#for (i in (1:K) {#}

# DT <- data.table(xK = sample(letters, N, TRUE),
#                  FxK = sample(letters, N, TRUE), key = c("xK","FxK"))
# print(object.size(DT), units = "Mb")
# key(DT)
# t1 <- system.time(ans1<-DT[.("g", 877L)])
# t1
# head(ans1)
# dim(ans1)

## Test if == to each other
identical(ans1$val, ans2$val)


# Setup --------------------------------------------------------------------
# log2.K <- log(K, base=exp(2))      # Binary Search max iterations = O(log2 k)


# Binary Search (Feng) ----------------------------------------------------
U = runif(N,0,1)
L = 0
R = K
for (i in 1:K) {
DT <- data.table(x = sample(letters, N, TRUE),
             y = sample(letters, N, TRUE), key = c("xi","F.xi"))
}

while (L < (R-1)) {
  m = floor((L+R)/2)
  if (F.xm < U) {
    L = m} else {
      R = m}  }

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


# Binary Search Using library(data.tables) ---- EXAMPLE --------------------------------
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

value=100*exp(-0.05*1)*fft(0.07)


value=100*exp(-0.05*0.5)*max(0,100/(100-exp(inverse fft Xt))

# -------------------------------------------------------------------------


