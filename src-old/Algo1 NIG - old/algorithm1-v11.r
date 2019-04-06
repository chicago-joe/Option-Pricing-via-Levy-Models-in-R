###########################################################################
#
#
###########################################################################


alpha= 15; beta= -5; delta= 0.5; 
r= 0.05; q= 0.02; s0= 100; T= 6/12; 
K = s0;
simulations= 64*10^3; 

N = T/delta;

mu = r-q + delta*(sqrt((alpha*alpha)-((beta+1)*(beta+1))) 
                        - sqrt((alpha*alpha)-(beta*beta)));
gamma = sqrt((alpha*alpha)-(beta*beta))


## Box-Muller RNG ##
BoxMuller <- function() {
  #set.seed(seed)
  n = 2*N; 
  z = numeric(n);
  
  u1 = runif(n/2,0,1); 
  u2 = runif(n/2,0,1)
  
  z1 = sqrt(-2*log(u1)) * cos(2*pi*u2)  # half of normal variates
  z2 = sqrt(-2*log(u1)) * sin(2*pi*u2)  # other half
  z[seq(1, n, by=2)] = z1   # interleave
  z[seq(2, n, by=2)] = z2   # two halves
  
  return(z1)
}
out <- rep(0,simulations)
P <- rep(0,simulations)
C <- rep(0,simulations)

for (j in 1:simulations) {
### Step 1: Compute zeta ###
G1 <- BoxMuller()                               # generate standard normal random variable G1
Z <- (G1^2)/gamma                              # compute Z
zeta <- rep(0,N)
for (i in 1:N) {
  zeta = (1/gamma)*((delta*i) + (0.5*Z[i]) - sqrt((delta*i*Z[i]) + (Z[i]^2)/4));   
}

### Step 2: Generate uniform random variable U on (0,1) ###
## Uniform RNG ##
Uniform <- function() {
  #set.seed(seed);
    U = runif(N,0,1)
    zt <- rep(0,N)
    for (i in 1:N) {
    if (U[i] < delta*i/(delta*i+gamma*zeta[i])) {
      zt[i] = zeta[i] 
      } else  {
      zt[i] = ((delta^2)*(i^2)) / (gamma^2*zeta[i]) }
    }
    return(zt) } 

zt = Uniform() 

### Step 3: Generate StdNorm variable G2
G2 <- BoxMuller()

### Now simulate Xt = mu*t + beta*zt + sqrt(zt)*G2
### then compute St = s0*e^Xt
St <- rep(0,N)
Xt <- rep(0,N)



  
  for (i in 1:N) {
    
    if (i==1)  {
      Xt[i] = mu*i + beta*zt[i] + sqrt(zt[i])*G2[i] 
      St[i] = s0*exp(Xt[i])
      } else {
        Xt[i] = mu*i + beta*zt[i] + sqrt(zt[i])*G2[i]
        St[i] = s0*exp(Xt[i]) 
      }
  }
  out[j]=St[N]
  P[j]=exp(-r*T)*max(0,K-out[j])
  C[j]=exp(-r*T)*max(0,out[j]-K)
}
sum(out)/simulations
sum(P)/simulations
sum(C)/simulations
# ERROR:  stock price rises towards infinity
St
plot(St)

Xt
plot(Xt)
tail(St)

plot(V)
V