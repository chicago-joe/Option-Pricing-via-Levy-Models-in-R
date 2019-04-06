###########################################################################
#
#
###########################################################################
#seed <- 1
#set.seed(seed)

alpha= 15; beta= -5.0; delta= 0.5; 
r= 0.05; q= 0.02; s0= 100; t= 1/12; 
K = 100; x = 0.01;

simulations = N = 64*10^3; 
dt = x/simulations;
#dt = (T*((1/12))/simulations)
time = seq(0,t,by=dt)
time = time[-1]


# old code:
# N = simulations * T
# dt = T/N or dt = 1/200000 works too


mu = r-q + delta*sqrt((alpha^2)-((beta+1)*(beta+1))) - sqrt((alpha^2) - (beta^2));
gamma = sqrt((alpha^2)-(beta^2))


## Box-Muller RNG ##
BoxMuller <- function() {
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

### Step 1: Compute zeta ###
G1 <- BoxMuller()                               # generate standard normal random variable G1
Z <- (G1^2)/gamma                              # compute Z
zeta = (1/gamma)*((dt) + (0.5*Z) - sqrt((dt*Z) + (Z^2)/4));   


### Step 2: Generate uniform random variable U on (0,1) ###
## Uniform RNG ##
Uniform <- function() {
  U = runif(N,0,1)
  zt <- rep(0,N)
    for (i in 1:N) {
    if (U[i] < (delta*time[i])/(delta*time[i]+gamma*zeta[i])) {
      zt[i] = zeta[i] 
      } else  {
      zt[i] = ((delta^2)*(time[i]^2)) / (gamma^2*zeta[i]) }
    }
    return(zt) 
} 

zt = Uniform() 

  
### Step 3: Generate StdNorm variable G2
G2 <- BoxMuller()


### Now simulate Xt = mu*t + beta*zt + sqrt(zt)*G2
St <- rep(0,N)
Xt <- rep(0,N)
for (i in 1:N) {
    Xt[i] = mu*time[i] + beta*zt[i] + sqrt(zt[i])*G2[i] 
    St[i] = s0*exp(Xt[i])       # compute St = s0*e^Xt
}


# St rises towards infinity?
# St
# Xt


# Option Value computation:
V = exp(-r*t)*max(0,K-St[N])
V.m = mean(sapply(St,FUN=function(x){exp(-r*t)*max(0,K-x)}))

V # value = $7.4971 with Put Strike = s0 = 100, T = 1/12
V.m # value = $2.5718 with Put Strike = 105, s0 = 100, T = 1/12 

