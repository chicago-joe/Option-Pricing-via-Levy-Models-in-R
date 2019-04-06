###########################################################################
#
#
###########################################################################
seed <- 1
set.seed(seed)

alpha= 15; beta= -5; delta= 0.5; 
r= 0.05; q= 0.02; s0= 100; T= 6/12; 
K = s0; x = 0.1;
simulations= 64; 

N = simulations*(10^3)
dt = T/N
time = seq(0,T,by=dt)

mu = r-q + delta*sqrt((alpha*alpha)-((beta+1)*(beta+1))) 
                        - sqrt((alpha*alpha)-(beta*beta));
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

### Step 1: Compute zeta ###
G1 <- BoxMuller()                               # generate standard normal random variable G1
Z <- (G1^2)/gamma                              # compute Z
zeta = (1/gamma)*((dt) + (0.5*Z) - sqrt((dt*Z) + (Z^2)/4));   


### Step 2: Generate uniform random variable U on (0,1) ###
## Uniform RNG ##
Uniform <- function() {
  #set.seed(seed);
    U = runif(N,0,1)
    zt <- rep(0,N)
    for (i in 1:N) {
    if (U[i] < dt/(dt+gamma*zeta[i])) {
      zt[i] = zeta[i] } else  {
      zt[i] = ((delta^2)*(time^2)) / (gamma^2*zeta[i]) }
    }
    return(zt) } 

zt = Uniform() 

  
### Step 3: Generate StdNorm variable G2
G2 <- rnorm(n=N,0,1)

### Now simulate Xt = mu*t + beta*zt + sqrt(zt)*G2
### then compute St = s0*e^Xt
St <- rep(0,N)
Xt <- rep(0,N)
for (i in 1:N) {
  if (i==1)  {
    Xt[i] = mu*time[i+1] + beta*zt[i] + sqrt(zt[i])*G2[i]
    St[i] = s0*exp(Xt[i]) } else {
      Xt[i] = mu*time[i+1] + beta*zt[i] + sqrt(zt[i])*G2[i]
      St[i] = St[i-1]*exp(Xt[i]) }
}

# ERROR:  stock price rises towards infinity
head(St)
tail(St)

mu=mu;
sigma=beta;
P0=100;
T=1/12;
nt=50;
n=2^(8)
dt=T/n;
t=seq(0,T,by=dt)
X=matrix(rep(0,length(t)*nt), nrow=nt)
for (i in 1:N) {
  
}