###########################################################################
#
#
###########################################################################
#seed <- 1
#set.seed(seed)

alpha= 15; beta= -5; delta= 0.5; 
r= 0.05; q= 0.02; s0= 100; T= 6/12; 
K = s0; x0= 0.1;
simulations= 64; 

N = simulations*(10)
dt = T/N
t = seq(0,T,by=dt)

x = c(x0,(alpha*dt)+1*sqrt(dt)*rnorm(n=N,0,1))

mu = r-q + delta*sqrt((alpha*alpha)-((beta+1)*(beta+1))) - sqrt((alpha*alpha)-(beta*beta));
gamma = sqrt((alpha*alpha)-(beta*beta))


for (i in 1:N) {
  seed<-1; 
  set.seed(seed);
  
### Box-Muller RNG ###
  BoxMuller <- function() {
    #set.seed(seed);
    n = 2*N; z = numeric(n);
    
    u1 = runif(n/2,0,1); 
    u2 = runif(n/2,0,1)
    z1 = sqrt(-2*log(u1)) * cos(2*pi*u2)    # half of normal variates
    z2 = sqrt(-2*log(u1)) * sin(2*pi*u2)    # other half
    
    z[seq(1, n, by=2)] = z1     # interleave
    z[seq(2, n, by=2)] = z2     # two halves
    return(z1)  }
  
### Step 1: Compute zeta
  G1 <- BoxMuller()                     # generate standard normal random variable G1
  Z <- (G1*G1)/gamma                    # compute Z
  zeta = (1/gamma)*((dt) + (0.5*Z) 
                    - sqrt((dt*Z) + (Z^2)/4));   

### Step 2: Generate uniform random variable U on (0,1)
  Uniform <- function() {             ## Uniform RNG ##
    #set.seed(seed); n=N;
    U = runif(n=N,0,1)
    if (U < (dt/(dt+gamma*zeta))) {
      z = zeta
      return(z) }
    else  {
      z = ((delta^2)*(t^2))/((gamma^2)*zeta)
      return(z) } }
  
  U <- Uniform() 
  
### Step 3: Generate G2
  G2 <- rnorm(n=N,0,1) 

  Xt = mu*t + beta*z + sqrt(z)*G2
  
  St = s0*exp(Xt)   ## Compute St = s0*e^(Xt) ##
  St                # ERROR:  stock price seems to slowly rise towards infinity
}