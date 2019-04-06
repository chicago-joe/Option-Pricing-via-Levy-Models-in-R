###########################################################################
# Pricing European Put Options Through Simulating Levy Processes 
# Pseudocode implemented from Professor Liming Feng's Paper (pg. 17)
#
# Created by Joseph Loss on 12/17/2018
# Additional Recognition: Yuchen Duan (UIUC MSFE) and Daniel Liberman (UIC Finance)
#
# Algorithm 1: Simulating a Normal Inverse Gaussian process through a Brownian subordination
#
###########################################################################
# seed = 1234               # for debugging
# set.seed(seed)            

# Input Parameters --------------------------------------------------------
alpha = 15;                       # parameters taken from Feng Paper pg. 22-23
beta = -5; 
delta = 0.5; 
r = 0.05; 
q = 0.02; 
s0 = 100; 
K = s0; 
T = 0.5; 
N = 1.0; 

no_of_simulations= 256*10^3;       # change number of iterations here 

# calculate mu and gamma
mu = r - q + delta*(sqrt(alpha^2 - (beta+1)^2) - sqrt(alpha^2 - beta^2));
gamma = sqrt(alpha^2 - beta^2)

stock_prc <- rep(0, no_of_simulations)      # create array of possible stock_prices (1 for each iteration)
put_prc <- rep(0, no_of_simulations)        # create array of possible put_prices (1 for each iteration)
call_prc <- rep(0, no_of_simulations)       # create array of possible call_prices (1 for each iteration)
euro_vanilla_put <- rep(0, no_of_simulations)    


# Box-Muller Function Template --------------------------------------------
BoxMuller <- function() {
  n = 2*N; 
  z = numeric(n);
  
  u1 = runif(n/2,0,1); 
  u2 = runif(n/2,0,1)
  
  z1 = sqrt(-2*log(u1)) * cos(2*pi*u2)      # half of normal variates
  z2 = sqrt(-2*log(u1)) * sin(2*pi*u2)      # other half
  z[seq(1, n, by=2)] = z1                   # interleave
  z[seq(2, n, by=2)] = z2                   # two halves
  
  return(z1)                                # return half
}


# Algorithm 1: Normal Inverse Gaussian ------------------------------------
for (j in 1:no_of_simulations) 
{
  # Step 1: Generate G1 and compute zeta ----------------------------------------------------
  G1 <- BoxMuller()           # generate standard normal random variable G1
  Z <- (G1^2)/gamma           # compute Z
  
  zeta <- rep(0, N)
  for (i in 1:N) 
    {
      zeta[i] = (1/gamma)*((delta*i) + (0.5*Z[i]) - sqrt((delta*i*Z[i]) + (Z[i]^2)/4))  
    }

  
  # Step 2: Generate uniform random variable U on (0,1) ----------------------------------
  Uniform <- function()       # function to generate uniform r.v. using Algorithm 1, Step 2, on pg. 17
  {
    U = runif(N,0,1)
    zt <- rep(0,N)
    for (i in 1:N)  {
      if (U[i] < (delta)/(delta+gamma*zeta[i])) {
        zt[i] = zeta[i] 
      } else  {
        zt[i] = (delta^2) / (gamma^2 * zeta[i]) 
      }
    }
    return(zt)
  } 
  
  zt = Uniform()          # generate zt using uniform function above  
  

  # Step 3: Generate StdNorm variable G2 -------------------------------------------------
  G2 <- BoxMuller()
  
  St <- rep(0,N)        # create empty vector for St
  Xt <- rep(0,N)        # create empty vector for Xt 
  
  for (i in 1:N)  
    {
    Xt[i] = mu*i + beta*zt[i] + sqrt(zt[i])*G2[i]     # compute Xt = mu*t + beta*zt + sqrt(zt)*G2
    St[i] = s0 * exp(Xt[i])                           # compute St = s0*e^Xt
    }
  
  stock_prc[j] = St[N]                                  # reassign variable, ie St = ST (stock value at maturity)
  put_prc[j] = exp(-r*T) * max(0, K - stock_prc[j])     # compute put price at maturity 
  call_prc[j] = exp(-r*T) * max(0, stock_prc[j] - K)    # compute call price at maturity


# Section 5.4, pg.19: European Vanilla Options ----------------------------
  euro_vanilla_put[j] = s0 * exp(-r * T) * max(0, K/s0 - exp(log(stock_prc[j] / s0)))

  # V[j] = exp(-r * T) * max(0, K - s0 * exp(log(stock_prc[j] / s0)))         # Section 5.4, pg. 19
  # V[j] = s0 * exp(-r * T) * max(0, K/s0 - exp(Xt))                          # Section 5.4, pg. 19

}
# ------------------------------------------------------------------------------
# END MC -----------------------------------------------------------------------


# Output: Average the calculated price: ---------------------------------------------------
"NIG Stock Price (t=T)"
sum(stock_prc) / no_of_simulations

nig.call.value <- sum(call_prc) / no_of_simulations
"NIG Call Value: "
nig.call.value

nig.put.value <- sum(put_prc) / no_of_simulations
"NIG Put Value: " 
nig.put.value

euro_vanilla_put.value <- sum(euro_vanilla_put) / no_of_simulations
"European Vanilla Put Value: " 
euro_vanilla_put.value

values.table <- cbind(nig.put.value,euro_vanilla_put.value)
values.table

# RQuantLib benchmark -----------------------------------------------------
#  recalculate using RQuantLib EuropeanOption (Black-Scholes) with same inputs (except volatility)
library(RQuantLib)

ql.put.value <- EuropeanOption(type = "put", underlying = 100, strike = 100,
                         dividendYield = 0.02, riskFreeRate = 0.05, maturity = 0.5, volatility = 0.2519)
put_prices <- cbind(nig.put.value,euro_vanilla_put.value,ql.put.value)
put_prices

# ql.call <- EuropeanOption(type = "call", underlying = 100, strike = 100,
#                                dividendYield = 0.02, riskFreeRate = 0.05, maturity = 0.5, volatility = 0.30)
# call_prices <- cbind(ql.call, nig.call)
# call_prices

put_prices<-cbind(nig.put.value,ql.put.value)
put_prices
