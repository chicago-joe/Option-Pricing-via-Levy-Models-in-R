###########################################################################
# Pricing European Put Options Through Simulated Levy Processes 
# Pseudocode implemented from Professor Liming Feng's Paper (pg. 08)
#
# Created by Joseph Loss on 1/01/2019
# Co-developers: Yuchen Duan (UIUC MSFE) and Daniel Liberman (UIC Finance)
#
# Algorithm 2: Inverse-Transform Algorithm
#
###########################################################################

# Inital Parameters  -------------------------------------------------------------
# taken from Feng, Table 1 pg. 23
chi.0 = -0.477
chi.K = 0
K = 22


# List Function (for initializing variable arrays) ------------------------
seqlist<-function(a, b, n)  {
  list1 <- rep(0, n)
  list1[1] = a
  list1[n] = b
  for (i in 2:(n-1)) {
    list1[i] = (list1[n] - list1[i-1])*.326 + list1[i-1]
  }
  return(list1)
}


# Initialize Variable Arrays ----------------------------------------------
chi.list = seqlist(chi.0, chi.K, K+1)
Fhat.list = seq(0, 1, by=(1)/K)
eta = (Fhat.list[K+1] - Fhat.list[1]) / K
Fhat.list

# Brute-Force-Search function: chi ------------------------------------
# input = -0.1
brute_force_search.chi <- function(input) {
  stock_prices.list = 0
  for (variable in chi.list) {
    if(variable > input)  {
      return(stock_prices.list)
    }
    stock_prices.list = stock_prices.list+1
  }
}

# Brute-Force-Search function: Fhat ------------------------------------
brute_force_search.Fhat <- function(input)  {
  stock_prices.list = 0
  for (variable in Fhat.list) {
    if(variable > input)  {
      return(stock_prices.list)
    }
    stock_prices.list = stock_prices.list+1
  }
}

# Fhat Distribution Function ----------------------------------------------
# piecewise function for Fhat(x) taken from Feng, pg. 8 (Equation 3.12)
Fhat.distribution <- function(input) 
{
  stock_prices.list = 0
  if(input < chi.list[1])                     # x < x0 
  {                  
    stock_prices.list = 0
    return(stock_prices.list)
  } 
  else if(input >= chi.list[K+1])             # x > 0 xK  
  {        
    stock_prices.list = 1
    return(stock_prices.list)
  } 
  else                                        # (xk-1) <= chi
  {
    xk_1 =  brute_force_search.chi(input)     # xk = bs(input); 
    
    # Fhat[k-1] location:
    stock_prices.list = Fhat.list[xk_1] + (Fhat.list[xk_1+1] - Fhat.list[xk_1])/eta * (input-chi.list[xk_1])  
  }
  return(stock_prices.list)
}

# Inverse Transform Function ------------------------------------------------
# Approximation to F-1(U) using brute-force search
inverse_transform_method <- function()
{
  U = runif(1,0,1)
  xk_1 = brute_force_search.Fhat(U)   
  if(xk_1 < (K-1))              # function returns K-1
  {                       
    return(chi.list[xk_1+1] + (chi.list[xk_1+2] - chi.list[xk_1+1]) / (Fhat.list[xk_1+2] - Fhat.list[xk_1+1])*(U-Fhat.list[xk_1+1]))
  } 
  else
  {
    return(0)
  }
}


# Inverse Transform Method: Parameters -------------------------------
# taken from Feng, Section 6.1 pg. 22
Xt = inverse_transform_method()
s0 = 100;
strike = 100;
T = 0.5;
r = 0.03;

no_of_simulations = 1024*10^3     # change number of MC iteriations here

# initialize price lists
stock_prices.list <- rep(0, no_of_simulations)
put_prices.list <- rep(0, no_of_simulations)


# Inverse Transform Method 1 -----------------------------------------------
# Calculate V using Section 5.4 of Feng's Paper (pg.19) 
V <- rep(0, no_of_simulations)
for (j in 1:no_of_simulations) 
{
  V[j] = s0 * exp(-r*T) * max(0, strike/s0 - exp(inverse_transform_method()))
}

# Method 1 Result
put.prc.fft1 <- sum(V) / no_of_simulations

# plot the cdf
plot(chi.list)

# END  --------------------------------------------------------------------


# OPTIONAL: Inverse Transform Method 2 ----------------------------------------------
# uses default calculation for option price: V = e^(-r*T) * max(0, K - St)
for (j in 1:no_of_simulations) 
{
  Xt = inverse_transform_method()
  stock_prices.list[j] = s0 * exp(Xt)
  put_prices.list[j] = exp(-r * T) * max(0, strike - stock_prices.list[j])
  
}

# Method 2 Result
sum(stock_prices.list) / no_of_simulations
put.prc.fft2 <- sum(put_prices.list) / no_of_simulations


# Output both methods to table:
# note both methods should output very similar results:
values.table <- cbind(put.prc.fft1,put.prc.fft2)
values.table


# RQuantLib benchmarking -----------------------------------------------------
# Use the EuropeanOption function of RQuantLib to calculate the Black-Scholes
# price of the European Put (all inputs are identical to the NIG inputs above, except volatility)
library(RQuantLib)
q = 0.0

# QuantLib Price of the European Put option using Black-Scholes
# to be compared with the FFT price and the PC-Parity Price
quantlib.bsm.put.prc <- EuropeanOption(type = "put", underlying = s0, strike = strike,
                               dividendYield = q, riskFreeRate = r, maturity = 0.5, volatility = 0.25)

put_prices<-cbind(put.prc.fft1, quantlib.bsm.put.prc)

print("NIG Put Price vs Black-Scholes:")
put_prices

## NMOF Library function for PCP:
# library(NMOF)
# putCallParity("call",put=put.prc.fft1,S=s0,X=strike,tau=T,r=r,q=q)
## result = $6.10

# QuantLib Price of the European Call option using Black-Scholes
quantlib.bsm.call.prc<-RQuantLib::EuropeanOption("c",s0,strike,q,r,T,0.25)
quantlib.bsm.call.prc

# Put-Call-Parity formula from scratch:
Put_Parity_Prc = quantlib.bsm.call.prc$value - (s0*exp((-q * T)) - (strike * exp(-r * T)))
Put_Parity_Prc

# Call_Prc_Parity_Table<-cbind(FFT_Call_Parity_Prc, quantlib.bsm.call.prc)
# Call_Prc_Parity_Table

Put_Prc_Parity_Table<-cbind(put.prc.fft1, Put_Parity_Prc, quantlib.bsm.put.prc)
Put_Prc_Parity_Table

