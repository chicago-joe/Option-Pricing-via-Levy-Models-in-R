###########################################################################
#       Pricing European Put Options Through Simulated Levy Processes
#       Research proposed by Liming Feng, Zisheng Chen, and Xiong Lin
#               University of Illinois, Urbana-Champaign
#
# Created by Joseph Loss on 1/01/2019
# Co-developers: Yuchen Duan (UIUC MSFE) and Daniel Liberman (UIC Finance)
#
# Algorithm 2: Inverse-Transform Algorithm
#
###########################################################################
#Setup Sims & Time
times <- matrix(rep(0,9), 3, 3)

no_of_simulations.1 <-  256 * 10 ^ 3     # change number of MC iteriations here
no_of_simulations.2 <- (256 * 4) * 10 ^ 3
no_of_simulations.3 <- (256 * 8) * 10 ^ 3
no_of_simulations.4 <- (256 * 16) * 10 ^ 3
no_of_simulations.5 <- (256 * 32) * 10 ^ 3
no_of_simulations.6 <- (256 * 64) * 10 ^ 3

no_of_simulations.list <- c(no_of_simulations.1, no_of_simulations.2, no_of_simulations.3, no_of_simulations.4, no_of_simulations.5, no_of_simulations.6)


# Brute-Force-Search function: chi ------------------------------------
# input = -0.1
brute_force_search.chi <- function(input) {
  stock_prices.list = 0
  for (variable in chi.list) {
    if (variable > input)  {
      return(stock_prices.list)
    }
    stock_prices.list = stock_prices.list + 1
  }
}

# Brute-Force-Search function: Fhat ------------------------------------
brute_force_search.Fhat <- function(input)  {
  stock_prices.list = 0
  for (variable in Fhat.list) {
    if (variable > input)  {
      return(stock_prices.list)
    }
    stock_prices.list = stock_prices.list + 1
  }
}

#Binary Search
binary_search <- function(values, input, start, end) {
  mid = floor((start + end) / 2)
  if (values[mid] <= input && values[mid + 1] >= input) {
    return(mid)
  }
  if (values[mid] > input) {
    end = mid
    return(binary_search(values, input, start, end))
  }
  if (values[mid] < input) {
    start = mid
    return(binary_search(values, input, start, end))
  }
}

# Box-Muller Function Template --------------------------------------------
BoxMuller <- function() 
{
  n = 2 * N
  z = numeric(n)
  
  u1 = runif(n / 2, 0, 1)
  u2 = runif(n / 2, 0, 1)
  
  z1 = sqrt(-2 * log(u1)) * cos(2 * pi * u2)      # half of normal variates
  z2 = sqrt(-2 * log(u1)) * sin(2 * pi * u2)      # other half
  z[seq(1, n, by = 2)] = z1                   # interleave
  z[seq(2, n, by = 2)] = z2                   # two halves
  
  return(z1)                                # return half
}

# Fhat Distribution Function ----------------------------------------------
# piecewise function for Fhat(x) taken from Feng, pg. 8 (Equation 3.12)
Fhat.distribution <- function(input)
{
  stock_prices.list = 0
  if (input < chi.list[1])
    # x < x0
  {
    stock_prices.list = 0
    return(stock_prices.list)
  }
  else if (input >= chi.list[K + 1])
    # x > 0 xK
  {
    stock_prices.list = 1
    return(stock_prices.list)
  }
  else
    # (xk-1) <= chi
  {
    xk_1 =  binary_search(chi.list, input, 1, length(chi.list))
    #print(binary_search(chi.list,input,1,length(chi.list)))# xk = bs(input);
    #xk_1 = brute_force_search.chi(input)
    #print(brute_force_search.chi(input))
    
    # Fhat[k-1] location:
    stock_prices.list = Fhat.list[xk_1] + (Fhat.list[xk_1 + 1] - Fhat.list[xk_1]) /
      eta * (input - chi.list[xk_1])
  }
  return(stock_prices.list)
}

# Inverse Transform Function ------------------------------------------------
# Approximation to F-1(U) using brute-force search
inverse_transform_method <- function() {
  U = runif(1,0,1)
  xk_1 = binary_search(Fhat.list, U, 1, length(Fhat.list))
  
  if (xk_1 < (K-1)) {                   # function returns K-1
    return(chi.list[xk_1 + 1] + (chi.list[xk_1 + 2] - chi.list[xk_1 + 1]) / 
             (Fhat.list[xk_1 + 2] - Fhat.list[xk_1 + 1]) * (U - Fhat.list[xk_1 + 1]))
  } else {
    return(0) 
  }
}


#For Loop
for(i in 3:3){
  no_of_simulations <- no_of_simulations.all[i]
  #NIG Method
  #Define Variables
  alpha = 15;                         # CORRECTED PARAMETER? 
  beta = -5;                          # Paper gave alpha=15, this generated an incorrect option price of $6.26  
  delta = 0.5; 
  r = 0.05; 
  q = 0.02; 
  s0 = 100; 
  K = s0; 
  T = 0.5; 
  N = 1.0; 
  
  # no_of_simulations= 4096 * 10^3;       # change number of iterations here 
  
  # calculate mu using the formula given at the top of pg. 19
  mu = r - q + delta*(sqrt(alpha^2 - (beta+1)^2) - sqrt(alpha^2 - beta^2));
  
  # calculate gamma using the formula given in the Algorithm 1 pseudocode
  gamma = sqrt(alpha^2 - beta^2)
  
  stock_prc <- rep(0, no_of_simulations)       # create array of possible stock_prices (1 for each iteration)
  put_prc <- rep(0, no_of_simulations)         # create array of possible put_prices (1 for each iteration)
  call_prc <- rep(0, no_of_simulations)        # create array of possible call_prices (1 for each iteration)
  euro_vanilla_put <- rep(0, no_of_simulations)    
  nigv <- rep(0, no_of_simulations)            # create NIG random vector
  
  #Start Iterations
  start.0 <- Sys.time()
  for (j in 1:no_of_simulations) 
  {
    # Step 1: Generate G1 and compute zeta --------------------------------------------------
    G1 <- BoxMuller()           # generate standard normal random variable G1
    Z <- (G1^2)/gamma              # compute Z
    
    zeta <- rep(0, N)
    for (i in 1:N) 
    {
      zeta[i] = (1/gamma)*((delta*T) + (0.5*Z[i]) - sqrt((delta*T*Z[i]) + (Z[i]^2)/4))  
    }
    
    
    # Step 2: Generate uniform random variable U on (0,1) -----------------------------------
    Uniform <- function()       # function to generate uniform r.v. using Algorithm 1, Step 2, on pg. 17
    {
      U = runif(N, 0, 1)
      zt <- rep(0, N)
      for (i in 1:N)  
      {
        if (U[i] < (delta*T) / ((delta*T) + gamma * zeta[i])) {
          zt[i] = zeta[i]
        } else  {
          zt[i] = (delta ^ 2 * T^2) / (gamma ^ 2 * zeta[i])
        }
      }
      return(zt)
    }
    
    zt = Uniform()          # generate zt using uniform function above
    
    
    # Step 3: Generate StdNorm variable G2 -------------------------------------------------
    G2 <- rnorm(1,0,1)
    
    St <- rep(0, N)        # create empty vector for St
    Xt <- rep(0, N)        # create empty vector for Xt 
    
    for (i in 1:N)
    {
      Xt[i] = mu * T + beta * zt[i] + sqrt(zt[i]) * G2[i]     # compute Xt = mu*t + beta*zt + sqrt(zt)*G2
      St[i] = s0 * exp(Xt[i])                         # from pg 16. compute St = S0 * e^(Xt)
    }
    
    stock_prc[j] = St[N]                                  # reassign variable, ie St = ST (stock value at maturity)
    put_prc[j] = exp(-r*T) * max(0, K - stock_prc[j])     # compute put price at maturity 
    call_prc[j] = exp(-r*T) * max(0, stock_prc[j] - K)    # compute call price at maturity
    nigv[j] = Xt[N]                                       
    
    
    # Section 5.4, pg.19: European Vanilla Options -------------------------------------------
    # We also calculate the option value "V" using the formula given in Section 5.4.
    
    # note that the resulting output is identical 
    # to the result we generated through St = S0 * e^(Xt) above. 
    euro_vanilla_put[j] = s0 * exp(-r * T) * max(0, K/s0 - exp(log(stock_prc[j] / s0)))
    
  }
  total.0 <- start.0 - Sys.time()
  #IFF Methods
  #Define Variables
  chi.0 = -0.477
  chi.K = 0
  K = 22
  times <- matrix(rep(0,9), 3, 3)
  
  # Inverse Transform Method: Parameters -------------------------------
  # taken from Feng, Section 6.1 pg. 22
  s0 = 100
  strike = 100
  T = 0.5
  r = 0.03
  
  # List Function (for initializing variable arrays) ------------------------
  seqlist <- function(a, b, n)  {
    list1 <- rep(0, n)
    list1[1] = a
    list1[n] = b
    for (i in 2:(n - 1)) {
      list1[i] = (list1[n] - list1[i - 1]) * .326 + list1[i - 1]
    }
    return(list1)
  }
  
  
  # Initialize Variable Arrays ----------------------------------------------
  chi.list = seqlist(chi.0, chi.K, K + 1)
  Fhat.list = seq(0, 1, by = (1) / K)
  eta = (Fhat.list[K + 1] - Fhat.list[1]) / K
  Fhat.list
  
  Xt = inverse_transform_method()
  
  #Setup empty Vectors
  stock_prices.list <- rep(0, no_of_simulations)
  put_prices.list <- rep(0, no_of_simulations)
  
  # Inverse Transform Method 1 -----------------------------------------------
  # Calculate V using Section 5.4 of Feng's Paper (pg.19)
  put_prcs <- rep(0, no_of_simulations)
  start.1 <- Sys.time()
  for (j in 1:no_of_simulations) {
    put_prcs[j] = s0 * exp(-r*T) * max(0, strike/s0 - exp(inverse_transform_method()))
  }
  total.1 <- Sys.time() - start.1
  
  # Method 1 Result
  put_prc.InvT <- sum(put_prcs) / no_of_simulations
  
  
  #plot(chi.list)
  # END  --------------------------------------------------------------------
  
  
  # OPTIONAL: Inverse Transform Method 2 ----------------------------------------------
  # uses default calculation for option price: V = e^(-r*T) * max(0, K - St)
  start.2 <- Sys.time()
  for (j in 1:no_of_simulations) 
  {
    Xt = inverse_transform_method()
    stock_prices.list[j] = s0 * exp(Xt)
    put_prices.list[j] = exp(-r * T) * max(0, strike - stock_prices.list[j])
  }
  total.2 <- Sys.time() - start.2
  
  # Method 2 Result
  sum(stock_prices.list) / no_of_simulations
  put.prc.fft2 <- sum(put_prices.list) / no_of_simulations
  
  #Time matrix
  
  times[i,1] <- total.0
  times[i,2] <- total.1
  times[i,3] <- total.2
  
}

colnames(times) <- c("NIG", "IFF 1", "IFF 2")
rownames(times) <- c("256x10^3", "1024x10^3", "2048x10^3")

# Output both methods to table:
# note both methods should output very similar results:
values.table <- cbind(put_prc.InvT, put.prc.fft2)
values.table


# RQuantLib benchmarking -----------------------------------------------------
# Use the EuropeanOption function of RQuantLib to calculate the Black-Scholes
# price of the European Put (all inputs are identical to the NIG inputs above, except volatility)
library(RQuantLib)
q = 0.0

# QuantLib Price of the European Put option using Black-Scholes
# to be compared with the FFT price and the PC-Parity Price
quantlib.bsm.put.prc <-
  EuropeanOption(
    type = "put",
    underlying = s0,
    strike = strike,
    dividendYield = q,
    riskFreeRate = r,
    maturity = 0.5,
    volatility = 0.19
  )

put_prices <- cbind(put.prc.fft1, quantlib.bsm.put.prc)

print("NIG Put Price vs Black-Scholes:")
put_prices