###########################################################################

#       Pricing European Put Options Through Simulated Levy Processes 

#       Research proposed by Liming Feng, Zisheng Chen, and Xiong Lin

#               University of Illinois, Urbana-Champaign

#

# Created by Joseph Loss on 12/17/2018

# Co-developers: Yuchen Duan (UIUC MSFE) and Daniel Liberman (UIC Finance)

#

# Algorithm 1: Simulating a Normal Inverse Gaussian process through a Brownian subordination

#

###########################################################################
# set.seed(2017)
library(tictoc)
no_of_simulations.1 <-  256 * 10 ^ 3     # change number of MC iteriations here
no_of_simulations.2 <- (256 * 4) * 10 ^ 3
no_of_simulations.3 <- (256 * 8) * 10 ^ 3
no_of_simulations.4 <- (256 * 16) * 10 ^ 3
no_of_simulations.5 <- (256 * 32) * 10 ^ 3
no_of_simulations.6 <- (256 * 64) * 10 ^ 3

no_of_simulations.list <- c(no_of_simulations.1, no_of_simulations.2, no_of_simulations.3, no_of_simulations.4, no_of_simulations.5, no_of_simulations.6)

for (i in 1:6)
{
    
    # Input Parameters --------------------------------------------------------
    
    # taken from Feng et al, pg.22-23
    
    alpha = 15;                         # CORRECTED PARAMETER? 
    
    beta = -5;                          # Paper gave alpha=15, this generated an incorrect option price of $6.26  
    
    delta = 0.5; 
    
    r = 0.05; 
    
    q = 0.02; 
    
    s0 = 100; 
    
    K = s0; 
    
    T = 0.5; 
    
    N = 1.0; 
    
    no_of_simulations = no_of_simulations.list[i];       # change number of iterations here 
    
    
    # Start the clock!
    # ptm <- proc.time()
    # tic("Monte Carlo Simulation")
    tic(no_of_simulations)
    # calculate mu using the formula given at the top of pg. 19
    
    mu = r - q + delta*(sqrt(alpha^2 - (beta+1)^2) - sqrt(alpha^2 - beta^2));
    
    
    
    # calculate gamma using the formula given in the Algorithm 1 pseudocode
    
    gamma = sqrt(alpha^2 - beta^2)
    
    
    
    stock_prc <- rep(0, no_of_simulations)       # create array of possible stock_prices (1 for each iteration)
    
    put_prc <- rep(0, no_of_simulations)         # create array of possible put_prices (1 for each iteration)
    
    call_prc <- rep(0, no_of_simulations)        # create array of possible call_prices (1 for each iteration)
    
    euro_vanilla_put <- rep(0, no_of_simulations)    
    
    nigv <- rep(0, no_of_simulations)            # create NIG random vector
    
    
    
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
    # Algorithm 1: Normal Inverse Gaussian, Monte-Carlo Simulation --------------------------------------------
    
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
    # Stop the clock
    # proc.time() - ptm
    # ctime = c(proc.time())
    # }
    
    # END MonteCarlo Simulation --------------------------------------------------------------
    # ----------------------------------------------------------------------------------------
    # ctime <- c(toc(no_of_simulations))  
    # }
    
    
    
    # Average the computed prices: -----------------------------------------------------------
    
    "NIG Stock Price (t=T)"
    
    sum(stock_prc) / no_of_simulations
    
    
    # 
    # NIG_Put_Prc <- sum(put_prc) / no_of_simulations
    # 
    # "NIG Put Value: " 
    # 
    # NIG_Put_Prc
    
    
    
    euro_vanilla_put.value <- sum(euro_vanilla_put) / no_of_simulations
    
    # "European Vanilla Put Value: " 
    
    # euro_vanilla_put.value
    toc(log = TRUE)
    # toc()
    tic.clearlog()
    tic.clear()
}



values.table <- cbind(NIG_Put_Prc, euro_vanilla_put.value)

values.table





# Model Benchmarking & Verification ------------------------------------------------------------------

# NIG Verification --------------------------------------------------------

library(GeneralizedHyperbolic)      # testing our NIG model with GeneralizedHyperbolic built-in "rnig"



GenHyp.NIG <- rnig(no_of_simulations, delta=delta, alpha=alpha, beta=beta)+mu

algo1.NIG <- nigv



par(mfrow=c(1,2))       # Plot NIG distributions



# GeneralizedHyperbolic NIG Histogram

GenHyp.NIG.hist <-
    
    hist(GenHyp.NIG,
         
         breaks = 40,
         
         col = "grey",
         
         main = "GeneralizedHyperbolic NIG")

xfit <- seq(min(GenHyp.NIG), max(GenHyp.NIG), length = 40)

yfit <- dnorm(xfit, mean = mean(GenHyp.NIG), sd = sd(GenHyp.NIG))

yfit <- yfit * diff(GenHyp.NIG.hist$mids[1:2]) * length(GenHyp.NIG)

lines(xfit, yfit, col = "blue", lwd = 3)



# Algorithm 1 NIG Histogram

algo1.NIG.hist <-
    
    hist(algo1.NIG,
         
         breaks = 40,
         
         col = "grey",
         
         main = "Algorithm 1 NIG")

xfit <- seq(min(algo1.NIG), max(algo1.NIG), length = 40)

yfit <- dnorm(xfit, mean = mean(algo1.NIG), sd = sd(algo1.NIG))

yfit <- yfit * diff(algo1.NIG.hist$mids[1:2]) * length(algo1.NIG)

lines(xfit, yfit, col = "blue", lwd = 3)





# NIG Mean Comparison

GenHyp.NIG.mu <- (mean(rnig(no_of_simulations, delta=delta, alpha=alpha, beta=beta)+mu))

GenHyp.NIG.mu

algo1.NIG.mu = mean(nigv)                            # our NIG 

algo1.NIG.mu



nig_comparison_table <- cbind(algo1.NIG.mu, GenHyp.NIG.mu)

nig_comparison_table



(algo1.NIG.mu - GenHyp.NIG.mu)/GenHyp.NIG.mu         # computed difference in model means



values.table



RQuantLib::EuropeanOption('put',s0,K,q,r,T,0.191)
# toc()
# }