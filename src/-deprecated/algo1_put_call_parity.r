# # BENCHMARKING PURPOSES ONLY ----------------------------------------------
#
# # Put-Call Parity Verification --------------------------------------------
# # PCP: Call Price
# Call_Parity_Prc = euro_vanilla_put.value + s0*exp((-q * T)) - (K * exp(-r * T))
# Call_Parity_Prc
#
# # PCP: Call Price
# Put_Parity_Prc = Call_Parity_Prc - (s0*exp((-q * T)) - (K * exp(-r * T)))
# Put_Parity_Prc
# Put_Parity_Table <- cbind(Put_Parity_Prc, NIG_Put_Prc)
#
# nig_comparison_table
# (algo1.NIG.mu - GenHyp.NIG.mu)/GenHyp.NIG.mu         # computed difference in models
# Put_Parity_Table
#
# # RQuantLib BlackScholes Benchmarking --------------------------------------------------
# # Used RQuantLib's built-in BlackScholesModel to compute the price of the European put option
#
# library(RQuantLib)
#
# # all inputs are identical to the NIG inputs above, except volatility)
# QuantLib_Put_Prc <- EuropeanOption(type = "put",
#                                    underlying = 100,
#                                    strike = 100,
#                                    dividendYield = 0.02,
#                                    riskFreeRate = 0.05,
#                                    maturity = 0.5,
#                                    volatility = 0.2519)
#
# put_prices<-cbind(NIG_Put_Prc, QuantLib_Put_Prc)
# print("NIG Put Price vs Black-Scholes:")
# put_prices
#
# # QuantLib Price of European Call option:
# QuantLib_Call_Prc <- EuropeanOption(type = "c",
#                                     underlying = s0,
#                                     strike = K,
#                                     dividendYield = q,
#                                     riskFreeRate = r,
#                                     maturity = T,
#                                     volatility = 0.25)
#
# Call_Parity_Table <- cbind(Call_Parity_Prc, QuantLib_Call_Prc)
# Call_Parity_Table
