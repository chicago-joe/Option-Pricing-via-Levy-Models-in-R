# *Options Pricing in 𝐋𝐞́𝐯𝐲 Models*
#####_Joseph Loss, Yuchen Duan, Daniel Liberman_
We consider options pricing in Lévy models, specifically the implementation of two algorithms
listed in Feng3 and compared our results with those published in the original paper. Like Feng, we
assume that the asset price follows a geometric Lévy process under a given and equal martingale
measure. 

This is expressed as:
$$ 𝑆_𝑡 = 𝑆_0 e ^ {X_t} $$

where $𝑋_𝑡$ is the Lévy process at 𝑡 = 0 and that $S_0$ is the asset price at 𝑡 = 0.
&nbsp;&nbsp;&nbsp;&nbsp; 

#### Algorithm 1: Normal Inverse Gaussian
Our first algorithm is an implementation of the normal inverse Gaussian (NIG) process, which we
simulated as a subordinated Brownian motion, to compute the price of a European put option contract.

The NIG process is characterized by:
$$X_t = \mu t + \beta z_t + B_{z_t}   $$ where $𝐵_𝑡$ is a standard Brownian motion, $𝑧_𝑡$ is an independent inverse Gaussian process and
$$ \mu = r - q + \delta(\sqrt{\alpha^2-(\beta+1)^2}  -\sqrt{\alpha^2-\beta^2} ) $$ Note that $\alpha, \beta, \delta, r, q $ are given by Feng as inputs; we will discuss these shortly.
&nbsp;&nbsp;&nbsp;&nbsp; 

For 𝑡 > 0, we simulate $𝑋_𝑡$ with the following steps:
1. Generate a standard normal random variable using the Box-Muller algorithm described on pg.66 of Glasserman4. With $ \lambda = \sqrt{\alpha^2-\beta^2}, Z = \frac{G^2_1}{\lambda}, $ we compute
$$ \zeta = \frac{1}{\lambda}(\delta t + \frac{1}{2}Z-\sqrt{\delta tZ + Z^2/4}). $$

2. Generate a uniform random variable 𝑈 on (0, 1). If
$$  U < \delta t / (\delta t +\gamma \zeta ) $$
3. 
4.  (𝑈 < 𝛿𝑡𝛿𝑡 +𝛾𝜁), then 𝑧𝑡 = 𝜁.
Otherwise, 𝑧𝑡 =
𝛿2𝑡2
𝛾2𝜁 .

3. Generate a standard normal random variable 𝐺2. Compute 𝑋𝑡 = 𝜇𝑡 + 𝛽𝑧𝑡 + √𝑧𝑡𝐺2.
1 University of Illinois at Urbana-Champaign, M.S. Financial Engineering
2 University of Illinois at Chicago, Finance
3 Feng, Liming, et al. “Simulating Lé vy Processes from Their Characteristic Functions and Financial Applications.”
University of Illinois, 30 July 2011.
4 Glasserman, Paul. “Monte Carlo Methods in Financial Engineering.” Springer-Verlag, 2003.To compute the price of a European vanilla put option, we use the inputs given in Section 6.1 (pg. 22) of
Feng: 𝛼 = 15, 𝛽 = −5, 𝛿 = 0.5, 𝑟 = 0.05, 𝑞 = 0.02, 𝑆0 = 𝐾 = 100, 𝑇 = 0.5.
The price of the option at 𝑡 = 0 can be calculated as:
𝑉 = 𝑆0𝑒−𝑟𝑇Ε[f(XT)]
where 𝑋𝑡 is a Lé vy process and 𝑓(𝑥) = max(0, 𝐾/𝑆0 − 𝑒𝑥).
Using the inputs above, we computed the value of the option, $6.25836. We compared this with
the Black-Scholes model for a European option in the “RQuantLib” package:


&nbsp;&nbsp;&nbsp;&nbsp; 
&nbsp;&nbsp;&nbsp;&nbsp; 
$$
Algorithm 2: Inverse Fourier Transform
The second option pricing model is an implementation of the inverse transform method from
tabulated probabilities which, depending on the desired accuracy required, could be multiple times
faster than the normal inverse gaussian process.
Note tInput parameters are taken from Section 6.1, pg. 22 of Feng’s paper.
First, we begin by initializing arrays for the lists of variables 𝜒 and 𝐹̂. Brute-Force-Search is used
in place of the binary search originally prescribed in Section 3.1. This can be seen in the function
templates for chi and Fhat.
There are two more function templates that build the foundation for this algorithm. The first,
“Fhat Distribution Function”, is the direct implementation of the distribution in equation 3.12 on pg. 8 of
Feng’s paper. The second, termed “Inverse Transform Function”, implements the approximation to
𝐹−1(𝑈) on pg. 8 of the paper:
This performs the inverse transform function for each generated U between 0 and 1, utilizing bruteforce search to find 0 ≤ 𝑘 ≤ 𝐾 − 1 so that 𝐹̂𝑘 ≤ 𝑈 < 𝐹̂𝑘+1:The Inverse Transform algorithm begins with its parameters on line 102; these parameters are
taken from Section 6.1 on pg. 22 of Feng’s paper. In computing the European put price, there are two
formulas used. The first method uses the formula for “European Vanilla Put Options” given in Section
5.4, pg. 19 of Feng’s paper (note the max function = max (0, 𝑠𝑡𝑟𝑖𝑘𝑒 𝑆0 − 𝑒𝑥). The second method uses the
more general, max(0, 𝑠𝑡𝑟𝑖𝑘𝑒 − 𝑆0𝑒𝑋𝑇).
The computed option prices are averaged over the number of MonteCarlo iterations and an example
output is shown below:
Note that both prices converge to a value of $4.58 for the European put option. This is extremely close
to the $4.589 value of Feng’s implementation (listed on pg. 24).
Final Thoughts
We are currently working on making both algorithms available in a package format, which can be
downloaded, used, and improved by the community.