# LeastSquareMonteCarlo-PricingAmericanPutOptions
Using Least Square Monte Carlo (LSMC) simulation with different polynomials basis to price American put option.

(a) Use the LSMC method with 100,000 paths simulations (50,000 + 50,000 antithetic) to price an American put option 
with strike price of X = $40, maturity of 0.5-years, 1-years, 2-years, and current stock prices of $36, $40, $44. 
Use the first k of the Laguerre polynomials for k = 2, 3, 4.

(b) Use the LSMC method with 100,000 paths simulations (50,000 + 50,000 antithetic) to price an American put option 
with strike price of X = $40, maturity of 0.5-years, 1-years, 2-years, and current stock prices of $36, $40, $44. 
Use the first k of the Hermite polynomials for k = 2, 3, 4.

(c) Use the LSMC method with 100,000 paths simulations (50,000 + 50,000 antithetic) to price an American put option 
with strike price of X = $40, maturity of 0.5-years, 1-years, 2-years, and current stock prices of $36, $40, $44. 
Use the first k of the Simple Monomials polynomials for k = 2, 3, 4.


Forward-start options: <br>
Forward-start options are path dependent options that have strike prices to be determined at a future date. For example,
 a forward start put option payoff at maturity is max(St - ST, 0) <br>
 where the strike price of the put option is St. Here 0 <= t <= T. <br>
 
 (a) Estimate the value of the forward-start European put option on a stock with these characteristics: S0 = $65, K = $60,
 sigma = 20% per annum, risk-free rate is r = 6% per annum, t = 0.2 and T = 1.
 
 (b) Estimate the value of the forward-start American put option on a stock with these characteristics: S0 = $65, K = $60,
 sigma = 20% per annum, risk-free rate is r = 6% per annum, t = 0.2 and T = 1.
