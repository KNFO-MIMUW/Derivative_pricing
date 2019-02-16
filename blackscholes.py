import numpy as np
import math
from scipy.stats import norm
from random import gauss
import matplotlib.pyplot as plt

"""
  Variables:
  S_0     underlying price at time 0
  S       stock price
  K       strike price
  r       risk free interest rate
  sigma   volatility
  t       time until option expiration
  T       maturity time of option
"""

def asset_price(S_0, r, sigma, t):

    return S_0 * math.exp(t * (r - 0.5 * sigma ** 2) + sigma * gauss(0,t))

def BSM(S_t, K, r, sigma, t):

        d1 = (np.log(S_t/K) + t * (r + (sigma**2)/2))/(sigma * math.sqrt(t))
        d2 = d1 - sigma * math.sqrt(t)

        return S_t * norm.cdf(d1) - math.exp(-r * t) * K * norm.cdf(d2)

def asset_price_T(S_0, r, sigma, T):

    return S_0 * math.exp(T * (r - 0.5 * sigma ** 2) + math.sqrt(T) * sigma * gauss(0,1))

def monte_carlo(S_0, K, r, sigma, T):

    num_simulations = 100000
    factor = math.exp(-r * T)
    payoffs = []

    for i in range(num_simulations):
        S_T = asset_price_T(S_0, r, sigma, T)
        payoffs.append(max(S_T-K,0))

    return factor * (sum(payoffs)/num_simulations)

def monte_carlo_full(S_0, K, r, sigma, t):

    num_days = int(t * 365)
    num_simulations = 100000
    dt = t / num_days

    S = np.zeros((num_days,num_simulations))
    S[0] = S_0

    for t in range(1,num_days):
        n = np.random.standard_normal(num_simulations)
        S[t] = (S[t-1] * np.exp((r - 0.5 * sigma ** 2) * dt + np.sqrt(dt) * sigma * n))

    # payoff = math.exp(-r * t) * np.sum(np.maximum(S[-1] - K, 0)) * 1 / num_simulations

    plt.plot(S[:,0:15]) #plot of 15 first paths
    plt.grid(True)
    plt.xlim([0,num_days])
    plt.show()



S_0 = float(input("underlying price at time 0: "))
K = float(input("strike price: "))
r = float(input("interest rate: "))
sigma = float(input("volatility: "))
t = float(input("time until option expiration in years: "))
S = asset_price(S_0,r,sigma,t)

print("the price of the option calculated by Black-Scholes formula: %.4f" %BSM(S,K,r,sigma,t))
print("the price of the option calculated by Monte Carlo method: %.4f" %monte_carlo(S,K,r,sigma,t))
monte_carlo_full(S_0,K,r,sigma,t)
