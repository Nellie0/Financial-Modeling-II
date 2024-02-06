import numpy as np
from scipy.stats import norm

def d1d2(S, K, r, q, sigma, T):
    d1 = (np.log(S/K) + ((r - q) + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    return d1, d2

def BS(S, K, r, q, sigma, T, isCall = True):
    d1 = d1d2(S, K, r, q, sigma, T)[0]
    d2 = d1d2(S, K, r, q, sigma, T)[1]
    if isCall == True:
        return S * norm.cdf(d1) - K * np.exp(-(r - q) * T) * norm.cdf(d2)
    else:
        return K * np.exp(-(r - q) * T) * norm.cdf(-d2) - S * norm.cdf(-d1)

def delta(S, K, r, q, sigma, T, isCall = True):
    d1 = d1d2(S, K, r, q, sigma, T)[0]
    if isCall == True:
        return norm.cdf(d1)
    else:
        return norm.cdf(d1) - 1
    
def gamma(S, K, r, q, sigma, T):
    d1 = d1d2(S, K, r, q, sigma, T)[0]
    return norm.pdf(d1) / (S * sigma * np.sqrt(T))

def vega(S, K, r, q, sigma, T):
    d1 = d1d2(S, K, r, q, sigma, T)[0]
    return S * norm.pdf(d1) * np.sqrt(T)

def theta(S, K, r, q, sigma, T, isCall = True):
    d1 = d1d2(S, K, r, q, sigma, T)[0]
    d2 = d1d2(S, K, r, q, sigma, T)[1]
    if isCall == True:
        return (-S * norm.pdf(d1) * sigma) / (2 * np.sqrt(T)) - (r-q) * K * np.exp(-(r-q) * T) * norm.cdf(d2)
    else:
        return (-S * norm.pdf(d1) * sigma) / (2 * np.sqrt(T)) + (r-q) * K * np.exp(-(r-q) * T) * norm.cdf(-d2)
    
def rho(S, K, r, q, sigma, T, isCall = True):
    d2 = d1d2(S, K, r, q, sigma, T)[1]
    if isCall == True:
        return K * T * np.exp(-(r-q) * T) * norm.cdf(d2)
    else:
        return -K * T * np.exp(-(r-q) * T) * norm.cdf(-d2)