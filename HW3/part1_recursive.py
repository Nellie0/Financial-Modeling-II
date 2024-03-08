import numpy as np
import HW1

# Bisection method
def bisection(S, K, V, r, q, T, upper_bound = 10, lower_bound = 0, count = 0, isCall = True):
    if count <= 1000:
        if upper_bound - lower_bound > 0.0001:
            if HW1.BS(S, K, r, q, (upper_bound - lower_bound) / 2, T, isCall) - V > 0:
                upper_bound = (upper_bound + lower_bound) / 2 # New middle, upper bound moves
            else:
                lower_bound = (upper_bound + lower_bound) / 2

        count += 1
        upper_bound = upper_bound
        lower_bound = lower_bound
        
         # Repeat unil count=1000
        return bisection(S, K, V, r, q, T, upper_bound, lower_bound, count, isCall = True)
    
    # Count = 1000
    else:
        return lower_bound

# Newton's method
def newton(S, K, V, r, q, T, sigma = 0.3, count = 0, isCall = True):
    vega = HW1.vega(S, K, r, q, sigma, T)
    BS_V = HW1.BS(S, K, r, q, sigma, T, isCall)

    if count <= 1000:
        if abs(BS_V - V) < 0.001: #tolerance level
            sigma = sigma

        else:
            # Using formula given to us
            sigma = (BS_V - V) / vega
    
        count += 1
        
        return newton(S, K, V, r, q, T, sigma, count, isCall) # Repeat until count=1000
    
    else:
        return sigma
