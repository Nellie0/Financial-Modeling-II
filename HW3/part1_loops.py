import numpy as np
import HW1

# This is my second attempt at completing part 1 when I wasn't sure if my recursive approach was working. I *think* it is working. I'd prefer you grade part1_recursive, but figured I'd still upload this attempt.

def bisection(S, K, V, r, q, T, isCall=True):
    sigma = 0.3
    up = 1
    down = 0.001
    count = 0
    diff = HW1.BS(S, K, r, q, sigma, T, isCall) - V

    # want difference to be small enough according to tolerance
    while abs(diff) > 0.0001 and count < 1000:
        if diff < 0:
            down = sigma
            sigma = (up + sigma) / 2
        else:
            up = sigma
            sigma = (sigma + down) / 2

        diff = HW1.BS(S, K, r, q, sigma, T, isCall) - V
        count += 1

    return sigma
    
def newton(S, K, V, r, q, T, isCall=True):
    sigma = 0.3
    for i in range(100):
        diff = HW1.BS(S, K, r, q, sigma, T, isCall) - V

        # End loop if the difference goes below tolerance level
        if abs(diff) < 0.0001:
            break

        sigma = sigma - (diff / HW1.vega(S, K, r, q, sigma, T))
    
    return sigma