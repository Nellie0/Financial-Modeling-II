import numpy as np
import HW1

# Option price function, CRR binomial model
def bm(S, K, r, q, sigma, T, n, t = 0, isCall = True, isEuropean = True):
    h = T / n
    u = np.exp(sigma * np.sqrt(h))
    d = 1 / u
    p = (np.exp((r - q) * h)  - d) / (u - d)

    # Greeks for later
    delta = HW1.delta(S, K, r, q, sigma, T, n)
    gamma = HW1.gamma(S, K, r, q, sigma, T)
    vega = HW1.vega(S, K, r, q, sigma, T)
    theta = HW1.theta(S, K, r, q, sigma, T)
    rho = HW1.rho(S, K, r, q, sigma, T)

    # Little t counts the number of time steps. It increases recursively.
    if t == n: # Expiration
        if isCall == True:
            return np.maximum(S - K, 0), delta, gamma, vega, theta, rho
        else:
             return np.maximum(K - S, 0), delta, gamma, vega, theta, rho
    elif t!=n:
        S_up = S * u
        S_down = S * d
        t += 1
        V_up = bm(S_up, K, r, q, sigma, T, n, t)[0]
        V_down = bm(S_down, K, r, q, sigma, T, n, t)[0]
        V_back = np.exp(-(r - q) * h) * (p * V_up + (1 - p) * V_down)
        # V_back goes back one step until it's at the present

        # Early value
        if isEuropean == True:
            return V_back, delta, gamma, vega, theta, rho
        else:
            early = 0
            if isCall == True:
                early = np.maximum(S - K, 0)
            else:
                early = np.maximum(K - S, 0)

            return np.maximum(V_back, early), delta, gamma, vega, theta, rho
    else:
        return ValueError

print(bm(100, 105, 0.05, 0, 0.4, 0.5, 4))