
import numpy as np
from numba import njit

@njit
def msoft_threshold( delta, lam, denom):
    
    if delta > lam:
        return (delta - lam) / denom
    
    elif delta < -lam:
        return (delta + lam) / denom
    
    else:
        return 0

@njit
def epoch_lasso_v2(X, beta, lam, r, c1, n, chosen_ps, s_new, s_diff):
    
    # O(9*np) = O(np)
    for j in chosen_ps:
    
        b_old = beta[j]

        delta = c1 * ( np.dot(X[:,j], r) + n*b_old) # O(n)
        denom = c1 * n
        beta[j] = msoft_threshold(delta, lam, denom)

        diff_bj = b_old - beta[j]
        if diff_bj != 0.0:    
            Xj_diff_bj = X[:,j] * diff_bj # O(2n)

            r      += Xj_diff_bj # O(2n)
            s_diff += Xj_diff_bj # O(2n)

        # rank 1 sums 
        if beta[j] != 0.0:
            s_new += X[:,j] * beta[j] # O(2n)

@njit
def update_beta_lasso_v2(X, beta, lam, r, c1, n, chosen_ps):
    
    # O(3*np) = O(np)
    for j in chosen_ps:
        b_old = beta[j]

        delta = c1 * ( np.dot(X[:,j], r) + n*b_old ) # O(n)
        denom = c1 * n
        beta[j] = msoft_threshold(delta, lam, denom)

        diff_bj = b_old - beta[j]
        if diff_bj != 0.0:
            r  += X[:,j] * diff_bj # O(3n)

# @jit
def calculate_dual_gap(R, X, y, lam, beta):
    n = len(y)

    p_obj = (1/n) * np.linalg.norm(R)**2  + lam * np.linalg.norm(beta, ord=1)

    theta = (2/n) * R
    norm_inf = np.linalg.norm(X.T @ theta, ord=np.inf)
    theta_ = theta/max(lam, norm_inf)

    d_obj = (1/n)*( np.linalg.norm(y)**2 - np.linalg.norm(y - (lam*n/2)*theta_)**2 )

    w_d_gap = (p_obj - d_obj)*n/np.linalg.norm(y)**2

    return w_d_gap


def dualpa(X, y, R, lam, beta, ny2):

    n = len(y)
    c1 = 1/n

    rt = 2*c1*R
    
    # O(np)
    norm_inf = np.linalg.norm(X.T @ rt, ord=np.inf)
    c2 = lam/max(norm_inf, lam)
    
    # O(n + p)
    d_gap = (1/n)*np.dot(R, (1 + c2**2)*R - 2*c2*y) + lam*np.linalg.norm(beta, ord=1)

    return d_gap*n/ny2

