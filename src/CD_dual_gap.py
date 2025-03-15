
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

# region epoch functions
@njit
def epoch_lasso(X, beta, lam, c1, Xj_T_y, Xj_T_Xj, X_T_X, chosen_ps, s_new, diff):

    for j in chosen_ps:
        # j = 1
        b_old = beta[j]
        
        A2 = np.where(beta != 0)[0]
        # take out j from the array A2
        A2 = A2[A2 != j]
        tmp_X_T_X = X_T_X[j,:]

        delta = c1 * ( Xj_T_y[j] - np.dot(tmp_X_T_X[A2], beta[A2]) )
        denom = c1 * Xj_T_Xj[j]
        beta[j] = msoft_threshold(delta, lam, denom)


        # beautiful rank 1 sums 
        if beta[j] != 0.0:
            s_new += X[:,j] * beta[j]

        diff_j = b_old - beta[j]
        if diff_j != 0.0:    
            diff += X[:,j] * diff_j


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
def epoch_enet(X, beta, c1, Xj_T_y, Xj_T_Xj, X_T_X, chosen_ps, lam_1_alpha, lam_alpha, s_new, diff):

    for j in chosen_ps:
        # j = 1
        b_old = beta[j]
        
        A2 = np.where(beta != 0)[0]
        # take out j from the array A2
        A2 = A2[A2 != j]
        tmp_X_T_X = X_T_X[j,:]

        delta = c1 * ( Xj_T_y[j] - np.dot(tmp_X_T_X[A2], beta[A2]) )
        denom = (c1 * Xj_T_Xj[j]) + lam_1_alpha

        beta[j] = msoft_threshold(delta, lam_alpha, denom)


        # beautiful rank 1 sums 
        if beta[j] != 0.0:
            s_new += X[:,j] * beta[j]

        diff_j = b_old - beta[j]
        if diff_j != 0.0:    
            diff += X[:,j] * diff_j



@njit
def epoch_enet_v2(X, beta, lam_1_alpha, lam_alpha, r, c1, n, chosen_ps, s_new, s_diff):

    for j in chosen_ps:
        b_old = beta[j]
        
        delta = c1 * ( np.dot(X[:,j], r) + n*b_old )
        denom = (c1 * n) + lam_1_alpha
        beta[j] = msoft_threshold(delta, lam_alpha, denom)

        diff_bj = b_old - beta[j]
        if diff_bj != 0.0:    
            Xj_diff_bj = X[:,j] * diff_bj # O(2n)

            r      += Xj_diff_bj # O(2n)
            s_diff += Xj_diff_bj # O(2n)

        # rank 1 sums 
        if beta[j] != 0.0:
            s_new += X[:,j] * beta[j] # O(2n)

# endregion

# region update functions
@njit
def update_beta_lasso(beta, lam, c1, Xj_T_y, Xj_T_Xj, X_T_X, chosen_ps):

    for j in chosen_ps:
        # j = 1
        # A2 = np.where(beta !=  0)[0]
        # # take out j from the array A2
        # A2 = A2[A2 != j]

        A2 = np.where(beta !=  beta[j])[0]

        tmp_X_T_X = X_T_X[j,:]

        delta = c1 * ( Xj_T_y[j] - np.dot(tmp_X_T_X[A2], beta[A2]) )
        denom = c1 * Xj_T_Xj[j]
        beta[j] = msoft_threshold(delta, lam, denom)



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


@njit
def update_beta_enet(beta, c1, Xj_T_y, Xj_T_Xj, X_T_X, 
                     chosen_ps, lam_1_alpha, lam_alpha):
    for j in chosen_ps:
        # j = 1
        A2 = np.where(beta != 0)[0]
        # take out j from the array A2
        A2 = A2[A2 != j]

        tmp_X_T_X = X_T_X[j,:]

        delta = c1 * ( Xj_T_y[j] - np.dot(tmp_X_T_X[A2], beta[A2]) )
        denom = (c1 * Xj_T_Xj[j]) + lam_1_alpha

        beta[j] = msoft_threshold(delta, lam_alpha, denom)


@njit
def update_enet_v2(X, beta, lam_1_alpha, lam_alpha, r, c1, n, chosen_ps):

    for j in chosen_ps:
        b_old = beta[j]
        
        delta = c1 * ( np.dot(X[:,j], r) + n*b_old )
        denom = (c1 * n) + lam_1_alpha
        beta[j] = msoft_threshold(delta, lam_alpha, denom)

        diff_bj = b_old - beta[j]
        if diff_bj != 0.0:
            r  += X[:,j] * diff_bj # O(3n)

# endregion

# region dual gap functions
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


def dualpa_elnet(R, X, y, beta, lam_1_alpha, lam_alpha, ny2):
    
    if lam_1_alpha == 0:
        return dualpa(X, y, R, lam_alpha, beta, ny2)

    n = len(y)
    b_norm_sq = np.dot(beta, beta)

    rt = (2/n) * R
    norm_inf = np.linalg.norm(X.T @ rt  - lam_1_alpha*beta, ord=np.inf)

    c2 = lam_alpha/max(norm_inf, lam_alpha)
    c3 = 1 + c2**2

    d_gap = (1/n)*np.dot(R, c3*R - 2*c2*y) +\
            lam_alpha*np.linalg.norm(beta, ord=1) +\
            c3*(lam_1_alpha/2)*b_norm_sq

    return d_gap*n/ny2

# endregion