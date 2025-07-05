import sys
import time
import random

import numpy as np

def max_lambda(X, y, alpha = None):
    # 1/n * max( |X^T * y| )
    # SLS book page 114
    n,_ = X.shape

    if alpha is not None:
        c1 = 2 / (alpha * n)
    else:
        c1 = 2 / n

    return c1 * np.max( np.abs( np.matmul(X.T, y) ) )


def get_lambdas(alpha, X_train, y_train, K, epsilon, verbose=False):

    # max_lambda = O(nT^4) operations, under the budget of
    # simulations, or ISLE ensemble, it is negligible
    max_lam = max_lambda(X_train, y_train, alpha=alpha)
    min_lam = max_lam * epsilon

    if verbose:
        print("Lambda range: ", min_lam, max_lam)

    lambdas =  np.logspace(np.log10(min_lam), np.log10(max_lam), K, endpoint=True)[::-1]
    return lambdas

def get_alpha_max_lam(X_train, y_train, alphas):
    
    all_max_lams = -np.inf
    all_max_lams_alpha_indx = -1
    for i,a in enumerate(alphas):
        tmp_a = max_lambda(X_train, y_train, alpha=a)
        if tmp_a > all_max_lams:
            all_max_lams = tmp_a
            all_max_lams_alpha_indx = i
            
    return alphas[all_max_lams_alpha_indx]


def rmse(y, X, B):
    return np.sqrt( np.mean( ( y - X.dot(B) )**2 ) )

def split_data(X,y,num_test, seed = 123):

    random.seed(seed)
    n,_ = X.shape

    test_idx  = random.sample(range(n), k = num_test)
    train_idx = list( set(range(n)) - set(test_idx) )

    X_train, X_test = X[train_idx,:], X[test_idx,:]
    y_train, y_test = y[train_idx]  , y[test_idx]

    return X_train, X_test, y_train, y_test

def write_test_errors(prefix, errors, lambdas):
    """
    Calculate the test errors for each lambda value in the path

    Parameters
    ----------
    prefix : str
        The prefix for the output file

    errors : list
        The list of RMSE values for each lambda

    lambdas : list
        The list of lambda values

    Returns
    -------
    numpy.ndarray
        The test errors, which is a 2D array with the first column
        indicating the lambda value and the second column indicating
        the RMSE value

    """
    test_errors = np.zeros((len(errors), 2))
    for j, rmse_j in enumerate(errors):
        test_errors[j, :] = [lambdas[j], rmse_j]

    np.savetxt(prefix + "_testErrors.csv",
                test_errors,
                delimiter=',',
                comments='')

    return test_errors

def progressbar(it, prefix="", size=60, out=sys.stdout): 
    """
    ref: https://stackoverflow.com/a/34482761
    Progress bar for loops

    """ 
    count = len(it)
    start = time.time() # time estimate start
    def show(j):
        x = int(size*j/count)
        # time estimate calculation and string
        remaining = ((time.time() - start) / j) * (count - j)        
        mins, sec = divmod(remaining, 60) # limited to minutes
        time_str = f"{int(mins):02}:{sec:03.1f}"
        print(f"{prefix}[{u'█'*x}{('.'*(size-x))}] {j}/{count} Est wait {time_str}", end='\r', file=out, flush=True)
    show(0.1) # avoid div/0 
    for i, item in enumerate(it):
        yield item
        show(i+1)
    print("\n", flush=True, file=out)


def standardize_Xy(X_train, y_train, X_test = None, y_test = None, nstdy = False):
    """
    Standardize the data
    """
    x_u = np.mean(X_train, axis = 0)
    x_sd = np.std(X_train, axis = 0)
    
    zero_sd_indices = x_sd == 0 # O(p)
    if np.any(zero_sd_indices): # O(p)
        x_sd[zero_sd_indices] = 1 # O(p)

    y_u = np.mean(y_train)
    y_sd = np.std(y_train)

    if y_sd == 0:
        y_sd = 1

    X_train = (X_train - x_u) / x_sd
    y_train = y_train - y_u if nstdy else (y_train - y_u) / y_sd

    if X_test is not None:
        X_test = (X_test - x_u) / x_sd
        y_test = y_test - y_u if nstdy else (y_test - y_u) / y_sd

    return X_train, X_test, y_train, y_test

def standardize_X(X_train,  X_test = None):
    """
    Standardize the data
    """
    x_u = np.mean(X_train, axis = 0)
    x_sd = np.std(X_train, axis = 0)
    
    x_sd_zero = x_sd == 0
    if np.any(x_sd_zero):
        x_sd[x_sd_zero] = 1

    X_train = (X_train - x_u) / x_sd

    if X_test is not None:
        X_test = (X_test - x_u) / x_sd

    return X_train, X_test
