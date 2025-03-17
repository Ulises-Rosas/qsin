import sys
import time

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


def _scaler(ref, dat, sted = True):
    u = np.mean(ref, axis=0)
    
    cred = dat - u

    if sted:
        sd = np.std(ref, axis=0)
        cred = cred / sd

    return cred

def rmse(y, X, B):
    return np.sqrt( np.mean( ( y - X.dot(B) )**2 ) )


def calculate_test_errors(args, path, params, X_test, y_test, write = True):
    """
    Calculate the test errors for each lambda value in the path

    Parameters
    ----------
    args : argparse.Namespace
        Arguments from the command line

    path : numpy.ndarray
        The path of coefficients

    params : dict
        Parameters for the path

    X_test : numpy.ndarray
        The test data

    y_test : numpy.ndarray
        The test response

    Returns
    -------
    numpy.ndarray
        The test errors, which is a 2D array with the first column
        indicating the lambda value and the second column indicating
        the RMSE value

    """
    if args.nwerror:
        test_errors = None
    else:
        test_errors = np.zeros((path.shape[1], 2))
        for j in range(path.shape[1]):
            beta_j = path[:, j]
            rmse_j = rmse(y_test, X_test, beta_j)
            test_errors[j, :] = [params['lam'][j], rmse_j]

        if write:
            np.savetxt(args.prefix + "_testErrors.csv",
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
