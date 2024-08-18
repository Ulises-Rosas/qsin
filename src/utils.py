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
