from collections import deque

import numpy as np

# path = np.random.randint(0, 2, size=(5, 3))
# non_zero = np.where(path[:,1] != 0)[0]
# sel_rows = list(non_zero + 1) 
# ",".join([str(i) for i in sel_rows])

def conv_errs(test_errors, tol = 1e-4):

    opt_j = np.inf
    for j in range(1, len(test_errors)):
        if np.abs(test_errors[j] - test_errors[j-1])/test_errors[j] < tol:
            opt_j = j
            break

    return opt_j

def constraint_test_error(test_errors, path):
    start_j = 0
    for j in range(len(test_errors)):
        if j == 0:
            # first j is always zero. If not, it is
            # floating point error
            continue
        sel = get_modelselection(path, j)
        if sel:
            start_j = j
            break

    return test_errors[start_j:], start_j


def data_driven_lambda_k(test_errors, path, verbose = False, tol = 1e-4, ignore_conv = False):
    """
    Choose the best j based on the test errors

    it  constraints the test errors to the first non-zero selection
    imagine path has the following set of non-zero selections:
    idx: 0    1   2  |  3    4  ...
        [ 0,  0,  0, | x1,  x2, ...] -> path_nz
        [e1, e2, e3, | e4, e_m, ...] -> test_errors

    then the start_j = 3. The opt_j =  start_j + opt_j_c.

    if the opt_j_c is at e_m, then opt_j_c = 1 and opt_j = 3 + 1 = 4


    Parameters
    ----------
    test_errors : numpy.ndarray
        The test errors with shape (k,2)
        where the first column is the lambda values and the second column is the RMSE values
    path : numpy.ndarray
        The path of coefficients with shape (p,k) or (M,k) if it is from ISLE
        where p is the number of features and k is the number of lambda values
    verbose : bool, optional
        If True, the function will print the chosen j
        The default is False.
    tol : float, optional
        The tolerance for the convergence of the test errors
        The default is 1e-4.
    ignore_conv : bool, optional
        If True, the function will ignore the convergence of the test errors
        The default is False.
    Returns
    -------
    int
        The index of the best j
        This is the index of the lambda value that minimizes the test error
        or the index of the lambda value that converges to the test error
        depending on the convergence of the test errors.
        If the convergence is ignored, then it returns the index of the minimum test error.
    """
    # make sure the opt_j is not an empty selection
    test_errors_c, start_j = constraint_test_error(test_errors, path)

    min_err_j = np.argmin(test_errors_c) #O(k) = O(1) for fixed k
    conv_err_j = conv_errs(test_errors_c, tol=tol) #O(k) = O(1) for fixed k

    if conv_err_j < min_err_j and not ignore_conv:
        if verbose:
            print("lambda_k choosen by test error convergence", end= "")

        opt_j_c = conv_err_j
    else:
        if verbose:
            print("lambda_k choosen by minimum test error", end= "")
        opt_j_c = min_err_j

    return opt_j_c + start_j

def get_modelselection(path, j):
    """
    Get the number of non-zero coefficients for the j-th lambda value
    """
    if isinstance(path, np.ndarray):
        # beta_j is the j-th column of path
        # which is obtained from the j-th lambda value
        beta_j = path[:,j]
        # selection by elastic net
        return np.sum(beta_j != 0)
    
    else:
        # it is comming from the ISLE new_path
        # which is a list of splitting variables
        # indeces
        return len(path[j])

def choose_j(path, test_errors = None, factor = -1,
              verbose = False, tol = 1e-4, p = None,
              ignore_conv = False):
    """
    Choose the best j based on the path and test errors
    Parameters
    ----------
    path : numpy.ndarray
        The path of coefficients with shape (p,k) or (M,k) if it is from ISLE
        where p is the number of features and k is the number of lambda values
    test_errors : numpy.ndarray, optional
        The test errors with shape (k,2)
        where the first column is the lambda values and the second column is the RMSE values
    factor : float, optional
        The factor to choose the best j
        if factor is -1, then the function will return the index of the minimum test error
        if factor is between 0 and 1, then the function will return the index of the best j
        based on the number of non-zero coefficients
        The default is 1/2.
    verbose : bool, optional
        If True, the function will print the chosen j
        The default is False.
    tol : float, optional
        The tolerance for the convergence of the test errors
        The default is 1e-4.

    Returns
    -------
    int
        The index of the best j
    
    """

    # path has (p,k) shape, where p is the number of features
    # and k is the number of lambda values
    if factor == -1:
        # tests_errors contains two columns
        # the first one is the lambda values
        # the second one is the RMSE values
        # check 'calculate_test_errors' function
        
        # O(np) for obtain test_errors
        return data_driven_lambda_k(test_errors[:,1], path, verbose, tol, ignore_conv) # O(k) = O(1) for fixed k
        # return np.argmin(test_errors[:,1]) # O(k) = O(1) for fixed k
    
    else:
        if factor < 0 or factor > 1:
            raise ValueError('Factor must be between 0 and 1 if factor is not -1.')            
        
        if verbose:
            print(f"lambda_k choosen in function of factor ({factor}).")

        # recall p is the number of features
        k = path.shape[1] if isinstance(path, np.ndarray) else len(path)
        user_selection = np.round(p*factor).astype(int)            

        best_dist = np.inf
        best_j = 0
        for j in range(k):

            model_selection = get_modelselection(path, j)
            # distance between of desired number of non-zero
            # coefficients and the current number of non-zero
            tmp_dist = np.abs(model_selection - user_selection)
            if tmp_dist < best_dist:
                best_dist = tmp_dist
                best_j = j

        return best_j

def get_nonzero(path, j):

    if isinstance(path, np.ndarray):
        # stills returns an np.array
        return np.where(path[:,j] != 0)[0] # O(T^4)
    
    else:
        return path[j] # O(1)

def add_offset(beta_j_nz):
    """
    Add offset to the path
    and if beta_j_nz is an array
    then convert it to a list
    """
    if isinstance(beta_j_nz, np.ndarray):
        return list(beta_j_nz + 1)
    
    else:
        return [i + 1 for i in beta_j_nz]

def read_CT(CT_file):
    CT = np.loadtxt(CT_file, delimiter=',', skiprows=1, dtype=str)
    CT_spps = CT[:, :4]
    n_spps = len(np.unique(CT_spps))
    return CT_spps, n_spps

def row_selection(path, test_errors,  p = None,
                  factor = 1/2, inbetween = 0, 
                  check_spps = False, CT_file = None,
                  verbose = False, tol = 1e-4, 
                  ignore_conv = False):
    
    if check_spps:
        CT_spps, n_spps = read_CT(CT_file) # O(T^4)
        

    # path has (p,k) shape, where p is the number of features
    j_opt = choose_j(path, test_errors, factor = factor, verbose=verbose, tol=tol, p = p, ignore_conv = ignore_conv) # O(k) = O(1) for fixed k
    chosen_j = np.linspace(0, j_opt, 2 + inbetween,
                            endpoint=True, dtype=int) 

    taken = set()
    new_batches = deque()
    for j_int in chosen_j:

        # once it is integer,
        # it might be the case
        # that there are repeated j's
        if j_int in taken:
            continue
        
        # in any of the paths: elastic net or ISLE
        # the first column is all zeros
        if j_int == 0:
            taken.add(j_int)
            continue
    
        beta_j_nz = get_nonzero(path, j_int) # O(1) for isle

        if check_spps:
            # check on the number of species
            if len(np.unique(CT_spps[beta_j_nz,:])) < n_spps:
                taken.add(j_int)
                continue

        # plus one as Julia starts from 1
        # O(\rho T^4) for isle
        new_batches.append(  add_offset(beta_j_nz) )
        taken.add(j_int)
    
    return list(new_batches), j_opt

def write_rows(outfile, rows):

    with open(outfile, 'w') as f:
        for row in rows:
            # convert the row to a string
            f.write(",".join([str(b) for b in row]) + "\n")
