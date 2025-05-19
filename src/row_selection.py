from collections import deque

import numpy as np

# path = np.random.randint(0, 2, size=(5, 3))
# non_zero = np.where(path[:,1] != 0)[0]
# sel_rows = list(non_zero + 1) 
# ",".join([str(i) for i in sel_rows])


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

def choose_j(path, test_errors = None, factor = 1/2):
    """
    Choose the best j based on the path and test errors
    Parameters
    ----------
    path : numpy.ndarray
        The path of coefficients with shape (p,k)
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

    Returns
    -------
    int
        The index of the best j
    
    """

    # path has (p,k) shape, where p is the number of features
    # and k is the number of lambda values
    if factor == -1 and test_errors is not None:
        # tests_errors contains two columns
        # the first one is the lambda values
        # the second one is the RMSE values
        # check 'calculate_test_errors' function
        
        # O(np) for obtain test_errors
        return np.argmin(test_errors[:,1]) # O(k) = O(1) for fixed k
    
    else:
        if factor < 0 or factor > 1:
            raise ValueError('Factor must be between 0 and 1 if nwerror is false and factor is not -1.')

        # recall p is the number of features
        p,k = path.shape
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

def row_selection(path, CT_spps, test_errors, n_spps = 15, 
                factor = 1/2, inbetween = 0, 
                check_spps = False):
    
    # path has (p,k) shape, where p is the number of features
    j_opt = choose_j(path, test_errors, factor = factor)
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

        if j_int == 0:
            taken.add(j_int)
            continue
    
        beta_j_nz = get_nonzero(path, j_int) # O(1) for isle

        if check_spps:
            # check on the number of species
            if len(np.unique(CT_spps[beta_j_nz,:])) < n_spps:
                taken.append(j_int)
                continue

        # plus one as Julia starts from 1
        # O(\rho T^4) for isle
        new_batches.append(  add_offset(beta_j_nz) )
        taken.append(j_int)
    
    return list(new_batches)

def write_rows(outfile, rows):

    with open(outfile, 'w') as f:
        for row in rows:
            # convert the row to a string
            f.write(",".join([str(b) for b in row]) + "\n")
