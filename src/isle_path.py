import time
import json
from copy import deepcopy
from collections import deque

from qsin.sparse_solutions import split_data
from qsin.utils import progressbar, standardize_Xy

import numpy as np
from sklearn.tree import DecisionTreeRegressor

def re_center_for_isle(T_test, T_train):
    """
    Center data for ISLE.  When using ISLE,
    the X matrix is a set of predictions from
    decision trees, which are not centered. 
    Since ISLE assumes there is an intercept term,
    and the lasso/elnet post-processing assumes the data
    is centered, we need to center the data.

    Parameters
    ----------
    T_test : numpy.ndarray
        The predictors for the test set
    
    T_train : numpy.ndarray
        The predictors for the training set

    Returns
    -------
    numpy.ndarray
        The rescaled predictors for the test set and the training set
    """

    u = np.mean(T_train, axis=0)
    sd = np.std(T_train, axis=0)
    
    sd_zero = sd == 0
    if np.any(sd_zero):
        sd[sd_zero] = 1

    return (T_test - u)/sd, (T_train - u)/sd


def Sm(X, y, f_m, sample_size, replace = False, rng = None):
    # f_m = f_0

    n,p = X.shape
    test_idx  = rng.choice(range(n), size = sample_size, replace = replace)
    # print("test_idx: ", test_idx)

    return X[test_idx,:], y[test_idx], f_m[test_idx]

def make_isle_ensemble(X_train, y_train, model, eta, nu,
                      M, rng = None, verbose = True):
    

    n_train = X_train.shape[0]
    train_sample_size = int(n_train*eta)
    if verbose:
        print("Starting ISLE ensemble")
        print("Random sample size for trees: ", train_sample_size)

    # initialize memory function
    f_m = np.repeat(np.mean(y_train), n_train)

    F_train = np.zeros((n_train,M))
    estimators = deque()
    # O(MnT^4log(n))
    for i in progressbar(range(M), "Computing trees: ", 40):

        model.set_params(random_state = rng)
        # random sample the data, including the memory function
        # O(nT^4)
        X_sm, y_sm, f_sm = Sm(X_train, y_train, f_m, train_sample_size, replace=False, rng=rng)
        # fit the model. O(p n log(n)), 
        # where p is considered number of features.
        # p <= T^4 -> O(p n log(n)) <= O(T^4 n log(n))
        model.fit(X_sm + f_sm.reshape(-1,1), y_sm )
        
        # update memory function
        f_m = f_m + nu*model.predict(X_train)
        
        F_train[:,i] = model.predict(X_train)  # O(n)
        # model largest object is the tree
        # and this has a size of O(d) = O(1) for
        # small d
        estimators.append(deepcopy(model)) # O(1)

    return F_train, list(estimators) # O(M)

def make_init_model(max_features = None, max_depth = 5, max_leaves = 6, param_file = None):

    if param_file is None:
        return DecisionTreeRegressor(max_features = max_features, 
                                     max_depth = max_depth, 
                                     max_leaf_nodes = max_leaves)
    
    else:
        # read json file
        # param_file = './tree_params.txt'
        with open(param_file) as f:
            params = json.load(f)

        return DecisionTreeRegressor(max_features = max_features,
                                     max_depth = max_depth,
                                     max_leaf_nodes = max_leaves
                                     **params)

def make_F_test(X_test, estimators):
    M = len(estimators)
    F_test = np.zeros((X_test.shape[0], M))
    
    # O(nM)
    for i,m in enumerate(estimators):
        F_test[:,i] = m.predict(X_test) # O(n)
    return F_test

# def split_data_isle(X, y, num_test, seed, 
#                     isle=True, 
#                     mx_p=1/2, max_depth=5, param_file=None, max_leaves=6,
#                     eta=0.5, nu=0.1, M=100,
#                     verbose=True, nstdy = True, args=None):

def split_data_isle(X, y, num_test, seed, 
                    isle=True, max_leaves=6, eta=0.5, nu=0.1,
                    nstdy = True, args=None):
    """
    Split data into training and testing sets, and apply ISLE if needed.

    Parameters
    ----------
    X : array-like
        Feature matrix.
    y : array-like
        Target vector.
    num_test : int
        Number of samples to use for testing.
    seed : int
        Random seed.
    isle : bool
        Whether to apply ISLE.
    mx_p : int or str
        Maximum number of features to use in each tree.'
        This follows the sklearn convention.
        If 'sqrt', then max_features = sqrt(p).
        If 'log2', then max_features = log2(p).
        If None, then max_features = p.
    max_depth : int
        Maximum depth of the decision tree.
    max_leaves : int
        Maximum number of leaves in the decision tree.
    param_file : str
        JSON file with parameters for the decision tree.
    eta : float
        Proportion of samples to use in each tree.
    nu : float
        Learning rate.
    M : int
        Number of trees in the ensemble.
    verbose : bool
        Whether to print progress.
    nstdy : bool
        Whether to standardize the target vector.

    Returns
    -------
    X_train : array-like
        Training feature matrix.
    X_test : array-like
        Testing feature matrix.
    y_train : array-like
        Training target vector.
    y_test : array-like
        Testing target vector.
    """
    rng = np.random.RandomState(seed)

    X_train, X_test, y_train, y_test = split_data(X, y, num_test=num_test, seed=seed)

    # improve stability of the algorithm
    X_train, X_test, y_train, y_test = standardize_Xy(X_train, y_train, X_test, y_test, nstdy)


    param_size = len(eta) * len(nu) * len(max_leaves)
    if isle and param_size == 1:
        # full certainity on the hyperaparameters used.
        # no need to do further testing via cross-validation
        eta = eta[0] if isinstance(eta, list) else eta
        nu = nu[0] if isinstance(nu, list) else nu
        max_leaves = max_leaves[0] if isinstance(max_leaves, list) else max_leaves

        # elastice net is applied directly on here
        (X_train, X_test, 
         estimators) = isle_ensemble_pipe(X_train, y_train, 
                                          X_test, max_leaves, 
                                          eta, nu, rng, args)
    else:
        estimators = None

    return X_train, X_test, y_train, y_test, estimators
    

def isle_ensemble_pipe(X_train, y_train, X_test, max_leaves, eta, nu, rng, args):
        """
        generates F_train and F_test and center them
        """
        # initialize model without random state
        # this is later set in the ensemble loop
        T_0 = make_init_model(max_leaves=max_leaves,
                              max_features=args.max_features, 
                              max_depth=args.max_depth,  
                              param_file=args.param_file)

        start = time.time()
        F_train, estimators = make_isle_ensemble(X_train, y_train, T_0, eta, nu, args.M, rng=rng)
        end_isle = time.time() - start

        if args.verbose:
            print("Isle ensemble done: ", end_isle, " seconds")

        F_test = make_F_test(X_test, estimators)

        # recenter
        if args.verbose:
            print("Re-standarize data for ISLE")
        F_test, F_train = re_center_for_isle(F_test, F_train)

        return F_train, F_test, estimators


def get_new_path(estimators, path):
    """
    this path is based on the feature importances that
    the selected estimators have. For each lambda, there 
    is an ensemble of estimators and the new_path contains
    the average feature importances of the ensemble

    The new path is a p x K matrix instead of M x K

    Parameters
    ----------
    estimators : list
        List of decision tree regressors.
    path : numpy.ndarray
        The path of coefficients with shape (M,k)
        where M is the number of estimators and k is the number of lambda values

    Returns
    -------
    list
        The new path, which is a list of lists of feature indices.
        Each set corresponds to a lambda value.
        The length of the list is equal to the number of lambda values.
        Each set contains the indices of the features selected by the ensemble
        of estimators for that lambda value.
    """

    estimators = np.array(estimators)
    new_path = deque()

    for j in range(path.shape[1]):
        # first iteration of the path
        # all the coefficients are 0: no selection
        if j == 0:
            # we add an empty set for backward compatibility
            # with the elastic net path
            new_path.append([])
            continue

        # j = 2
        coeffs = path[:,j]
        coeffs_logic = coeffs != 0

        tmp_ensemble = estimators[coeffs_logic]

        I_k = set()
        for m in tmp_ensemble:
            I_k_m = set(m.tree_.feature[m.tree_.feature != -2])
            I_k |= I_k_m

        I_k = list(I_k)
        new_path.append(I_k)

    return list(new_path)


def create_param_grid(eta, nu, leaves, cv_sample = 100, rng = None):
    """
    return an array with combinations of hyperparameters
    for the ISLE ensemble. If cv_sample is lower than all possible
    combinations of hyperparameters, a random sample of size cv_sample
    is taken
    """
    # eta = args.eta
    # nu = args.nu
    # leaves = args.max_leaf_nodes
    # cv_sample = args.cv_sample

    param_size = len(eta) * len(nu) * len(leaves)
    
    if param_size <= cv_sample:
        grid = deque()
        for nj in eta:
            for vj in nu:
                for lj in leaves:
                    grid.append((nj, vj, lj))           
        grid = np.array(grid)
    else:
        grid = np.zeros((cv_sample, 3))
        grid[:,0] = rng.choice(eta, cv_sample)
        grid[:,1] = rng.choice(nu, cv_sample)
        # needs to be changed into integers
        grid[:,2] = rng.choice(leaves, cv_sample)

    return grid

# from argparse import Namespace
# args = Namespace()

# args.alpha = [1,2,3]
# args.eta = [0.1,0.2,0.3]
# args.nu = [0.4,0.5,0.6]
# args.max_leaf_nodes = [2,3,4,5]
# args.cv_sample = 100 # currently only randomize ISLE ensemble params only
# # alpha of elastice net are untouched
# rng = np.random.RandomState(12038)



# a = np.zeros((2,2))
# def mysterious_fun(a):
#     a[0,0] = 100
#     return 1
# b = mysterious_fun(a)