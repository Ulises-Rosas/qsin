import time
import json
import numpy as np
from copy import deepcopy


from qsin.sparse_solutions import split_data
from sklearn.tree import DecisionTreeRegressor

def Sm(X, y, f_m, sample_size, replace = False, seed = 12038):
    np.random.seed(seed)
    # f_m = f_0

    n,p = X.shape
    test_idx  = np.random.choice(range(n), size = sample_size, replace = replace)
    return X[test_idx,:], y[test_idx], f_m[test_idx]

def make_isle_ensemble(X_train, y_train, model, eta, nu, M, seed = 12038, verbose = True):
    np.random.seed(seed)


    n_train = X_train.shape[0]
    # n_test = X_test.shape[0]
    f_m = np.repeat(np.mean(y_train), n_train)

    F_train = np.zeros((n_train,M))
    # F_test = np.zeros((n_test,M))

    rs = [np.random.choice(range(seed)) for _ in range(M)]
    estimators = []

    train_sample_size = int(n_train*eta)

    if verbose:
        print("Starting ISLE ensemble")
        print("Random sample size for trees: ", train_sample_size)

    for i in range(M):
        model.set_params(random_state =  rs[i])
        X_sm, y_sm, f_sm = Sm(X_train, y_train, f_m, train_sample_size, replace=False)


        model.fit(X_sm + f_sm.reshape(-1,1), y_sm )
        f_m = f_m + nu*model.predict(X_train)
        
        # X_sm.shape
        F_train[:,i] = model.predict(X_train)
        # F_test[:,i] = model.predict(X_test)
        estimators.append(deepcopy(model))

    return F_train, estimators

def make_init_model(max_features = None, max_depth = 5, param_file = None):

    if param_file is None:
        return DecisionTreeRegressor(max_features = max_features, 
                                     max_depth = max_depth,)
    
    else:
        # read json file
        # param_file = './tree_params.txt'
        with open(param_file) as f:
            params = json.load(f)

        return DecisionTreeRegressor(max_features = max_features,
                                     max_depth = max_depth,
                                     **params)

def make_F_test(X_test, estimators):
    M = len(estimators)
    F_test = np.zeros((X_test.shape[0], M))
    for i,m in enumerate(estimators):
        F_test[:,i] = m.predict(X_test)
    return F_test


def split_data_isle(X, y, num_test, seed, 
                    isle=True, nwerror=False, 
                    mx_p=1/2, max_depth=5, param_file=None, 
                    eta=0.5, nu=0.1, M=100,
                    verbose=True):
    
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
    nwerror : bool
        Whether to use the entire dataset for training.
    mx_p : float
        Maximum proportion of features to use in each tree.
    max_depth : int
        Maximum depth of the decision tree.
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

    if nwerror:
        X_train, y_train = X, y
        X_test, y_test = None, None
    else:
        X_train, X_test, y_train, y_test = split_data(X, y, num_test=num_test, seed=seed)

    if isle:
        n, p = X_train.shape
        model = make_init_model(max_features=round(p * mx_p), max_depth=max_depth, param_file=param_file)

        start = time.time()
        F_train, estimators = make_isle_ensemble(X_train, y_train, model, eta, nu, M, seed=seed)
        end_isle = time.time() - start

        if verbose:
            print("Isle ensemble done: ", end_isle, " seconds")

        X_train = F_train

        if not nwerror:
            F_test = make_F_test(X_test, estimators)
            X_test = F_test

    else:
        estimators = None

    return X_train, X_test, y_train, y_test, estimators
    

def get_new_path(estimators, path,p):
    """
    this path is based on the feature importances that
    the selected estimators have. For each lambda, there 
    is an ensemble of estimators and the new_path contains
    the average feature importances of the ensemble

    The new path is a p x K matrix instead of M x K
    """

    estimators = np.array(estimators)
    new_path = np.zeros((p, path.shape[1]))

    for j in range(path.shape[1]):

        if j == 0:
            new_path[:,j] = np.repeat(0, p)
            continue
        # j = 2
        coeffs = path[:,j]
        coeffs_logic = coeffs != 0

        # tmp_ensemble_w = coeffs[coeffs_logic]
        tmp_ensemble = estimators[coeffs_logic]
        tmp_ensemble_fi = np.zeros(( p, len(tmp_ensemble)))

        for i,m in enumerate(tmp_ensemble):
            # tmp_ensemble_fi[:,i] = m.feature_importances_*tmp_ensemble_w[i]
            tmp_ensemble_fi[:,i] = m.feature_importances_

        # inf norm avoids numerical instabilities associated
        # with the sum of the feature importances or mean
        new_path[:,j] = np.linalg.norm(tmp_ensemble_fi, axis = 1, ord = np.inf)

    return new_path



# from comparing_solutions import _scaler, ElasticNet, lasso_path, max_lambda, rmse
# file = './test_sims/test_n15_qll.csv'
# CT_file = './1_seqgen.CFs_n15.csv'

# data = np.loadtxt(file, delimiter=',',skiprows=1, )
# X,y = data[:,:-1], data[:,-1]
# n,p = X.shape



# nwerror = False
# isle = True

# # elastic net parameters
# # alpha = 0.999
# alpha = 0.999
# e =  0.001
# K = 100

# # decision tree parameters
# mx_p = 1/2
# max_depth = 5
# param_file = None

# # isle parameters
# # eta = 0.5
# eta = 0.1
# nu = 0.1
# M = 500

# seed = 12038
# verbose = True



# X = _scaler(X, X, sted = True)
# y = _scaler(y, y, sted = True)


# num_test = int(n*0.35)

# (X_train,X_test,
#  y_train,y_test,
#  estimators) = split_data_isle(X, y, 
#                     num_test=num_test, seed=seed,
#                     isle=True, nwerror=False, 
#                     mx_p=mx_p, max_depth=max_depth, param_file=None, 
#                     eta=eta, nu=nu, M=M,
#                     verbose=True)


# max_lam = max_lambda(X_train, y_train, alpha=alpha)
# min_lam = max_lam * e
# params = {'lam': np.logspace(np.log10(min_lam), np.log10(max_lam), K, endpoint=True)[::-1]}

# model = ElasticNet(fit_intercept=False, 
#                     max_iter=1000,
#                     init_iter=1, 
#                     copyX=True, 
#                     alpha=alpha, 
#                     tol=0.00001)

# start = time.time()
# path = lasso_path(X_train, y_train, params, model, sequential_screening=False)
# end_lasso = time.time() - start
# print("Elastic net path done: ", end_lasso, " seconds")


# lam_path = np.concatenate((params['lam'].reshape(-1, 1), path.T), axis=1)
# test_errors = np.zeros((path.shape[1], 2))
# for j in range(path.shape[1]):
#     beta_j = path[:, j]
#     rmse_j = rmse( y_test, X_test, beta_j)
#     test_errors[j, :] = [ params['lam'][j] , rmse_j]




# from matplotlib import pyplot as plt
# fs = 18
# fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

# # Upper panel - Coefficient plot

# ax1.plot(lam_path[:,0], lam_path[:,1:], marker='o', alpha=0.8)
# ax1.set_xscale('log')
# ax1.set_ylabel('Coefficient', fontsize=fs)
# ax1.axhline(0, color='black', lw=2)
# ax1.set_title(f'Lasso Path ($\\alpha = {alpha}$)', fontsize=fs)


# ax2.plot(test_errors[:,0], test_errors[:,1], marker='o', alpha=0.8, label='Weighted average error of selected trees')
# ax2.set_xscale('log')
# ax2.set_yscale('log')
# ax2.set_xlabel('$\lambda$', fontsize=fs)
# ax2.set_ylabel('RMSE', fontsize=fs)
# ax2.legend()



# new_path = get_new_path(estimators, path, p)
# j_min = np.argmin(test_errors[:,0])
# new_path_j = new_path[:,j_min]
# np.sum(new_path_j != 0)


