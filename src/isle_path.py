import time
import json
from copy import deepcopy
from collections import deque

from qsin.sparse_solutions_hd import ElasticNet, lasso_path
from qsin.utils import progressbar, get_lambdas

import numpy as np
from sklearn.tree import DecisionTreeRegressor

def Sm(X, y, f_m, sample_size, replace = False, rng = None):
    # f_m = f_0

    n,p = X.shape
    test_idx  = rng.choice(range(n), size = sample_size, replace = replace)
    # print("test_idx: ", test_idx)

    return X[test_idx,:], y[test_idx], f_m[test_idx]

def make_isle_ensemble(X_train, y_train, model, eta, nu,
                      M, rng = None, verbose = True):
    

    n_train = X_train.shape[0]
    train_sample_size = np.clip( int(n_train*eta), 1, n_train )
    if verbose:
        print("Starting ISLE ensemble")
        print("Random sample size for trees: ", train_sample_size)

    # initialize memory function
    f_m = np.repeat(np.mean(y_train), n_train)

    F_train = np.zeros((n_train,M))
    estimators = deque()

    pb = progressbar(range(M), "Computing trees: ", 40) if verbose else range(M)
    # O(MnT^4log(n))
    for i in pb:
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
        coeffs_logic = coeffs != 0.0

        tmp_ensemble = estimators[coeffs_logic]

        I_k = set()
        for m in tmp_ensemble:
            I_k_m = set(m.tree_.feature[m.tree_.feature != -2])
            I_k |= I_k_m

        I_k = list(I_k)
        new_path.append(I_k)

    return list(new_path)

class ISLEPath:
    """
    Path solution for the ISLE algorithm
    """
    def __init__(self,
                 # ISLE Ensemble params
                 eta = 0.5,
                 nu = 0.1,
                 max_leaves = 2,
                 M = 1000,
                 max_features = None,
                 max_depth = None,
                 param_file = None,
                 rng = None,
                 make_ensemble = True,
                 
                 # Elastic Net params
                 fit_intercept = True,
                 max_iter = 1000,
                 alpha = 0.5,
                 tol = 0.0001,
                 zero_thresh = 1e-15,  # threshold for zero coefficients
                 # Path params
                 epsilon  = 0.0001,
                 K = 100,
                 # other params
                 verbose = True
                 ):
        
        # ISLE ensemble params
        self.eta = eta
        self.nu = nu
        self.max_leaves = max_leaves
        self.M = M
        self.max_features = max_features
        self.max_depth = max_depth
        self.param_file = param_file
        self.rng = rng
        self.make_ensemble = make_ensemble

        # Elastic Net params
        self.fit_intercept = fit_intercept
        self.max_iter = max_iter
        self.alpha = alpha
        self.zero_thresh = zero_thresh # internal parameter for the elastic net

        self.tol = tol

        # Path params
        self.epsilon = epsilon
        self.K = K

        # other params
        self.verbose = verbose

        # internal variables
        self.estimators = None
        self.path = None
        self.intercepts = None

        self.mu = None
        self.sigma = None
        self.lambdas = []  # will be set in fit method (pruning)

    def set_params(self, **params):
        """
        Set parameters for the ISLEPath instance.
        """
        for key, value in params.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                raise ValueError(f"Parameter {key} is not recognized.")

    def generate_ensemble(self, X, y):
        base_tree = make_init_model(
            max_leaves=self.max_leaves,
            max_features=self.max_features,
            max_depth=self.max_depth,
            param_file=self.param_file
        )

        start = time.time()
        T, self.estimators = make_isle_ensemble(X, y, base_tree, 
                                                self.eta, self.nu, self.M, 
                                                rng=self.rng, 
                                                verbose=self.verbose)
        end_isle = time.time() - start

        if self.verbose:
            print("Isle ensemble done: ", end_isle, " seconds")

        # centering T for Elastic Net
        self.set_mu_sigma(T)
        T = (T - self.mu)/self.sigma

        return T

    def elnet_pruning(self, T, y):

        extra_str = f"Ensemble(eta={self.eta}, nu={self.nu}, max_leaves={self.max_leaves})"

        self.set_lambdas(T, y)
        elnet = ElasticNet(fit_intercept = self.fit_intercept,
                           max_iter = self.max_iter,
                           init_iter = 1,
                           copyX = True,
                           alpha = self.alpha,
                           tol = self.tol,)
        
        elnet.zero_thresh = self.zero_thresh

        (self.path, 
         self.intercepts) = lasso_path(T, y, self.lambdas, elnet,
                                       print_progress = self.verbose,
                                       extra_str = extra_str)
        # the effect of the intercept when (X,y) are centered
        # is negligible.
        # print("Intercepts: ", self.intercepts)
    
    def set_mu_sigma(self, T):
        self.mu = np.mean(T, axis=0)
        self.sigma = np.std(T, axis=0)
        
        sd_zero = self.sigma == 0
        if np.any(sd_zero):
            self.sigma[sd_zero] = 1

    def set_lambdas(self, X, y):
        self.lambdas = get_lambdas(self.alpha, X, y, self.K, self.epsilon, verbose=self.verbose)
        # max_lam = 1.7468125572517783
        # self.lambdas =  np.logspace(np.log10(max_lam*self.epsilon), np.log10(max_lam), self.K, endpoint=True)[::-1]

    def fit(self, X, y):

        if self.make_ensemble:
            T = self.generate_ensemble(X, y)
            # no need for centering in y.
        else:
            # if we are not making an ensemble,
            # then T is just X
            T = X

        # Lasso will re-center (X,y) and then 
        # optimize for \beta. The \beta will account for the
        # uncentered by getting \beta_0 and 
        # \hat \beta = \beta/\sigma, where \sigma is the
        # original variance. \sigma = ones if X is standarized.

        # notice that if T and y are already centered,
        # then the intercept (\beta_0) should be close to 0.
        self.elnet_pruning(T, y)

    def predict(self, X):
        n,_ = X.shape
        K = len(self.lambdas)

        if self.make_ensemble:
            assert self.estimators is not None, "Model fitting is necessary"
            T_test = make_F_test(X, self.estimators)
            T_test = (T_test - self.mu)/self.sigma

        else:
            # if we are not making an ensemble,
            # then F_test is just X
            T_test = X

        y_pred = np.zeros((n,K))
        for j in range(K):
            y_pred[:,j] = self.intercepts[j] + T_test @ self.path[:, j]
        
        return y_pred

    def score(self, X, y, metric = 'rmse'):
        
        y_pred = self.predict(X)
        # y_pred \in R^{n,K}, e.g.,

        # y_pred = np.ones((3,5)) + np.random.normal(size=(3,5))
        # y = np.array([1,2,3]).reshape(-1,1)

        if metric == 'rmse':
            return np.sqrt(np.mean((y_pred - y.reshape(-1,1))**2, axis=0))
        
        elif metric == 'sse':
            # sum of squared errors 
            return np.sum((y_pred - y.reshape(-1,1))**2, axis=0)
        
        else:
            raise ValueError(f"Unknown metric: {metric}")
    