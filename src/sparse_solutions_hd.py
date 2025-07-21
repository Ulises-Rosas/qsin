# high-dimensional sparse solutions

from collections import deque


import numpy as np
from qsin.utils import progressbar
from qsin.CD_dual_gap import (epoch_lasso_v2, update_beta_lasso_v2, dualpa,
                              epoch_enet_v2 , update_enet_v2, dualpa_elnet)


class Lasso:
    """
    lasso model.
    It assumes that the input data has been standardized. 
    If not, lasso does it nonetheless on (X,y)

    init beta is a zero vector
    """
    def __init__(self, 
                 max_iter=300, 
                 lam=0.1,
                 prev_lam = None,
                 fit_intercept  = True,
                 warm_start = False,
                #  beta = None,
                 tol = 0.001,
                 init_iter = 1,
                 copyX = False,):
        
        self.max_iter = max_iter
        self.lam = lam
        self.prev_lam = prev_lam
        self.tol = tol
        self.warm_start = warm_start

        self.fit_intercept = fit_intercept
        self.beta = np.array([])


        self.intercept = 0
        
        self.init_iter = init_iter
        self.copyX = copyX

        self.X = np.array([])
        self.ny2 = np.array([])
        self.r = np.array([])

        self.y = np.array([])

        self.x_bar = np.array([])
        self.y_bar = 0.0

        self.x_sd = np.array([])
        self.y_sd = 1.0

        self._verbose = False

        self.zero_thresh = 1e-15
        self.message_conv = ""

    def __str__(self):
        return f"Lasso(tol={self.tol}, lam={self.lam}, max_iter={self.max_iter})"
    
    def update_beta(self, c1, n, all_p):
        update_beta_lasso_v2(
            self.X, self.beta, self.lam, self.r,
            c1, n, all_p, self.zero_thresh
        )

    def initial_iterations(self, c1, all_p):

        n,_ = self.X.shape
        for _ in range(self.init_iter):
            # O(np)
            self.update_beta(c1, n, all_p)

    def sorted_init_active_set(self):
        """
        Get the active set A: non-zero coefficients

        Returns
        -------

        Aq: the active set as a deque sorted with respect to the indices\\
        A: the active set as a set, not necessarily sorted.
        """

        Aq = deque()
        A = set()
        # O(2p)
        for j, b in enumerate(self.beta):
            # if abs(b) >= self.zero_thresh:
            if b != 0.0:
                # O(2)
                A.add(j)
                Aq.append(j)

        return Aq, A

    def initial_active_set(self, c1, all_p):

        # few iterations of coordinate descent
        # O(np*T), where T is the number of initial iterations
        self.initial_iterations(c1, all_p)

        # we define an active set A as the set of indices
        # O(2p)
        return self.sorted_init_active_set()
    
    def set_Xy(self, X,y):
        """
        set X

        This is used under the 
        following logic

        if copy and x == 0:
            set_X
        elif not copy and x == 0:
            set_X (overwrite)
        elif copy and x != 0:
            do nothing
        elif not copy and x != 0:
            set_X (overwrite)
        """

        if self.copyX and len(self.X) != 0:
            pass

        else:
            # as contiguous array
            # makes the accessing of columns
            # faster. The same for residuals.
            # this fast block access is important
            # for multiple dot products in the
            # coordinate descent
            # XXX
            # asfortranarray makes a copy of the array
            self.X = np.asfortranarray(X) # O(np)
            
            # ascontiguousarray does not make a copy
            # but ensures that the array is contiguous
            self.y = np.ascontiguousarray(y.copy(), dtype=y.dtype) # O(n)
            self.center_Xy() # O(np)

            self.ny2 = np.dot(self.y, self.y) # O(n)
            self.r = self.y - self.X @ self.beta # O(np)
            self.r = np.ascontiguousarray(self.r, dtype=self.X.dtype) # O(n)
        
    def center_Xy(self):
        """
        center X and y
        """
        self.x_bar = np.mean(self.X, axis=0, dtype=self.X.dtype) # O(np)
        self.x_sd  = np.std(self.X , axis=0, dtype=self.X.dtype) # O(np)
        
        # avoid division by zero
        zero_sd_indices = self.x_sd == 0 # O(p)
        if np.any(zero_sd_indices): # O(p)
            self.x_sd[zero_sd_indices] = 1 # O(p)

        self.X = (self.X - self.x_bar) / self.x_sd # O(np)

        self.y_bar = np.mean(self.y) # O(n)
        self.y -= self.y_bar # O(n)

    def dual_gap(self, y):
        return dualpa(self.X, y, self.r, self.lam, self.beta, self.ny2)

    def get_sorted_complement(self, A, all_p):
        """
        get the sorted complement of the active set

        A is the active set and all_p is the set of all
        indices. It returns ordered list of indices
        """
        Ac_q = deque()
        # O(2p)  in average
        for j in all_p:
            if j not in A:
                Ac_q.append(j)

        return np.array(Ac_q) # O(p)
    
    def update_sorted_active_set(self, A, Ac_f, all_p):
        """
        update the active set based on the exclusion test
        it returns the new active set ordered with respect
        to the indices
        """
        Anew = deque()
        # O(3p) in average
        for j in all_p:
            if (j in A) or (j in Ac_f):
                Anew.append(j)

        return Anew
        
    def fit(self, X, y):
        # X = X_train
        # y = y_train
        # self.max_iter = 10000
        # self.set_params(max_iter = 100, lam = 0.1)
        
        n, p = X.shape
        c1 = 2 / n
        
        if not self.warm_start:            
            # zero initialization
            self.beta = np.zeros(p, dtype=X.dtype, order='F')

        # XXX, O(np)
        self.set_Xy(X, y)

        # few iterations of coordinate descent
        all_p = np.array(range(p))

        # O(2p + np*T_init) = O(np)
        Aq, A = self.initial_active_set(c1, all_p)

        
        if len(A) == 0:
            # if the active set is empty
            # then the model is converged
            if self.fit_intercept:
                self.set_intercept()
            return

        A_rr = np.array(Aq, dtype=np.int64)
        # print("Active set: ", A_rr)

        left_iter = self.max_iter - self.init_iter
        
        # O(T*c1*np + c2*p + c3*n) = O(np) for p >> n and T small,
        # where T is the number of left iterations 
        # and c1, c2, c3 are constants
        # Small T are practially possible with the warm starts
        # and thorough path of the lasso
        for i in range(left_iter):

            # O(n)
            xb_diff = np.zeros(n, dtype=np.float64)
            xb_new = np.zeros(n, dtype=np.float64)
            # O(np)
            self.cd_epoch(c1, n, A_rr, xb_new, xb_diff)

            # O(n)
            # checking convergence on the beta values only
            # might be not correct as in high dimensions (p > n) there 
            # does not exist a unique solution (i.e., X'X is singular). 
            # However, there is a unique solution for XB. so, we can check
            # the convergence on the XB values
            # ref: Sparsity, the Lasso, and Friends 
            # (section 3.2), Ryan Tibshirani, 2017.
            max_updt = np.max(np.abs(xb_diff))
            w_max = np.max(np.abs(xb_new))

            if self._verbose:
                print("iteration:", i,"Max update: ", max_updt, "Max weight: ", w_max)

            if w_max == 0 or max_updt/w_max < self.tol:
                
                # O(3p)
                Ac_arr = self.get_sorted_complement(A, all_p)

                #TODO: streamline updt. beta and exclusion test
                # O(np)
                self.update_beta(c1, n, all_p)

                # O(np + p) = O(np)
                # A_c that failed the exclusion test
                Ac_f_arr = self.exclusion_test(self.X, c1, Ac_arr)
                Ac_f = set(Ac_f_arr)
                

                if len(Ac_f_arr) == 0:
                    # it means that all
                    # coefficients follow the
                    # KKT conditions
                    if self._verbose:
                        print('kkt, finished at iteration: ', i)
                    break

                else:
                    # O(np + n + p) = O(np)
                    w_d_gap = self.dual_gap(self.y)
                    if w_d_gap < self.tol:
                        if self._verbose:
                            print('dual, finished at iteration: ', i)
                        break

                    else:
                        # O(3p)
                        Anew = self.update_sorted_active_set(A, Ac_f, all_p)
                        # O(2p)
                        A = set(Anew)
                        A_rr = np.array(Anew)

        if i == left_iter - 1:
            # if the iterations reach this point,
            # it means that there is still an active set.
            # then, the model did not converge
            self.message_conv = f"{self.__str__()} did not converge. Try to increase either max_iter or tol params."

        if self.fit_intercept:
            self.set_intercept()

    def set_intercept(self):
        # this scaling will consider the unscaled
        # data. The optimization was done assuming the
        # scaled data. 
        # If x_sd ever had a zero value, then
        # it was replaced by 1.
        self.beta /= self.x_sd
        self.intercept = self.y_bar - np.dot(self.x_bar, self.beta) # O(p)

        if abs(self.intercept) < self.zero_thresh:
            # if the intercept is close to zero, then
            # we can set it to zero
            self.intercept = 0.0

    def cd_epoch(self, c1, n, chosen_ps, s_new, s_diff):
        epoch_lasso_v2(
            self.X, self.beta, self.lam, self.r,
            c1, n, chosen_ps, s_new, s_diff, self.zero_thresh
        )

    def exclusion_test(self, X, c1, A_c_arr):
        """
        Exclusion test (SLS book, page 114)

        The optimization problem dictates that:

        X^T * r = lam * gamma
        
        gamma being partial of ||B||_1 and whose
        subgradient is given by:

        gamma_j = + 1     if beta_j > 0 (sign of beta_j)\\
        gamma_j = - 1     if beta_j < 0 (sign of beta_j)\\
        gamma_j = (-1, 1) if beta_j = 0 (sign of 0, undefined)\\

        Since the range of gamma_j is between -1 and 1,\\
        then the range of lam * gamma_j is between -lam and lam.\\
        So, the range of Xj^T * r is between -lam and lam.

        Leading the following set of inequalities:
        -lam <= Xj^T * r <= lam

        which from:\\
        lam >= Xj^T * r\\
        lam >= -Xj^T * r

        Implies: 
        lam >= |Xj^T * r|

        Now, we can define the exclusion test as:
        lam > |Xj^T * r|
        whose strict inequiality test for beta_j = 0
        as +1 and -1 are not in the range of the subgradient

        If we pass this test over the omited  variables
        that were supposed to be zero and fails, it means
        they were actually non-zero and should be included
        in the active set.
        """

        if len(A_c_arr) == 0:
            return A_c_arr
        
        # exclusion test (SLS book, page 114) based on
        # KTT conditions

        # O(np)
        e_test = c1 * np.abs( np.matmul( X.T[A_c_arr, :], self.r) ) >= self.lam
        # those that are in the new active set
        # are those who did not pass the KTT conditions.
        # true means that they failed the exclusion
        # test and should be included in the active set
        
        return A_c_arr[e_test] # O(p)

    def predict(self, X):
        if self.beta is None:
            raise Exception("Model has not been fitted yet.")
        
        return X @ self.beta + self.intercept
    
    def score(self, X, y):
        y_pred = self.predict(X)
        return np.sqrt( np.mean( (y - y_pred)**2 ) )
    
    def set_params(self, **params):
        if "max_iter" in params:
            self.max_iter = int(params["max_iter"])

        if "lam" in params:
            self.lam = params["lam"]

        if "intercept" in params:
            self.intercept = params["intercept"]

        if "tol" in params:
            self.tol = params["tol"]

        if "warm_start" in params:
            self.warm_start = params["warm_start"]

        if "beta" in params:
            self.beta = params["beta"]

        if "prev_lam" in params:
            self.prev_lam = params["prev_lam"]

    def get_params(self):
        return {'max_iter': self.max_iter, 
                'lam': self.lam, 
                'intercept': self.intercept, 
                'tol': self.tol}
    
    def set_beta(self, beta):
        self.beta = beta

# region: ElasticNet
class ElasticNet(Lasso):
    def __init__(self, 
                 max_iter=300, 
                 alpha = 0.5, 
                 lam=0.1, 
                 prev_lam=None, 
                 fit_intercept=True, 
                 warm_start=False,
                 tol=1e-4,
                 init_iter=1,
                 copyX=False):
        # passes arguments to the Lasso class
        super().__init__(max_iter, lam, prev_lam, fit_intercept, warm_start,  tol, init_iter, copyX)

        self.alpha = alpha
        self.lam = lam

        # when lam or alpha
        # change, these values are
        # updated at set_params
        self.lam_alpha = alpha * lam
        self.lam_1_alpha = (1 - alpha) * lam

    def __str__(self):
        return f"ElasticNet(tol={self.tol}, lam={self.lam}, alpha={self.alpha}, max_iter={self.max_iter})"

    def set_lam_alpha(self, lam, alpha):
        self.lam_alpha = alpha * lam
        self.lam_1_alpha = (1 - alpha) * lam

    def set_params(self, **params):
        if "max_iter" in params:
            self.max_iter = int(params["max_iter"])

        if "lam" in params:
            self.lam = params["lam"]

            if "alpha" in params:
                self.alpha = params["alpha"]
            
            self.set_lam_alpha(self.lam, self.alpha)
            # self.lam_alpha = self.alpha * self.lam
            # self.lam_1_alpha = (1 - self.alpha) * self.lam

        if 'alpha' in params:
            self.alpha = params['alpha']

            if "lam" in params:
                self.lam = params["lam"]

            self.set_lam_alpha(self.lam, self.alpha)
            # self.lam_alpha = self.alpha * self.lam
            # self.lam_1_alpha = (1 - self.alpha) * self.lam

        if "intercept" in params:
            self.intercept = params["intercept"]

        if "tol" in params:
            self.tol = params["tol"]

        if "warm_start" in params:
            self.warm_start = params["warm_start"]

        if "beta" in params:
            self.beta = params["beta"]

        if "prev_lam" in params:
            self.prev_lam = params["prev_lam"]

    def cd_epoch(self, c1, n, chosen_ps, xb_new, xb_diff):
        epoch_enet_v2(self.X, self.beta, self.lam_1_alpha, self.lam_alpha, self.r,
                      c1, n, chosen_ps, xb_new, xb_diff, self.zero_thresh)

    def update_beta(self, c1, n, chosen_ps):
        update_enet_v2(self.X, self.beta, self.lam_1_alpha, self.lam_alpha, self.r,
                       c1, n, chosen_ps, self.zero_thresh )

    def dual_gap(self, y):
        return dualpa_elnet(self.r, self.X, y, self.beta, self.lam_1_alpha, self.lam_alpha, self.ny2)
    
    def exclusion_test(self, X, c1, A_c_arr):

        if len(A_c_arr) == 0:
            return A_c_arr

        # exclusion test (SLS book, page 114) based on
        # KTT conditions
        if self.lam_1_alpha == 0:
            elastic_term = 0
        else:
            elastic_term = (1/c1) * self.lam_1_alpha * self.beta[A_c_arr]

        elastic_test = c1 * np.abs( np.matmul( X.T[A_c_arr, :], self.r ) - elastic_term ) >= self.lam_alpha
        # those that are in the new active set
        # are those who did not pass the KTT conditions
        return A_c_arr[elastic_test]
# endregion

# region: Cross-validation

def k_fold_cv(X, y, model, num_folds):
    n, p = X.shape
    fold_size = n // num_folds
    mse_sum = 0

    for i in range(num_folds):

        test_idx = list(range(i * fold_size, (i + 1) * fold_size))
        train_idx = list(set(range(n)) - set(test_idx))

        X_train, X_test = X[train_idx, :], X[test_idx, :]
        y_train, y_test = y[train_idx], y[test_idx]

        model.fit(X_train, y_train)
        # y_pred = lasso.predict(X_test)
        mse_sum += model.score( X_test, y_test )

    return mse_sum / num_folds

def k_fold_cv_random(X, y, 
                     model, 
                     params,
                     num_folds = 3, 
                     sample = 100,
                     verbose = False,
                     seed = 123,
                     warm_starts = False
                     ):
    
    # model = Lasso()
    # X = X_train
    # y = y_train
    
    np.random.seed(seed=seed)
    
    all_params = params.keys()
    tested_params = np.ones((sample, len(params)))

    for n,k in enumerate(all_params):
        tested_params[:,n] = np.random.choice(params[k], sample)
    
    if warm_starts:
        # check index where 'lam' is in 
        # all_params
        idx = list(all_params).index('lam')
        # sort the tested_params by using the idx in the decreasing
        # order
        tested_params = tested_params[tested_params[:,idx].argsort()[::-1]]
        # tested_params = tested_params[tested_params[:,idx].argsort()]
    
    all_errors = []
    for vec in tested_params:
        # vec = tested_params[1]
        tmp_params = dict(zip(all_params, vec))

        if warm_starts and len(model.beta):
            model.set_params(**tmp_params, 
                             warm_start=True, 
                             beta=model.beta)
        else:
            model.set_params(**tmp_params)


        tmp_err = k_fold_cv(X, y, model, num_folds)
        all_errors.append([tmp_params, tmp_err])

        if verbose:
            print('Error: %s, tested params: %s' % (tmp_err, vec))

    best_ = sorted(all_errors, key=lambda kv: kv[1], reverse=False)[0]
    if verbose:
        print("CV score: ", best_[1])

    return best_[0]

def lasso_path(X_train, y_train, lams, model, print_progress = True,
               extra_str = ""):
    """
    compute the lasso path based on the training set
    and  with errors based on the test set

    it also returns the intercepts
    for each lambda value in the path
    so that the model can be used for prediction
    with the intercepts
    """

    _,p = X_train.shape
    K  = len(lams)

    path = np.ones((p, K))
    intercepts = np.zeros(K)

    model.set_params(lam=lams[0])
    model.fit(X_train, y_train)
    # all values at the beginning of the path are zeros.
    # convergence warning are due to floating point errors
    # and it should be ignored
    if model.message_conv:
        model.message_conv = ""

    path[:,0] = model.beta
    intercepts[0] = model.intercept

    if print_progress:
        index_set = progressbar(range(1, len(lams)), "Computing lasso path: ", 40)
        
    else:
        index_set = range(1, len(lams))
    
    for i in index_set:

        model.set_params(lam = lams[i],
                          warm_start = True, 
                          prev_lam = None)
        model.fit(X_train, y_train)

        if model.message_conv:
            print(f"Warning: {extra_str}+{model.message_conv}")
            model.message_conv = ""

        path[:,i] = model.beta
        intercepts[i] = model.intercept

    return path, intercepts

# endregion

