import random
import numpy as np

# import matplotlib.pyplot as plt
import copy
from numba import njit


# @njit
def sparse_XB(X, B):
    # return sparse.csr_matrix(X).multiply(sparse.csr_matrix(B))
    return np.matmul(X, B)

# sparse_XB = lambda X, B: np.matmul(X, B)


# np.random.seed(123)
# random.seed(123)

# standarize data
def std_data(X):
    return (X - np.mean(X, axis=0)) / np.std(X, axis=0)

def mse(y, X, B):
    return np.mean( ( y - X.dot(B) )**2 ) 


@njit
def msoft_threshold( delta, lam, denom):
    if delta > lam:
        return (delta - lam) / denom
    elif delta < -lam:
        return (delta + lam) / denom
    else:
        return 0

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

    # return beta


# @jit
def calculate_dual_gap(R, X, y, lam, beta):
    n = len(y)

    p_obj = (1/n) * np.linalg.norm(R)**2  + lam * np.linalg.norm(beta, ord=1)

    theta = (2/n) * R
    norm_inf = np.linalg.norm(X.T @ theta, ord=np.inf)
    # max([1,2])
    # deno = norm_inf if norm_inf > lam else lam
    theta_ = theta/np.max([lam, norm_inf])

    d_obj = (1/n)*( np.linalg.norm(y)**2 - np.linalg.norm(y - (lam*n/2)*theta_)**2 )

    w_d_gap = (p_obj - d_obj)*n/np.linalg.norm(y)**2

    return w_d_gap


def dualpa_max(norm_inf, lam):
    return lam if norm_inf < lam else norm_inf

def dualpa(X, y, lam, beta, ny2):
    R = y - X @ beta
    

    n= len(y)
    c1 = 1/n

    rt = 2*c1*R
    norm_inf = np.linalg.norm(X.T @ rt, ord=np.inf)
    f = dualpa_max(norm_inf, lam)

    c2 = lam/f
    
    d_gap = (1/n)*np.dot(R, (1 + c2**2)*R - 2*c2*y) + lam*np.linalg.norm(beta, ord=1)
    # d_gap = dualpa_gap(R, y, lam, beta, c2, n)

    return d_gap*n/ny2

# @jit
def duality_gap_elnet(R, X, y, beta, lam_1_alpha, lam_alpha):

    n = len(y)

    R_norm_sq = np.linalg.norm(R)**2
    b_norm_sq = np.linalg.norm(beta)**2

    c1 = 1/n
    p_obj = c1 * R_norm_sq  + lam_alpha * np.linalg.norm(beta, ord=1) + (lam_1_alpha/2)*b_norm_sq
    rt = (2/n) * R

    norm_inf = np.linalg.norm(X.T @ rt  - lam_1_alpha*beta, ord=np.inf)

    # f = dualpa_max(norm_inf, lam_alpha)
    f = norm_inf if norm_inf < lam_alpha else lam_alpha
    c2 = lam_alpha/f

    d_obj   = (2/n) * c2 * np.dot(R, y) - (c2**2) * (c1 * R_norm_sq + (lam_1_alpha/2)*b_norm_sq)
    w_d_gap = (p_obj - d_obj)*n/np.linalg.norm(y)**2

    return w_d_gap

# @njit
def dualpa_elnet(R, X, y, beta, lam_1_alpha, lam_alpha, ny2):
    
    # R_norm_sq = np.linalg.norm(R)**2
    if lam_1_alpha == 0:
        return dualpa(X, y, lam_alpha, beta, ny2)
    
    # R = y - X @ beta

    n = len(y)
    # b_norm_sq = np.linalg.norm(beta)**2
    b_norm_sq = (beta**2).sum()

    rt = (2/n) * R
    norm_inf = np.linalg.norm(X.T @ rt  - lam_1_alpha*beta, ord=np.inf)

    # f = dualpa_max(norm_inf, lam_alpha)
    f = norm_inf if norm_inf < lam_alpha else lam_alpha
    c2 = lam_alpha/f
    c3 = 1 + c2**2

    d_gap = (1/n)*np.dot(R, c3*R - 2*c2*y) +\
            lam_alpha*np.linalg.norm(beta, ord=1) +\
            c3*(lam_1_alpha/2)*b_norm_sq

    return d_gap*n/ny2


class Lasso:
    """
    lasso model.
    It assumes that the input data has been standardized
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
                 copyX = False,
                 seed = 123,
                 **kwargs):
        
        self.max_iter = max_iter
        self.lam = lam
        self.prev_lam = prev_lam
        self.tol = tol
        self.warm_start = warm_start
        self.seed = seed

        self.fit_intercept = fit_intercept
        


        # if self.warm_start:
        #     self.beta = beta
        # else:
        self.beta = np.array([])


        self.intercept = 0
        
        self.init_iter = init_iter
        self.copyX = copyX

        self.X = np.array([])
        self.X_T_X = np.array([])
        self.Xj_T_Xj = np.array([])
        self.Xj_T_y = np.array([])
        self.ny2 = np.array([])

    # @njit
    def soft_threshold(self, delta, lam, denom):
        return msoft_threshold(delta, lam, denom)
        
    def initial_iterations(self, c1, Xj_T_y, Xj_T_Xj, X_T_X, all_p):

        for _ in range(self.init_iter):
            self.update_beta(c1, Xj_T_y, Xj_T_Xj, X_T_X, all_p)

    def initial_active_set(self, X, y, c1, Xj_T_y, Xj_T_Xj, X_T_X, all_p):

        if self.prev_lam is None:
            # few iterations of coordinate descent
            self.initial_iterations(c1, Xj_T_y, Xj_T_Xj, X_T_X, all_p)
            # we define an active set A as the set of indices
            A = np.where(self.beta != 0)[0]

        else:
            A = self.sequential_strong_rule(X, y, c1, all_p)

        return A
    
    def sparse_XB(self, X):
        """
        This function returns the sparse matrix
        X[:, A] * B[A], where A is the index set of
        non-zero coefficients
        """
        # A = np.where(self.beta != 0)[0]
        # return np.matmul(X[:, A], self.beta[A])

        return np.matmul(X, self.beta)
        # return sparse_XB(X, self.beta)
    
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

        if self.copyX and len(self.X) == 0:
            self.X = copy.deepcopy(X)
            self.X_T_X = np.matmul(X.T, X)
            self.Xj_T_Xj = np.diag(self.X_T_X)
            self.Xj_T_y = np.matmul(X.T, y)
            self.ny2 = np.linalg.norm(y)**2

        elif not self.copyX and len(self.X) == 0:
            self.X = X
            self.X_T_X = np.matmul(X.T, X)
            self.Xj_T_Xj = np.diag(self.X_T_X)
            self.Xj_T_y = np.matmul(X.T, y)
            self.ny2 = np.linalg.norm(y)**2

        elif self.copyX and len(self.X) != 0:
            pass

        elif not self.copyX and len(self.X) != 0:
            self.X = X
            self.X_T_X = np.matmul(X.T, X)
            self.Xj_T_Xj = np.diag(self.X_T_X)
            self.Xj_T_y = np.matmul(X.T, y)
            self.ny2 = np.linalg.norm(y)**2

    def dual_gap(self, y):
        return dualpa(self.X, y, self.lam, self.beta, self.ny2)
        # R = y - self.sparse_XB(self.X)
        # return dualpa(R, self.X, y, self.lam, self.beta)
    
    def fit(self, X, y):
        # X = X_train
        # y = y_train
        # self.max_iter = 10000
        # self.lam = 10
        # self.set_params(max_iter = 1000, lam = 300)
        
        n, p = X.shape
        c1 = 2 / n

        self.set_Xy(X,y)

        # X_T_X = self.X_T_X
        # # precompute Xj^T * Xj
        # Xj_T_Xj = self.Xj_T_Xj
        # # precompute Xj^T * y
        # Xj_T_y = np.matmul(self.X.T, y)
        
        if not self.warm_start:
            # He-styled 
            # initialization 
            # of the coefficients
            np.random.seed(self.seed)
            self.beta = np.random.normal(0, np.sqrt(2/p), size=p) 
         
        # few iterations of coordinate descent
        all_p = np.array(range(p))
        A = self.initial_active_set(self.X, y, c1, 
                                    self.Xj_T_y, 
                                    self.Xj_T_Xj, 
                                    self.X_T_X, 
                                    all_p)
        

        left_iter = self.max_iter - self.init_iter
        for i in range(left_iter):
            # i = 1
            diff = np.zeros(n)
            
            # s_old = np.matmul( self.X[:, A], self.beta[A] )
            s_new = np.zeros(n)
            self.cd_epoch(c1, A, s_new, diff)

            # self.update_beta(c1, self.Xj_T_y, self.Xj_T_Xj, self.X_T_X, A)
            # s_new = np.matmul( self.X[:, A], self.beta[A] )
            # s_new =  self.beta

            # diff = s_new - s_old
            max_updt = np.max(np.abs(diff))
            w_max = np.max(np.abs(s_new))

            # max_updt = np.linalg.norm(s_new - s_old, ord=np.inf)
            # w_max = np.linalg.norm(s_new, ord=np.inf)
            

            if w_max == 0 or max_updt/w_max < self.tol:
            # if np.linalg.norm(s_new - s_old) < self.tol:
                # print("Max updates: ", max_updt)

                # A_ = set(all_p) - set(A)
                A_ = np.setdiff1d(all_p, A)

                #TODO: streamline updt. beta and exclusion test
                self.update_beta(c1, 
                                 self.Xj_T_y, 
                                 self.Xj_T_Xj,
                                 self.X_T_X, A_)
                A_new = self.exclusion_test(self.X, y, c1, A, all_p)
                

                if len(A_new) == 0:
                    # it means that all
                    # coefficients follow the
                    # KKT conditions
                    # print('kkt, finished at iteration: ', i)
                    break

                else:

                    # R = y - self.sparse_XB(self.X)
                    # w_d_gap = dualpaR(R, self.X, y, self.lam, self.beta)
                    w_d_gap = self.dual_gap(y)
                    # w_d_gap = dualpa(self.X, y, self.lam, self.beta)

                    if w_d_gap < self.tol:
                        # print('dual, finished at iteration: ', i)
                        break

                    else:
                        A = np.concatenate((A, A_new))
                        A.sort(kind = 'mergesort')

        if i == left_iter - 1:
            # if the iterations reach this point,
            # it means that there is still an active set.
            # then, the model did not converge
            print("Model did not converge")

        if self.fit_intercept:
            self.intercept = np.mean( y - self.sparse_XB(X) )
            # print("Intercept: ", self.intercept)

    def cd_epoch(self, c1, chosen_ps, s_new, diff):
        epoch_lasso(self.X, self.beta, self.lam, c1,
                     self.Xj_T_y, self.Xj_T_Xj, self.X_T_X, 
                     chosen_ps, s_new, diff)

    def update_beta(self, c1, Xj_T_y, Xj_T_Xj, X_T_X, chosen_ps):
        update_beta_lasso(self.beta, self.lam, c1, Xj_T_y, Xj_T_Xj, X_T_X, chosen_ps)

        # for j in chosen_ps:
        #     # j = 1
        #     A2 = np.where(self.beta != 0)[0]
        #     # take out j from the array A2
        #     A2 = A2[A2 != j]

        #     delta = c1 * ( Xj_T_y[j] - np.dot(X_T_X[j, A2], self.beta[A2]) )
        #     denom = c1 * Xj_T_Xj[j]
        #     self.beta[j] = self.soft_threshold(delta, self.lam, denom)

    def exclusion_test(self, X, y, c1, A, all_p):
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
        # A_ = set(all_p) - set(A)
        A_ = np.setdiff1d(all_p, A)

        if len(A_) == 0:
            return A_

        r = y - self.sparse_XB(X)

        # A_ = np.array(list(A_))
        # exclusion test (SLS book, page 114) based on
        # KTT conditions
        e_test = c1 * np.abs( np.matmul( X.T[A_, :], r) ) >= self.lam
        # those that are in the new active set
        # are those who did pass the KTT conditions
        return A_[e_test]

    def sequential_strong_rule(self, X, y, c1, all_p):
        """
        rule based on the dual polytope projection 
        and futher modified by the strong rule
        showed at SLS book page 128
        """
        A_ = np.array(all_p)

        r = y - self.sparse_XB(X)

        # sequential strong rule (SLS book, page 128)
        e_test = c1 * np.abs( np.matmul( X.T[A_, :], r) ) >= 2*self.lam - self.prev_lam
        return A_[e_test]
    
    def predict(self, X):
        if self.beta is None:
            raise Exception("Model has not been fitted yet.")
        
        return self.sparse_XB(X) + self.intercept
    
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

class ElasticNet(Lasso):
    def __init__(self, 
                 max_iter=300, 
                 alpha = 0.5, 
                 lam=0.1, 
                 prev_lam=None, 
                 fit_intercept=True, 
                 warm_start=False,
                #  beta=None,
                 tol=0.001, 
                 **kwargs):
        super().__init__(max_iter, lam, prev_lam, fit_intercept, warm_start,  tol, **kwargs)

        self.alpha = alpha
        self.lam = lam

        # when lam or alpha
        # change, these values are
        # updated at set_params
        self.lam_alpha = alpha * lam
        self.lam_1_alpha = (1 - alpha) * lam

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

    def cd_epoch(self, c1, chosen_ps, s_new, diff):
        epoch_enet(self.X, self.beta, c1,
                   self.Xj_T_y, self.Xj_T_Xj, self.X_T_X, 
                   chosen_ps, self.lam_1_alpha, self.lam_alpha, s_new, diff)


    def update_beta(self, c1, Xj_T_y, Xj_T_Xj, X_T_X, chosen_ps):
        update_beta_enet(self.beta, c1, Xj_T_y, Xj_T_Xj,
                                     X_T_X, chosen_ps, self.lam_1_alpha,
                                     self.lam_alpha)

    def dual_gap(self, y):
        
        R = y - self.sparse_XB(self.X)
        return dualpa_elnet(R, self.X, y, self.beta, self.lam_1_alpha, self.lam_alpha, self.ny2)
        # return duality_gap_elnet(R, self.X, y, 
        #                           self.beta, self.lam_1_alpha,
        #                             self.lam_alpha)
    
    def exclusion_test(self, X, y, c1, A, all_p):
        # A_ = set(all_p) - set(A)
        A_ = np.setdiff1d(all_p, A)

        if len(A_) == 0:
            return A_

        r = y - self.sparse_XB(X)

        # A_ = np.array(list(A_))
        # exclusion test (SLS book, page 114) based on
        # KTT conditions

        # Xj_T_r = np.matmul( X.T[A_, :], r)

        if self.lam_1_alpha == 0:
            elastic_term = 0
        else:
            elastic_term = (1/c1) * self.lam_1_alpha * self.beta[A_]

        elastic_test = c1 * np.abs(np.matmul( X.T[A_, :], r) - elastic_term) >= self.lam_alpha
        # those that are in the new active set
        # are those who did pass the KTT conditions
        return A_[elastic_test]

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


def lasso_path(X_train, y_train, params, model, 
               sequential_screening = True):
    """
    compute the lasso path based on the training set
    and  with errors based on the test set
    """
    # model = Lasso()
    # X = X_train
    # y = y_train
    # params = {'lam': np.logspace(-2, max_lambda(X,y, alpha), 3)}

    # if X_test is None and y_test is None:
    #     X_test = X_train
    #     y_test = y_train


    _,p = X_train.shape
    lams = params['lam']

    # errors = np.zeros(len(lams))
    path = np.ones((p, len(params['lam'])))

    model.set_params(lam=lams[0])
    model.fit(X_train, y_train)

    path[:,0] = model.beta
    
    for i in range(1, len(lams)):

        model.set_params(lam = lams[i],
                          warm_start = True, 
                        #   prev_lam = lams[i-1],
                          prev_lam = None
                          )
            
        model.fit(X_train, y_train)            
        # errors[i] = None
        path[:,i] = model.beta

    return path

def myenet_path(X_train, y_train, X_test, y_test, params, model, 
               sequential_screening = True,
               alpha = 0.5):
    """
    compute the lasso path based on the training set
    and  with errors based on the test set
    """
    # model = Lasso()
    # X = X_train
    # y = y_train
    # params = {'lam': np.logspace(-2, max_lambda(X,y), 3)}

    if X_test is None and y_test is None:
        X_test = X_train
        y_test = y_train


    _,p = X_train.shape
    lams = params['lam']

    errors = np.zeros(len(lams))
    path = np.ones((p, len(params['lam'])))
    y_train_enet = np.concatenate((y_train, np.zeros(p)))

    all_gammas = np.zeros(len(lams))
    for i,lam in enumerate(lams):
        # i = 0
        # lam = lams[i]
        alfa = lam * alpha
        beta = lam * (1 - alpha)
        sf = (1 + beta) ** (-1/2)
        gamma = alfa * sf
        all_gammas[i] = gamma

        if i > 1:
            if sequential_screening:
                prev_lam = all_gammas[i-1]

            else:
                prev_lam = None

            model.set_params(lam = gamma,
                             warm_start = True, 
                             beta = model.beta,
                             prev_lam = prev_lam)
            
        else:

            model.set_params(lam=gamma)
            # model.warm_start

        X_train_enet =  np.vstack([X_train, np.sqrt(beta) * np.eye(p)])
        model.fit(X_train_enet, y_train_enet)
        errors[i] = None
        path[:,i] = model.beta*sf

    return path, all_gammas

def theta(path_i):

    path = np.copy(path_i)
    theta_path = np.zeros(path.shape)

    for i in range(path.shape[1]):
        tmp_path = path[:,i]
        logic_gate = tmp_path != 0

        tmp_path[logic_gate] = 1
        theta_path[:,i] = tmp_path

    return theta_path

def split_data(X,y,num_test, seed = 123):

    random.seed(seed)
    n,_ = X.shape

    test_idx  = random.sample(range(n), k = num_test)
    train_idx = list( set(range(n)) - set(test_idx) )

    X_train, X_test = X[train_idx,:], X[test_idx,:]
    y_train, y_test = y[train_idx]  , y[test_idx]

    return X_train, X_test, y_train, y_test

def get_ZO(path, spps_array, n_spps):

    Z = np.zeros((n_spps, path.shape[0]))
    for i in range(path.shape[0]):
        Z[:,i][spps_array[i,:] - 1] = 1

    theta_path = theta(path)
    return Z @ theta_path

def get_non_zero_coeffs(path, ZO, thresh = 0.5):
    n_features = path.shape[0]
    string_version = []
    non_zero_coeffs = []
    for j in range(path.shape[1]):
        # j = 10
        path_j = path[:,j]
        ZO_j = ZO[:,j]

        # checking if all species are
        # covered by the selected coefficients
        if np.any(ZO_j == 0):
            continue

        non_zero = path_j != 0

        if np.all(non_zero):
            break

        if np.sum(non_zero) <= thresh*n_features:
            non_zero_idx = list(np.where(non_zero)[0])
            non_zero_idx_str = str(set(np.sort(non_zero_idx)))

            if non_zero_idx_str not in string_version:
                string_version.append(non_zero_idx_str)
                non_zero_coeffs.append(non_zero_idx)

    return non_zero_coeffs


