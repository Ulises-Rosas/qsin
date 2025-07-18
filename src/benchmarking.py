from sparse_solutions_hd import ElasticNet, Lasso
from qsin.utils import split_data

import numpy as np
import time
seed = 12038
np.random.seed(seed)
n = 600
p = 100000
lam=0.0001
alpha = 0.5
tol = 1e-5
X = np.random.randn(n, p)
y = np.random.randn(n)
n = len(y)

X_train, X_test, y_train, y_test = split_data(X, y, 200)
X_train = (X_train - np.mean(X_train, axis=0)) / np.std(X_train, axis=0)
X_test = (X_test - np.mean(X_test, axis=0)) / np.std(X_test, axis=0)


# self = Lasso(max_iter=1000, lam=lam, seed=seed, tol=tol, fit_intercept=True)
self = ElasticNet(max_iter=1000, lam=lam, seed=seed, tol=tol, alpha=alpha, fit_intercept=True)
self._verbose = False
start = time.time()
self.fit(X_train, y_train)
print("Time: ", time.time() - start)
print(np.sum(self.beta != 0))
test_error = self.score(X_test, y_test)
print("Test error: ", test_error)

from sklearn.linear_model import ElasticNet as skElasticNet
from sklearn.linear_model import Lasso as skLasso

sklasso = skElasticNet(alpha=lam/2, l1_ratio=alpha, max_iter=1000, tol=tol, fit_intercept=True)
# sklasso = skLasso(alpha=lam/2, max_iter=10000, tol=tol, fit_intercept=True)
start = time.time()
sklasso.fit(X_train, y_train)
print("\nSklearn Time: ", time.time() - start)
sktest_error = np.sqrt(np.mean((y_test - sklasso.predict(X_test))**2))
print("Sklearn Test error: ", sktest_error)
print(sklasso.coef_[sklasso.coef_ != 0].shape)
