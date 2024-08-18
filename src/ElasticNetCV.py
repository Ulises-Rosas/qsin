
import random
import numpy as np
import multiprocessing as mp

from qsin.utils import calculate_test_errors
from qsin.sparse_solutions import ElasticNet, lasso_path

def error_fn(tmp_alpha,
             args = None, i = None, 
             X_train_t = None, X_test_t = None, 
             y_train_t = None, y_test_t = None,
             params = None):
    # tmp_alpha = alphas[1]

    if args.verbose:
        print("Fold: ", i, " alpha: ", tmp_alpha)

    model = ElasticNet(fit_intercept = False,
                       max_iter = args.max_iter,
                       init_iter = 1,
                       copyX = True,
                       alpha = tmp_alpha,
                       tol = args.tol,
                       seed = args.seed)
    
    path = lasso_path(X_train_t, y_train_t, 
                      params, model, 
                      sequential_screening = False)
    
    tmp_errors = calculate_test_errors(args, path, params,
                                       X_test_t, y_test_t,
                                       write = False)
        
    # import matplotlib.pyplot as plt
    # plt.plot(tmp_errors[:,0], tmp_errors[:,1])
    # plt.xscale('log')

    return (tmp_errors[:, 1], tmp_alpha)

def get_parallel_errors(args, X_train, y_train, alphas, params, num_folds, ncores):
    """
    let X_t be a fold in the training set 
    and a_j be an alpha in alphas. Then 
    each thread takes the pair (X_t, a_j)
    for all j and t and computes the RMSE for 
    the path of all Lambda values in params.

    For a given alpha j and f folds,
    the RMSE for all Lambda values:

    [ [ RMSE_1,j \in R^{1 x K} ]   -> (X_1, a_j)
       ...
      [ RMSE_f,j \in R^{1 x K} ] ]  -> (X_f, a_j)

    Where  K is the number of Lambda values.

    If the average of the column j is taken,
    then it will effectively be the CV_error 
    for the pair (alpha_j, lambda_i) hyperparameters.
    """
    # X = X_train
    # y = y_train
    # num_folds = 5

    n = X_train.shape[0]
    fold_size = n // num_folds

    #shuffle the data
    random.seed(args.seed)
    all_index = list(range(n))
    random.shuffle(all_index)

    X_train = X_train[all_index, :] # check to shuffle X
    y_train = y_train[all_index] # check to shuffle y!


    out = []
    with mp.Pool( processes = ncores ) as pool:

        preout = []
        for i in range(num_folds):

            test_idx = list(range(i * fold_size, (i + 1) * fold_size))
            train_idx = list(set(range(n)) - set(test_idx))

            X_train_t, X_test_t = X_train[train_idx, :], X_train[test_idx, :]
            y_train_t, y_test_t = y_train[train_idx], y_train[test_idx]

            for al in alphas:

                errors = pool.apply_async(
                    error_fn, 
                    (al, args, i, 
                     X_train_t, X_test_t, 
                     y_train_t, y_test_t, 
                     params)
                )
                preout.append(errors)

        for errors in preout:
            out.append(errors.get())

    # from matplotlib import pyplot as plt
    # for i in range(len(out)):
    #     # i = 0
    #     plt.plot(params['lam'], out[i][0])
    #     plt.xscale('log')

    return out

def get_best_params(all_errors, alphas, params):
    # all_errors = out

    fold_alpha = np.array([b for _,b in all_errors])
    fold_error = np.array([a for a,_ in all_errors])

    best_alpha = 0
    best_lam = 0
    min_rmse = np.inf

    for al in alphas:
        # al
        idx = np.where(fold_alpha == al)[0]
        CV_error =  np.mean(fold_error[idx,:], 0)
        tmp_cv_err = np.min(CV_error)
        if tmp_cv_err < min_rmse:
            min_rmse = tmp_cv_err
            best_alpha = al
            best_lam = params['lam'][np.argmin(CV_error)]
            # print(best_alpha, best_lam, min_rmse)

    return best_alpha, best_lam, min_rmse

def ElasticNetCV_alpha(args, X_train, y_train, alphas, 
                 params, folds, ncores):
    
    """
    Find the best set of hyperparameters for the ElasticNet
    using a cross-validation. The function returns
    the best alpha hyperparameter.

    Higltights: it parallelizes over all the folds and 
    alpha values. 
    """

    if args.verbose:
        print("alphas: ", alphas)
        print("Performing CV with ", folds, " folds")
        
    all_errors = get_parallel_errors(
        args, X_train, y_train, alphas, 
        params, folds, ncores
    )

    (best_alpha, 
     best_lam,
     min_rmse) = get_best_params(all_errors, alphas, params)

    if args.verbose:
        print("best alpha: ", best_alpha)
        print("best lambda: ", best_lam)
        print("min RMSE: ", min_rmse)

    return best_alpha

# alphas = [ 0.5, 0.9 ]
# folds = 2
# ncores = 6
# args.verbose = True
# ElasticNetCV_alpha(args, X_train, y_train, [0.5, 0.9], params, 2, 6)




# move all this below upwards to test the function
# """
# sim_nets.R 15 --max_iter 500 --prefix test --out_path ./test_data/test_sims --ncores 5
# # 10 min
# infer_qlls.jl ./test_data/1_seqgen.CFs_n15.csv\
#               ./test_data/test_sims/test*.txt\
#               --outfile ./test_data/test_qll.csv\
#               --ncores 5 
# """
# from phylokrr.utils import k_folds

# import numpy as np
# import random

# from sparse_solutions import ElasticNet, lasso_path
# from utils import max_lambda, _scaler, calculate_test_errors
# from isle_path import split_data_isle, get_new_path
# from path_subsampling import calculate_test_errors

# from argparse import Namespace
# # print(args)
# Xy_file = "/Users/ulises/Desktop/ABL/software/qsin/test_data/test_qll.csv"
# args = {'p_test': 0.2, 'seed': 0, 
#         'isle': False, 'nwerror': False,
#         'max_features': 0.5, 'max_depth': 5, 
#         'param_file': None, 'eta': 0.1, 
#         'nu': 0.1, 'M': 100, 'verbose': False, 
#         'alpha': 0.5, 'e': 0.01, 'K': 100, 
#         'max_iter': 1000, 'tol': 1e-06, 
#         'wpath': False, 'prefix': 'test', 
#         'CT_file': '../test_data/1_seqgen.CFs_n15.csv', 
#         'factor': 0.5, 'inbetween': 0.5, 'window': 1, 
#         'nstdy': False,}
# args = Namespace(**args)
# args.Xy_file = Xy_file




# data = np.loadtxt(args.Xy_file, delimiter=',', skiprows=1)
# X, y = data[:, :-1], data[:, -1]
# n,p = X.shape

# X = _scaler(X, X, sted = True)
# y = _scaler(y, y, sted = False if args.nstdy else True)

# num_test = int(n*args.p_test)

# (X_train,X_test,
#     y_train,y_test,
#     estimators # None if isle is False
#     ) = split_data_isle(X, y,
#         num_test=num_test, seed=args.seed,
#         isle=args.isle, nwerror=args.nwerror, 
#         mx_p=args.max_features, max_depth=args.max_depth, param_file=args.param_file, 
#         eta=args.eta, nu=args.nu, M=args.M,
#         verbose=args.verbose)

# max_lam = max_lambda(X_train, y_train, alpha=args.alpha)
# min_lam = max_lam * args.e
# params = {'lam': np.logspace(np.log10(min_lam), np.log10(max_lam), args.K, endpoint=True)[::-1]}
