#!/usr/bin/env python3

import time
import argparse


import numpy as np
# from qsin.sparse_solutions import ElasticNet, lasso_path
from qsin.sparse_solutions_hd import ElasticNet, lasso_path
from qsin.utils import  calculate_test_errors, max_lambda
from qsin.isle_path import split_data_isle, get_new_path
from qsin.ElasticNetCV import ElasticNetCV_alpha
from qsin.row_selection import row_selection, write_rows


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

def max_features_type(value):
    try:
        return float(value) if value not in ['sqrt', 'log2'] else value
    
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid value for --max_features: {value}. Must be a float, 'sqrt', or 'log2'.")

def main():

    # region: parse arguments
    parser = argparse.ArgumentParser(description="""
    Generate batches from ElasticNet path.
    """, 
    add_help=False,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    req_args = parser.add_argument_group("Required arguments")
    req_args.add_argument("Xy_file", help="Path to the Xy input data file.")
    req_args.add_argument("CT_file", help="Path to the CT file.")

    opt_args = parser.add_argument_group("Optional arguments")
    opt_args.add_argument("-h","--help", action="help", help="Show this help message and exit.")
    opt_args.add_argument("--verbose", action="store_true", help="Whether to print verbose output.")
    opt_args.add_argument("--isle", action="store_true", help="Whether to use path from decision tree-based ISLE (i.e., ensemble learning).")
    opt_args.add_argument("--p_test", type=float, default=0.35, metavar="", help="Proportion of observations to use for testing.")
    opt_args.add_argument("--seed", type=int, default=12038, metavar="", help="Random seed.")
    


    isle_args = parser.add_argument_group("ISLE parameters (if --isle is used)")
    isle_args.add_argument("--eta", type=float, default=0.5, metavar="", help="Proportion of observations to use in each tree.")
    isle_args.add_argument("--nu", type=float, default=0.1, metavar="", help="Learning rate.")
    isle_args.add_argument("--M", type=int, default=500, metavar="", help="Number of trees in the ensemble.")
    isle_args.add_argument("--max_depth", type=int, default=5, metavar="", help="Maximum depth of the decision tree.")
    isle_args.add_argument("--max_leaf_nodes", type=int, default=6, metavar="", help="Maximum number of leaf nodes in the decision tree.")
    # default according Hastie et al. 2009, pp. 592
    isle_args.add_argument("--max_features", type=max_features_type, default=0.3, metavar="", help="Maximum proportion of features it is considered to grow nodes in a regression tree. It can also be 'sqrt' or 'log2'.")
    isle_args.add_argument("--param_file", type=str, default=None, metavar="", help="""JSON file with parameters for the decision tree
                           different from max_depth and mx_p. The decision trees are made using
                           sklearn's DecisionTreeRegressor. Then a complete list of parameters can be found
                           at https://scikit-learn.org/stable/modules/generated/sklearn.tree.DecisionTreeRegressor.html.""")

    elnet_args = parser.add_argument_group("Elastic Net parameters")
    elnet_args.add_argument("--alpha", type=float, nargs='+', default=[0.999], metavar="", help="Alpha value controling l1 and l2 norm balance in ElasticNet.")
    elnet_args.add_argument("--e", type=float, default=1e-4,metavar="", help="Epsilon value, which is used to calculate the minimum lambda (i.e., min_lambda =  max_lambda * e).")
    elnet_args.add_argument("--K", type=int, default=100, metavar="",help="Number of lambda values to test between max_lambda and min_lambda. ")
    elnet_args.add_argument("--tol", type=float, default=0.00001,metavar="", help="Tolerance for convergence.")
    elnet_args.add_argument("--max_iter", type=int, default=1000,metavar="", help="Maximum number of iterations.")
    elnet_args.add_argument("--nwerror", action="store_true",  help="Not write test RMSE error for every lambda used in the path.")
    elnet_args.add_argument("--wpath", action="store_true",  help="Write ElasticNet path. Warning: This can be large.")
    elnet_args.add_argument("--nstdy", action="store_true", help="Not standardize y. Standarizing y helps to numerical stability. ")
    elnet_args.add_argument("--cv", action="store_true", help="Use cross-validation to select the best lambda and alpha value.")
    elnet_args.add_argument("--folds", type=int, default=5, metavar="", help="Number of folds for cross-validation when cv is True.")
    elnet_args.add_argument("--ncores", type=int, default=1, metavar="", help="Number of cores to use for cross-validation.")


    batch_args = parser.add_argument_group("Batch selection parameters")
    batch_args.add_argument("--prefix", type=str, default='batches', metavar="", help="Prefix for output files.")
    batch_args.add_argument("--factor", type=float, default=1/2, metavar="", help="Reduction factor for selecting overlapped batches. If -1 and nwerror is False, then the batch with the minimum RMSE is selected.")
    batch_args.add_argument("--inbetween", type=int, default=5, metavar="",help="Number of in-between batches for selecting overlapped batches.")
    batch_args.add_argument("--window", type=int, default=1,metavar="", help="Window size for agglomerating disjoint batches. With the current there is no agglomeration.")
    batch_args.add_argument("--check_spps", action="store_true", help="Check the number of unique species in the batches.")
    

    args = parser.parse_args()
    # print(args)

    # assert args.factor >= -1 and args.factor <= 1, "Factor must be between 0 and 1."
    assert args.inbetween >= 0, "Inbetween must be greater or equal to 0."
    assert args.window >= 1, "Window must be greater or equal to 1."
    assert args.p_test > 0 and args.p_test < 1, "Proportion of test samples must be between 0 and 1."
    assert args.e > 0, "Epsilon must be greater than 0."
    assert args.K > 0, "K must be greater than 0."
    assert args.tol > 0, "Tolerance must be greater than 0."
    assert args.max_iter > 0, "Maximum number of iterations must be greater than 0."
    assert args.seed >= 0, "Seed must be greater or equal to 0."
    assert args.M > 0, "Number of trees in the ensemble must be greater than 0."
    assert args.max_depth > 0, "Maximum depth of the decision tree must be greater than 0."
    assert (isinstance(args.max_features, float) and args.max_features > 0 and args.max_features <= 1) or args.max_features in ['sqrt', 'log2'], "Maximum proportion of features must be between 0 and 1 or 'sqrt' or 'log2'."
    assert args.eta > 0 and args.eta <= 1, "Proportion of observations to use in each tree must be between 0 and 1."
    assert args.nu >= 0 and args.nu <= 1, "Learning rate must be between 0 and 1."
    assert args.prefix != "", "Prefix must not be empty."
    # endregion: parse arguments

    # O(nT^4) for reading only
    data = np.loadtxt(args.Xy_file, delimiter=',', skiprows=1)
    X, y = data[:, :-1], data[:, -1]
    n,p = X.shape


    num_test = int(n*args.p_test)
    start = time.time()
    # O(nT^4) for centering or O(M T^4 nlog(n)) for ISLE
    (X_train,X_test,
     y_train,y_test,
     estimators # None if isle is False
     ) = split_data_isle(X, y,
            num_test=num_test, seed=args.seed,
            isle=args.isle, nwerror=args.nwerror, 
            mx_p=args.max_features, max_depth=args.max_depth, max_leaves=args.max_leaf_nodes,
            param_file=args.param_file, eta=args.eta, nu=args.nu, M=args.M,
            verbose=args.verbose, nstdy = args.nstdy)
    
    if args.isle:
        # standarize for ISLE post-processing uses lasso/elnet
        # and lasso/elnet assumes the data is standarized
        if args.verbose:
            print("Re-standarize data for ISLE")
        X_test, X_train = re_center_for_isle(X_test, X_train)


    all_max_lams = [max_lambda(X_train, y_train, alpha=a) for a in args.alpha]
    if args.verbose:
        print("Max lambda(s): ", all_max_lams)
    
    max_lam = np.max(all_max_lams)
    min_lam = max_lam * args.e
    params = {'lam': np.logspace(np.log10(min_lam), np.log10(max_lam), args.K, endpoint=True)[::-1]}

    assert all([a > 0 and a <= 1 for a in args.alpha]), "Alpha values must be between 0 and 1."
    if len(args.alpha) > 1:
        assert args.cv, "If alpha is a list, then cv must be True."

        # do cross-validation
        args.alpha = ElasticNetCV_alpha(args, X_train, y_train,
                                        args.alpha, params,
                                        args.folds, args.ncores)
        
        max_lam = max_lambda(X_train, y_train, alpha=args.alpha)
        min_lam = max_lam * args.e

        if args.verbose:
            print("Lambda range: ", min_lam, max_lam)

        params = {'lam': np.logspace(np.log10(min_lam), np.log10(max_lam), args.K, endpoint=True)[::-1]}

    else:
        assert not args.cv, "If alpha is a single value, then cv must be False."
        args.alpha = args.alpha[0]

    # O(npk) = O(nT^4) for fixed k 
    # or O(nM) if isle is True
    model = ElasticNet(fit_intercept=True, 
                        max_iter=args.max_iter,
                        init_iter=1, 
                        copyX=True, 
                        alpha=args.alpha, 
                        tol=args.tol, seed = args.seed)

    # O(npk) = O(nT^4) for fixed k
    # or O(Mk) if isle is True
    path = lasso_path(X_train, y_train, params, model)
    end_lasso = time.time() - start

    if args.verbose:
        print("Elastic net path done: ", end_lasso, " seconds")


    if args.wpath:
        lam_path = np.concatenate((params['lam'].reshape(-1, 1), path.T), axis=1)
        np.savetxt(args.prefix + "_elnetPath.csv", 
                   lam_path, 
                   delimiter=',',
                   comments='')
    # O(n)
    test_errors = calculate_test_errors(args, path, params, X_test, y_test)

    #O(T^4)
    CT = np.loadtxt(args.CT_file, delimiter=',', skiprows=1)
    CT_spps = CT[:, :4]
    n_spps = len(np.unique(CT_spps))

    if args.isle:
        picked_file = args.prefix + "_overlappedBatches_isle.txt"
        # batch_file = args.prefix + "_disjointBatches_isle.txt"
        
        # transform the path
        # O(kp) = O(kT^4) = O(T^4) for fixed k
        # it returns a list of sets
        path = get_new_path(estimators, path) 
        
    else:
        picked_file = args.prefix + "_overlappedBatches.txt"
        # batch_file = args.prefix + "_disjointBatches.txt"


    # overlapping batches
    # O(\rho T^4) for isle
    rows_selected = row_selection(path, CT_spps, test_errors, n_spps, 
                                   args.factor, args.inbetween, args.check_spps) 
    
    # O(\rho T^4) where \rho is the fraction of non-zero coefficients
    write_rows(picked_file, rows_selected)


    if not args.nwerror and args.verbose:
        # print(test_errors)
        min_err_sel = len(rows_selected[-1]) # the best is the last one
        print("Number of rows selected at min error: ", min_err_sel)

if __name__ == "__main__":
    main()

