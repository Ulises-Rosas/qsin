#!/usr/bin/env python3

import time
import argparse

import numpy as np

from qsin.isle_path import get_new_path, ISLEPath
from qsin.ISLEPathCV import create_full_grid, ISLEPathCV
from qsin.row_selection import row_selection, write_rows
from qsin.utils import (standardize_Xy, write_test_errors, 
                        split_data)

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
    isle_args.add_argument("--eta", type=float, nargs='+', default=[0.5], metavar="", help="Proportion of observations to use in each tree.")
    isle_args.add_argument("--nu", type=float, nargs='+', default=[0.1], metavar="", help="Learning rate.")
    isle_args.add_argument("--max_leaf_nodes", type=int, nargs='+' , default=[6], metavar="", help="Maximum number of leaf nodes in the decision tree.")
    isle_args.add_argument("--cv_sample", type=int, default=100, metavar="", help="Number of randomly choose ISLE hyperparameter combinations")
    isle_args.add_argument("--max_depth", type=int, default=5, metavar="", help="Maximum depth of the decision tree.")
    isle_args.add_argument("--M", type=int, default=500, metavar="", help="Number of trees in the ensemble.")

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
    elnet_args.add_argument("--wpath", action="store_true",  help="Write ElasticNet path. Warning: This can be large.")
    elnet_args.add_argument("--nstdy", action="store_true", help="Not standardize y. Standarizing y helps to numerical stability. ")
    elnet_args.add_argument("--cv", action="store_true", help="Use cross-validation to select the best lambda and alpha value.")
    elnet_args.add_argument("--folds", type=int, default=5, metavar="", help="Number of folds for cross-validation when cv is True.")
    elnet_args.add_argument("--ncores", type=int, default=1, metavar="", help="Number of cores to use for cross-validation.")


    batch_args = parser.add_argument_group("Row selection parameters (After the model is trained)")
    batch_args.add_argument("--prefix", type=str, default='batches', metavar="", help="Prefix for output files.")
    batch_args.add_argument("--factor", type=float, default=-1, metavar="", help="Proportion of row selection. If factor is -1, the proportion of row selected is based on test error.")
    batch_args.add_argument("--inbetween", type=int, default=5, metavar="",help="Number of in-between sets from the solution of the beginning of the path to the final row selection.")
    batch_args.add_argument("--check_spps", action="store_true", help="Check the number of unique species in the selected rows.")
    batch_args.add_argument("--ignore_conv", action="store_true", help="Ignore convergence of test errors. If True, the function will ignore the convergence of the test errors and return the index of the minimum test error.")
    batch_args.add_argument("--tol_test", type=float, default=1e-4, metavar="", help="""
                            Tolerance for convergence of test errors for all k solutions in the ElNet path to pick the best lambda_k.
                            If the test errors are not converging, then lambda_k with the minimum test error is picked.""")

    args = parser.parse_args()
    # print(args)

    # assert args.factor >= -1 and args.factor <= 1, "Factor must be between 0 and 1."
    assert args.inbetween >= 0, "Inbetween must be greater or equal to 0."
    assert args.p_test > 0 and args.p_test < 1, "Proportion of test samples must be between 0 and 1."
    assert args.e > 0, "Epsilon must be greater than 0."
    assert args.K > 0, "K must be greater than 0."
    assert args.tol > 0, "Tolerance must be greater than 0."
    assert args.max_iter > 0, "Maximum number of iterations must be greater than 0."
    assert args.seed >= 0, "Seed must be greater or equal to 0."
    assert args.M > 0, "Number of trees in the ensemble must be greater than 0."
    assert args.max_depth > 0, "Maximum depth of the decision tree must be greater than 0."
    assert (isinstance(args.max_features, float) and args.max_features > 0 and args.max_features <= 1) or args.max_features in ['sqrt', 'log2'], "Maximum proportion of features must be between 0 and 1 or 'sqrt' or 'log2'."
    assert args.prefix != "", "Prefix must not be empty."
    # endregion: parse arguments

    # O(nT^4) for reading only
    data = np.loadtxt(args.Xy_file, delimiter=',', skiprows=1)
    X, y = data[:, :-1], data[:, -1]
    n,p = X.shape # p (\approx T^4) used for row selection only
    rng = np.random.RandomState(seed=args.seed)


    num_test = int(n*args.p_test)
    start = time.time()
    # O(M T^4 nlog(n)) for ISLE ensemble

    # O(nT^4) for centering
    X_train, X_test, y_train, y_test = split_data(X, y, num_test=num_test, seed=args.seed)

    # improve stability of the algorithm
    X_train, X_test, y_train, y_test = standardize_Xy(X_train, y_train, X_test, y_test, args.nstdy)

    base_model = ISLEPath(
            # elastic net parameters
            fit_intercept = True,
            max_iter = args.max_iter,
            tol = args.tol,
            zero_thresh = 1e-12,  # threshold for zero coefficients
            # path parameters
            epsilon = args.e,
            K = args.K,
            # tree parameters
            M = args.M,
            max_features = args.max_features,
            max_depth = args.max_depth,
            param_file = args.param_file,
            make_ensemble = args.isle,
            verbose = False,
        )
    
    full_grid = create_full_grid(isle=args.isle, eta=args.eta, nu=args.nu,
                                  leaves=args.max_leaf_nodes, alphas=args.alpha, 
                                  cv_sample=args.cv_sample, rng=rng)


    (nj, vj, lj, alpha) = ISLEPathCV(base_model, X_train, y_train, full_grid, 
                                     args.folds, args.ncores, 
                                     verbose = args.verbose, rng=rng, tmp_seed=args.seed)

    rng = np.random.RandomState(seed=args.seed) # reinitialize rng for the ensemble
    base_model.set_params(eta=nj, nu=vj, max_leaves=lj, alpha=alpha, 
                          verbose=args.verbose, rng=rng)

    # O(npk) = O(nT^4) for fixed k 
    # or O(nM) if isle is True
    base_model.fit(X_train, y_train)
    end_lasso = time.time() - start
    
    path = base_model.path
    lambdas = base_model.lambdas

    if args.verbose:
        print("ISLE path done: ", end_lasso, " seconds")


    if args.wpath:
        # the first column is the lambda values
        # which adds one extra column into the path
        lam_path = np.concatenate((lambdas.reshape(-1, 1), path.T), axis=1)
        np.savetxt(args.prefix + "_elnetPath.csv", 
                   lam_path, 
                   delimiter=',',
                   comments='')
    # O(n)
    # write the test errors as well.
    test_rmse = base_model.score(X_test, y_test)
    if args.verbose:
        print("Min test RMSE: ", np.min(test_rmse), "at lambda: ",
              lambdas[np.argmin(test_rmse)])

    # print(base_model.intercepts)
    # sel = base_model.path[:,np.argmin(test_rmse)]
    # print(sel)
    # print(np.sum(sel != 0), "non-zero coefficients at the best lambda_k solution.")

    # test errors contain the first column with lambda values
    # and the second column with RMSE values.
    # this is compatible with the row selection. 
    # TODO: streamline from test_rmse
    test_errors = write_test_errors(args.prefix, test_rmse, lambdas)

    if args.isle:
        picked_file = args.prefix + "_overlappedBatches_isle.txt"
        
        # transform the path
        # O(kp) = O(kT^4) = O(T^4) for fixed k
        # it returns a list of sets
        path = get_new_path(base_model.estimators, path) 
        
    else:
        picked_file = args.prefix + "_overlappedBatches.txt"

    # overlapping batches
    # O(\rho T^4) for isle. If check_spps is True, it is O(T^4)
    rows_selected, j_opt = row_selection(path, test_errors, p,
                                  args.factor, args.inbetween, args.check_spps,
                                  args.CT_file, args.verbose, args.tol_test, 
                                  args.ignore_conv)
    
    if args.verbose:
        print(f" at lambda: {lambdas[j_opt]}.")
    # O(\rho T^4) where \rho is the fraction of non-zero coefficients
    write_rows(picked_file, rows_selected)


    if args.verbose:
        # print(test_errors)
        min_err_sel = len(rows_selected[-1]) # the best is the last one
        print("Number of rows selected at lambda_k solution: ", min_err_sel)

if __name__ == "__main__":
    main()

