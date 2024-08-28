#!/usr/bin/env python3

# placeholder for the model
import time
import numpy as np
import argparse


from qsin.sparse_solutions import ElasticNet, lasso_path
from qsin.utils import  calculate_test_errors, max_lambda, _scaler
from qsin.isle_path import split_data_isle, get_new_path
from qsin.ElasticNetCV import ElasticNetCV_alpha

def write_batches(outfile, batches):

    with open(outfile, 'w') as f:
        for batch in batches:
            f.write(",".join([str(b) for b in batch]) + "\n")

def make_batches(path, CT_spps, n_spps = 15):
    batches = []
    picked_idx = []
    for k in range(path.shape[1]):

        beta_k = path[:,k]
        idx_k = np.where(beta_k != 0)[0]

        if len(idx_k) == 0:
            continue

        if len(batches) == 0:
            batches.append(idx_k)
            picked_idx.extend(idx_k)
            continue


        tilde_idx_k = np.array(list(set(idx_k) - set(picked_idx)))
        if len(tilde_idx_k) == 0:
            continue

        idx_k_1 = batches[-1]
        spps_k_1 = len(np.unique(CT_spps[idx_k_1,:]))

        # tilde_spps_k = len(np.unique(CT_spps[tilde_idx_k,:]))
        
        if spps_k_1 < n_spps:
            batches[-1] = np.concatenate((batches[-1], tilde_idx_k))
            picked_idx.extend(tilde_idx_k)

        else:
            batches.append(tilde_idx_k)
            picked_idx.extend(tilde_idx_k)

    return [b+1 for b in batches]

# agglomerate batches
def agglomerate_batches(batches, window = 2):
    new_batches = []
    for i in range(0, len(batches), window):
        agglomerated = [b for b in batches[i:i+window] if len(b) > 0]
        agglomerated = np.concatenate(agglomerated)
        new_batches.append(agglomerated)
    return new_batches


def choose_j(path, test_errors = None, factor = 1/2):
    # factor = -1
    # test_errors = 1
    if factor == -1 and test_errors is not None:
        # tests_errors contains two columns
        # the first one is the lambda values
        # the second one is the RMSE values
        # check 'calculate_test_errors' function
        return np.argmin(test_errors[:,1])
    
    else:
        if factor < 0 or factor > 1:
            raise ValueError('Factor must be between 0 and 1 if nwerror is false and factor is not -1.')

        k = np.round(path.shape[0]*factor)
        non_zeros_len = []
        for j in range(path.shape[1]):
            beta_j = path[:,j]
            non_zeros_len.append( np.sum(beta_j != 0) )

        return np.argmin( np.abs(non_zeros_len - k ))


def select_path(path, CT_spps, test_errors, n_spps = 15, factor = 1/2, inbetween = 0):

    # k = np.round(path.shape[0]*factor)

    # non_zeros = []
    # for j in range(path.shape[1]):
    #     beta_j = path[:,j]
    #     non_zeros.append( np.sum(beta_j != 0) )
    
    # inbetween =  10
    # j_opt = np.argmin( np.abs(non_zeros - k ))

    j_opt = choose_j(path, test_errors, factor = factor)

    chosen_j = np.linspace(0, j_opt, 2 + inbetween, endpoint=True )
    chosen_j = chosen_j[chosen_j != 0]


    taken = []
    new_batches = []
    for j in chosen_j:
        # make sure it is the integer
        j_int = np.int64(np.round(j))


        # once it is integer,
        # it might be the case
        # that there are repeated j's
        if j_int in taken:
            continue
        
    
        beta_j_nz = np.where(path[:,j_int] != 0)[0]

        # check on the number of species
        if len(np.unique(CT_spps[beta_j_nz,:])) < n_spps:
            taken.append(j_int)
            continue

        new_batches.append(  beta_j_nz + 1)
        taken.append(j_int)
    
    return new_batches


def main():
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
    opt_args.add_argument("--version", action="version", version="%(prog)s 1.0")
    opt_args.add_argument("--verbose", action="store_true", help="Whether to print verbose output.")
    opt_args.add_argument("--isle", action="store_true", help="Whether to use path from decision tree-based ISLE (i.e., ensemble learning).")
    opt_args.add_argument("--p_test", type=float, default=0.35, metavar="", help="Proportion of observations to use for testing.")
    opt_args.add_argument("--seed", type=int, default=12038, metavar="", help="Random seed.")
    


    isle_args = parser.add_argument_group("ISLE parameters (if --isle is used)")
    isle_args.add_argument("--eta", type=float, default=0.5, metavar="", help="Proportion of observations to use in each tree.")
    isle_args.add_argument("--nu", type=float, default=0.1, metavar="", help="Learning rate.")
    isle_args.add_argument("--M", type=int, default=500, metavar="", help="Number of trees in the ensemble.")
    isle_args.add_argument("--max_depth", type=int, default=2, metavar="", help="Maximum depth of the decision tree.")
    isle_args.add_argument("--max_features", type=float, default=1/2, metavar="", help="Maximum proportion of features to use in each tree.")
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
    assert args.max_features > 0 and args.max_features <= 1, "Maximum proportion of features must be between 0 and 1."
    assert args.eta > 0 and args.eta <= 1, "Proportion of observations to use in each tree must be between 0 and 1."
    assert args.nu >= 0 and args.nu <= 1, "Learning rate must be between 0 and 1."
    assert args.prefix != "", "Prefix must not be empty."


    data = np.loadtxt(args.Xy_file, delimiter=',', skiprows=1)
    X, y = data[:, :-1], data[:, -1]
    n,p = X.shape

    X = _scaler(X, X, sted = True)
    y = _scaler(y, y, sted = False if args.nstdy else True)

    num_test = int(n*args.p_test)

    start = time.time()

    (X_train,X_test,
     y_train,y_test,
     estimators # None if isle is False
     ) = split_data_isle(X, y,
            num_test=num_test, seed=args.seed,
            isle=args.isle, nwerror=args.nwerror, 
            mx_p=args.max_features, max_depth=args.max_depth, param_file=args.param_file, 
            eta=args.eta, nu=args.nu, M=args.M,
            verbose=args.verbose)

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

    else:
        assert not args.cv, "If alpha is a single value, then cv must be False."
        args.alpha = args.alpha[0]
            
    # else:
    #     assert not args.cv, "If cv is true, then alpha must be a list."
    #     assert args.alpha > 0 and args.alpha <= 1, "Alpha must be between 0 and 1."

    model = ElasticNet(fit_intercept=False, 
                        max_iter=args.max_iter,
                        init_iter=1, 
                        copyX=True, 
                        alpha=args.alpha, 
                        tol=args.tol, seed = args.seed)


    path = lasso_path(X_train, y_train, params, model, sequential_screening=False)
    end_lasso = time.time() - start

    if args.verbose:
        print("Elastic net path done: ", end_lasso, " seconds")


    if args.wpath:
        lam_path = np.concatenate((params['lam'].reshape(-1, 1), path.T), axis=1)
        np.savetxt(args.prefix + "_elnetPath.csv", 
                   lam_path, 
                   delimiter=',',
                   comments='')

    test_errors = calculate_test_errors(args, path, params, X_test, y_test)


    CT = np.loadtxt(args.CT_file, delimiter=',', skiprows=1)
    CT_spps = CT[:, :4]
    n_spps = len(np.unique(CT_spps))

    if args.isle:
        picked_file = args.prefix + "_overlappedBatches_isle.txt"
        batch_file = args.prefix + "_disjointBatches_isle.txt"
        # transform the path
        path = get_new_path(estimators, path, p)
        
    else:
        picked_file = args.prefix + "_overlappedBatches.txt"
        batch_file = args.prefix + "_disjointBatches.txt"


    # overlapping batches
    picked_batches = select_path(path, CT_spps, test_errors, n_spps, args.factor, args.inbetween)
    write_batches(picked_file, picked_batches)

    # disjoint batches creation
    batches = make_batches(path, CT_spps, n_spps)

    if args.window <= 1:
        write_batches(batch_file, batches)

    else:
        new_batches = agglomerate_batches(batches, window=args.window)
        write_batches(batch_file, new_batches)


    if not args.nwerror and args.verbose:
        # print(test_errors)
        j_min = np.argmin(test_errors[:,1])
        new_path_j = path[:,j_min]
        min_err_sel = np.sum(new_path_j != 0)
        print("Number of rows selected at min error: ", min_err_sel)

if __name__ == "__main__":
    main()

