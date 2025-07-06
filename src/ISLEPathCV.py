import numpy as np
import multiprocessing as mp
from collections import deque

def error_fn(theta_j,
             base_model = None, i = None, 
             X_train_t = None, X_test_t = None, 
             y_train_t = None, y_test_t = None,
             verbose = True, seed = None):
    
    rng = np.random.RandomState(seed)
    nj, vj, lj, alpha = theta_j
    
    base_model.set_params(
        eta = nj,
        nu = vj,
        max_leaves = lj,
        alpha = alpha,
        rng = rng,
        verbose = False
    )
    base_model.fit(X_train_t, y_train_t)
    tmp_errors = base_model.score(X_test_t, y_test_t)

    if verbose:
        min_err = np.round(np.min(tmp_errors), 3)
        fold_message  = f"Fold: {i}, alpha: {alpha}"
        fold_message += f", eta: {nj}" if nj else ""
        fold_message += f", nu: {vj}" if vj else ""
        fold_message += f", leaves: {lj}" if lj else ""
        fold_message += f", err: {min_err}, seed: {seed}"
        print(fold_message)

    return (tmp_errors, theta_j)


def get_parallel_errors(base_model, X_train, y_train, full_grid, 
                           rng, num_folds, ncores, verbose = False):
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
    all_index = list(range(n))
    rng.shuffle(all_index)

    X_train = X_train[all_index, :] # check to shuffle X
    y_train = y_train[all_index] # check to shuffle y!
    

    out = deque([])
    with mp.Pool( processes = ncores ) as pool:

        preout = deque([])
        for i in range(num_folds):

            test_idx = list(range(i * fold_size, (i + 1) * fold_size))
            train_idx = list(set(range(n)) - set(test_idx))

            X_train_t, X_test_t = X_train[train_idx, :], X_train[test_idx, :]
            y_train_t, y_test_t = y_train[train_idx], y_train[test_idx]

            for theta_j in full_grid:
                tmp_seed = rng.randint(0, 2**31 - 1)

                errors = pool.apply_async(
                    error_fn, 
                    (theta_j, base_model, i, 
                     X_train_t, X_test_t, 
                     y_train_t, y_test_t, 
                     verbose, tmp_seed)
                )
                preout.append(errors)

        for errors in preout:
            out.append(errors.get())

    return list(out)

def get_best_params(all_errors, lambdas, folds = 5):
    # all_errors = out
    """
    Fold_alpha has the following structure:
    [  theta_{1,j}, ..., theta_{f,j}  ] 

    Fold error has the following structure:
    [ [ RMSE_1,j \in R^{1 x K} ]   -> theta_{1,j}
       ...
      [ RMSE_f,j \in R^{1 x K} ] ]  -> theta_{f,j}

    where f is the fold index and j is the hyperparameter index.
    """
    # getting row-wise average of the errors
    dict_errs = {}
    for (e_fj, theta_j) in all_errors:
        if isinstance(theta_j, list):
            # make it hashable if it is a list
            theta_j = tuple(theta_j)  

        if theta_j not in dict_errs:
            dict_errs[theta_j] = e_fj/folds

        else:
            dict_errs[theta_j] += e_fj/folds

    best_theta_j = 0
    best_lam = 0
    min_rmse = np.inf

    for theta_j, ave_e_j in dict_errs.items():
        tmp_cv_err = np.min(ave_e_j)

        if tmp_cv_err < min_rmse:
            min_rmse = tmp_cv_err
            best_theta_j = theta_j
            best_lam = lambdas[np.argmin(ave_e_j)]
            # print(best_theta_j, best_lam, min_rmse)

    return best_theta_j, best_lam, min_rmse

def ISLEPathCV(base_model, X_train, y_train, full_grid, 
               lambdas, folds, ncores, verbose = True, rng = None):
    """
    Find the best set of hyperparameters for the ISLEPath
    using a cross-validation. The function returns
    the best hyperparameter set.

    Higltights: it parallelizes over all the folds and 
    hyperparameter values.
    """

    if len(full_grid) == 1:
        return full_grid[0]

    if verbose:
        print("Hyperparameter grid size: ", len(full_grid))
        print("Performing CV with ", folds, " folds")
        
    all_errors = get_parallel_errors(
        base_model, X_train, y_train, full_grid,
        rng, folds, ncores, verbose
    )

    (best_theta_j,
     best_lam,
     min_rmse) = get_best_params(all_errors, lambdas, folds=folds)
    
    if verbose:
        print("CV best (eta, nu, leaves, alpha): ", best_theta_j)
        print("CV best lambda: ", best_lam)
        print("CV min average RMSE: ", min_rmse)

    return best_theta_j


def create_tree_param_grid(eta, nu, leaves, cv_sample = 100, rng = None):
    """
    return an array with combinations of hyperparameters
    for the ISLE ensemble. If cv_sample is lower than all possible
    combinations of hyperparameters, a random sample of size cv_sample
    is taken
    """
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


def create_full_grid(isle, eta, nu, leaves, alphas, cv_sample, rng):

    assert all([a > 0 and a <= 1 for a in alphas]), "Alpha values must be between 0 and 1."

    if not isle:
        return [[None, None, None, a] for a in alphas]
    
    grid = create_tree_param_grid(eta, nu, leaves, cv_sample, rng)
    
    out_grid  = deque([])
    for alpha_j in alphas:
        for tree_parms_j in grid:
            nj, vj, lj = list(tree_parms_j)
            out_grid.append( [nj, vj, int(lj), alpha_j] )

    return out_grid
