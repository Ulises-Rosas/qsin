import argparse
from copy import deepcopy
from collections import deque
import math


from netio import is_leaf, get_cycles, is_HybLeaf, writeNewick, parseNewickNet

# true_net = "((B,#H1),((D)#H1,C),A);"
# infered_net = "((C,#H1),((D)#H1,B),A);"


# region: get hybrid cycle with subtree leaves
def merge_queues(p_split, n_lr):

    # case when an internal node does not have
    # a left or right child. This happens (left case)
    # in child hybrid nodes.
    if not n_lr:
        return deque()
    
    # this not a real taxon
    if is_HybLeaf(n_lr):
        return deque()

    # append the leaf to the queue
    if is_leaf(n_lr):
        return deque([n_lr.name])
    
    # if not a leaf, then this node
    # has already been traversed
    # and the leaves are already at p_split
    l_leaves, r_leaves = deepcopy(p_split[n_lr])
    l_leaves.extend(r_leaves) # O(n)
    return l_leaves

def dfs_s(n, p_split):
    """
    DFS traversal of the tree to find the leaves in each subtree.

    For each node, we traverse its left and right children
    and store their leaves in a dictionary p_split.
    For each node we merge the leaves of its left and right children
    and store them in p_split.
    In the worst case, the merging costs the total number of leaves
    on the subtree. Lets call this number x_i. If let n_i be the number
    of nodes in the subtree, then the total cost is O(x_i*n_i). Since
    n_i is a function of x_i, we can say that the cost is O(n_i^2).

    If the tree is balanced, then the cost is O(n_i*log(n_i)).
    """
    # meomization
    # this ensure every node is traversed only once
    # this is needed as in each cycle it is possible
    # to look for the same node
    if n in p_split: 
        return None

    if not n:
        return None
    
    # O(n)
    if is_HybLeaf(n) or is_leaf(n):
        return None
    
    # stops at (hyb) leaf, it returns None
    # and goes into the previous call in the stack
    # this is the base case. The previous call
    # is an internal node.
    dfs_s(n.left , p_split) 
    dfs_s(n.right, p_split)

    # merge the queues
    # n here is an internal node
    left_leaves  = merge_queues(p_split, n.left ) # O(n)
    right_leaves = merge_queues(p_split, n.right) # O(n)

    p_split[n] = [left_leaves, right_leaves]

def sub_leaves(n, p_split):
    """
    Traversal of the tree to find the leaves in each subtree.
    given a node n
    """
    if n in p_split:
        return None

    if is_leaf(n):
        p_split[n] = [deque([n.name]), deque()]
    
    else:    
        # this only works for internal node
        # then if n is leaf nothing it is done
        # to p_split
        dfs_s(n, p_split)

def get_cycles_taxa(nodes, hyb_nodes, root):
    """
    Get the taxa of each cycle in the tree.
    The taxa are stored in a dictionary where the key is the
    hybrid node and the value is a dictionary with the taxa
    of the cycle and the node where they are located.
    The taxa are sorted and stored as a tuple of strings.

    All nodes are traversed once, and the leaves of each
    subtree are computed once. On each node visit
    leaves are appended to the queue of the parent node.
    The cost of this is O(n) for each node.
    The total cost is O(n^2) in the worst case.
    If the tree is balanced, including the cycles!,
    then this is O(n*log(n)).

    Parameters
    ----------
    hyb_nodes : dict
        dict of hybrid nodes

    root : Node
        Root of the tree.
    Returns
    -------

    all_cycle_taxa : dict
        Dictionary with the taxa of each cycle. For each cycle 
        nodes it is presented a list of sorted taxa names that has
        this node as ancestor.
        e.g. {hyb_node: [[taxa1, taxa2, ...], n]}, where n is the 
        ancestral node of the taxa1, taxa2, ...
        The order of the cycle nodes starts from either left or right from the cycle root and 
        then goes all the way to the root of the cycle.
    """
    all_taxa = set([n.name for n in nodes if is_leaf(n)])

    blobs = get_cycles(hyb_nodes, root)
    all_cycle_taxa = {}
    # shared for all nodes
    p_split = {}
    for hyb_label, (V_q, cr) in blobs.items():
        
        # (V_q, cr) = blobs["#H18"]
        
        # set of cycle nodes, order is lost
        # v_q the order is kept, anticlockwise.
        V_c = set(V_q)

        # taxa of the cycle per node
        cycle_taxa_n = deque()
        
        # meomization is used to avoid
        # recomputing the leaves of each subtree

        # in the worst case scenario, the first
        # cycle is close to the root, and 
        # all the tree is traversed when
        # looking for the leaves of the cycle
        # this is O(n^2) in the worst case. But this 
        # is done only once as results are stored.
        # if the tree is balanced, including the cycles!,
        # then this is O(n*log(n))

        non_root_taxa = set()

        for n in V_q:
            # print(n)
            if n == cr:
                cycle_taxa_n.append([None, n])
                continue
                # continue
        
            # hybrid do not have
            # subtrees.
            if is_HybLeaf(n):
                cycle_taxa_n.append([None, n])
                continue
            
            # get subtrees leaves
            outward_n = n.right if n.left in V_c else n.left
            sub_leaves(outward_n, p_split) # O(1) to O(n^2)
            # print(p_split[outward_n])
            left, right = deepcopy(p_split[outward_n])
            left.extend(right) # right node can be empty

            cycle_taxa_n.append([left, n])
            non_root_taxa |= set(left)

        
        cycle_taxa_n = list(cycle_taxa_n)

        for i,(taxa,n) in enumerate(cycle_taxa_n):

            if (not taxa) and (not n.hybrid):
                root_taxa = all_taxa - non_root_taxa
                cycle_taxa_n[i] = [root_taxa, n]

        all_cycle_taxa[hyb_label] = cycle_taxa_n


    for hyb_label, info in all_cycle_taxa.items():

        tmp_taxa_n = deque()
        for taxa, n in info:
            # print(taxa)
            if taxa is not None:
                taxa = tuple(sorted(taxa)) # O(n_i*log(n_i))

            tmp_taxa_n.append([taxa, n])

        all_cycle_taxa[hyb_label] = list(tmp_taxa_n)

    return all_cycle_taxa
# endregion: get hybrid cycle with subtree leaves


# region: change direction of hybrid node
def flip_pointers(n, V_c):
    """
    Flip the pointers of the node n and its children
    in the cycle V_c.
    The pointers are flipped in the following way:
    - If the left/right child of n is in V_c, then
      the left/right child of n is set to the ancestor of n.
    - Ancestor of n are set to the left/right child of n.
    """

    n_anc = n.ancestor
    n_lr = None

    # pointers change should be 
    # inside the cycle
    if n.left in V_c:
        n_lr   = n.left
        n.left = n_anc

    else:
        # check which right
        if isinstance(n.right, list):

            if n.right[0] in V_c:
                n_lr = n.right[0]
                n.right[0] = n_anc

            else:
                n_lr = n.right[1]
                n.right[1] = n_anc
        
        else:
            n_lr = n.right
            n.right = n_anc

    n.ancestor = n_lr


def node_replacement(u, v_old, v_new):

    # pointers change should be 
    # inside the cycle
    if u.left == v_old:
        u.left = v_new
        
    else:
        if isinstance(u.right, list):
            # check which right is the 
            # the v_old
            k = 0 if u.right[0] == v_old else 1
            u.right[k] = v_new

        else:
            u.right = v_new

def dismount_hybrid(hc):
    gamma = hc.gamma
    label = hc.label

    # hc.ancestor = hc.ancestor[0]
    hc.hybrid = False
    hc.isChild = False
    hc.gamma = 1
    hc.label = None
    return gamma, label

# find descendant in cycle
def incycle_descendant(u, cycle):
    """
    find descendant in cycle V_c
    """
    if u.left in cycle:
        return u.left
    else:
        if isinstance(u.right, list):
            k = 0 if u.right[0] in cycle else 1
            return u.right[k]
        
        else:
            return u.right
        
# find descendant outside the cycle
def outcycle_descendant(u, cycle):
    """
    find descendant not in cycle V_c
    """
    if u.left not in cycle:
        return u.left
    else:
        if isinstance(u.right, list):
            k = 1 if u.right[0] in cycle else 0
            return u.right[k]
        
        else:
            return u.right

# add attributes to hybrid node
def mount_hybrid(n, gamma, label):
    n.hybrid = True
    n.isChild = True
    n.gamma = gamma
    n.label = label


def changeHybrid(hc_old, hc_new, V_c, is_anc):
    """
    Change the direction of the hybrid node hc
    to the target node n_target.
    The direction is changed by flipping the pointers.

    parameters
    ----------
    hc_old : myNode
        hybrid child node. That is the node that
        has two ancesotrs and a left descendant

    hc_new : myNode
        target node. 
        This is the node that will be
        the new child hybrid node.

    V_c : set
        set of nodes in the cycle. This is only
        used to assess presence or absence of a node in the cycle.

    is_anc : bool
        True if the target node is an ancestor of the hybrid child node
        False if the target node is a descendant of the hybrid child node
        in the cycle.

    Returns
    -------
    None
    """
    
    # this is true for level 1 networks
    # the second ancestor is the hybrid
    # leaf node. The first one is the tree
    # based
    hl = hc_old.ancestor[1]

    # dismount child hybrid node
    gamma, label = dismount_hybrid(hc_old)


    # ancestor of leaf hybrid node
    X = hl.ancestor
    # get new rights and 
    # ancestors. Recall that child hybrid nodes
    # do not have right and have two ancestors
    # print("is_anc: ", is_anc)
    # print("hc_old: ", hc_old)
    # print("hc_old.ancestor: ", hc_old.ancestor)

    if is_anc:
        new_right = hc_old.ancestor[0]
        new_anc =  X
        
    else:
        new_right = X
        new_anc = hc_old.ancestor[0]

    node_replacement(u = X, v_old=hl, v_new=hc_old)
    hc_old.right = new_right
    hc_old.ancestor = new_anc

    n_curr = new_right
    # print()
    while n_curr != hc_new:
        # print("n_curr: ", n_curr)
        n_curr_anc = n_curr.ancestor
        flip_pointers(n_curr, V_c)
        n_curr = n_curr_anc

    # find descendant that belong 
    # to cycle node set
    n_prev = incycle_descendant(u = hc_new, cycle = V_c)

    # make hybrid leaf to be a descendant
    # of n_prev
    node_replacement(u = n_prev, v_old = hc_new, v_new = hl)
    hl.ancestor = n_prev

    # update attributes of n_target
    mount_hybrid(hc_new, gamma, label)

    # step 1: replace left of hybrid leaf 
    # by targe node
    hl.left = hc_new 
    # step 2: create new ancestors
    hc_new.ancestor = [hc_new.ancestor, hl]

    # step3: make left descendant only to the target node
    new_left = outcycle_descendant(u = hc_new, cycle = V_c)
    hc_new.left = new_left
    hc_new.right = None

def flatten_taxa_info(cycles_info):
    """
    Get taxa hash from the cycles info (flat it out).
    This essentially creates data structures that are
    useful for the changeHybrid function.
    O(n) algorithm.

    parameters
    ----------
    cycles_info : dict
        dictionary of cycles info. More description on 
        the structure of this dictionary can be found
        in the docstring of the `get_cycles_taxa` function.

    Returns
    -------
    
    taxa_hash : dict
        dictionary of taxa hash. The keys are the 
        hashable taxa and the values [H, n]. Where H is the
        hybrid identifier and n is a node in the H cycle.
    """
    taxa_hash = {}
    
    # O(n_1 + n_2 + ... + n_H) = O(n),
    # where n is total number of nodes in all cycles
    for H, (taxa_nodes) in cycles_info.items():
        for (taxa,n) in taxa_nodes: # O(n_i)
            if taxa:
                taxa_hash[taxa] = [H, n]

    return taxa_hash

# find child hybrid nodes in cycles
def get_Child_taxa(cycles_info):
    """
    find sub leaves of child hybrid node in cycles.
    """
    child_nodes = deque()
    # O(n_1 + n_2 + ... + n_H) = O(n),
    for H, (taxa_nodes) in cycles_info.items():
        for (taxa,n) in taxa_nodes: # O(n_i)
            if n.isChild:
                child_nodes.append(taxa)
    return child_nodes

def is_target_ancestor(ordered_cycle, hc, target_n):
    """
    Check whether target node is ancestor of hc in the ordered cycle.
    O(n) algorithm.

    parameters
    ----------
    ordered_cycle : deque
        ordered cycle of nodes. This ordering starts
        from an adjacent node of the root node and 
        goes all the way to the hybrid child node always
        and then back to the root node.

    hc : node
        hybrid child node.

    target_n : node
        target node.

    Returns
    -------
    bool
        True if target node is ancestor of hc, False otherwise.
    """
    cycle_list = list(ordered_cycle)
    is_anc = False

    # find whether target node is ancestor of hc
    k = 0
    n = cycle_list[k]
    # go from n_x to hc
    while n != hc:
        # if in between there is a target node
        # then it is an ancestor of hc
        if n == target_n:
            is_anc = True
            break

        k += 1
        n = cycle_list[k]

    return is_anc

def is_target_ancestor_v2(hc_new, hc_old, cr):
    """
    Check whether target node is ancestor of hc in the ordered cycle.

    parameters
    ----------
    hc_new : node
        hybrid child node.
    hc_old : node   
        hybrid child node.

    cr : node
        root of the cycle.

    Returns
    -------
    bool
        True if target node is ancestor of hc, False otherwise.
    """


    # tree-based ancestor
    n = hc_old.ancestor[0]
    # go back until the root of the cycle
    # and if you did not find hc_new
    # it means that it is not an ancestor
    is_anc = False
    while n != cr:
        if n == hc_new:
            is_anc = True
            break

        n = n.ancestor

    return is_anc



# endregion: change direction of hybrid node

def split_symdiff(true_cycle_taxa, tmp_cycle_taxa):
    # then check if they have the same splits
    true_splits = set()
    for splits_nodes in true_cycle_taxa.values():
        for (split, _) in splits_nodes:
            if split is None:
                continue
            true_splits.add(split)
        
    mod_splits = set()
    for splits_nodes in tmp_cycle_taxa.values():
        for (split, _) in splits_nodes:
            if split is None:
                continue
            mod_splits.add(split)

    return (true_splits | mod_splits) - (true_splits & mod_splits) # symmetric difference

def changeDirection(true_cycle_taxa, true_child_taxa, tmp_cycle_taxa, 
                    check_all_node_splits = True,
                    verbose=False):
    """
    Change the direction of the hybridization events in the tmp network
    to match the true network. O(n) algorithm

    Parameters
    ----------

    true_cycle_taxa : dict
        dictionary containing information about the cycles in the true network.
        This dictionary comes from the function get_cycles_taxa.


    true_child_taxa : list
        list of taxa that are children of hybrid child node in the true network.
        This list comes from the function get_Child_taxa.

    tmp_cycle_taxa : dict
        dictionary containing information about the cycles in the tmp network.
        This dictionary comes from the function get_cycles_taxa.

    check_all_node_splits : bool, optional
        if True, check all node splits in cycles.
        This will allow some nodes in the cycle
        have different decendants, including the hybrid node.
        When the latter happens,
        the corresponding cycle does not change direction.
        The default is True.

    verbose : bool, optional
        if True, print the differences between the two networks.
        The default is False.
    
    Returns
    -------
    bool
        True if the direction of the hybridization events in the tmp network
    """

    # true_net = get_cycles_taxa(true_hyb_nodes, true_root)
    # # taxa of target nodes
    # hcs = get_Child_taxa(true_net)
    # tmp_cycle_taxa = get_cycles_taxa(tmp_hyb_nodes, tmp_root)

    # check if they both have the same number of hybridization events
    if len(true_cycle_taxa) != len(tmp_cycle_taxa):
        if verbose:
            print("The networks have different number of hybridization events")

        return False

    # then check if they cycles nodes of true net and
    # tmp net have the same set descendants
    if check_all_node_splits:
        new_set = split_symdiff(true_cycle_taxa, tmp_cycle_taxa)

        if len(new_set) != 0:
            if verbose:
                print("The networks have different splits")
                print("The splits are: ", new_set)
            return False

    # if they have do, then they have the same nodes in
    # the undirected cycles
    # we change the direction of the cycle to match
    # the same as the true network

    # taxa hash from tmp network
    tmp_flat_taxa_hash = flatten_taxa_info(tmp_cycle_taxa)
    at_least_one_changed = False
    for hc_new_taxa in true_child_taxa:
        # hc_taxa
        # get the target node
        if hc_new_taxa not in tmp_flat_taxa_hash:
            continue
        
        # this is always correct if both cycles 
        # have same descedants
        H, hc_new = tmp_flat_taxa_hash[hc_new_taxa]
        # if H in updated_cycles:
        #     # this means that the cycle has already been updated.
        #     continue
        cycle_taxa_n = tmp_cycle_taxa[H]

        cycle_q = deque()
        hc_old = None
        hc_old_taxa = None
        # O(n_i)
        for taxa,n in cycle_taxa_n:
            if n.isChild:
                hc_old = n
                hc_old_taxa = taxa

            cycle_q.append(n)

        if hc_old_taxa == hc_new_taxa:
            # this means that it already
            # has the right subleaves.
            # no need to change direction
            continue

        if verbose:
            print(f"Changing direction of {hc_new_taxa}")
            print("cycle_taxa_n: ", cycle_taxa_n)
            print("hc_old: ", hc_old)
            print("hc_new: ", hc_new)

        # root cycle
        # this is the last node in the cycle
        cr = cycle_taxa_n[-1][1]
        if cr == hc_new:
            if verbose:
                print("Attempting to make the root of the cycle to be the new hybrid child node, we skip this")
            continue

        is_anc = is_target_ancestor_v2(hc_new, hc_old, cr)
        changeHybrid(hc_old = hc_old, # from
                     hc_new = hc_new,  # to
                     V_c    = set(cycle_q), # Set of a cycle. Used to fast find only
                     is_anc = is_anc)
        
        at_least_one_changed = True
        
        if verbose: 
            print("new cycle: ", cycle_q)
            print("new hc: ", hc_new)
            print("new hc.ancestor: ", hc_new.ancestor)
            print()
    
    return at_least_one_changed # success in changing the direction, at least to 1 hyb node


# region: changeDirection iterator
def changeDirectionIterator(true_nets_file, infered_nets_files,
                            check_all_node_splits = True,
                            verbose = False,
                             suffix = '_checked'):

    all_true_nets_newick = []
    with open(true_nets_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            all_true_nets_newick.append(line.strip())


    infered_nets_newick = []
    infered_nets_name = []
    for file in infered_nets_files:
        with open(file, 'r') as f:
            infered_nets_newick.append(f.readline().strip())
            infered_nets_name.append(file)


    all_true_cycle_child_taxa = deque()
    for true_net_newick in all_true_nets_newick:

        true_nodes, true_hyb_nodes, true_root = parseNewickNet(true_net_newick)
        # cycle taxa
        true_cycle_taxa = get_cycles_taxa(true_nodes, true_hyb_nodes, true_root) # O(n^2)
        # child taxa
        true_child_taxa = get_Child_taxa(true_cycle_taxa) # O(n)

        all_true_cycle_child_taxa.append(
            (true_cycle_taxa,true_child_taxa)
        )


    for i, tmp_newick in enumerate(infered_nets_newick):
        file_name = infered_nets_name[i]

        if verbose:
            print(f"Processing file {file_name}")

        best_str = ''
        best_hyb_num = -math.inf
        
        for j,(true_cycle_taxa,true_child_taxa) in enumerate(all_true_cycle_child_taxa): # O(1)
            # (true_cycle_taxa,true_child_taxa)
            tmp_nodes, tmp_hyb_nodes, tmp_root = parseNewickNet(tmp_newick)
            tmp_cycle_taxa = get_cycles_taxa(tmp_nodes, tmp_hyb_nodes, tmp_root) # O(n^2)

            if verbose:
                print(f"true net net {j+1}")


            is_changed = changeDirection(true_cycle_taxa, true_child_taxa, tmp_cycle_taxa,
                                         check_all_node_splits=check_all_node_splits, 
                                         verbose=verbose)
            
            if is_changed:
                if len(true_cycle_taxa) > best_hyb_num:
                    
                    try:
                        best_str = writeNewick(tmp_nodes, tmp_root.index)
                        best_hyb_num = len(true_cycle_taxa)
                    except RecursionError:
                        # this can happen when there are loops between
                        # pointers. One reason for this is the
                        # new hybrid node pointed to what was previously
                        # a root of the cycle, creating a loop.
                        print("Error on node directions")
                        pass


        if best_str:
            if verbose:
                print(f"Best hybridization number: {best_hyb_num}")
                print()

            with open(file_name.replace('.txt', f'{suffix}.txt'), 'w') as f:
                f.write(best_str)

        else:
            with open(file_name.replace('.txt', f'{suffix}.txt'), 'w') as f:
                f.write(tmp_newick)
# endregion: changeDirection iterator



def arg_parse():
    parser = argparse.ArgumentParser(description='Change direction of hybridization events in the tmp network to match the true network.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('true_nets_file', type=str, help='File with the true networks in newick format')
    parser.add_argument('infered_nets_files', type=str, nargs='+', help='Files with the inferred networks in newick format')
    parser.add_argument('-c', '--check_all_node_splits', action='store_true', help='''
                        Check all node splits in cycles. 
                        Not using this allow some nodes in the cycle 
                        have different decendants, including the hybrid node.
                        When the latter happens, the corresponding cycle does not change direction.
                        A tmp net passing this checking makes it more likely that
                        to be closer to the true network.
                        Otherwise, the tmp net can be very different and 
                        more likely to wierd things to happen. For instance, the new hybrid node can be 
                        located into the root of the cycle, creating a pointer loop when transversing the network/tree.
                        ''')
    parser.add_argument('-s','--suffix', type=str, default='_checked', help='Suffix to add to the output files')
    parser.add_argument('-v','--verbose', action='store_true', help='Verbose output')
    return parser.parse_args()


if __name__ == "__main__":
    args = arg_parse()
    true_nets_file = args.true_nets_file
    infered_nets_files = args.infered_nets_files
    suffix = args.suffix
    check_all_node_splits = args.check_all_node_splits
    verbose = args.verbose

    # import glob
    # true_nets_file = "/Users/ulisesrosas/Desktop/qsin/test_data/n15_hyb3.txt"
    # infered_nets_files_str = "/Users/ulisesrosas/Desktop/experiments_qsin/eta0.3_nu0.1_leaves6_models4000_hd_optBL/n15_from_*_uniq_row6_boot*_nets.txt"
    # infered_nets_files = glob.glob(infered_nets_files_str)
    # suffix = "_checked"
    # check_all_node_splits = False
    # verbose = False

    changeDirectionIterator(true_nets_file, infered_nets_files, 
                            check_all_node_splits, verbose, suffix)
