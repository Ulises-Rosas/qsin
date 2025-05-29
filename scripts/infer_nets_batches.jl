#!/usr/bin/env julia


# make default values for the arguments
startfile = "";
buckyCFfile = "";
batches = "";

h_max = 1;
n_epochs = 1;
tk = Inf;
seed = 12038;
nruns = 10;
Nfail = 75;
prefix = "NetInference";
verbose = true;
ncores = 1;



function help_func()
    # make a help message
    help_message = """

        Infer phylogenetic networks from batches of the Concordance Factors (CFs) 
        using a simulated annealing algorithm. This algorithm uses the SNaQ algorithm
        and warm starts.

        Notice that if you give a single batch and 1 epoch, then you will have the 
        phylogenetic networks for that batch using the starting tree when the
        default parameters are set.

    Usage: $(PROGRAM_FILE) startfile CFfile batches 
            --h_max h_max --n_epochs n_epochs --tk tk --seed seed
            --nruns nruns --Nfail Nfail --prefix prefix 

    Required arguments:
        startfile: str; path to the file with the starting network.
        CFfile: str; path to the file with the CFs
        batches: str; path to the file with the batches

    Optional arguments:
        --h_max h_max: int; maximum number of hybridizations. (default: $h_max)
        --n_epochs n_epochs: int; number of epochs. (default: $n_epochs)
        --tk tk: float; temperature. Inf temperature
            accepts all suboptimal moves. Lower than Inf
            a probability of accepting a suboptimal move is
            calculated. (default: $tk)
        --seed seed: int; seed for the random number generator. (default: $seed)
        --nruns nruns: int; number of runs. (default: $nruns)
        --Nfail Nfail: int; number of failures. (default: $Nfail)
        --prefix prefix: str; prefix for the output files. (default: $prefix)
        --ncores: int; number of cores for running SNaQ (default: $ncores)    
    """;
    println(help_message);
    exit(0);    
end

if length(ARGS) < 3
    help_func();
end


for i in eachindex(ARGS)
    if i == 1 && !startswith( ARGS[i], "--" )
        global startfile = ARGS[i];
    elseif i == 2  && !startswith( ARGS[i], "--" )
        global buckyCFfile = ARGS[i];
    elseif i == 3  && !startswith( ARGS[i], "--" )
        global batches = ARGS[i];
    elseif ARGS[i] == "--h_max"
        global h_max = parse(Int, ARGS[i+1]);
    elseif ARGS[i] == "--n_epochs"
        global n_epochs = parse(Int, ARGS[i+1]);
    elseif ARGS[i] == "--tk"
        global tk = parse(Float64, ARGS[i+1]);
    elseif ARGS[i] == "--seed"
        global seed = parse(Int, ARGS[i+1]);
    elseif ARGS[i] == "--nruns"
        global nruns = parse(Int, ARGS[i+1]);
    elseif ARGS[i] == "--Nfail"
        global Nfail = parse(Int, ARGS[i+1]);
    elseif ARGS[i] == "--prefix"
        global prefix = ARGS[i+1];
    elseif ARGS[i] == "--ncores"
        global ncores = parse(Int, ARGS[i+1]);
    elseif ARGS[i] == "--help" || ARGS[i] == "-h"
        help_func();
    end
end

if startfile == "" || buckyCFfile == "" || batches == ""
    help_func();
end


using Suppressor;
using Distributed;
using CSV;
@suppress using DataFrames;
using Random;
using Distributed;

addprocs(ncores)
@suppress @everywhere using PhyloNetworks;


function LookConvergence(all_nets)
    # look for relative convergence
    if length(all_nets) >= 3

        curr_lik = all_nets[end].loglik;
        prev_lik = all_nets[end-1].loglik;
        diff = abs(curr_lik - prev_lik)/abs(prev_lik);
    
        println("\nrel diff: ", diff);
    end
end

"""
CT: DataFrame
    DataFrame with the CFs\\
batch: str
    string with the indices of the rows to subsample\\
"""
function subsampleCF(CT, batch)
    idx = [parse(Int, j) for j in split(batch, ",")];
    return readTableCF(CT[idx, :]);
end

function get_net(N_k, CT_k, h_max, nruns, Nfail, seed)

    try
        oldstd = stdout
        redirect_stdout(devnull)
        return snaq!(N_k, CT_k,
                      hmax=h_max,
                      filename="", 
                      runs=nruns, 
                      verbose=false, 
                      Nfail=Nfail,
                      seed=seed
                      );
        redirect_stdout(oldstd) # recover original stdout
    catch
        return nothing;
    end
end

function writeNets_Liks(all_nets, prefix)
    lik_file = prefix * "_liks.txt";
    net_file = prefix * "_nets.txt";

    all_nets = all_nets[2:end];

    if length(all_nets) == 0
        return;
    end

    open(lik_file, "w") do io
        for i in eachindex(all_nets);
            write(io, string(all_nets[i].loglik), "\n");
        end
    end

    open(net_file, "w") do io
        for i in eachindex(all_nets);
            write(io, writeTopology(all_nets[i]), "\n");
        end
    end

end

"""
get the networks for each batch in a lasso-like path (warm starts).
It infers the phylogenetic networks for each batch and considers the best network as a 
starting solution for the next batch.

when tk = Inf, 
    it always takes the previous solution as a starting solution.
when tk < Inf, 
    i) takes the best network so far as a starting solution or
    ii) takes the previous solution as a starting solution with probability p.

tk: float
    temperature\\
e: int
    epoch\\
l_best: float
    best pseudo-deviance so far\\

all_batches: Array{String}
    array with the batches\\
all_nets: Array{PhyloNetwork}
    array with the networks\\
CT: DataFrame
    DataFrame with the CFs\\
h_max: int
    maximum number of hybridizations\\
nruns: int
    number of runs\\
Nfail: int
    number of failures\\
seed: int
    seed for the random number generator\\
verbose: bool
    print messages\\
"""
function batches_path(tk, e, all_batches, l_best, all_nets, CT, h_max,  nruns, Nfail, seed, verbose)
    
    net_k = nothing;
    for (i,batch) in enumerate( all_batches )
        # i = 5
        # batch = all_batches[i]
        if verbose
            println("Processing batch ", i, " epoch ", e);
        end

        # subsample rows using batch's indices
        CT_k = subsampleCF(CT, batch);
        # get the network by using the previous solution as a warm start
        # even if this solution is suboptimal, accepted with probability p
        try
            oldstd = stdout
            redirect_stdout(devnull)
            net_k = snaq!(all_nets[end], CT_k, hmax=h_max,
                          runs=nruns, Nfail=Nfail, seed=seed, 
                          verbose=false, filename="",);
            # TODO: control for prob of NNI.
            redirect_stdout(oldstd)
        catch
            println("\nFailed to find a network for batch ", i, " epoch ", batch, "\n");
            continue;        
        end

        # if net_k == nothing
        #     println("\nFailed to find a network for batch ", i, " epoch ", batch, "\n");
        #     continue;
        # end
        
        dE = net_k.loglik - l_best;
        if dE <= 0
            # The quartet pseudo-deviance is such that a perfect fit corresponds 
            # to a deviance of 0.0. You want to minimize the deviance. 
            # If there was an improvement, then we added to the list of networks
            # and update the best pseudo-deviance. 
            if verbose
                println("\nOptimal move: Accepted with dE ", dE, "\n");
            end

            push!(all_nets, deepcopy(net_k));
            l_best = net_k.loglik;
        else
            # if there was not an improvement, then we take the new network
            # with probability p. 

            # notice that if tk = Inf, then p = 1 and we always take the new network
            # as a starting solution, which it is like a lasso path with warm starts.
            p = exp(-dE/tk);
            if rand() <= p
                if verbose
                    println("\nSuboptimal move: Accepted with probability ", p, "\n");
                end

                push!(all_nets, deepcopy(net_k));
            end
        end

        if verbose
            LookConvergence(all_nets);
        end
    end

    return l_best, all_nets;
end


"""
startfile: str
    path to the file with the starting network\\
buckyCFfile: str
    path to the file with the CFs\\
batches: str
    path to the file with the batches\\
h_max: int
    maximum number of hybridizations\\
n_epochs: int
    number of epochs\\
tk: float
    temperature\\
seed: int
    seed for the random number generator\\
nruns: int
    number of runs. SNaQ has 10 runs. This one has 1.\\
Nfail: int
    number of failures. SNaQ has 75. This one has 75.\\
prefix: str
    prefix for the output files\\
verbose: bool
    print messages\\
rate: int
    rate for decreasing the temperature\\
"""
function main(startfile, buckyCFfile, batches, 
    h_max = 1, n_epochs = 1, tk = Inf, seed = 12038,
    nruns = 10, Nfail = 75, 
    prefix = "./test_sims/disjointInference",
    verbose = true, 
    rate = 10
    )
    
    # startfile    = "./test_data/1_seqgen.QMC.tre";
    # buckyCFfile  = "./test_data/1_seqgen.CFs.csv";
    # batches      = "./test_data/linear_batches_overlappedBatches.txt";
    # h_max = 1;
    

    # read batches file
    all_batches = readlines(batches);
    # read intial network
    netstart = readTopology(startfile);
    # read csv file buckyCFfile
    CT = CSV.read(buckyCFfile, DataFrame);
        
    # best pseudo deviance
    l_best = Inf;
    all_nets = [deepcopy(netstart)];

    for e in 1:n_epochs
        # e = 1
        if (e > 1) & (tk < Inf) & (e % rate == 0)
            tk = 0.98 * tk;
        end
        # fit all batches in a lasso-like path (warm starts).
        # It infers the phylogenetic networks for each batch
        # and considers the best network as a starting solution
        # for the next batch.
        l_best, new_nets = batches_path(tk, e, all_batches, l_best, all_nets, CT, h_max, 
                                        nruns, Nfail, seed, verbose);

        all_nets = vcat(all_nets, new_nets);

        println("\nEpoch ", e, " done\n");
        println("Best pseudo-deviance: ", l_best, "\n");

    end
    # write nets and its corresponding pseudo-deviance
    # from index 2 to the end
    writeNets_Liks(all_nets, prefix);

end

@time main(startfile, buckyCFfile, batches, h_max, n_epochs, tk, seed, nruns, Nfail, prefix);
