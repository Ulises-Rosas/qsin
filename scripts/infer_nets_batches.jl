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

# println("startfile: ", startfile);
# println("buckyCFfile: ", buckyCFfile);
# println("batches: ", batches);
# println("h_max: ", h_max);
# println("n_epochs: ", n_epochs);
# println("tk: ", tk);
# println("seed: ", seed);
# println("nruns: ", nruns);
# println("Nfail: ", Nfail);
# println("prefix: ", prefix);
# println("ncores: ", ncores)

using Suppressor;
using Distributed;
using CSV;
@suppress using DataFrames;
using Random;
using Distributed;

addprocs(ncores)
@suppress @everywhere using PhyloNetworks;


function checkconvergence(all_liks, l_k)
    l_best_abs = maximum(abs.(all_liks));
    diff0 = l_k - all_liks[end];
    diff = abs(diff0)/l_best_abs;
    # take the maximum all_liks
    println("\nDiff: ", diff, "; lik = ", l_k, "\n");
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
    number of failures. SNaQ has 75. This one has 6.\\
"""
function main(startfile, buckyCFfile, batches, 
    h_max = 2, n_epochs = 1, tk = 1000, seed = 120,
    nruns = 10, Nfail = 6, prefix = "./test_sims/disjointInference")
    
    # startfile = "/Users/ulises/Desktop/SLL/SparseQuartets copy/1_seqgen.QMC_n15.tre";
    # buckyCFfile = "/Users/ulises/Desktop/SLL/SparseQuartets copy/1_seqgen.CFs_n15.csv";
    # # batches = "/Users/ulises/Desktop/SLL/SparseQuartets copy/test_sims/batches_n15_qll_g2.txt";
    # # batches = "/Users/ulises/Desktop/SLL/SparseQuartets copy/test_sims/linear_batches_overlappedBatches.txt";
    # batches = "/Users/ulises/Desktop/SLL/SparseQuartets copy/test_sims/linear_batches_disjointBatches.txt";
    # h_max = 1;
    # read batches file

    all_batches = readlines(batches);
    # split csv the first line
    # ARGS = 
    # batch_files[1]
    
    
    netstart    = readTopology(startfile);
    all_buckyCF = readTableCF(buckyCFfile);
    # h_max = hmax;
    
    # read csv file buckyCFfile
    CT = CSV.read(buckyCFfile, DataFrame);
    
    
    
    N_prev = deepcopy(netstart);
    all_liks = [];
    all_nets = [];
    # tk = 1000;
    
    i = 0;
    # n_epochs = 1;
    error_at = [];
    net_k = deepcopy(netstart);
    for epoch in 1:n_epochs
        i = 0
        if epoch > 1
            tk = 0.98 * tk;
        end
        for batch in all_batches
            i += 1
            println("Processing batch ", i, " epoch ", epoch)

            # batch = all_batches[i]
    
            idx = [parse(Int, j) for j in split(batch, ",")];
            CT_k = readTableCF(CT[idx, :]);
            try
                oldstd = stdout
                redirect_stdout(devnull)
                net_k = snaq!(N_prev, CT_k,
                    hmax=h_max,
                    filename="", 
                    runs=nruns, 
                    verbose=false, 
                    Nfail=Nfail,
                    seed=seed, 
                    );
                redirect_stdout(oldstd) # recover original stdout
            catch
                println("Error in ", batch);
                push!(error_at, N_prev);
                if error_at[end] == N_prev
                    println("Error in ", batch);
                    N_prev = all_nets[end-1];
                end
                continue
            end
            
            if i == 1 && epoch == 1
                l_k = topologyQPseudolik!(net_k, all_buckyCF);

                push!(all_liks, l_k);
                push!(all_nets, net_k);
                N_prev = deepcopy(net_k);
                continue;
            end
    
            l_k = topologyQPseudolik!(net_k, all_buckyCF);
            l_best = minimum(all_liks);

            dE = l_k - l_best;
                
            if dE < 0
                """
                The quartet pseudo-deviance is the
                    **negative** log pseudo-likelihood, up
                    to an additive constant, such that
                    a perfect fit corresponds to a deviance
                    of 0.0. 
    
                => it is a minimization problem              
                """
                println("\nOptimal move: Accepted with likelihood diff. ", abs(dE), "\n");
                N_prev = deepcopy(net_k);
                checkconvergence(all_liks, l_k);

                push!(all_liks, l_k);
                push!(all_nets, net_k);
            
            else
                # make a random uniform number
                # from 0 to 1
                xi = rand();
                p = exp(-dE/tk);
                println("\np =  ", p, "");
                if xi <= p
                    println("\nSuboptimal move: Accepted with probability ", p, "\n");
                    N_prev = deepcopy(net_k);
                    checkconvergence(all_liks, l_k);

                    push!(all_liks, l_k);
                    push!(all_nets, net_k);
                end
            end
        end
    end
    
    lik_file = prefix * "_liks.txt";
    net_file = prefix * "_nets.txt";

    open(lik_file, "w") do io
        for i in eachindex(all_liks);
            write(io, string(all_liks[i]), "\n");
        end
    end
    
    open(net_file, "w") do io
        for i in eachindex(all_nets);
            write(io, writeTopology(all_nets[i]), "\n");
        end
    end
end

@time main(startfile, buckyCFfile, batches, h_max, n_epochs, tk, seed, nruns, Nfail, prefix);


# # using PhyloPlots;
# # using RCall;
# # netstart = readInputTrees("./test_sims/all_nets.txt");
# # plot(netstart[2], showgamma=true);

