#!/usr/bin/env julia


# make default values for the arguments
startfile = ""
buckyCFfile = ""
batches = ""

h_max = 1
warm_start = false
seed = 12038
nruns = 15
Nfail = 75
prefix = "NetInference"
ncores = 1



function help_func()
    # make a help message
    help_message = """

        Infer phylogenetic networks from index subsamples of the Concordance Factors (CFs).
        This algorithm uses the SNaQ algorithm. It can also use previous network as warm start.

        Notice that if you give a single batch and 1 epoch, then you will have the
        phylogenetic networks for that batch using the starting tree when the
        default parameters are set.

    Usage: $(PROGRAM_FILE) startfile CFfile batches 
            --h_max h_max  
            [--seed seed] [--nruns nruns] [--Nfail Nfail] 
            [--prefix prefix]
            [--warm_start] [--ncores ncores] [--help]

    Required arguments:
        startfile: str; path to the file with the starting network.
        CFfile: str; path to the file with the CFs
        batches: str; path to the file with the batches

    Optional arguments:
        --h_max h_max: int; maximum number of hybridizations. (default: $h_max)
        --warm_start: bool; use warm start from previous network. (default: $warm_start)
        --seed seed: int; seed for the random number generator. (default: $seed)
        --nruns nruns: int; number of runs. (default: $nruns)
        --Nfail Nfail: int; number of failures. (default: $Nfail)
        --prefix prefix: str; prefix for the output files. (default: $prefix)
        --ncores: int; number of cores for running SNaQ (default: $ncores)
    """
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
    elseif ARGS[i] == "--warm_start"
        global warm_start = true;
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
seed: int
    seed for the random number generator\\
nruns: int
    number of runs. SNaQ has 10 runs. This one has 1.\\
Nfail: int
    number of failures. SNaQ has 75. This one has 75.\\
"""
function main(startfile, buckyCFfile, batches, 
    h_max = 2, seed = 120,
    nruns = 10, Nfail = 75, 
    prefix = "./test_sims/disjointInference"; warm_start = true)

    all_batches = readlines(batches);
    netstart  = readTopology(startfile);
    CT = CSV.read(buckyCFfile, DataFrame);
        
    
    N_prev = deepcopy(netstart);
    all_liks = []
    all_nets = []

    i = 1
    net_k = nothing
    for batch in all_batches

        println("Processing batch ", i)
        idx = [parse(Int, j) for j in split(batch, ",")]
        CT_k = readTableCF(CT[idx, :])
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
                )
            redirect_stdout(oldstd) # recover original stdout
        catch
            println("Error in ", batch)
            continue
        end

        if warm_start
            N_prev = deepcopy(net_k)
        end

        push!(all_liks, net_k.loglik)
        push!(all_nets, net_k)
        i += 1
    end
    
    lik_file = prefix * "_liks.txt"
    net_file = prefix * "_nets.txt"

    open(lik_file, "w") do io
        for i in eachindex(all_liks)
            write(io, string(all_liks[i]), "\n")
        end
    end
    
    open(net_file, "w") do io
        for i in eachindex(all_nets)
            write(io, writeTopology(all_nets[i]), "\n")
        end
    end
end

@time main(startfile, buckyCFfile, batches, h_max, seed, nruns, Nfail, prefix;
    warm_start = warm_start);
