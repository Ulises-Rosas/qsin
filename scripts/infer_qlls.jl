#!/usr/bin/env julia

using Distributed;


CFfile = "";
nets = [];
outfile = "qlls.csv";
ncores = 1;

ftolRel = 1e-6
ftolAbs = 1e-6
xtolAbs = 1e-3
xtolRel = 1e-2

function help_func()
    # make a help message
    help_message = """

    Calculate expected CF and overall pseudolikelihood score from
    a set of defined phylogenetic networks

    Usage: $(PROGRAM_FILE) CFfile [network files]
            --outfile outfile
            --ncores ncores

    Required arguments:
        CFfile: str; file with the CFs
        [network files]: [str]; a set of phylogenetic network files

    Optional arguments:
        --outfile outfile: str; output file name. (default: $outfile)
        --ftolRel: float; relative tolerance for the objective function (default: $ftolRel)
        --ftolAbs: float; absolute tolerance for the objective function (default: $ftolAbs)
        --xtolRel: float; relative tolerance for parameter changes (default: $xtolRel)
        --xtolAbs: float; absolute tolerance for parameter changes  (default: $xtolAbs)
        --ncores: int; number of cores (default: $ncores)        
        --help: display this help message
""";
    println(help_message);
    exit(0);    
end

if length(ARGS) < 2
    help_func();
end

# touched another argument?
toa = false

for i in eachindex(ARGS)

    if i == 1 && !startswith( ARGS[i], "--" )
        global CFfile = ARGS[i];
        continue
    end
        
    if !startswith( ARGS[i], "--" ) && !toa
        push!(nets, ARGS[i]);

    else
        global toa = true

        if ARGS[i] == "--ncores"
        global ncores = parse(Int, ARGS[i+1]);

        elseif ARGS[i] == "--outfile"
            global outfile = ARGS[i+1];

        elseif ARGS[i] == "--ftolRel"
            global ftolRel = parse(Float64, ARGS[i+1]);
        
        elseif ARGS[i] == "--ftolAbs"
            global ftolAbs = parse(Float64, ARGS[i+1]);
        
        elseif ARGS[i] == "--xtolRel"
            global xtolRel = parse(Float64, ARGS[i+1]);
        
        elseif ARGS[i] == "--xtolAbs"
            global xtolAbs = parse(Float64, ARGS[i+1]);

        elseif ARGS[i] == "--help" || ARGS[i] == "-h"
            help_func();
        end
    end

end

if CFfile == "" || length(nets) == 0 
    help_func();
end

# println("CFfile: ", CFfile);
# println("nets : ", length(nets));
# println("outfile: ", outfile);
# println("ncores: ", ncores);

using Suppressor;

addprocs(ncores)

@everywhere using CSV;
@suppress @everywhere using DataFrames;
@everywhere using PhyloNetworks;

@everywhere function QuartetCounts(ngenes, df_long)
    """
    ngenes: number of genes
    df_long: dataframe after using fittedQuartetCF with :long
    df_long[:,6] is the observed probability of the quartet
    """
    return repeat(ngenes, inner = 3) .* df_long[:,6]
end

@everywhere function std_loglik(ngenes, df_long)
    """
    standard log-likelihood

    ngenes: number of genes
    df_long: dataframe after using fittedQuartetCF with :long
    df_long[:,7] is the expected probability of the quartet

    From the documentation:
    "if long, the output has one row per quartet,
    i.e. 3 rows per 4-taxon sets, and *7 columns*:
    4 columns for the taxon names, one column to give 
    the quartet resolution, one column for the 
    observed CF and the *last column for 
    the expected CF."

    """
    QC = QuartetCounts(ngenes, df_long)
    return sum( QC .* log.( df_long[:,7] ) )
end

@everywhere function spps_code(df)
    # make rows to collapse in a string in a
    code_spps = []
    quartets = unique(df)
    for i in 1:size(quartets,1)
        # collapse all columns in a string
        tmp_code = join(quartets[i,:], ".")
        push!(code_spps, tmp_code)
    end
    
    return code_spps
end



@everywhere function QLL(ngenes, df_long)
    """
    quartet log-likelihood

    ngenes: number of genes
    df_long: dataframe after using fittedQuartetCF with :long
    """
    QC = QuartetCounts(ngenes, df_long)
    all_qlls = QC .* log.( df_long[:,7] )
    
    # loop that takes 3 rows at a time of all_qlls
    qlls = []
    for i in 1:3:size(all_qlls,1)
        push!(qlls, sum(all_qlls[i:i+2]))
    end
    
    return qlls
end

@everywhere function iter_df(ngenes, df_long)
    
    qll = QLL(ngenes, df_long)
    spps = spps_code(df_long[:, 1:4])

    push!(qll, sum(qll))
    push!(spps, "sum")

    return DataFrame(qll', spps)
end


function evaluate_sims(networks, buckyCFfile, outputfile, ftolRel, ftolAbs, xtolRel, xtolAbs)

    # buckyCFfile = "./test_data/1_seqgen.CFs_n15.csv"
    # netfile = "./test_data/n15_sim_netv2/n15_sim_netv2_sim2796.txt"

    @everywhere function process_network(netfile, all_buckyCF, dat, ftolRel, ftolAbs, xtolRel, xtolAbs)
        # O(1)

        netstart = readTopology(netfile) # O(n)

        try
            # branch lengths from simulation are clock time-based
            # and the time snaq considers branch lengths on
            # coalescent units. For that reason we use 
            # topologyMaxQPseudolik!

            # it returns a new network with updated branch lengths
            # and gamma values, which is not stored in any variable here.
            # The original network is not modified.
            # However, what is modified is all_buckyCF
            topologyMaxQPseudolik!(netstart, all_buckyCF, 
                                    ftolRel = ftolRel, 
                                    ftolAbs = ftolAbs, 
                                    xtolRel = xtolRel, 
                                    xtolAbs = xtolAbs)
            # println(new_net)
            
            # we transform all_buckyCF into a proper format
            df_long = fittedQuartetCF(all_buckyCF, :long)
            result = iter_df(dat.ngenes, df_long)
            return result
        catch e
            println("Error in ", netfile, ": ", e)
            return DataFrame()
        end
    end

    all_buckyCF = readTableCF(buckyCFfile)
    dat = DataFrame(CSV.File(buckyCFfile); copycols=false)

    main_df = @distributed (vcat) for netfile in networks
        println(netfile)

        dat_tmp = deepcopy(dat)
        all_buckyCF_tmp = deepcopy(all_buckyCF)
        process_network(netfile, all_buckyCF_tmp, dat_tmp, ftolRel, ftolAbs, xtolRel, xtolAbs)
    end
    CSV.write(outputfile, main_df)

end

@time evaluate_sims(nets, CFfile, outfile, ftolRel, ftolAbs, xtolRel, xtolAbs)
