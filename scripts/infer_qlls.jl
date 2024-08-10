#!/usr/bin/env julia

using Distributed;


CFfile = "";
nets = [];
outfile = "qlls.csv";
ncores = 1;


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
        --ncores: int; number of cores (default: $ncores)        
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


function simlated_QLL(networks, buckyCFfile, outputfile)

    # buckyCFfile = "/Users/ulises/Desktop/SLL/SparseQuartets/1_seqgen.CFs_n15.csv"
    # netfile = "/Users/ulises/Desktop/ABL/comps/claudia/UnderstandingNetworks/n15_sim_v2/test_500.txt"
    # netfile2 = "/Users/ulises/Desktop/ABL/comps/claudia/UnderstandingNetworks/n15_sim_v2/test_499.txt"
    # networks = [ netfile2, netfile]

    @everywhere function process_network(netfile, all_buckyCF, dat)
        netstart = readTopology(netfile)
        try
            topologyQPseudolik!(netstart, all_buckyCF)
            df_long = fittedQuartetCF(all_buckyCF, :long)
            return iter_df(dat.ngenes, df_long)
        catch
            println("Error in ", netfile)
            return DataFrame()
        end
    end

    function simlated_QLL(networks, buckyCFfile, outputfile)
        # buckyCFfile = "/Users/ulises/Desktop/SLL/SparseQuartets/1_seqgen.CFs_n15.csv"
        all_buckyCF = readTableCF(buckyCFfile)
        dat = DataFrame(CSV.File(buckyCFfile); copycols=false)

        main_df = @distributed (vcat) for netfile in networks
            println(netfile)

            dat_tmp = deepcopy(dat)
            all_buckyCF_tmp = deepcopy(all_buckyCF)
            process_network(netfile, all_buckyCF_tmp, dat_tmp)

            # process_network(netfile, all_buckyCF, dat)
        end
        CSV.write(outputfile, main_df)
    end

    simlated_QLL(networks, buckyCFfile, outputfile)

end


@time simlated_QLL(nets, CFfile, outfile)
