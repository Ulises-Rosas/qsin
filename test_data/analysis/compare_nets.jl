#!/usr/bin/env julia


using  Suppressor;
@suppress using PhyloNetworks;

main_net = "";
nets = [];
outfile = "diffs.txt";
root = "1";


function help_func()
    # make a help message
    help_message = """

    Calculate expected CF and overall pseudolikelihood score from
    a set of defined phylogenetic networks

    Usage: $(PROGRAM_FILE) main_net [network files]
            --outfile outfile

    Required arguments:
        main_net: str; main net
        [network files]: [str]; a set of phylogenetic network files

    Optional arguments:
        --outfile outfile: str; output file name. (default: $outfile)
        --root root: int; root for all nets. (default: $root)
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
        global main_net = ARGS[i];
        continue
    end
        
    if !startswith( ARGS[i], "--" ) && !toa
        push!(nets, ARGS[i]);

    else
        global toa = true

        if ARGS[i] == "--root"
            global root =  ARGS[i+1];

        elseif  ARGS[i] == "--outputfile"
            global outputfile = ARGS[i+1];

        elseif ARGS[i] == "--help" || ARGS[i] == "-h"
            help_func();
        end
    end

end

if main_net == "" || length(nets) == 0 
    help_func();
end

println("main_net: ", main_net);
println("nets : ", nets);
println("outfile: ", outfile);
println("root: ", root);

function root_exist(net, net_file, root, verbose = false)

    if root in [i.name for i in net.leaf]
        return true
    
    else
        if verbose
            println("$root taxa is not at $net_file")
        end

        return false
    end
    
end



function main(main_net, other_nets, root, outfile);

    target_net = readInputTrees(main_net)[1];
    if ~root_exist(target_net, main_net, root, true)
        # println("next time it will be better")
        return 1
    end

    rootatnode!(target_net, root);

    distances = [];
    for net_file in other_nets

        tmp_net = readInputTrees(net_file)[1];
        if ~root_exist(tmp_net, net_file, root, true)
            continue
        end
        rootatnode!(tmp_net, root)
        
        dist = hardwiredClusterDistance(target_net, tmp_net, true);
        push!(distances, dist);
    end

    println(distances)
    return distances
end


@time main(main_net, nets, root, outfile);

"""
./compare_nets.jl ../full_data_net6.txt\
                  ../n6/n6_lin_boot1_row*_nets.txt\
                  --root 1
"""


# if length(ARGS) >= 2

#     # net_file = "/Users/ulises/Desktop/SLL/SparseQuartets/just_networks.txt";
#     # root = "1";

#     net_file = ARGS[1];
#     root = ARGS[2];

#     outfile = replace(net_file, ".txt" => "_distances_root_$root.csv");


#     # read target_net
#     networks = readInputTrees(net_file);
#     # get the first element of the array

#     target_net = networks[1];


#     rootatnode!(target_net, root);
    


#     # R"par"(mar=[0,0,0,0]); # to reduce margins (no margins at all here)
#     # plot(target_net, showgamma=true);


#     distances = [];
#     for net in networks[2:end]
#         # println(tipLabels(net));
        
#         rootatnode!(net, root);
#         dist = hardwiredClusterDistance(target_net, net, true);
#         push!(distances, dist);
#     end

#     CSV.write(outfile, DataFrame(distances=distances), writeheader=false);

# else
#     println("Usage: estimate_network.jl network_file root_node");

# end





# using Suppressor;

# addprocs(ncores)

# @everywhere using CSV;
# @suppress @everywhere using DataFrames;
# @everywhere using PhyloNetworks;

# @everywhere function QuartetCounts(ngenes, df_long)
#     """
#     ngenes: number of genes
#     df_long: dataframe after using fittedQuartetCF with :long
#     df_long[:,6] is the observed probability of the quartet
#     """
#     return repeat(ngenes, inner = 3) .* df_long[:,6]
# end

# @everywhere function std_loglik(ngenes, df_long)
#     """
#     standard log-likelihood

#     ngenes: number of genes
#     df_long: dataframe after using fittedQuartetCF with :long
#     df_long[:,7] is the expected probability of the quartet

#     From the documentation:
#     "if long, the output has one row per quartet,
#     i.e. 3 rows per 4-taxon sets, and *7 columns*:
#     4 columns for the taxon names, one column to give 
#     the quartet resolution, one column for the 
#     observed CF and the *last column for 
#     the expected CF."

#     """
#     QC = QuartetCounts(ngenes, df_long)
#     return sum( QC .* log.( df_long[:,7] ) )
# end

# @everywhere function spps_code(df)
#     # make rows to collapse in a string in a
#     code_spps = []
#     quartets = unique(df)
#     for i in 1:size(quartets,1)
#         # collapse all columns in a string
#         tmp_code = join(quartets[i,:], ".")
#         push!(code_spps, tmp_code)
#     end
    
#     return code_spps
# end



# @everywhere function QLL(ngenes, df_long)
#     """
#     quartet log-likelihood

#     ngenes: number of genes
#     df_long: dataframe after using fittedQuartetCF with :long
#     """
#     QC = QuartetCounts(ngenes, df_long)
#     all_qlls = QC .* log.( df_long[:,7] )
    
#     # loop that takes 3 rows at a time of all_qlls
#     qlls = []
#     for i in 1:3:size(all_qlls,1)
#         push!(qlls, sum(all_qlls[i:i+2]))
#     end
    
#     return qlls
# end

# @everywhere function iter_df(ngenes, df_long)
    
#     qll = QLL(ngenes, df_long)
#     spps = spps_code(df_long[:, 1:4])

#     push!(qll, sum(qll))
#     push!(spps, "sum")

#     return DataFrame(qll', spps)
# end


# function simlated_QLL(networks, buckyCFfile, outputfile)

#     # buckyCFfile = "/Users/ulises/Desktop/SLL/SparseQuartets/1_seqgen.CFs_n15.csv"
#     # netfile = "/Users/ulises/Desktop/ABL/comps/claudia/UnderstandingNetworks/n15_sim_v2/test_500.txt"
#     # netfile2 = "/Users/ulises/Desktop/ABL/comps/claudia/UnderstandingNetworks/n15_sim_v2/test_499.txt"
#     # networks = [ netfile2, netfile]

#     @everywhere function process_network(netfile, all_buckyCF, dat)
#         netstart = readTopology(netfile)
#         try
#             topologyQPseudolik!(netstart, all_buckyCF)
#             df_long = fittedQuartetCF(all_buckyCF, :long)
#             return iter_df(dat.ngenes, df_long)
#         catch
#             println("Error in ", netfile)
#             return DataFrame()
#         end
#     end

#     function simlated_QLL(networks, buckyCFfile, outputfile)
#         # buckyCFfile = "/Users/ulises/Desktop/SLL/SparseQuartets/1_seqgen.CFs_n15.csv"
#         all_buckyCF = readTableCF(buckyCFfile)
#         dat = DataFrame(CSV.File(buckyCFfile); copycols=false)

#         main_df = @distributed (vcat) for netfile in networks
#             println(netfile)

#             dat_tmp = deepcopy(dat)
#             all_buckyCF_tmp = deepcopy(all_buckyCF)
#             process_network(netfile, all_buckyCF_tmp, dat_tmp)

#             # process_network(netfile, all_buckyCF, dat)
#         end
#         CSV.write(outputfile, main_df)
#     end

#     simlated_QLL(networks, buckyCFfile, outputfile)

# end


# @time simlated_QLL(nets, CFfile, outfile)

