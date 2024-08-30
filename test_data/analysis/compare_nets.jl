#!/usr/bin/env julia

using  DelimitedFiles;
using  Suppressor;
@suppress using PhyloNetworks;

main_net = "";
nets = [];
outfile = "diffs.txt";
root = "1";


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


function get_dist(target_net, target_file, net_file, root);

    if ~root_exist(target_net, target_file, root, true)
        # println("next time it will be better")
        return error("Root taxa is not at $main_net")
    end

    rootatnode!(target_net, root);

    tmp_net = readInputTrees(net_file)[1];    
    if ~root_exist(tmp_net, net_file, root, true)
        return error("Root taxa is not at $net_file")
    end

    rootatnode!(tmp_net, root)
    # dist = hardwiredClusterDistance(target_net, tmp_net, true);
    return hardwiredClusterDistance(target_net, tmp_net, true);
end

function root_and_dist(all_taxa, target_net, target_file,net_file, root, verbose)

    while true
        try
            my_dist = get_dist(target_net, target_file, net_file, root);

            if verbose
                println("Root: ", root);
                println("Distances: ", my_dist);
            end
            
            return my_dist, root;

        catch e
            # choose randomdly a new root
            root = all_taxa[rand(1:end)];

            if verbose
                println(e);
                println("New root (randomdly chosen): ", root);
            end

            filter!(e->e ≠ root, all_taxa);
    
            if length(all_taxa) == 0
                return error("No more taxa to choose from");
            end

        end

    end
end


function main(target_file, root, net_files, outfile)

    target_net = readInputTrees(target_file)[1];

    all_taxa = [i.name for i in target_net.leaf];
    filter!(e->e ≠ root, all_taxa);
    
    all_dists = [];
    for net_file in net_files
        copied_taxa = copy(all_taxa);
        tmp_dist, tmp_root = root_and_dist(copied_taxa, target_net, target_file, net_file, root, true);
        
        # # get base name
        out_name = basename(net_file);
        row = match(r".*_row([0-9]+)_.*", out_name)[1];
        boot =  match(r".*_boot([0-9]+)_.*", out_name)[1];

        push!(all_dists, [row, boot, tmp_dist, tmp_root]);
    end
    println(all_dists)
    writedlm(outfile, all_dists, ',');
    # CSV.write(outfile, DataFrame(all_dists), writeheader=false);
end


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

# main_net = "./test_data/full_data_net6.txt";
# root = "3";
# # net_file = "./test_data/n6/n6_lin_boot2_row8_nets.txt";
# net_files = ["./test_data/n6/n6_lin_boot1_row1_nets.txt", "./test_data/n6/n6_lin_boot1_row2_nets.txt", "./test_data/n6/n6_lin_boot1_row3_nets.txt", "./test_data/n6/n6_lin_boot1_row4_nets.txt", "./test_data/n6/n6_lin_boot1_row5_nets.txt", "./test_data/n6/n6_lin_boot1_row6_nets.txt", "./test_data/n6/n6_lin_boot1_row7_nets.txt", "./test_data/n6/n6_lin_boot1_row8_nets.txt"];
# nets = net_files;
@time main(main_net, root, nets, outfile);

# all_dists = [];
# for i in eachindex(nets)
#     base_name = basename(nets[i]);
#     row = match(r".*_row([0-9]+)_.*", base_name);
#     boot =  match(r".*_boot([0-9]+)_.*", base_name);
#     push!(all_dists, [boot[1], row[1]]);
# end
# writedlm(outfile, all_dists, ',');

#TODO: test it out. It should work with the following command
"""
./test_data/analysis/compare_nets.jl ./test_data/full_data_net6.txt\
                  ./test_data/n6/n6_lin_boot1_row*_nets.txt\
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

