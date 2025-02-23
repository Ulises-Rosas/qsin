#!/usr/bin/env julia

using  DelimitedFiles;
using  Suppressor;
using Statistics;
@suppress using PhyloNetworks;

main_net = "";
nets = [];
outfile = "diffs.txt";
# root = "1";
thresh = 0.0;


"""
net: PhyloNetwork; phylogenetic network
net_file: str; network file. It is just for printing purposes
root: str; root taxa
verbose: bool; verbose mode
"""
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


function get_dist(true_net, net_file, root, thresh);
    
    # target_file it is just for printing purposes
    if ~root_exist(true_net, "true_net_file", root, true)
        return error("Root taxa is not at $main_net")
    end

    # thresh = 0.1;
    # deleteHybridThreshold!(target_net, thresh);
    rootatnode!(true_net, root);

    tmp_net = readInputTrees(net_file)[1];    
    if ~root_exist(tmp_net, net_file, root, true)
        return error("Root taxa is not at $net_file")
    end

    deleteHybridThreshold!(tmp_net, thresh);
    rootatnode!(tmp_net, root)

    return hardwiredClusterDistance(true_net, tmp_net, true);
end

function root_and_dist(all_taxa, true_net, net_file, thresh, verbose)

    # all_taxa = deepcopy(all_taxa);
    # curr_root = deepcopy(root);
    # thresh = 0.0 ;
    # verbose = false;
    
    best_dist = Inf;
    best_root = NaN;

    for curr_root in all_taxa
        # curr_root=all_taxa[2];
        try
            curr_dist = get_dist(true_net, net_file, curr_root, thresh);

            if verbose
                println("Root: ", curr_root);
                println("Distances: ", curr_dist);
            end

            if curr_dist < best_dist
                best_dist = curr_dist;
                best_root = curr_root;
            end

        catch e
            if verbose
                println(e);
            end

            continue
        end
    end

    if best_dist == Inf
        best_dist = NaN;
    end

    return best_dist, best_root;
end

function report_results(all_dists)

    table = Dict();
    for i in eachindex(all_dists)
        row = all_dists[i][1];
        dist = all_dists[i][3];

        if isnan(dist)
            continue
        end

        if !haskey(table, row)
            table[row] = [dist];
        else
            push!(table[row], dist);
        end
    end

    all_row_names = sort(collect(keys(table)));

    for k in all_row_names
        println(
            "Row: ", k, 
            " Mean: ", round( mean(table[k]), digits=3), 
            " Std: ",  round( std(table[k]) , digits=3)
        );
    end

end


function main(true_nets_file, thresh, net_files, outfile)

    # true_nets_file = "/Users/ulisesrosas/Desktop/qsin/test_data/n15_hyb3.txt";
    # net_files = [
    #     "/Users/ulisesrosas/Desktop/experiments_qsin/n15_h3_eta0.5_nu0.1_d3_Tc_progress/n15_from_0_uniq_row2_boot4_nets.txt",
    #     "/Users/ulisesrosas/Desktop/experiments_qsin/n15_h3_eta0.5_nu0.1_d3_Tc_progress/n15_from_52_uniq_row5_boot53_nets.txt",
    #     "/Users/ulisesrosas/Desktop/experiments_qsin/n15_h3_eta0.5_nu0.1_d3_Tc_progress/n15_from_28_uniq_row1_boot29_nets.txt",
    #     "/Users/ulisesrosas/Desktop/experiments_qsin/n15_h3_eta0.5_nu0.1_d3_Tc_progress/n15_from_28_uniq_row2_boot30_nets.txt ",
    #     "/Users/ulisesrosas/Desktop/experiments_qsin/n15_h3_eta0.5_nu0.1_d3_Tc_progress/n15_from_8_uniq_row1_boot10_nets.txt",
    # ];

    true_nets = readInputTrees(true_nets_file);
    all_taxa = [i.name for i in true_nets[1].leaf];
    
    all_dists = [];
    for tmp_net_file in net_files

        best_net_tmp_dist = Inf;
        best_net_tmp_root = NaN;

        for true_net in true_nets

            # target_net = true_nets[3];
            # println("target_net: ", target_net);

            tmp_dist, tmp_root = root_and_dist(all_taxa, true_net, tmp_net_file, thresh, false);

            if !isnan(tmp_dist) && tmp_dist < best_net_tmp_dist
                best_net_tmp_dist = tmp_dist;
                best_net_tmp_root = tmp_root;

                if best_net_tmp_dist == 0
                    break
                end
                
            end
        end

        if best_net_tmp_dist == Inf
            best_net_tmp_dist = NaN;
        end

        # # get base name
        tmp_name_base = basename(tmp_net_file);

        row = match(r".*_row([0-9]+)_.*", tmp_name_base)[1];
        boot =  match(r".*_boot([0-9]+)_.*", tmp_name_base)[1];

        push!(all_dists, [row, boot, best_net_tmp_dist, best_net_tmp_root]);

        println( 
            "tmp_net_file: ", tmp_name_base,
            " row: ", row, 
            " boot: ", boot, 
            " best_net_tmp_dist: ", best_net_tmp_dist, 
            " best_net_tmp_root: ", best_net_tmp_root);

    end

    # println(all_dists);
    println("");
    report_results(all_dists);
    writedlm(outfile, all_dists, ',');
end


function help_func()
    # make a help message
    help_message = """

    Calculate distance between a main network and
    a set of defined phylogenetic networks

    Usage: $(PROGRAM_FILE) main_net [network files]
            --outfile outfile

    Required arguments:
        main_net: str; main net
        [network files]: [str]; a set of phylogenetic network files

    Optional arguments:
        --outfile outfile: str; output file name. (default: $outfile)
        --thresh thresh: float; threshold for hybrid nodes. (default: $thresh)
        --help, -h: show this help message
""";
    println(help_message);
    exit(0);    
end

# --root root: int; root for all nets. (default: $root)

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

        # if ARGS[i] == "--root"
        #     global root =  ARGS[i+1];

        if ARGS[i] == "--thresh"
            global thresh =  parse(Float64, ARGS[i+1]);

        elseif  ARGS[i] == "--outfile"
            global outfile = ARGS[i+1];

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
# println("root: ", root);
println("thresh: ", thresh);


@time main(main_net, thresh, nets, outfile);

#TODO: test it out. It should work with the following command
"""
./test_data/analysis/compare_nets.jl ./test_data/full_data_net6.txt\
                  ./test_data/n6/n6_lin_boot1_row*_nets.txt
"""
