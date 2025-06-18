#!/usr/bin/env julia

using Distributed;

main_net = "";
nets = [];
outfile = "diffs.txt";
thresh = 0.0;
ncores = 1;
let_root = false;

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
        --ncores: int; number of cores (default: $ncores)
        --let_root: bool; let the distance be obtained after rooting (default: $let_root)
        --help, -h: show this help message
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
        global toa = true;

        if ARGS[i] == "--thresh"
            global thresh =  parse(Float64, ARGS[i+1]);

        elseif  ARGS[i] == "--outfile"
            global outfile = ARGS[i+1];

        elseif ARGS[i] == "--let_root"
            global let_root = true;

        elseif ARGS[i] == "--ncores"
            global ncores = parse(Int, ARGS[i+1]);

        elseif ARGS[i] == "--help" || ARGS[i] == "-h"
            help_func();
        end
    end

end

if main_net == "" || length(nets) == 0 
    help_func();
end


using  Suppressor;

addprocs(ncores);


using  DelimitedFiles;
@suppress using Statistics;
@suppress @everywhere using DataFrames;
@suppress @everywhere using PhyloNetworks;

@everywhere function get_dist(true_net, net_file, root, thresh; let_root = true);

    # true_net = true_nets_tmp[end]
    # net_file = net_files[1]
    
    tmp_net = readInputTrees(net_file)[1]
    deleteHybridThreshold!(tmp_net, thresh)

    if let_root
        rootatnode!(tmp_net, root)
        rootatnode!(true_net, root)
    end


    return hardwiredClusterDistance(true_net, tmp_net, let_root);
end

@everywhere function root_and_dist(all_taxa, true_net, net_file, thresh; let_root = true)

    # all_taxa = deepcopy(all_taxa);
    # curr_root = deepcopy(root);
    # thresh = 0.0 ;
    # verbose = false;
    # println(all_taxa);
    best_dist = Inf;
    best_root = NaN;

    if !let_root
        try
            best_dist = get_dist(true_net, net_file, nothing, thresh; let_root = false);
        catch e
            nothing
        end

        return best_dist, NaN;
    end

    for curr_root in all_taxa
        # curr_root=all_taxa[2];
        try
            curr_dist = get_dist(true_net, net_file, curr_root, thresh; let_root = true);

            # if verbose
            # println("Root: ", curr_root);
            # println("Distances: ", curr_dist);
            # end

            if curr_dist < best_dist
                best_dist = curr_dist;
                best_root = curr_root;
            end

        catch e
            # println(e);
            continue
        end
    end

    if best_dist == Inf
        best_dist = NaN;
    end

    return best_dist, best_root;
end

function report_results(all_dists_mat)

    table = Dict();
    for i in eachrow(all_dists_mat)
        row  = i[1];
        dist = i[3];

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

@everywhere function evaluate_true_net( true_nets, all_taxa, tmp_net_file, thresh; let_root = true)

    best_net_tmp_dist = Inf;
    best_net_tmp_root = NaN;

    for true_net in true_nets
        # println("true_net: ", true_net);

        all_taxa_tmp = deepcopy(all_taxa);
        # println("all_taxa_tmp: ", all_taxa_tmp);
        tmp_dist, tmp_root = root_and_dist(all_taxa_tmp, true_net, tmp_net_file, thresh; let_root = let_root);

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
    
    row = match(r".*_row(.+)_[bn].*", tmp_name_base)[1];
    boot =  match(r".*_boot([0-9]+)_.*", tmp_name_base)[1];
    # boot =  match(r".*_boot([0-9]+).net", tmp_name_base)[1];

    println( 
        "tmp_net_file: ", tmp_name_base,
        " row: ", row, 
        " boot: ", boot, 
        " best_net_tmp_dist: ", best_net_tmp_dist, 
        " best_net_tmp_root: ", best_net_tmp_root
    );
    # println([row, boot, best_net_tmp_dist, best_net_tmp_root])
    return [row, boot, best_net_tmp_dist, best_net_tmp_root]
end

function main(true_nets_file, thresh, net_files, outfile; let_root = true)

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

    all_dists = @distributed (vcat) for tmp_net_file in net_files
        true_nets_tmp = deepcopy(true_nets);
        all_taxa_tmp = deepcopy(all_taxa);
        tmp = evaluate_true_net(true_nets_tmp, all_taxa_tmp, tmp_net_file, thresh; let_root = let_root)

        [tmp]
    end

    println("\n");
    # arr = Matrix(all_dists);
    # report_results(arr);
    # println("all_dists: ", all_dists);
    writedlm(outfile, all_dists, ',');

end


println("main_net: ", main_net);
println("len of nets : ", length(nets));
println("outfile: ", outfile);
println("thresh: ", thresh);


@time main(main_net, thresh, nets, outfile; let_root = let_root);

#TODO: test it out. It should work with the following command
"""
./test_data/analysis/compare_nets.jl ./test_data/full_data_net6.txt\
                  ./test_data/n6/n6_lin_boot1_row*_nets.txt
"""
