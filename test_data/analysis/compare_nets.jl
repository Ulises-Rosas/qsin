#!/usr/bin/env julia

using  DelimitedFiles;
using  Suppressor;
@suppress using PhyloNetworks;

main_net = "";
nets = [];
outfile = "diffs.txt";
root = "1";
thresh = 0.0;


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


function get_dist(target_net, target_file, net_file, root, thresh);

    if ~root_exist(target_net, target_file, root, true)
        # println("next time it will be better")
        return error("Root taxa is not at $main_net")
    end

    # thresh = 0.1;
    # deleteHybridThreshold!(target_net, thresh);
    rootatnode!(target_net, root);

    tmp_net = readInputTrees(net_file)[1];    
    if ~root_exist(tmp_net, net_file, root, true)
        return error("Root taxa is not at $net_file")
    end

    deleteHybridThreshold!(tmp_net, thresh);
    rootatnode!(tmp_net, root)

    return hardwiredClusterDistance(target_net, tmp_net, true);
end

function root_and_dist(all_taxa, target_net, target_file,net_file, root, thresh, verbose)

    while true
        try
            my_dist = get_dist(target_net, target_file, net_file, root, thresh);

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
                println("No more taxa to choose from");
                return NaN, NaN;
                # return error("No more taxa to choose from");
            end

        end

    end
end


function main(target_file, root, thresh, net_files, outfile)

    # target_net = readInputTrees(target_file)[1];
    true_nets = readInputTrees(target_file);

    all_taxa = [i.name for i in true_nets[1].leaf];
    # println("All taxa: ", all_taxa);
    filter!(e->e ≠ root, all_taxa);
    
    all_dists = [];
    for net_file in net_files

        best_tmp_dist = Inf;
        best_tmp_root = NaN;
        for target_net in true_nets
            copied_taxa = copy(all_taxa);
            # target_net = true_nets[i];
            # println("target_net: ", target_net);

            tmp_dist, tmp_root = root_and_dist(copied_taxa, target_net, target_file, net_file, root, thresh, false);
            if !isnan(tmp_dist) && tmp_dist < best_tmp_dist
                best_tmp_dist = tmp_dist;
                best_tmp_root = tmp_root;
            end
        end

        if best_tmp_dist == Inf
            best_tmp_dist = NaN;
        end
                
        # tmp_dist, tmp_root = root_and_dist(copied_taxa, target_net, target_file, net_file, root, thresh, true);
        
        # # get base name
        out_name = basename(net_file);
        row = match(r".*_row([0-9]+)_.*", out_name)[1];
        boot =  match(r".*_boot([0-9]+)_.*", out_name)[1];

        push!(all_dists, [row, boot, best_tmp_dist, best_tmp_root]);
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
        --thresh thresh: float; threshold for hybrid nodes. (default: $thresh)
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
        global toa = true

        if ARGS[i] == "--root"
            global root =  ARGS[i+1];

        elseif ARGS[i] == "--thresh"
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
println("root: ", root);
println("thresh: ", thresh);


@time main(main_net, root, thresh, nets, outfile);

#TODO: test it out. It should work with the following command
"""
./test_data/analysis/compare_nets.jl ./test_data/full_data_net6.txt\
                  ./test_data/n6/n6_lin_boot1_row*_nets.txt\
                  --root 1
"""
