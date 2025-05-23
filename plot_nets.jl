
# you can compare the obtained networks with the original network, which is available in the test_data folder.
# You  might need to install the PhyloPlots and RCall packages to run this code.
# The first plot is the full data network, and the second plot is the half data network.

using PhyloNetworks;
using PhyloPlots; # you may need to install PhyloPlots package
using RCall;      # you may need to install RCall package

net_files = [
    "/Users/ulisesrosas/Desktop/experiments_qsin/cui_eta0.5_nu0.1_leaves4_progress/n15_from_0_uniq_row1_boot1_nets.txt",
    "/Users/ulisesrosas/Desktop/experiments_qsin/cui_eta0.5_nu0.1_leaves4_progress/n15_from_0_uniq_row2_boot1_nets.txt",
    "/Users/ulisesrosas/Desktop/experiments_qsin/cui_eta0.5_nu0.1_leaves4_progress/n15_from_0_uniq_row3_boot1_nets.txt",
]

all_rows_file = "/Users/ulisesrosas/Desktop/experiments_qsin/cui_eta0.5_nu0.1_leaves4_progress/n15_from_0_uniq.txt"
total_number_rows = 10626

# read all_rows_file csv file
all_rows = open(all_rows_file) do f
    readlines(f)
end



row_sizes = length.(split.(all_rows, ","))


# make sorted list of files
net_files = sort(net_files)
all_networks = []

for f in net_files
    net = readInputTrees(f)[1];
    push!(all_networks, net)
end


function left_rights(n, root)
    if n.leaf
        return nothing, nothing
    end
    if n == root
        left = n.edge[1].node[1]
        right = [ n.edge[2].node[1], n.edge[3].node[1] ]
    else
        left = n.edge[1].node[1]
        right = n.edge[2].node[1]
    end

    return left, right
end

function  get_ancestor(n)
    n_anc  = nothing
    for e in n.edge
        if e.isMajor
            u, u_anc  = e.node
            # ancestry direction goes like: v <- u
            # println(u.number," ",v.number)
            if u == n
                n_anc = u_anc
            end
        end
    end
    return n_anc
end


function dfs(n, path, paths, root)
    
    if n !== nothing

        push!(path, n)

        if n.leaf
            push!(paths, path)
        end

        left, rights = left_rights(n, root)

        dfs(left, path, paths, root)

        if n == root
            dfs(rights[1], path, paths, root)
            dfs(rights[2], path, paths, root)
        else
            dfs(rights, path, paths, root)
        end

    end
end

tmp_net = all_networks[1]
plot(tmp_net,shownodenumber=true,)




for n in tmp_net.node
    if n == root
        continue
    end
    println(n.number)
    n_anc = get_ancestor(n)
    println(n_anc.number)
    println()
end



root.edge

n = hyb_nodes[1]



R"layout(matrix(1:3, 1, 3))"; 
R"par"(mar=[0,0,2,0], bg = "white"); 
# R"mtext"("rho")
for (i,net) in enumerate(all_networks)
    size = row_sizes[i]
    rho = round( size*100 / total_number_rows, digits=2)
    # R"mtext"("$(size) rows, $(rho) rho")
    # add title to the plot
    # R"mtext"("$(size) rows, $(rho) rho", side=3, line=0.5)
    # title!("My Title")
    plot(net, showgamma=true,);
    R"title"("$(size) rows, $(rho) perc.", line=0.5)
    println("$(size) rows, $(rho) perc.")
end
# R"png"("./imgs/plot_net_cui.png", width=10, height=4, res=340, units="in");


# R"layout(matrix(1:2, 1, 2))"; 
# R"par"(mar=[0,0,1,0], bg = "white"); 
# plot(all_nets[1], showgamma=true);
# R"mtext"("Full data network (15 rows)")
# plot(all_nets[2], showgamma=true);
# R"mtext"("Half data network (8 rows)")
# # R"png"("./imgs/plot_net_cui.png", width=10, height=4, res=340, units="in");
