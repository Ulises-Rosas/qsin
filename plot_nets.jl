
using PhyloNetworks;
using PhyloPlots; # you may need to install PhyloPlots package
qsin_net = readInputTrees("./linear_overlapped_nets.txt");
ori_net = readInputTrees("./test_data/full_data_net.txt");
# make the two plots in the same window

using RCall;  # you may need to install RCall package
R"layout(matrix(1:2, 1, 2))"; 
R"par"(mar=[0,0,1,0]) 
plot(ori_net[1], showgamma=true);
R"mtext"("Full data network")
plot(qsin_net[1], showgamma=true);
R"mtext"("Half data network")
# save plot
# R"dev.copy(png, 'plot_nets.png')";
