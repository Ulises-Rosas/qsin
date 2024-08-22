
# you can compare the obtained networks with the original network, which is available in the test_data folder.
# You  might need to install the PhyloPlots and RCall packages to run this code.
# The first plot is the full data network, and the second plot is the half data network.

using PhyloNetworks;
using PhyloPlots; # you may need to install PhyloPlots package
using RCall;      # you may need to install RCall package

qsin_net = readInputTrees("./linear_overlapped_nets.txt");
ori_net = readInputTrees("./test_data/full_data_net6.txt");


R"png"("./imgs/plot_nets.png", width=7, height=4, res=340, units="in");
R"layout(matrix(1:2, 1, 2))"; 
R"par"(mar=[0,0,1,0], bg = "white"); 
plot(ori_net[1], showgamma=true);
R"mtext"("Full data network (15 rows)")
plot(qsin_net[1], showgamma=true);
R"mtext"("Half data network (8 rows)")
R"dev.off()";
