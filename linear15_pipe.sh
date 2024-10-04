# !/bin/bash

# This script is an example of how to run the pipeline for the analysis of a single dataset.
# The script is divided into three main parts:
# 1. Simulation of networks
# 2. Inference of networks
# 3. Bootstrapping of networks
# The script uses the test data provided in the repository.

set -e

ncores=15
nfolds=5
h_max=2
ntaxa=15
nsims=1000
prefix="n15_linear"
prefix2="linear_5e-4"

CT_file="./test_data/1_seqgen.CFs_n15.csv"
init_tree="./test_data/1_seqgen.QMC_n15.tre"


# sim_nets.R $ntaxa --out_path ./test_data/$prefix\
#             --prefix $prefix"_sim" --ncores $ncores --max_iter $nsims

# infer_qlls.jl $CT_file\
#               ./test_data/$prefix/$prefix"_sim"*.txt\
#               --ncores $ncores\
#               --outfile ./test_data/$prefix/test_qll_$prefix.csv


path_subsampling.py ./test_data/$prefix/test_qll_n15.csv\
        $CT_file\
        --wpath --verbose --e 5e-4 --factor -1 --inbetween 5\
        --prefix ./test_data/$prefix/$prefix2\
        --cv --alpha 0.5 0.95 0.99 1\
        --ncores $ncores --folds $nfolds

cat ./test_data/$prefix/$prefix2"_overlappedBatches.txt" | uniq > ./test_data/$prefix/$prefix2"_overlappedBatches_uniq.txt" 
nrows=$(cat ./test_data/$prefix/$prefix2"_overlappedBatches_uniq.txt" | wc -l)

boostraps=15
for i in $(seq 1 $boostraps)
do     
     for j in $(seq 1 $nrows)
     do   
          seed=$RANDOM
          
          echo "Bootstrap $i, row $j"
          echo "Seed: $seed"

          awk "NR==$j" ./test_data/$prefix/$prefix2"_overlappedBatches_uniq.txt" > ./test_data/$prefix/$prefix2"_row"$j".txt"
          
          infer_nets_batches.jl $init_tree\
                    $CT_file\
                    ./test_data/$prefix/$prefix2"_row"$j".txt"\
                    --h_max $h_max\
                    --prefix ./test_data/$prefix/$prefix2"_boot"$i"_row"$j \
                    --ncores $ncores\
                    --seed $seed
     done
done > ./test_data/$prefix/$prefix2"_nets.log"


