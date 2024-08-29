# !/bin/bash

# This script is an example of how to run the pipeline for the analysis of a single dataset.
# The script is divided into three main parts:
# 1. Simulation of networks
# 2. Inference of networks
# 3. Bootstrapping of networks
# The script uses the test data provided in the repository.


ncores=10
nfolds=2
h_max=1
CT_file="./test_data/1_seqgen.CFs.csv"
init_tree="./test_data/1_seqgen.QMC.tre"


sim_nets.R 6 --out_path ./test_data/n6\
            --prefix n6_sim --ncores $ncores --max_iter 500

infer_qlls.jl $CT_file\
              ./test_data/n6/n6_sim*.txt\
              --ncores $ncores\
              --outfile ./test_data/n6/test_qll_n6.csv

path_subsampling.py ./test_data/n6/test_qll_n6.csv\
        $CT_file\
        --wpath --verbose --e 1e-2 --factor 0.90 --inbetween 10\
        --prefix ./test_data/n6/linear_batches\
        --cv --alpha 0.5 0.95 0.99 1 --ncores $ncores --folds $nfolds


cat ./test_data/n6/linear_batches_overlappedBatches.txt | uniq > ./test_data/n6/linear_batches_overlappedBatches_uniq.txt 
nrows=$(cat ./test_data/n6/linear_batches_overlappedBatches_uniq.txt | wc -l)
boostraps=2
for i in $(seq 1 $boostraps)
do     
     for j in $(seq 1 $nrows)
     do   
          seed=$RANDOM
          
          echo "Bootstrap $i, row $j"
          echo "Seed: $seed"

          awk "NR==$j" ./test_data/n6/linear_batches_overlappedBatches_uniq.txt > "./test_data/n6/lbo_$j.txt"
          
          ./scripts/infer_nets_batches.jl $init_tree\
                    $CT_file\
                    "./test_data/n6/lbo_$j.txt"\
                    --h_max $h_max\
                    --prefix "./test_data/n6/n6_lin_boot"$i"_row"$j \
                    --ncores $ncores\
                    --seed $seed
     done
done > ./test_data/n6/n6_lin_nets.log

