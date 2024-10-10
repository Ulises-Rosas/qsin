# !/bin/bash

set -e


from=$1
# set path
cd $2
ncores=6


echo "From_count: $from"
prefix="n15_depth3_v0.9_n1_from_"$from
mkdir -p ./test_data/$prefix




nfolds=5
h_max=2
ntaxa=15
nsims=1000
prefix2="n15_depth3_v0.9_n1"

CT_file="./test_data/1_seqgen.CFs_n15.csv"
init_tree="./test_data/1_seqgen.QMC_n15.tre"

./src/path_subsampling.py ./test_data/test_qll_n15.csv\
        $CT_file\
        --wpath --verbose --e 1e-4\
        --isle --eta 0.9 --nu 0.9 --max_depth 3\
        --factor -1 --inbetween 5\
        --prefix ./test_data/$prefix/$prefix2\
        --cv --alpha 0.5 0.95 0.99 1\
        --ncores $ncores --folds $nfolds

cat ./test_data/$prefix/$prefix2"_overlappedBatches_isle.txt" | uniq > ./test_data/$prefix/$prefix2"_overlappedBatches_uniq.txt" 
nrows=$(cat ./test_data/$prefix/$prefix2"_overlappedBatches_uniq.txt" | wc -l)


# boostraps=10
# for i in $(seq 1 $boostraps)
# do     
#      for j in $(seq 1 $nrows)
#      do   
#           seed=$RANDOM
#           boot=$(( $i + $from ))
          
#           echo "Bootstrap $boot, row $j"
#           echo "Seed: $seed"
#           awk "NR==$j" ./test_data/$prefix/$prefix2"_overlappedBatches_uniq.txt" > ./test_data/$prefix/$prefix2"_row"$j".txt"
          
#           ./scripts/infer_nets_batches.jl $init_tree\
#                     $CT_file\
#                     ./test_data/$prefix/$prefix2"_row"$j".txt"\
#                     --h_max $h_max\
#                     --prefix ./test_data/$prefix/$prefix2"_boot"$boot"_row"$j \
#                     --ncores $ncores\
#                     --seed $seed
#      done
# done > ./test_data/$prefix/$prefix2"_nets.log"
# echo "Done"


# sim_nets.R $ntaxa --out_path ./test_data/$prefix\
#             --prefix $prefix"_sim" --ncores $ncores --max_iter $nsims

# infer_qlls.jl $CT_file\
#               ./test_data/$prefix/$prefix"_sim"*.txt\
#               --ncores $ncores\
#               --outfile ./test_data/$prefix/test_qll_$prefix.csv
