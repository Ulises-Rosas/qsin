
./scripts/sim_nets.R 15 --out_path ./test_data/test_sims --prefix test --ncores 2 --max_iter 50

./scripts/infer_qlls.jl ./test_data/1_seqgen.CFs_n15.csv ./test_data/test_sims/test2*.txt\
                                  --ncores 2 --outfile ./test_data/test_n15_qll.csv

./src/path_subsampling.py ./test_data/test_n15_qll.csv ./test_data/1_seqgen.CFs_n15.csv\
        --wpath --verbose --e 1e-3 --factor -1\
        --prefix ./test_data/linear_batches

./src/path_subsampling.py ./test_data/test_n15_qll.csv ./test_data/1_seqgen.CFs_n15.csv\
    --wpath --verbose --isle --e 1e-3 --factor -1\
    --prefix ./test_data/non_linear_batches




./scripts/infer_nets_batches.jl ./test_data/1_seqgen.QMC_n15.tre\
                        ./test_data/1_seqgen.CFs_n15.csv\
                        ./test_data/non_linear_batches_overlappedBatches_isle.txt\
                        --h_max 2 --nruns 2 --Nfail 15\
                        --prefix non_linear_overlapped \
                        --ncores 2 

# ./infer_nets_batches.jl 1_seqgen.QMC_n15.tre\
#                         1_seqgen.CFs_n15.csv\
#                         linear_batches_overlappedBatches.txt\
#                         --h_max 2 --nruns 10 --Nfail 75\
#                         --prefix linear_overlapped


# ./infer_nets_batches.jl 1_seqgen.QMC_n15.tre\
#                         1_seqgen.CFs_n15.csv\
#                         non_linear_batches_disjointBatches_isle.txt\
#                         --h_max 2 --nruns 10 --Nfail 1\
#                         --prefix non_linear_overlapped

