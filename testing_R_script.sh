# print(N)
# print(out_path)
# print(prefix)
# print(ncores)
# print(max_iter)


./scripts/sim_networks_v2.R 15 --out_path ./test_data/test_sims --prefix test --ncores 10 --max_iter 500

# quartets outfile [sim. networks]
# julia -p 15 ./scripts/estimate_qlls_v2.jl ./test_data/1_seqgen.CFs_n15.csv ./test_data/test_n15_qll.csv ./test_data/test_sims/test2*.txt

./scripts/estimate_qlls_v2.jl ./test_data/1_seqgen.CFs_n15.csv ./test_data/test_sims/test2*.txt\
                                  --ncores 10 --outfile ./test_data/test_n15_qll.csv


./path_subsampling.py test_sims/test_n15_qll.csv 1_seqgen.CFs_n15.csv\
        --wpath --verbose --e 1e-3 --factor -1\
        --prefix ./test_sims/linear_batches



./path_subsampling.py test_sims/test_n15_qll.csv 1_seqgen.CFs_n15.csv\
    --wpath --verbose --isle --e 1e-3 --factor -1\
    --prefix ./test_sims/non_linear_batches



./infer_nets_batches.jl 1_seqgen.QMC_n15.tre\
                        1_seqgen.CFs_n15.csv\
                        non_linear_batches_overlappedBatches_isle.txt\
                        --h_max 2 --nruns 10 --Nfail 75\
                        --prefix non_linear_overlapped



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

