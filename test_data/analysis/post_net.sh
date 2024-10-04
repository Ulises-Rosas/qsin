
# full_data="./test_data/full_data_net15.txt"
# inferd_nets_loc="/Users/ulisesrosas/Desktop/qsin/test_data/n15/n15_linear"


# full_data="./test_data/full_data_net15.txt"
# inferd_nets_loc="/Users/ulisesrosas/Desktop/qsin/test_data/n15/n15_linear"
# log_file=$inferd_nets_loc/linear_nets.log

# ./test_data/analysis/compare_nets.jl $full_data\
#                   $inferd_nets_loc/*_nets.txt\
#                   --outfile $inferd_nets_loc/compared_nets_linear.csv\
#                   --root 15

# ./test_data/analysis/process_time.py $log_file\
#                   -o $inferd_nets_loc/processed_time_linear.csv


# conda activate  /Users/ulisesrosas/miniconda3/envs/qsin 

full_data="./test_data/full_data_net15.txt"
inferd_nets_loc="/Users/ulisesrosas/Desktop/qsin/test_data/n15/n15_depth4_tmp"
log_file=$inferd_nets_loc/nonlinear_depth4_nets_15more.log

./test_data/analysis/compare_nets.jl $full_data\
                  $inferd_nets_loc/*_nets.txt\
                  --outfile $inferd_nets_loc/compared_nets_depth3_30boots.csv\
                  --root 15

./test_data/analysis/process_time.py $log_file\
                  -o $inferd_nets_loc/processed_time_nonlinear.csv


echo $inferd_nets_loc/compared_nets_depth3_30boots.csv
ls $inferd_nets_loc/*uniq.txt
echo $inferd_nets_loc/processed_time_nonlinear.csv



ls $inferd_nets_loc/*testErrors.csv
ls $inferd_nets_loc/*elnetPath.csv