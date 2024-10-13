
full_data="./test_data/full_data_net15.txt"
inferd_nets_loc="/Users/ulises/Desktop/ABL/software/qsin/test_data/n15_non_linear_oscer"

log_file=$(ls $inferd_nets_loc/*.log)

# get base name of $inferd_nets_loc
inferd_nets_loc_base=$(basename $inferd_nets_loc)

./test_data/analysis/compare_nets.jl $full_data\
                  $inferd_nets_loc/*_nets.txt\
                  --outfile $inferd_nets_loc/compared_nets_$inferd_nets_loc_base.csv\
                  --root 12
# 12, 9 
./test_data/analysis/process_time.py $log_file\
                  -o $inferd_nets_loc/processed_time_$inferd_nets_loc_base.csv

./test_data/analysis/process_pseudodeviances.py $inferd_nets_loc/*liks.txt\
                  -o $inferd_nets_loc/processed_pseudo_$inferd_nets_loc_base.csv


echo
echo

echo file = "'$inferd_nets_loc/compared_nets_$inferd_nets_loc_base.csv'"
echo combs_file = "'$(ls $inferd_nets_loc/*uniq.txt)'"
echo time_file ="'$inferd_nets_loc/processed_time_$inferd_nets_loc_base.csv'"
echo pseudo_file ="'$inferd_nets_loc/processed_pseudo_$inferd_nets_loc_base.csv'"

echo
echo

error_file=$(ls $inferd_nets_loc/*testErrors.csv)
path_file=$(ls $inferd_nets_loc/*elnetPath.csv)


echo 'errors_path = np.loadtxt("'$error_file'", delimiter="," )'
echo 'elnet_path  = np.loadtxt("'$path_file'", delimiter="," )'
