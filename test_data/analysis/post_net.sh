#!/bin/bash

full_data=$1
inferd_nets_loc=$2
compare_nets=$3
process_time=$4
process_pseudo=$5
ncores=$6


log_file=$(ls $inferd_nets_loc/*.log)

# get base name of $inferd_nets_loc
inferd_nets_loc_base=$(basename $inferd_nets_loc)

$compare_nets $full_data\
                  $inferd_nets_loc/*_nets.txt\
                  --outfile $inferd_nets_loc/compared_nets_$inferd_nets_loc_base.csv\
                  --thresh 0.0 --ncores $ncores

$process_time $log_file\
                  -o $inferd_nets_loc/processed_time_$inferd_nets_loc_base.csv

$process_pseudo $inferd_nets_loc/*liks.txt\
                  -o $inferd_nets_loc/processed_pseudo_$inferd_nets_loc_base.csv


echo
echo

echo file = "'$( readlink -f $inferd_nets_loc/compared_nets_$inferd_nets_loc_base.csv)'"
echo combs_file = "'$(readlink -f $(ls $inferd_nets_loc/*uniq.txt))'"
echo time_file ="'$(readlink -f $inferd_nets_loc/processed_time_$inferd_nets_loc_base.csv)'"
echo pseudo_file ="'$(readlink -f $inferd_nets_loc/processed_pseudo_$inferd_nets_loc_base.csv)'"

echo
echo

error_file=$(ls $inferd_nets_loc/*testErrors.csv)
path_file=$(ls $inferd_nets_loc/*elnetPath.csv)


echo 'errors_path = np.loadtxt("'$error_file'", delimiter="," )'
echo 'elnet_path  = np.loadtxt("'$path_file'", delimiter="," )'
