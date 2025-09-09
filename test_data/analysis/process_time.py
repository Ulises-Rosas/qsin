#!/usr/bin/env python3

import re
import argparse

# log_file = "/Users/ulises/Desktop/ABL/software/experiments_qsin/CV_5foldLeaf2_batches_progress/CV_5foldLeaf2_batches_progress.log"
# out_file = "/Users/ulisesrosas/Desktop/qsin/test_data/n15/n15_linear/linear_nets.csv"


def main(log_file, out_file):
    D = {}
    with open(log_file, 'r') as f:
        lines = f.readlines()
        curr_key = None

        for line in lines:
            line = line.strip()

            if 'Bootstrap' in line:
                # line = "      From worker 6:    [121] signal (15): TerminatedBootstrap 110, row 1"
                boot,row = line.split(', ')
                boot = re.sub(".*Bootstrap ", "", boot)
                row = row.replace('row ', '')
                
                assert boot.isdigit(), f"Boot: {boot}"
                assert row.isdigit(), f"Row: {row}"
                curr_key = (int(boot), int(row))
                D[curr_key] = None
                continue

            if 'seconds' in line and curr_key is not None:
                time = re.findall(r'^\d+\.\d+ seconds', line).pop()
                time = float(time.split(' ')[0])
                D[curr_key] = time


    with open(out_file, 'w') as f:
        f.write("bootstrap,row,time\n")
        for k,v in D.items():
            f.write(f"{k[0]},{k[1]},{v}\n")


def parse_args():
    parser = argparse.ArgumentParser(description='Process time from log file')
    parser.add_argument('log_file', type=str, help='Path to log file')
    parser.add_argument('-o','--out_file', type=str, help='Path to output file', default='out.csv')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    main(args.log_file, args.out_file)




# line = "3258.509218 seconds (24.03 M allocations: 1.645 GiB, 0.03% gc time, 0.45% compilation time: 4% of which was recompilation)"
# line.isdigit()
# 'secondsads' in line
