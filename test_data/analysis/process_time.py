#!/usr/bin/env python3

import re
import argparse

# log_file = "/Users/ulisesrosas/Desktop/qsin/test_data/n15/n15_linear/linear_nets.log"
# out_file = "/Users/ulisesrosas/Desktop/qsin/test_data/n15/n15_linear/linear_nets.csv"


def main(log_file, out_file):

    all_lines = []
    with open(log_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()

            if 'Bootstrap' in line:
                boot,row = line.split(', ')
                boot = boot.replace('Bootstrap ', '')
                row = row.replace('row ', '')
                
                assert boot.isdigit(), f"Boot: {boot}"
                assert row.isdigit(), f"Row: {row}"

                all_lines.append(int(boot))
                all_lines.append(int(row))

            if 'seconds' in line:
                time = re.findall(r'^\d+\.\d+ seconds', line).pop()
                time = time.split(' ')[0]
                all_lines.append(float(time))


    with open(out_file, 'w') as f:
        f.write("bootstrap,row,time\n")
        for i in range(0, len(all_lines), 3):
            f.write(f"{all_lines[i]},{all_lines[i+1]},{all_lines[i+2]}\n")


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
