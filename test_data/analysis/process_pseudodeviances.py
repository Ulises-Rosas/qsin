#!/usr/bin/env python3

import re
import os
import argparse
# import glob

def parse_files(files):
    out = [['row','boot' ,'likdev']]*(len(files) + 1)

    for i, file in enumerate(files):
        base_name = os.path.basename(file)
        row = re.findall(r"row(\d+)_", base_name)[0]
        boot = re.findall(r"boot(\d+)_", base_name)[0]
        with open(file) as f:
            likdev = f.readlines()[0].strip()
        out[i+1] = [int(row), int(boot), float(likdev)]

    return out

def main(files, out_file):
    # out_file = 'n15_depth3_v0.9_n0.5_liks.csv'
    out = parse_files(files)
    with open(out_file, 'w') as f:
        for i in out:
            f.write(f"{i[0]},{i[1]},{i[2]}\n")

def parse_args():
    parser = argparse.ArgumentParser(description='Process liks from log file')
    parser.add_argument('files', type=str, nargs='+', help='Path to log file')
    parser.add_argument('-o','--out_file', type=str, help='Path to output file', default='out.csv')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    main(args.files, args.out_file)


