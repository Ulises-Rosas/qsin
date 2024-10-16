#!/usr/bin/env python

import numpy as np
np.random.seed(12038)

# uniq = '/Users/ulises/Desktop/ABL/software/qsin/test_data/n15_depth3_v0.9_n0.5/nonlinear_depth3_nu0.9_overlappedBatches_uniq.txt'
# max_size = 1365

def main(uniq, max_size):
    new_subset = []
    with open(uniq, 'r') as f:
        for line in f.readlines():        
            tmp_len = len(line.strip().split(','))
            picked = np.random.choice(range(1, max_size + 1), tmp_len, replace=False)
            new_subset.append(",".join([str(p) for p in sorted(picked)]))
            
    with open(uniq.replace('uniq', 'uniq_random'), 'w') as f:
        for line in new_subset:
            f.write(line + '\n')

def arg_parser():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('uniq', type=str, help='uniq file')
    parser.add_argument('max_size', type=int, help='max size')
    parser.add_argument('--seed', type=int, default=12038, help='random seed')
    return parser.parse_args()

if __name__ == '__main__':
    args = arg_parser()
    np.random.seed(args.seed)
    main(args.uniq, args.max_size)
