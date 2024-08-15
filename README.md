# Quartet subsampling for phylogenetic network inference via sparse machine learning

## Installation

The installation currently works for Linux and MacOS only. The two currently available options for installing `qsin` use the files `environment.yml` and `build.sh`, which are located among this repository files.
#### Option 1: Using conda

```bash
# construct the environment
conda env create -f environment.yml
conda activate qsin
# install julia dependencies at qsin
./build.sh 
```

#### Option 2: Using Mamba

```bash
# construct the environment
mamba env create -f environment.yml
conda activate qsin
# install julia dependencies at qsin
./build.sh
```


## Overview: A minimal example

#### 1. Simulate networks 

```bash
# create dir where sim nets will be stored
mkdir -p ./test_data/test_sims

# simulate networks
sim_nets.R 6 --max_iter 1000 --prefix test --out_path ./test_data/test_sims
```

#### 2. Expected concordance factors for simulated networks

```bash
infer_qlls.jl ./test_data/1_seqgen.CFs.csv\
              ./test_data/test_sims/test*.txt\
              --outfile ./test_data/test_qll.csv
```

#### 3. Create subsamples

```bash
path_subsampling.py ./test_data/test_qll.csv\
        ./test_data/1_seqgen.CFs.csv\
        --wpath     \
        --factor 0.5\
        --e 1e-2    \
        --alpha 0.5 \
        --verbose   \
        --prefix ./test_data/linear_batches
```

`--wpath` specify whether the elastic net path information is output it.  if `--factor` is  set to -1 then,  get the row combination with the lowest RMSE.

You can check the elastic net path produced for creating the subsample by running the following code:
```python
import numpy as np
import matplotlib.pyplot as plt

errors_path = np.loadtxt("./test_data/linear_batches_testErrors.csv", delimiter=',')
elnet_path  = np.loadtxt("./test_data/linear_batches_elnetPath.csv", delimiter=',')

fs = 18
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
ax1.plot(elnet_path[:,0], elnet_path[:,1:], marker='o', alpha=0.8)
ax1.set_xscale('log')
ax1.set_ylabel('Coefficient', fontsize=fs)
ax1.axhline(0, color='black', lw=2)
ax1.set_title('Elastic Net Path', fontsize=fs)
ax2.plot(errors_path[:,0], errors_path[:,1], marker='o', alpha=0.8)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel('$\\lambda$', fontsize=fs)
ax2.set_ylabel('RMSE', fontsize=fs)
```

![here](https://github.com/Ulises-Rosas/qsin/blob/main/imgs/linear_batches_elnetPath.png)

#### 4.  Phylogenetic network inference

```bash
# lets take the last
tail -n1 test_data/linear_batches_overlappedBatches.txt > ./test_data/linear_batches_last.txt

infer_nets_batches.jl ./test_data/1_seqgen.QMC.tre\
        ./test_data/1_seqgen.CFs.csv\
        ./test_data/linear_batches_last.txt\
        --h_max 1\
        --prefix linear_overlapped\
        --ncores 10
```

