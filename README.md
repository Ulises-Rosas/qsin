# Quartet subsampling for inferring phylogenetic networks via sparse learning

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

#### Option 3: Manual Installation

For this it requires that you have `julia`, `R` and `python` installed in your system. You can install the dependencies for `julia` by running the following command in the Julia console:

```julia
using Pkg
Pkg.add("CSV"); Pkg.add("DataFrames"); Pkg.add("PhyloNetworks"); Pkg.add("Suppressor")
```
You can install the dependencies for `R` by running the following command in the R console:

```R
install.packages("SiPhyNetwork")
```

You can install the dependencies for `python` by running the following command in the terminal:

```bash
pip install qsin
```

## Overview: A minimal example

### 1. Simulate networks 

```bash
# create dir where sim nets will be stored
mkdir -p ./test_data/test_sims

# simulate networks
sim_nets.R 6 --max_iter 500 --prefix test --out_path ./test_data/test_sims
```

### 2. Expected concordance factors for simulated networks

```bash
infer_qlls.jl ./test_data/1_seqgen.CFs.csv\
              ./test_data/test_sims/test*.txt\
              --outfile ./test_data/test_qll.csv
```

### 3. Create subsamples

```bash
path_subsampling.py ./test_data/test_qll.csv\
        ./test_data/1_seqgen.CFs.csv\
        --wpath     \
        --factor 0.5\
        --e 1e-2    \
        --cv\
        --folds 2\
        --alpha 0.5 0.95 0.99 1\
        --ncores 4\
        --verbose\
        --prefix ./test_data/linear_batches
```

* `--wpath` specify whether the elastic net path information is output it.
* `--factor` is the subsampling factor. 0.5 means that we write the selected rows of the input file until around half of the rows are selected.
* `--e` is the constant that multiplies the lambda_max to get the lambda_min value.
* `--cv` specifies that cross-validation is to be used.
* `--folds` is the number of folds for cross-validation.
* `--alpha` is the alpha values to be used.
* `--ncores` is the number of cores to be used in cross-validation.
* `--verbose` specifies that the output should be verbose.
* `--prefix` is the prefix for the output files.


You can check the elastic net path produced for creating the subsample by running the following `python` code:

```python
import numpy as np
import matplotlib.pyplot as plt

errors_path = np.loadtxt("./test_data/linear_batches_testErrors.csv", delimiter=',')
elnet_path  = np.loadtxt("./test_data/linear_batches_elnetPath.csv", delimiter=',')

fs = 14
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 5))

ax1.plot(elnet_path[:,0], elnet_path[:,1:], marker='o', alpha=0.8, markersize=5)
ax1.set_xscale('log')
ax1.set_ylabel('Coefficient', fontsize=fs)
ax1.axhline(0, color='black', lw=2)
ax1.set_title('Elastic Net Path', fontsize=fs)

ax2.plot(errors_path[:,0], errors_path[:,1], marker='o', alpha=0.8, markersize=5)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel('$\\lambda$', fontsize=fs)
ax2.set_ylabel('RMSE', fontsize=fs)
ax2.set_title('Test Error', fontsize=fs)
plt.tight_layout()
```

![here](https://github.com/Ulises-Rosas/qsin/blob/main/imgs/linear_batches_elnetPath.png)

### 4.  Phylogenetic network inference

```bash
# lets take the last row from the subsamples obtained above
tail -n1 test_data/linear_batches_overlappedBatches.txt > ./test_data/linear_batches_last.txt

# lets infer a network using that subsample
infer_nets_batches.jl ./test_data/1_seqgen.QMC.tre\
        ./test_data/1_seqgen.CFs.csv\
        ./test_data/linear_batches_last.txt\
        --h_max 1\
        --prefix linear_overlapped\
        --ncores 10
```

We can compare the obtained networks with the original network, which is available in the test_data folder.  You  might need to install the PhyloPlots and RCall packages to run this code in `julia`. The first plot is the full data network, and the second plot is the half data network.

```julia
using PhyloNetworks;
using PhyloPlots; # you may need to install PhyloPlots package
using RCall;      # you may need to install RCall package

qsin_net = readInputTrees("./linear_overlapped_nets.txt");
ori_net = readInputTrees("./test_data/full_data_net6.txt");

R"layout(matrix(1:2, 1, 2))"; 
R"par"(mar=[0,0,1,0], bg = "white"); 
plot(ori_net[1], showgamma=true);
R"mtext"("Full data network (15 rows)")
plot(qsin_net[1], showgamma=true);
R"mtext"("Half data network (8 rows)")
```
![here](https://github.com/Ulises-Rosas/qsin/blob/main/imgs/plot_nets.png)
