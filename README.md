# Quartet subsampling for phylogenetic network inference via sparse machine learning

## Installation

These two currently available options for installing `qsin` use the files `environment.yml` and `build.sh`, which are located among this repository files

### Option 1: Using conda

```bash
# construct the environment
conda env create -f environment.yml
conda activate qsin
# install julia and r dependencies at qsin
./build.sh 
```

### Option 2: Using Mamba

```bash
# construct the environment
mamba env create -f environment.yml
conda activate qsin
# install julia and r dependencies at qsin
./build.sh
```

