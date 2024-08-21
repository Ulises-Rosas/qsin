#!/bin/bash

set -e

julia -e 'using Pkg; Pkg.add("Suppressor")'
# DataFrames and CSV are required for PhyloNetworks
julia -e 'using Pkg, Suppressor; println("Installing DataFrames"); @suppress Pkg.add("DataFrames"); Pkg.add("CSV"); println("Installing PhyloNetworks"); @suppress Pkg.add("PhyloNetworks")'
# Rscript -e 'install.packages("SiPhyNetwork", repos="https://cloud.r-project.org")'

# julia -e 'using Pkg, Suppressor; println("Installing PhyloNetworks"); @suppress Pkg.add("PhyloNetworks")'

pip install qsin -U

echo
echo -e "\033[7;32""m qsin installation complete \033[0m"
echo