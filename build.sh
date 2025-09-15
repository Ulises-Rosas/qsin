#!/bin/bash

set -e

julia -e 'using Pkg; Pkg.add("Suppressor")'
# DataFrames and CSV are required for PhyloNetworks
julia -e 'using Pkg, Suppressor; println("Installing DataFrames"); @suppress Pkg.add("DataFrames"); Pkg.add("CSV"); println("Installing PhyloNetworks"); @suppress Pkg.add(Pkg.PackageSpec(;name="PhyloNetworks", version="0.16.4"))'

pip install qsin -U

echo
echo -e "\033[7;32""m qsin installation complete \033[0m"
echo