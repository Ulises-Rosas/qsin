#!/bin/bash

julia -e 'using Pkg; Pkg.add("Suppressor")'
julia -e 'using Pkg, Suppressor; println("Installing DataFrames"); @suppress Pkg.add("DataFrames"); Pkg.add("CSV"); Pkg.add("PhyloNetworks")'
Rscript -e 'install.packages("SiPhyNetwork", repos="https://cloud.r-project.org")'


echo
echo -e "\033[7;32""m qsin installation complete \033[0m"
echo