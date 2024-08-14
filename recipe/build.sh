#!/bin/bash

# 'Autobrew' is being used by more and more packages these days
# to grab static libraries from Homebrew bottles. These bottles
# are fetched via Homebrew's --force-bottle option which grabs
# a bottle for the build machine which may not be macOS 10.9.
# Also, we want to use conda packages (and shared libraries) for
# these 'system' dependencies. See:
# https://github.com/jeroen/autobrew/issues/3
echo $PWD
ls -lhtr py_dep
ls -lhtr r_dep

cd r_dep

export DISABLE_AUTOBREW=1

# R refuses to build packages that mark themselves as Priority: Recommended
mv DESCRIPTION DESCRIPTION.old
grep -va '^Priority: ' DESCRIPTION.old > DESCRIPTION
# shellcheck disable=SC2086
${R} CMD INSTALL --build . ${R_ARGS}
ls -lhtr
# Add more build steps here, if they are necessary.
cd ../py_dep


ls -lhtr 
# See
# https://docs.conda.io/projects/conda-build
# for a list of environment variables that are set during the build process.
${PYTHON} -m pip install . -vv --no-deps --no-build-isolation
ls -lhtr






# #!/bin/bash


# export JULIA_DEPOT_PATH_BACKUP=${JULIA_DEPOT_PATH:-}
# export JULIA_PROJECT_BACKUP=${JULIA_PROJECT:-}
# export JULIA_LOAD_PATH_BACKUP=${JULIA_LOAD_PATH:-}

# export JULIA_DEPOT_PATH="$CONDA_PREFIX/share/julia:$JULIA_DEPOT_PATH"
# # Ensure the correct RPATH is set
# export JULIA_PROJECT="@${CONDA_PREFIX##*/}"
# # Modify load path so that projects stack on the conda-named environment
# export JULIA_LOAD_PATH="@:$JULIA_PROJECT:@stdlib"

# export JULIA_SSL_CA_ROOTS_PATH_BACKUP=${JULIA_SSL_CA_ROOTS_PATH:-}
# if [[ $OSTYPE != 'darwin'* ]]; then
#   export JULIA_SSL_CA_ROOTS_PATH=$CONDA_PREFIX/ssl/cacert.pem
# fi

# export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"

# julia -e 'using Pkg; Pkg.add("Suppressor")'
# # Rscript -e 'install.packages("SiPhyNetwork", repos="https://cloud.r-project.org")'


