# Use a stable Ubuntu LTS version as the base image
FROM ubuntu:22.04

# Set environment variables to allow non-interactive installation
ENV DEBIAN_FRONTEND=noninteractive
ENV JULIA_VERSION=1.10.4

# --- System Setup ---
# Update package lists and install essential build tools and utilities
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    software-properties-common \
    ca-certificates \
    curl \
    wget \
    gnupg \
    && rm -rf /var/lib/apt/lists/*

# --- Python Installation ---
# Install Python 3, pip, and create a symlink for 'python'
RUN apt-get update && \
    apt-get install -y python3 python3-pip python3-venv && \
    ln -s /usr/bin/python3 /usr/bin/python && \
    rm -rf /var/lib/apt/lists/*

# --- R Installation ---
# Add the CRAN repository to get the latest version of R
RUN apt-get update && \
    apt-get install -y --no-install-recommends software-properties-common dirmngr && \
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" && \
    apt-get install -y --no-install-recommends r-base r-base-dev && \
    rm -rf /var/lib/apt/lists/*

# --- Julia Installation ---
# Download and install the official Julia binary
RUN curl -fsSL "https://julialang-s3.julialang.org/bin/linux/x64/${JULIA_VERSION%.*}/julia-${JULIA_VERSION}-linux-x86_64.tar.gz" | \
    tar -xz -C /opt && \
    ln -s /opt/julia-*/bin/julia /usr/local/bin/julia

# Julia packages
RUN julia -e 'using Pkg; Pkg.add(["DataFrames", "CSV", "Suppressor"]); Pkg.add(Pkg.PackageSpec(;name="PhyloNetworks", version="0.16.4")); Pkg.precompile()'


# R packages
RUN Rscript -e 'install.packages("SiPhyNetwork", repos="https://cloud.r-project.org")'


# Python packages
RUN pip install qsin -U


# --- Final Configuration ---
# Set a working directory
WORKDIR /work

# Start a bash shell by default when the container runs
CMD ["/bin/bash"]