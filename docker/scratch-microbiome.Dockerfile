# Use a specific version of Ubuntu as the base image
FROM --platform=linux/x86_64 rocker/verse:latest

# Set the working directory inside the container
WORKDIR /opt

# Timezone settings
ENV TZ=US/Central
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    echo $TZ > /etc/timezone

# Install system dependencies
RUN apt-get update && apt-get install -y \
    software-properties-common \
    dirmngr \
    gnupg \
    apt-transport-https \
    ca-certificates \
    wget \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    build-essential \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libpng-dev \
    libtiff-dev \
    libjpeg-dev

# Updating Quarto to v1.4.553
RUN wget https://github.com/quarto-dev/quarto-cli/releases/download/v1.4.553/quarto-1.4.553-linux-amd64.deb -O quarto.deb
RUN dpkg -i quarto.deb && rm quarto.deb

# Install Python3 and essential tools
RUN apt-get install -y \
    python3 \
    python3-pip

# Install bioinformatics tools
RUN apt-get install -y \
    samtools \
    bedtools \
    bowtie2 \
    kraken \
    kraken2 \
    minimap2 \
    hdf5-tools

# Install KrakenUniq and Bracken
RUN git clone https://github.com/fbreitwieser/krakenuniq.git /opt/krakenuniq && \
    cd /opt/krakenuniq && \
    make && \
    ln -s /opt/krakenuniq/krakenuniq /usr/local/bin/krakenuniq

RUN git clone https://github.com/jenniferlu717/Bracken.git /opt/bracken && \
    cd /opt/bracken && \
    make && \
    ln -s /opt/bracken/bracken /usr/local/bin/bracken

# Install MetaPhlAn
RUN pip3 install --no-cache-dir metaphlan

# Install HUMAnN for functional profiling
RUN pip3 install --no-cache-dir humann

# Install UMI-tools
RUN pip3 install --no-cache-dir umi_tools

# Install fundamental R packages
ARG R_DEPS="c(\
    'tidyverse', \
    'devtools', \
    'rmarkdown', \
    'BiocManager', \
    'remotes', \
    'optparse', \
    'R.utils', \
    'here' \
    )"

ARG R_BIOC_DEPS="c(\
    'Biobase', \
    'BiocGenerics', \
    'DelayedArray', \
    'DelayedMatrixStats', \
    'S4Vectors',\
    'SingleCellExperiment', \
    'SummarizedExperiment', \
    'HDF5Array', \ 
    'limma', \
    'lme4', \
    'ggrastr', \
    'Rsamtools' \
    )"

# Setting repository URL
ARG R_REPO="http://cran.us.r-project.org"

# Install R dependencies
RUN Rscript -e "install.packages(${R_DEPS}, Ncpus = 8, repos = '${R_REPO}', clean = TRUE)"
RUN Rscript -e "BiocManager::install(${R_BIOC_DEPS})"

# Install Python data science libraries
RUN python3 -m pip install --no-cache-dir numpy pandas scikit-learn matplotlib seaborn jupyter
RUN python3 -m pip install --no-cache-dir jupyter-cache
RUN python3 -m pip install --no-cache-dir papermill

# Cleaning up
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

# Command to run on container start
CMD ["bash"]
