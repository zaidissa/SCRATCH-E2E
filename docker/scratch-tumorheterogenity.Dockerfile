# Use a specific version of Ubuntu as the base image
FROM --platform=linux/x86_64 rocker/verse:latest

# Set the working directory inside the container
WORKDIR /opt

# Timezone settings
ENV TZ=US/Central
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    echo $TZ > /etc/timezone
# pass your PAT at build time so remotes::install_github can auth
ARG GITHUB_PAT
ENV GITHUB_PAT=${GITHUB_PAT}

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
    default-jre \
    libgfortran5 \
    liblapack-dev \
    libopenblas-dev \
    libjpeg-dev \
    libpng-dev \
    libtiff5-dev \
    zlib1g-dev \
    libxt-dev


# Updating quarto to Quarto v1.4.553
RUN wget https://github.com/quarto-dev/quarto-cli/releases/download/v1.4.553/quarto-1.4.553-linux-amd64.deb -O quarto-1.4.553-linux-amd64.deb
RUN dpkg -i quarto-1.4.553-linux-amd64.deb


RUN apt-get update && apt-get install -y python3 python3-pip python3-venv python3-dev build-essential

# Install core R packages
RUN Rscript -e "install.packages(c('R.utils','rmarkdown','devtools','tidyverse','readr', 'dplyr', 'ggplot2', 'cowplot', 'remotes', 'BiocManager','reticulate', 'HGNChelper'), repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages(c('leiden', 'RColorBrewer', 'reshape2', 'scales', 'NMF', 'colorspace', 'tibble', 'data.table', 'stringr', 'Matrix', 'bigmemory', 'patchwork', 'pheatmap'), repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('viridis',dependencies = TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('bigmemory',dependencies = TRUE)"
RUN Rscript -e "install.packages('doMC', dependencies = TRUE, repos='http://R-Forge.R-project.org')"
RUN Rscript -e "install.packages('optparse', dependencies = TRUE)"
RUN Rscript -e "install.packages('pheatmap', dependencies = TRUE)"


RUN Rscript -e "BiocManager::install(c('S4Vectors','DelayedMatrixStats','BiocGenerics','Biobase', 'SummarizedExperiment', 'AnnotationDbi', 'org.Hs.eg.db'), ask=FALSE, update=TRUE)"
RUN Rscript -e "BiocManager::install(c( \
    'HDF5Array','rhdf5','rhdf5lib', \
    'SingleCellExperiment', \
    'GOSemSim','MatrixGenerics','treeio','DOSE','ggtree','enrichplot', \
    'clusterProfiler','DirichletMultinomial','rtracklayer','GenomicFeatures', \
    'BSgenome','ensembldb','TFBSTools', \
    'BSgenome.Hsapiens.UCSC.hg38','EnsDb.Hsapiens.v86', 'Homo.sapiens'), \
  ask=FALSE, update=FALSE )"

# Setting repository URL
ARG R_REPO="http://cran.us.r-project.org"



# # Install BiocManager
RUN sed -i 's/-Werror=format-security//g' /usr/local/lib/R/etc/Makeconf
RUN Rscript -e "remotes::install_github('jlaffy/scalop', dependencies=TRUE)"
RUN Rscript -e "install.packages('NMF')"


RUN Rscript -e 'BiocManager::install("readr", dependencies = TRUE)'


# Install Seurat Wrappers
RUN wget https://github.com/satijalab/seurat/archive/refs/heads/seurat5.zip -O /opt/seurat-v5.zip
RUN wget https://github.com/satijalab/seurat-data/archive/refs/heads/seurat5.zip -O /opt/seurat-data.zip
RUN wget https://github.com/satijalab/seurat-wrappers/archive/refs/heads/seurat5.zip -O /opt/seurat-wrappers.zip

RUN Rscript -e "devtools::install_local('/opt/seurat-v5.zip')"
RUN Rscript -e "devtools::install_local('/opt/seurat-data.zip')"
RUN Rscript -e "devtools::install_local('/opt/seurat-wrappers.zip')"
RUN Rscript -e "devtools::install_github('BlakeRMills/MetBrewer')"
RUN Rscript -e "devtools::install_github('sjmgarnier/viridis')"





# Build tools so pip can compile wheels if needed (arm64, etc.)
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential python3-dev \
  && rm -rf /var/lib/apt/lists/*

# Create venv and install Python deps (including python-igraph + leidenalg)
RUN python3 -m venv /opt/venv \
  && /opt/venv/bin/pip install --upgrade pip setuptools wheel \
  && /opt/venv/bin/pip install --no-cache-dir \
       numpy pandas scipy scikit-learn matplotlib seaborn \
       jupyter jupyter-cache papermill \
       anndata scanpy session_info scSpectra metatime celltypist \
       python-igraph==0.11.* leidenalg==0.10.* \
  && /opt/venv/bin/python - <<'PY'
import importlib, sys
for m in ("igraph","leidenalg"):
    importlib.import_module(m)
print("OK: igraph/leidenalg present in", sys.executable)
PY

# Make the venv the default for PATH and for R/reticulate
ENV PATH="/opt/venv/bin:${PATH}"
ENV RETICULATE_PYTHON=/opt/venv/bin/python



# Additional packages
RUN apt-get update && apt-get install -y libhdf5-dev && rm -rf /var/lib/apt/lists/*
RUN Rscript -e "install.packages('hdf5r')"


# Java + Fortran 
RUN apt-get update && apt-get install -y default-jre libgfortran5

# JAGS
RUN apt-get install -y jags



RUN apt-get update && \
    apt-get install -y --no-install-recommends \
       libgsl-dev \
    && rm -rf /var/lib/apt/lists/*  


# Cleaning apt-get cache
RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/*

# 5) validate loads at build time
# RUN Rscript -e "library(SingleCellExperiment); library(SCP)"

# Command to run on container start
CMD ["bash"]

