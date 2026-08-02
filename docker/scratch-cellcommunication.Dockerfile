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
RUN Rscript -e "install.packages(c('R.utils','rmarkdown','devtools','tidyverse','magrittr','readr', 'dplyr', 'ggplot2', 'cowplot', 'remotes', 'BiocManager','reticulate', 'HGNChelper'), repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages(c('optparse', 'leiden', 'RColorBrewer', 'viridis', 'reshape2', 'scales', 'NMF', 'MetBrewer', 'colorspace', 'tibble', 'data.table'), repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages(c('stringr', 'Matrix', 'bigmemory', 'doMC', 'patchwork', 'pheatmap'), repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('BiocManager', repos = 'https://cloud.r-project.org')" && \
    R -e "BiocManager::install(c('purrr', 'tidyr', 'forcats','ComplexHeatmap', 'circlize'), ask = FALSE, update = TRUE)"

RUN Rscript -e "BiocManager::install(c('S4Vectors','DelayedMatrixStats','BiocGenerics','Biobase', 'SummarizedExperiment', 'AnnotationDbi', 'org.Hs.eg.db'), ask=FALSE, update=TRUE)"

RUN Rscript -e "BiocManager::install(c( \
    'HDF5Array','rhdf5','rhdf5lib', \
    'SingleCellExperiment', \
    'GOSemSim','MatrixGenerics','treeio','DOSE','ggtree','enrichplot', \
    'clusterProfiler','DirichletMultinomial','rtracklayer','GenomicFeatures', \
    'BSgenome','ensembldb','TFBSTools', 'glmGamPoi',  \
    'BSgenome.Hsapiens.UCSC.hg38','EnsDb.Hsapiens.v86', 'Homo.sapiens'), \
  ask=FALSE, update=FALSE )"

# Setting repository URL
ARG R_REPO="http://cran.us.r-project.org"

# # Install BiocManager
RUN sed -i 's/-Werror=format-security//g' /usr/local/lib/R/etc/Makeconf
RUN Rscript -e "remotes::install_github('jlaffy/scalop', dependencies=TRUE)"
RUN Rscript -e "install.packages('NMF')"


RUN R -q -e "install.packages( \
    c('harmony','future','future.apply','knitr'), \
    repos='https://cloud.r-project.org')"

RUN Rscript -e "remotes::install_github('theislab/kBET')"
RUN Rscript -e "remotes::install_github('saezlab/liana')"

# Install Bioconductor packages
RUN R -q -e "BiocManager::install(c('batchelor','BiocParallel','BiocSingular'), ask=FALSE, update=FALSE)"

# SeuratWrappers is GitHub-only
# RUN R -q -e 'remotes::install_github("satijalab/seurat-wrappers", upgrade="never", dependencies=TRUE)'
RUN R -q -e "remotes::install_github('futureverse/future', ref='develop')"


# RUN R -q -e "install.packages('https://cran.r-project.org/src/contrib/Archive/SeuratObject/SeuratObject_5.2.0.tar.gz', repos=NULL, type='source')"

# Install Seurat Wrappers
# RUN wget https://github.com/satijalab/seurat/archive/refs/heads/seurat5.zip -O /opt/seurat-v5.zip
RUN wget https://github.com/satijalab/seurat/archive/refs/tags/v5.3.0.zip -O /opt/seurat-v5.zip
RUN wget https://github.com/satijalab/seurat-data/archive/refs/heads/seurat5.zip -O /opt/seurat-data.zip
RUN wget https://github.com/satijalab/seurat-wrappers/archive/refs/heads/seurat5.zip -O /opt/seurat-wrappers.zip


RUN Rscript -e "install.packages('Seurat')"

# RUN Rscript -e "devtools::install_local('/opt/seurat-v5.zip')"
RUN Rscript -e "devtools::install_local('/opt/seurat-data.zip')"
RUN Rscript -e "devtools::install_local('/opt/seurat-wrappers.zip')"

# RUN R -q -e "install.packages('https://cran.r-project.org/src/contrib/Archive/SeuratObject/SeuratObject_5.2.0.tar.gz', repos=NULL, type='source')"
# RUN R -q -e "remotes::install_version('SeuratObject', version='5.2.0', repos='https://cloud.r-project.org', force = TRUE)"
RUN R -q -e "install.packages('SeuratObject', version = '4.1.4')"


RUN Rscript -e "devtools::install_github('saeyslab/nichenetr')"
RUN Rscript -e "devtools::install_github('jinworks/CellChat')"
RUN Rscript -e "devtools::install_github('immunogenomics/presto')"

# https://github.com/satijalab/seurat-object/archive/refs/tags/v5.2.0.zip


# # Create and activate virtual environment
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:${PATH}"
# Ensure both reticulate and Quarto/Jupyter use THIS Python (stops "Downloading cpython...")
ENV RETICULATE_PYTHON="/opt/venv/bin/python" \
    QUARTO_PYTHON="/opt/venv/bin/python"


# --- Upgrade installers and install the Python stack you need ---
# NOTE: 'session-info' is the PyPI name; it imports as session_info in code.
#       Install igraph + leidenalg for Seurat's Leiden via reticulate.
RUN pip install --no-cache-dir --upgrade pip setuptools wheel \
 && pip install --no-cache-dir \
      numpy pandas scipy scikit-learn \
      matplotlib seaborn \
      umap-learn \
      jupyter jupyter-cache papermill ipykernel \
      anndata scanpy scSpectra metatime \
      igraph leidenalg \
      session-info


# Additional packages
RUN apt-get install -y libhdf5-dev
RUN Rscript -e "install.packages('hdf5r')"
RUN Rscript -e "install.packages('ggraph')"

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

