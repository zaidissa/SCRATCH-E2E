# Use a specific version of Ubuntu as the base image
FROM --platform=linux/x86_64 rocker/verse:latest

# Set the working directory inside the container
WORKDIR /opt

# Timezone settings
ENV TZ=US/Central
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
  && echo $TZ > /etc/timezone


ENV R_REPOS="https://packagemanager.posit.co/cran/latest"

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
    libxt-dev \

    jags \
    python3 \
    python3-pip \
    python3-venv \
    libhdf5-dev \
    libgsl-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


# Updating quarto to Quarto v1.4.553
RUN wget https://github.com/quarto-dev/quarto-cli/releases/download/v1.4.553/quarto-1.4.553-linux-amd64.deb -O quarto-1.4.553-linux-amd64.deb
RUN dpkg -i quarto-1.4.553-linux-amd64.deb


# Install core R packages
RUN Rscript -e "install.packages(c( \
  'R.utils', \
  'rmarkdown', \
  'devtools', \
  'tidyverse', \
  'readr', \
  'dplyr', \
  'ggplot2', \
  'cowplot', \
  'remotes', \
  'BiocManager',\
  'reticulate', \
  'HGNChelper', \
  'hdf5r', \
  'optparse' \
  ), repos='${R_REPOS}')"

# Bioconductor packages
RUN Rscript -e "BiocManager::install(c( \
  'S4Vectors', \
  'DelayedMatrixStats', \
  'BiocGenerics',\
  'Biobase', \
  'SummarizedExperiment', \
  'AnnotationDbi', \
  'org.Hs.eg.db' \
  ), ask=FALSE, update=TRUE)" \
  && Rscript -e "BiocManager::install(c( \
    'HDF5Array', \
    'rhdf5', \
    'rhdf5lib', \
    'SingleCellExperiment', \
    'GOSemSim', \
    'MatrixGenerics', \
    'treeio', \
    'DOSE',\
    'ggtree',\
    'enrichplot', \
    'clusterProfiler',\
    'DirichletMultinomial',\
    'rtracklayer',\
    'GenomicFeatures', \
    'BSgenome',\
    'ensembldb',\
    'TFBSTools', \
    'BSgenome.Hsapiens.UCSC.hg38', \
    'EnsDb.Hsapiens.v86' \
  ), ask=FALSE, update=FALSE )"


RUN Rscript -e 'remotes::install_github("ctlab/fgsea")'

# Seurat and wrappers
RUN wget https://github.com/satijalab/seurat/archive/refs/tags/v5.3.1.zip -O /opt/seurat-v5.zip \
  && wget https://github.com/satijalab/seurat-data/archive/refs/heads/seurat5.zip -O /opt/seurat-data.zip \
  && wget https://github.com/satijalab/seurat-wrappers/archive/refs/heads/seurat5.zip -O /opt/seurat-wrappers.zip \
  && Rscript -e "devtools::install_local('/opt/seurat-v5.zip')" \
  && Rscript -e "devtools::install_local('/opt/seurat-data.zip')" \
  && Rscript -e "devtools::install_local('/opt/seurat-wrappers.zip')" \
  && rm /opt/seurat-v5.zip /opt/seurat-data.zip /opt/seurat-wrappers.zip \
  && Rscript -e "devtools::install_github('satijalab/azimuth', ref = 'master', dependencies=TRUE, upgrade='never')"

# Install SCP package from GitHub
RUN Rscript -e "devtools::install_github('PaulingLiu/ROGUE', dependencies = TRUE, upgrade = 'never', force = TRUE)" \
  && Rscript -e "devtools::install_github('zhanghao-njmu/SCP', dependencies = TRUE, upgrade = 'never', force = TRUE)" \
  && Rscript -e "devtools::install_github('cellgeni/sceasy', dependencies = TRUE, upgrade = 'never', force = TRUE)" \
  # uninstall reactome.db that SCP pulls (1.7Gb)
  && Rscript -e "remove.packages('reactome.db')"

# Install annotables
RUN Rscript -e "devtools::install_github('stephenturner/annotables')" \
  && Rscript -e "devtools::install_github('miccec/yaGST', dependencies = TRUE, upgrade = 'never')"

# Create and activate a virtual environment before installing Python packages 
# Note the no-gpu version of torch to avoid ~5Gb dependencies and keep image size smaller
RUN python3 -m venv /opt/venv \
    && /opt/venv/bin/pip install --no-cache-dir --upgrade pip \
    && /opt/venv/bin/pip install --no-cache-dir torch --index-url https://download.pytorch.org/whl/cpu \
    && /opt/venv/bin/pip install --no-cache-dir numpy pandas scikit-learn \
     matplotlib seaborn jupyter jupyter-cache papermill anndata scanpy scipy \
     session_info scSpectra metatime celltypist \
    && ln -s /opt/venv/bin/python /usr/local/bin/python \
    && ln -s /opt/venv/bin/pip /usr/local/bin/pip

ENV PATH="/opt/venv/bin:$PATH"
ENV CELLTYPIST_FOLDER=/opt/celltypist

# Command to run on container start
CMD ["bash"]

