# Use a specific version of Ubuntu as the base image
FROM --platform=linux/amd64 rocker/verse:4.4.1

WORKDIR /opt

ENV TZ=US/Central \
    DEBIAN_FRONTEND=noninteractive \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    CELLTYPIST_FOLDER=/opt/celltypist \
    PATH=/opt/venv/bin:$PATH

RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# ---------------------------------------------------------
# System dependencies
# ---------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    pkg-config \
    make \
    g++ \
    gfortran \
    git \
    wget \
    curl \
    unzip \
    ca-certificates \
    bzip2 \
    perl \
    python3 \
    python3-pip \
    python3-venv \
    python3-dev \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libgit2-dev \
    libjpeg-dev \
    libpng-dev \
    libtiff5-dev \
    zlib1g-dev \
    libxt-dev \
    libopenblas-dev \
    liblapack-dev \
    libgfortran5 \
    libhdf5-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    imagemagick \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------
# Quarto
# ---------------------------------------------------------
# Quarto (install + cleanup)
RUN wget -q https://github.com/quarto-dev/quarto-cli/releases/download/v1.4.553/quarto-1.4.553-linux-amd64.deb -O /tmp/quarto.deb \
 && dpkg -i /tmp/quarto.deb \
 && rm -f /tmp/quarto.deb

# ---- R packages (ordered) ----
ARG R_REPO="http://cran.us.r-project.org"

RUN apt-get update && apt-get install -y python3 python3-pip python3-venv python3-dev build-essential

# ---------------------------------------------------------
# R bootstrap
# ---------------------------------------------------------
ARG R_DEPS="c('tidyverse','devtools','rmarkdown','patchwork','BiocManager','remotes','optparse','R.utils','here','HGNChelper','reticulate')"
ARG WEB_DEPS="c('shiny','DT','kableExtra','flexdashboard','plotly')"

RUN Rscript -e "remotes::install_version('xfun', version = '0.55', repos = 'https://cran.rstudio.com', dependencies = TRUE)"
RUN Rscript -e "install.packages(c('knitr','rmarkdown','htmlTable','Hmisc','kableExtra'), dependencies = TRUE) \
"

# ---------------------------------------------------------
# CRAN packages actually used in the TCR QMDs
# ---------------------------------------------------------
RUN Rscript -e "\
  options(repos = c(CRAN = 'https://cran.rstudio.com')); \
  ncpus <- max(1L, parallel::detectCores() - 1L); \
  pk <- c( \
    'Seurat','SeuratObject','Matrix', \
    'dplyr','data.table','tibble','stringr','readr','tidyr','purrr','magrittr','glue','forcats', \
    'ggplot2','patchwork','scales','RColorBrewer','ggalluvial', \
    'knitr','rmarkdown','kableExtra','htmltools','htmlwidgets','jsonlite', \
    'circlize' \
  ); \
  inst <- rownames(installed.packages()); \
  need <- setdiff(pk, inst); \
  if (length(need)) install.packages(need, Ncpus = ncpus, dependencies = TRUE) \
"

# ---------------------------------------------------------
# Bioconductor packages actually used
# ---------------------------------------------------------
RUN Rscript -e "\
  ncpus <- max(1L, parallel::detectCores() - 1L); \
  pk <- c('ComplexHeatmap','scRepertoire'); \
  inst <- rownames(installed.packages()); \
  need <- setdiff(pk, inst); \
  if (length(need)) BiocManager::install(need, ask = FALSE, update = TRUE, dependencies = TRUE) \
"

# ---------------------------------------------------------
# Python via micromamba
# ---------------------------------------------------------
ENV CONDA_DIR=/opt/conda
ENV PATH=${CONDA_DIR}/bin:${PATH}

RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p ${CONDA_DIR} && \
    rm -f /tmp/miniconda.sh && \
    ${CONDA_DIR}/bin/conda clean -afy

# remove defaults completely and use conda-forge only
RUN conda config --system --remove-key channels || true && \
    conda config --system --add channels conda-forge && \
    conda config --system --set channel_priority strict && \
    conda config --system --set auto_update_conda false && \
    conda config --system --set show_channel_urls true

RUN conda create -y -n tcrenv --override-channels -c conda-forge \
      python=3.10 \
      pip \
      numpy \
      pandas \
      scipy \
      matplotlib \
      seaborn \
      networkx \
      umap-learn \
      scikit-learn \
      statsmodels \
      numba \
      pytables \
      scanpy \
      python-igraph \
      leidenalg \
      louvain \
      pyyaml \
      biopython \
      faiss-cpu \
      h5py \
      cairosvg && \
    conda clean -afy
    

ENV PATH=/opt/conda/envs/tcrenv/bin:${PATH}
ENV RETICULATE_PYTHON=/opt/conda/envs/tcrenv/bin/python

# ---------------------------------------------------------
# Python-side TCR tools
# ---------------------------------------------------------
# RUN /opt/conda/envs/tcrenv/bin/pip install --no-cache-dir tcrdist3


RUN Rscript -e "install.packages('gsl', repos='https://cloud.r-project.org')"


# CoNGA from source repo
RUN git clone https://github.com/phbradley/conga.git /opt/tools/conga && \
    /opt/conda/envs/tcrenv/bin/pip install --no-cache-dir -e /opt/tools/conga

# GIANA from source repo
RUN git clone https://github.com/s175573/GIANA.git /opt/tools/GIANA

# GLIPH from source repo (command-line scripts in bin/)
RUN git clone https://github.com/immunoengineer/gliph.git /opt/tools/gliph

ENV PATH=/opt/tools/gliph/bin:${PATH}
ENV PATH=/opt/tools/conga/scripts:${PATH}

RUN apt-get update && apt-get install -y --no-install-recommends \
    libgsl-dev \
    && rm -rf /var/lib/apt/lists/*


# RUN /opt/conda/envs/tcrenv/bin/pip install --no-cache-dir \
#     tcrdist3 \
#     torch \
#     pyro-ppl \
#     gseapy \
#     scvi-tools
RUN /opt/conda/envs/tcrenv/bin/pip install --no-cache-dir \
    "torch>=2.4.1" \
    --index-url https://download.pytorch.org/whl/cpu


RUN /opt/conda/envs/tcrenv/bin/pip install --no-cache-dir \
    tcrdist3 \
    "pyro-ppl>=1.9.1" \
    "scvi-tools>=1.3.0" \
    "mpltern>=1.0.4" \
    "gseapy>=1.1.4" \
    "tqdm>=4.66.5" \
    "torch-geometric>=2.6.1" \
    "torchmetrics>=1.6.1" \
    "daft"

# Force Conda to repair llvmlite and numba so the C++ compiler links correctly
# RUN /opt/conda/bin/conda install -n tcrenv -y -c conda-forge numba llvmlite
RUN /opt/conda/envs/tcrenv/bin/pip uninstall -y numba llvmlite

# 2. FORCE Conda to physically reinstall the working binaries
RUN /opt/conda/bin/conda install -n tcrenv -y -c conda-forge --force-reinstall numba llvmlite
# --------------------


RUN git clone https://github.com/nceglia/tcri.git /opt/tools/tcri && \
    cd /opt/tools/tcri && \
    /opt/conda/envs/tcrenv/bin/pip install --no-cache-dir .
RUN Rscript -e "BiocManager::install('scRepertoire', ask = FALSE, update = FALSE)"

# RUN R -e "BiocManager::install(c('BiocFileCache'), ask=FALSE)" \
#     && R -e "remotes::install_github('BorchLab/immGLIPH', upgrade='never')"

# ---------------------------------------------------------
# Validation
# ---------------------------------------------------------
# RUN Rscript -e "\
#   suppressPackageStartupMessages({ \
#     library(Seurat); \
#     library(SeuratObject); \
#     library(scRepertoire); \
#     library(ComplexHeatmap); \
#     library(circlize); \
#     library(ggplot2); \
#     library(dplyr); \
#     library(data.table); \
#     library(tidyr); \
#     library(stringr); \
#     library(readr); \
#     library(purrr); \
#     library(forcats); \
#     library(glue); \
#     library(knitr); \
#     library(kableExtra); \
#     library(ggalluvial); \
#   }); \
#   cat('R package validation OK\\n') \
# "

# RUN /opt/micromamba/envs/tcrenv/bin/python - <<'PY'
# import importlib

# mods = {
#     "numpy": "numpy",
#     "pandas": "pandas",
#     "scipy": "scipy",
#     "networkx": "networkx",
#     "umap": "umap",
#     "scanpy": "scanpy",
#     "tcrdist": "tcrdist",
#     "sklearn": "sklearn",
#     "Bio": "Bio",
#     "faiss": "faiss",
#     "conga": "conga",
# }
# for label, mod in mods.items():
#     importlib.import_module(mod)
#     print(f"{label} OK")

# print("Python validation OK")
# PY

# RUN test -f /opt/tools/GIANA/GIANA.py
# RUN test -f /opt/tools/gliph/bin/gliph-group-discovery.pl
# RUN test -f /opt/tools/gliph/bin/gliph-group-scoring.pl

RUN quarto --version

# ---------------------------------------------------------
# Cleanup
# ---------------------------------------------------------
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

CMD ["/bin/bash"]

