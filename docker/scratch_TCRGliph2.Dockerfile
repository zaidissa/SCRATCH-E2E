FROM --platform=linux/amd64 rocker/verse:4.5.1

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
    libgsl-dev \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------
# Quarto
# ---------------------------------------------------------
RUN wget -q https://github.com/quarto-dev/quarto-cli/releases/download/v1.4.553/quarto-1.4.553-linux-amd64.deb -O /tmp/quarto.deb \
 && dpkg -i /tmp/quarto.deb \
 && rm -f /tmp/quarto.deb

# ---------------------------------------------------------
# R bootstrap
# ---------------------------------------------------------
RUN Rscript -e "options(repos = c(CRAN='https://cran.rstudio.com')); install.packages(c('remotes','devtools','BiocManager'))"

# ---------------------------------------------------------
# Core R packages needed for GLIPH2 report / QMD
# ---------------------------------------------------------
RUN Rscript -e "\
  options(repos = c(CRAN = 'https://cran.rstudio.com')); \
  ncpus <- max(1L, parallel::detectCores() - 1L); \
  pk <- c( \
    'Seurat','SeuratObject', \
    'data.table','dplyr','tidyr','stringr','ggplot2', \
    'forcats','scales','glue','knitr','kableExtra','patchwork','gsl' \
  ); \
  inst <- rownames(installed.packages()); \
  need <- setdiff(pk, inst); \
  if (length(need)) install.packages(need, Ncpus = ncpus, dependencies = TRUE) \
"

# ---------------------------------------------------------
# Bioconductor packages needed
# ---------------------------------------------------------
RUN Rscript -e "\
  ncpus <- max(1L, parallel::detectCores() - 1L); \
  pk <- c('ComplexHeatmap','circlize','BiocFileCache'); \
  inst <- rownames(installed.packages()); \
  need <- setdiff(pk, inst); \
  if (length(need)) BiocManager::install(need, ask = FALSE, update = FALSE, dependencies = TRUE) \
"

# ---------------------------------------------------------
# immGLIPH
# ---------------------------------------------------------
RUN Rscript -e "remotes::install_github('BorchLab/immGLIPH', upgrade = 'never')"

# ---------------------------------------------------------
# Optional compatibility path placeholders
# Keeps path conventions similar without installing conda
# ---------------------------------------------------------
ENV CONDA_DIR=/opt/conda
RUN mkdir -p /opt/conda/envs/tcrenv/bin && \
    ln -sf /usr/bin/python3 /opt/conda/envs/tcrenv/bin/python && \
    ln -sf /usr/bin/python3 /opt/conda/envs/tcrenv/bin/python3 && \
    ln -sf /usr/bin/pip3 /opt/conda/envs/tcrenv/bin/pip

ENV PATH=/opt/conda/envs/tcrenv/bin:${PATH}
ENV RETICULATE_PYTHON=/opt/conda/envs/tcrenv/bin/python

# ---------------------------------------------------------
# Validation
# ---------------------------------------------------------
RUN Rscript -e "suppressPackageStartupMessages({ \
  library(Seurat); \
  library(SeuratObject); \
  library(data.table); \
  library(dplyr); \
  library(tidyr); \
  library(stringr); \
  library(ggplot2); \
  library(forcats); \
  library(scales); \
  library(glue); \
  library(knitr); \
  library(kableExtra); \
  library(ComplexHeatmap); \
  library(circlize); \
  library(patchwork); \
  library(immGLIPH); \
}); cat('R package validation OK\\n')"

RUN quarto --version
RUN python3 --version

# ---------------------------------------------------------
# Cleanup
# ---------------------------------------------------------
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

CMD ["/bin/bash"]