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
    libxml2-dev


# Updating quarto to Quarto v1.4.553
RUN wget https://github.com/quarto-dev/quarto-cli/releases/download/v1.4.553/quarto-1.4.553-linux-amd64.deb -O quarto-1.4.553-linux-amd64.deb
RUN dpkg -i quarto-1.4.553-linux-amd64.deb

# # Install remotes package
# RUN R -e "install.packages('remotes')"

# Install Python3
# RUN apt-get install -y \
#     python3 \
#     python3-pip
RUN apt-get update && apt-get install -y python3 python3-pip python3-venv

# Install fundamental R packages
ARG R_DEPS="c(\
    'tidyverse', \
    'devtools', \
    'rmarkdown', \
    'patchwork', \
    'BiocManager', \
    'remotes', \
    'optparse', \
    'R.utils', \
    'here', \
    'HGNChelper', \
    'reticulate' \
    )"

ARG DEV_DEPS="c(\
    'bnprks/BPCells', \
    'cellgeni/sceasy', \
    'zhanghao-njmu/SCP', \
    'immunogenomics/presto' \
    )"
    
ARG WEB_DEPS="c(\
    'shiny', \
    'DT', \
    'kable', \
    'kableExtra', \
    'flexdashboard', \
    'plotly' \
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
    'terra', \ 
    'ggrastr', \
    'Rsamtools', \
    'UCell', \
    'DropletUtils', \
    'MAST', \
    'DESeq2', \
    'batchelor', \
    'scDblFinder' \
    )"

# Setting repository URL
ARG R_REPO="http://cran.us.r-project.org"

# Caching R-lib on the building process
RUN Rscript -e "install.packages(${R_DEPS}, Ncpus = 8, repos = '${R_REPO}', clean = TRUE)"
RUN Rscript -e "install.packages(${WEB_DEPS}, Ncpus = 8, repos = '${R_REPO}', clean = TRUE)"

# Install BiocManager
RUN Rscript -e "BiocManager::install(${R_BIOC_DEPS})"
# RUN Rscript -e 'BiocManager::install("readr", dependencies = TRUE)'
# RUN Rscript -e 'BiocManager::install("dplyr", dependencies = TRUE)'
# RUN Rscript -e 'BiocManager::install("ggplot2", dependencies = TRUE)'
# RUN Rscript -e 'BiocManager::install("Seurat", dependencies = TRUE)'
# RUN Rscript -e 'BiocManager::install("DT", dependencies = TRUE)'
# RUN Rscript -e 'BiocManager::install("SingleCellExperiment", dependencies = TRUE)'
# RUN Rscript -e 'BiocManager::install("scDblFinder", dependencies = TRUE, force = TRUE)'
# RUN Rscript -e 'BiocManager::install("lpsymphony", dependencies = TRUE, force = TRUE)'
# RUN Rscript -e 'BiocManager::install("IHW", dependencies = TRUE, force = TRUE)'
# RUN Rscript -e 'BiocManager::install("scp", dependencies = TRUE, force = TRUE)'
RUN Rscript -e 'BiocManager::install(c("DOSE", "enrichplot", "clusterProfiler"), force = TRUE)'

# Install Seurat Wrappers
RUN wget https://github.com/satijalab/seurat/archive/refs/heads/seurat5.zip -O /opt/seurat-v5.zip
RUN wget https://github.com/satijalab/seurat-data/archive/refs/heads/seurat5.zip -O /opt/seurat-data.zip
RUN wget https://github.com/satijalab/seurat-wrappers/archive/refs/heads/seurat5.zip -O /opt/seurat-wrappers.zip
# Install SCP package from GitHub
# RUN R -e "remotes::install_github('zhanghao-njmu/SCP', upgrade = 'always', force = TRUE, dependencies = TRUE)"

# # Download the Miniconda installer
# RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
#     chmod +x /tmp/miniconda.sh && \
#     /tmp/miniconda.sh -b -p /opt/miniconda && \
#     rm /tmp/miniconda.sh

# # Update PATH environment variable
# ENV PATH=/opt/miniconda/bin:$PATH


RUN Rscript -e 'remotes::install_github("ctlab/fgsea")'

# # Install R packages
# RUN Rscript -e 'install.packages("remotes")' && \
#     Rscript -e 'remotes::install_github("zhanghao-njmu/SCP", upgrade = "always", force = TRUE, quiet = TRUE)' \
#     Rscript -e 'SCP::PrepareEnv( \
#             miniconda_repo = "https://mirrors.bfsu.edu.cn/anaconda/miniconda", \
#             pip_options = "-i https://pypi.tuna.tsinghua.edu.cn/simple")'

# Set the conda binary path and prepare the SCP environment
# RUN Rscript -e 'options(reticulate.conda_binary = "/opt/miniconda/bin/conda"); SCP::PrepareEnv(force = TRUE)'


# RUN Rscript -e 'renv::activate()'
# RUN wget https://github.com/zhanghao-njmu/SCP/archive/refs/heads/main.zip -O /opt/SCP.zip
# RUN unzip -o /opt/SCP.zip -d /opt/SCP
# RUN Rscript -e "devtools::install_local('/opt/SCP/SCP-main')"


# RUN Rscript -e "devtools::install_local('/opt/SCP.zip')"
RUN Rscript -e "devtools::install_local('/opt/seurat-v5.zip')"
RUN Rscript -e "devtools::install_local('/opt/seurat-data.zip')"
RUN Rscript -e "devtools::install_local('/opt/seurat-wrappers.zip')"
# RUN Rscript -e 'devtools::install_github("zhanghao-njmu/SCP")'
# RUN Rscript -e 'remotes::install_github("zhanghao-njmu/SCP", dependencies = TRUE, force = TRUE)'


# Install packages on Github
RUN Rscript -e "devtools::install_github(${DEV_DEPS})"


# Create and activate virtual environment
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Upgrade pip and install Python packages in venv
RUN pip install --upgrade pip && \
    pip install numpy pandas scikit-learn matplotlib seaborn jupyter jupyter-cache papermill

# # Download and install SCP manually
# RUN wget https://github.com/zhanghao-njmu/SCP/archive/refs/heads/main.zip -O /opt/SCP.zip
# RUN unzip /opt/SCP.zip -d /opt/SCP

# # Force SCP to use system Python (modify reticulate before installing SCP)
# RUN Rscript -e "install.packages('reticulate', repos='http://cran.us.r-project.org')" && \
#     Rscript -e "library(reticulate); use_python('/usr/bin/python3', required=TRUE); options(reticulate.conda_binary=NULL, SCP_env_name=NULL)" && \
#     Rscript -e "devtools::install_local('/opt/SCP/SCP-main')"
    
# Install Python packages for data science
# RUN python3 -m pip install --no-cache-dir numpy pandas scikit-learn matplotlib seaborn jupyter
# RUN python3 -m pip install --no-cache-dir jupyter-cache
# RUN python3 -m pip install --no-cache-dir papermill

# Install Python packages for data science
RUN python3 -m pip install --no-cache-dir numpy pandas scikit-learn matplotlib seaborn jupyter
RUN python3 -m pip install --no-cache-dir jupyter-cache
RUN python3 -m pip install --no-cache-dir papermill


# Additional packages
RUN apt-get install -y libhdf5-dev
RUN Rscript -e "install.packages('hdf5r')"

#RUN Rscript -e "remotes::install_version('Matrix', version = '1.6-1')"
#RUN Rscript -e "install.packages('SeuratObject')"
#RUN Rscript -e "install.packages('scCustomize')"

# Cleaning apt-get cache
RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/*

# Command to run on container start
CMD ["bash"]
