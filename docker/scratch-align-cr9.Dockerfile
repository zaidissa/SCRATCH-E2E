FROM --platform=linux/x86_64 rocker/verse:latest

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

# Install CellRanger 9
ARG CELLRANGER_VERSION='9.0.1'

ENV PATH=/opt/cellranger-${CELLRANGER_VERSION}:$PATH

COPY cellranger-${CELLRANGER_VERSION}.tar.gz /opt/cellranger-${CELLRANGER_VERSION}.tar.gz
RUN tar -zxvf /opt/cellranger-${CELLRANGER_VERSION}.tar.gz -C /opt && \
    rm -f /opt/cellranger-${CELLRANGER_VERSION}.tar.gz

# Install Python3
RUN apt-get update && apt-get install -y python3 python3-pip python3-venv

# Create and activate virtual environment
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Install Python packages in venv
RUN pip install --upgrade pip && \
    pip install --no-cache-dir numpy pandas scikit-learn matplotlib seaborn jupyter jupyter-cache papermill

# Cleaning apt-get cache
RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/*

# Command to run on container start
CMD ["cellranger"]

