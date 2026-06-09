FROM python:3.13-slim-bookworm

LABEL maintainer="Waldeyr Mendes Cordeiro da Silva"

SHELL ["/bin/bash", "-c"]

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/usr/local/bin:/opt/ncbi-igblast-1.22.0/bin:${PATH}"

# Install system dependencies
# fastqc: quality control (automatically installs Java dependency)
# ea-utils: provides fastq-join directly
# perl: required for running prinseq-lite
# wget, tar, ca-certificates: tools for downloading and unpacking components
# libgomp1: required by the compiled igblastp binary
RUN apt-get update && apt-get install -y --no-install-recommends \
    fastqc \
    ea-utils \
    perl \
    wget \
    tar \
    libgomp1 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install prinseq-lite (standalone perl script)
RUN wget --no-check-certificate https://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz -O /tmp/prinseq.tar.gz \
    && tar -xzf /tmp/prinseq.tar.gz -C /tmp \
    && mv /tmp/prinseq-lite-0.20.4/prinseq-lite.pl /usr/local/bin/prinseq-lite \
    && chmod +x /usr/local/bin/prinseq-lite \
    && rm -rf /tmp/prinseq*

# Install NCBI IgBlast v1.22.0
RUN mkdir -p /opt \
    && wget ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/1.22.0/ncbi-igblast-1.22.0-x64-linux.tar.gz -O /tmp/igblast.tar.gz \
    && tar -xzf /tmp/igblast.tar.gz -C /opt \
    && chmod +x /opt/ncbi-igblast-1.22.0/bin/* \
    && rm -f /tmp/igblast.tar.gz

# Set up working directory inside the container
WORKDIR /app

# Copy requirements and install python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy all project files into the container
COPY . .
RUN chmod +x run.sh install.sh programs/*.py

# Validate installation: runs pre-flight checks at build time
# Build fails here if any required tool is missing from the image
RUN python3 programs/test_attila.py TestInstallation -v

# Command to run bash by default
CMD ["/bin/bash"]
