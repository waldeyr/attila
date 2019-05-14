FROM fedora:29
MAINTAINER "Waldeyr Mendes Cordeiro da Silva"
SHELL ["/bin/bash", "-c"]
RUN mkdir attila 
WORKDIR "attila"
RUN whoami
RUN pwd
RUN dnf update -y && dnf install which nano git autoconf automake libtool gcc-c++ wget java-1.8.0-openjdk-devel R libnsl -y
RUN Rscript -e 'install.packages("ggplot2", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("scales",  repos = "http://cran.us.r-project.org")'
RUN wget https://github.com/s-andrews/FastQC/archive/v0.11.8.tar.gz -O fastqc.tar.gz && tar -xf fastqc.tar.gz && chmod +x FastQC-0.11.8/fastqc && ln -sf FastQC-0.11.8/fastqc /usr/local/bin/fastqc && rm fastqc.tar.gz -f
RUN wget https://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz -O prinseq-lite.tar.gz && tar -xf prinseq-lite.tar.gz && chmod +x prinseq-lite-0.20.4/prinseq* && ln -sf prinseq-lite-0.20.4/prinseq* /usr/local/bin/ && rm prinseq-lite.tar.gz -f
RUN wget ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/1.14.0/ncbi-igblast-1.14.0-x64-linux.tar.gz -O igblast.tar.gz && tar -xf igblast.tar.gz && chmod +x ncbi-igblast-1.14.0/bin/* && ln -sf ncbi-igblast-1.14.0/bin/* /usr/local/bin/ && rm igblast.tar.gz -f
RUN git clone https://github.com/waldeyr/attila.git attila-v1.0 && cd attila-v1.0 && chmod +x programs/*
RUN cd /attila/attila-v1.0 && ln -sf programs/* . && ln -sf /attila/prinseq-lite-0.20.4/prinseq* . && ln -sf /attila/FastQC-0.11.8/fastqc . && ln -sf /attila/ncbi-igblast-1.14.0/bin/* .
RUN echo 'To run it, type: docker run --memory="2048m" -ti waldeyr/attilav1.0 bash'
