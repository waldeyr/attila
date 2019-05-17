FROM fedora:29
MAINTAINER Waldeyr Mendes Cordeiro da Silva
SHELL ["/bin/bash", "-c"]
RUN mkdir attila 
WORKDIR "attila"

RUN dnf update -y && dnf install which nano git autoconf automake libtool gcc-c++ wget java-1.8.0-openjdk-devel R libnsl gsl-devel 'perl(App::Prove)' 'perl(Test::Builder)' -y

RUN Rscript -e 'install.packages("ggplot2", repos = "http://cran.us.r-project.org")'

RUN Rscript -e 'install.packages("scales",  repos = "http://cran.us.r-project.org")'

RUN wget https://github.com/s-andrews/FastQC/archive/v0.11.8.tar.gz -O fastqc.tar.gz && tar -xf fastqc.tar.gz && chmod +x FastQC-0.11.8/fastqc && rm fastqc.tar.gz -f

RUN wget https://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz -O prinseq-lite.tar.gz && tar -xf prinseq-lite.tar.gz && mv prinseq-lite-0.20.4/prinseq-lite.pl prinseq-lite-0.20.4/prinseq-lite && chmod +x prinseq-lite-0.20.4/prinseq* && rm prinseq-lite.tar.gz -f

RUN git clone https://github.com/ExpressionAnalysis/ea-utils.git && cd /attila/ea-utils/clipper && CC=g++ PREFIX=/usr/local make install

RUN wget ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/1.14.0/ncbi-igblast-1.14.0-x64-linux.tar.gz -O igblast.tar.gz && tar -xf igblast.tar.gz && chmod +x ncbi-igblast-1.14.0/bin/* && rm igblast.tar.gz -f

RUN git clone https://github.com/waldeyr/attila.git attila-v1.0 && cd attila-v1.0 && chmod +x programs/*

RUN ln -sf /attila/prinseq-lite-0.20.4/prinseq* . && ln -sf /attila/FastQC-0.11.8/fastqc . && ln -sf /attila/ncbi-igblast-1.14.0/bin/* .

RUN cd /usr/local/bin/ && ln -sf /attila/FastQC-0.11.8/fastqc .
RUN cd /usr/local/bin/ && ln -sf /attila/ncbi-igblast-1.14.0/bin/* .
RUN cd /attila/attila-v1.0/prgrams/ && ln -sf /attila/ncbi-igblast-1.14.0/bin/* .
RUN cd /usr/local/bin/ && ln -sf /attila/prinseq-lite-0.20.4/prinseq-lite .

RUN echo 'To run it, type: docker run -v $(whoami)/$(pwd):/attila/shared --memory="2048m" -ti waldeyr/attila:v1.0 bash'
