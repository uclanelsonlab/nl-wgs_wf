FROM ubuntu:20.04

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    build-essential \
    libbz2-dev \
    libcurl4-openssl-dev \
    lib32ncurses5-dev \
    liblzma-dev \
    libssl-dev \
    git \
    openjdk-17-jdk \
    python3 \
    python3-pip \
    perl \
    r-base \
    r-base-dev \
    tabix \
    unzip \
    zlib1g \
    wget \ 
    zlib1g-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2 && \
    tar -xf bcftools-1.16.tar.bz2 && \
    cd bcftools-1.16 && \
    ./configure && \
    make && \
    make install

ENV PATH=${PATH}:/bcftools-1.16/

RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static -P /usr/local/bin/ && \
    mv /usr/local/bin/bedtools.static /usr/local/bin/bedtools && \
    chmod a+x /usr/local/bin/bedtools

RUN git clone https://github.com/mquinodo/AutoMap.git && \
    chmod +x /AutoMap/AutoMap_v1.3.sh

ENV PATH=${PATH}:/AutoMap/