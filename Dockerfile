# x86_64
FROM debian:stable-slim

ENV PATH="/pipeline/bin:$PATH"
ENV SAMTOOLS_VERSION="1.21"

# install system dependencies
RUN DEBIAN_FRONTEND=noninteractive apt-get update && DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y --no-install-recommends install ca-certificates tzdata apt-utils wget curl bzip2 make xsltproc gcc g++ pkg-config zlib1g-dev libxml2-dev python3 python3-pip python3-distutils unzip && apt-get clean && rm -rf /var/lib/apt/lists/*
# isntall samtools
RUN mkdir -p /pipeline/samtools && wget -qO- https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 | tar -xjvf - -C /pipeline/samtools --strip-components 1 && cd /pipeline/samtools/ && ./configure --without-curses --disable-bz2 --disable-lzma --prefix=/pipeline/ && make -j
# install hisat2-3n
RUN wget -qO- wget https://github.com/DaehwanKimLab/hisat2/archive/refs/heads/hisat-3n.tar.gz | tar -xvz -C /pipeline/ && cd /pipeline/hisat2-hisat-3n/ && make -j
# install package by python-pip
RUN python3 -m pip install --break-system-packages --no-cache-dir snakemake cutseq polars
# clean up and reduce size
RUN apt-get purge -y wget bzip2 make xsltproc gcc g++ pkg-config unzip && apt-get clean && rm -rf /var/lib/apt/lists/*

COPY ./bin /pipeline/bin
COPY ./VERSION ./Snakefile ./default.yaml ./entrypoint /pipeline/
RUN chmod +x /pipeline/entrypoint
ENTRYPOINT ["/pipeline/entrypoint"]
