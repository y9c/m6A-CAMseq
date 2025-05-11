# x86_64
FROM debian:stable-slim

ENV PATH="/pipeline/bin:$PATH"
ENV SAMTOOLS_VERSION="1.21"
ENV GATK_VERSION="4.6.2.0"

# install system dependencies
RUN DEBIAN_FRONTEND=noninteractive apt-get update && DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y --no-install-recommends install ca-certificates tzdata apt-utils wget curl bzip2 unzip make gcc g++ pkg-config xsltproc zlib1g-dev libxml2-dev python3 python3-pip python3-distutils python-is-python3 && apt-get clean && rm -rf /var/lib/apt/lists/*
# isntall samtools
RUN mkdir -p /pipeline/samtools && wget -qO- https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 | tar -xjvf - -C /pipeline/samtools --strip-components 1 && cd /pipeline/samtools/ && ./configure --without-curses --disable-bz2 --disable-lzma --prefix=/pipeline/ && make -j
# install hisat2-3n
RUN wget -qO- wget https://github.com/DaehwanKimLab/hisat2/archive/refs/heads/hisat-3n.tar.gz | tar -xvz -C /pipeline/ && cd /pipeline/hisat2-hisat-3n/ && make -j
# install umicollapse
RUN mkdir -p /pipeline/UMICollapse/ && wget -q -P /pipeline/UMICollapse/ https://github.com/Daniel-Liu-c0deb0t/UMICollapse/raw/refs/heads/master/umicollapse.jar && wget -q -P /pipeline/UMICollapse/lib https://repo1.maven.org/maven2/com/github/samtools/htsjdk/2.19.0/htsjdk-2.19.0.jar && wget -q -P /pipeline/UMICollapse/lib https://repo1.maven.org/maven2/org/xerial/snappy/snappy-java/1.1.7.3/snappy-java-1.1.7.3.jar
# install gatk
RUN mkdir -p /pipeline/gatk && wget -q -P /tmp/ https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip && unzip /tmp/gatk-${GATK_VERSION}.zip -d /tmp/ && cp /tmp/gatk-${GATK_VERSION}/gatk-package-${GATK_VERSION}-local.jar /pipeline/gatk/gatk.jar && rm -rf /tmp/gatk-*
# install package by python-pip
RUN python3 -m pip install --break-system-packages --no-cache-dir snakemake cutseq polars
# clean up and reduce size
RUN apt-get purge -y wget bzip2 make xsltproc gcc g++ pkg-config unzip && apt-get clean && rm -rf /var/lib/apt/lists/*

COPY ./bin /pipeline/bin
COPY ./VERSION ./Snakefile ./default.yaml ./entrypoint /pipeline/
COPY ./external/trichromat/Snakefile ./external/trichromat/workflow_utils.py ./external/trichromat/default.yaml ./external/trichromat/config.schema.yaml /pipeline/external/trichromat/
RUN chmod +x /pipeline/entrypoint
ENTRYPOINT ["/pipeline/entrypoint"]
