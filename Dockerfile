# x86_64
FROM debian:bullseye-slim

ENV PATH="/pipeline/bin:/bin:/pipeline/micromamba/bin:$PATH"
ENV MAMBA_ROOT_PREFIX="/pipeline/micromamba"
ENV MAMBA_EXE="/bin/micromamba"

# install system dependencies
RUN DEBIAN_FRONTEND=noninteractive apt-get update && DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y --no-install-recommends install tzdata apt-utils wget curl bzip2 make xsltproc gcc g++ pkg-config zlib1g-dev libxml2-dev python3 python3-distutils default-jre unzip && apt-get clean && rm -rf /var/lib/apt/lists/*
# install package by micromamba
RUN wget -qO- https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C /tmp bin/micromamba && cp /tmp/bin/micromamba $MAMBA_EXE  && $MAMBA_EXE shell init -s bash -r $MAMBA_ROOT_PREFIX && eval "$($MAMBA_EXE shell hook --shell=posix)" && micromamba install -n base -y -c conda-forge -c bioconda libzlib falco bowtie2 star samtools bedtools fastp && micromamba clean -afy
# install hisat2-3n
RUN wget -qO- wget https://github.com/DaehwanKimLab/hisat2/archive/refs/heads/hisat-3n.tar.gz | tar -xvz -C /pipeline/ && cd /pipeline/hisat2-hisat-3n/ && make -j
# install umicollapse
RUN wget -q -P /bin/ https://github.com/Daniel-Liu-c0deb0t/UMICollapse/raw/refs/heads/master/umicollapse.jar && \
  mkdir -p /bin/lib && wget -q -P /bin/lib https://repo1.maven.org/maven2/com/github/samtools/htsjdk/2.19.0/htsjdk-2.19.0.jar && wget -q -P /bin/lib https://repo1.maven.org/maven2/org/xerial/snappy/snappy-java/1.1.7.3/snappy-java-1.1.7.3.jar
# install gatk
RUN wget -q -P /tmp/ https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip && unzip /tmp/gatk-4.3.0.0.zip -d /tmp/ && cp /tmp/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar /bin/gatk.jar && rm -rf /tmp/gatk-*
# install package by python-pip
RUN pip install --no-cache-dir snakemake pysam cutseq
# clean up and reduce size
RUN apt-get purge -y wget bzip2 make xsltproc gcc g++ pkg-config unzip && apt-get clean && rm -rf /var/lib/apt/lists/*

COPY ./bin /pipeline/bin
COPY ./VERSION ./Snakefile ./default.yaml ./entrypoint ./validator /pipeline/
RUN chmod +x /pipeline/entrypoint && chmod +x /pipeline/validator
ENTRYPOINT ["/pipeline/entrypoint"]
