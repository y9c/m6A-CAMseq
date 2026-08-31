# trichromat base image
# Includes optimized HISAT-3N, samtools, and core pipeline tools

ARG SAMTOOLS_VERSION="1.21"
ARG PICARD_VERSION="3.4.0"
ARG PYTHON_VERSION_FOR_APP="3.13"
ARG UV_BASE_IMAGE_TAG="python${PYTHON_VERSION_FOR_APP}-bookworm-slim"

# ----------- Builder Stage -----------
FROM ghcr.io/astral-sh/uv:${UV_BASE_IMAGE_TAG} AS builder

ARG SAMTOOLS_VERSION
ARG PICARD_VERSION
ARG PYTHON_VERSION_FOR_APP

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

RUN apt-get update && \
    apt-get -y --no-install-recommends install \
    ca-certificates tzdata apt-utils wget curl bzip2 unzip make cmake gcc g++ pkg-config xsltproc \
    zlib1g-dev libxml2-dev git \
    gfortran libopenblas-dev liblapack-dev \
    default-jre coreutils procps libjemalloc-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# --- Create a virtual environment for Python libraries ---
ENV APP_VENV_PATH=/opt/app_venv
RUN python${PYTHON_VERSION_FOR_APP} -m venv ${APP_VENV_PATH}

# Install Python packages
# polars>=1.43 required: sink_ipc streams join_pileup group_by to disk (~64GB->~19GB)
# and provides pl.map_batches used by the vectorized binomial p-value in filter_sites.
ENV PYTHON_PACKAGES="scipy polars>=1.43 cutseq snakemake==9.9.0 pyyaml"
RUN uv pip install --python ${APP_VENV_PATH}/bin/python --no-cache ${PYTHON_PACKAGES}

# --- Build samtools ---
WORKDIR /build/samtools_src
RUN wget -qO- https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 | tar -xjvf - --strip-components 1 && \
    ./configure --without-curses --disable-bz2 --disable-lzma && \
    make -j$(nproc) samtools

# --- Build HISAT-3N (official, hisat-3n branch) ---
WORKDIR /build/hisat2_build
RUN git clone --depth 1 --branch hisat-3n https://github.com/DaehwanKimLab/hisat2.git . && \
    make -j$(nproc)
# --- Prepare UMICollapse ---
WORKDIR /build/umicollapse_build
RUN wget -q -P ./ https://github.com/y9c/UMICollapse/releases/download/latest-prerelease/umicollapse-release.zip && \
    unzip umicollapse-release.zip && \
    rm umicollapse-release.zip

# --- Prepare picard ---
WORKDIR /build/picard_build
RUN wget -qO picard.jar https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar

# ----------- Final Stage -----------
FROM ghcr.io/astral-sh/uv:${UV_BASE_IMAGE_TAG} AS final

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC
ENV PIPELINE_HOME=/pipeline
ENV APP_VENV_PATH=/opt/app_venv

# PATH order: venv python, then other tool locations
ENV PATH="${APP_VENV_PATH}/bin:${PIPELINE_HOME}/bin:${PIPELINE_HOME}/samtools:${PIPELINE_HOME}/hisat2-hisat-3n:/usr/local/bin:$PATH"

RUN apt-get update && \
    apt-get -y --no-install-recommends install \
    ca-certificates tzdata default-jre zlib1g libxml2 libjemalloc2 libgomp1 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p ${PIPELINE_HOME}/bin \
             ${PIPELINE_HOME}/UMICollapse \
             ${PIPELINE_HOME}/picard \
             ${PIPELINE_HOME}/samtools \
             ${PIPELINE_HOME}/hisat2-hisat-3n \
             ${APP_VENV_PATH}

# Copy the application's Python virtual environment
COPY --from=builder ${APP_VENV_PATH} ${APP_VENV_PATH}

# Copy tools
COPY --from=builder /build/samtools_src/samtools ${PIPELINE_HOME}/samtools/samtools
COPY --from=builder /build/hisat2_build/hisat-3n ${PIPELINE_HOME}/hisat2-hisat-3n/hisat-3n
COPY --from=builder /build/hisat2_build/hisat2-align-s ${PIPELINE_HOME}/hisat2-hisat-3n/hisat2-align-s
COPY --from=builder /build/hisat2_build/hisat2-align-l ${PIPELINE_HOME}/hisat2-hisat-3n/hisat2-align-l
COPY --from=builder /build/hisat2_build/hisat-3n-build ${PIPELINE_HOME}/hisat2-hisat-3n/hisat-3n-build
COPY --from=builder /build/hisat2_build/hisat2-build-s ${PIPELINE_HOME}/hisat2-hisat-3n/hisat2-build-s
COPY --from=builder /build/hisat2_build/hisat2-build-l ${PIPELINE_HOME}/hisat2-hisat-3n/hisat2-build-l
COPY --from=builder /build/hisat2_build/hisat-3n-table ${PIPELINE_HOME}/hisat2-hisat-3n/hisat-3n-table
COPY --from=builder /build/umicollapse_build/ ${PIPELINE_HOME}/UMICollapse/
COPY --from=builder /build/picard_build/picard.jar ${PIPELINE_HOME}/picard/picard.jar

# Set library path for jemalloc
ENV LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so.2

WORKDIR /workspace
