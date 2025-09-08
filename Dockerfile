# x86_64
# Use ARGs for versions for easier management and clarity
ARG SAMTOOLS_VERSION="1.21"
ARG PICARD_VERSION="3.4.0"
ARG PYTHON_VERSION_FOR_APP="3.13" # Python version for the app and base image
# Construct the correct tag for the uv base image
ARG UV_BASE_IMAGE_TAG="python${PYTHON_VERSION_FOR_APP}-bookworm-slim"

# ----------- Builder Stage -----------
FROM ghcr.io/astral-sh/uv:${UV_BASE_IMAGE_TAG} AS builder

ARG SAMTOOLS_VERSION
ARG PICARD_VERSION
ARG PYTHON_VERSION_FOR_APP # Make it available in this stage too
# uv and Python ${PYTHON_VERSION_FOR_APP} are pre-installed in this base image.

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

# Install build-time dependencies for system tools.
RUN apt-get update && \
    apt-get -y --no-install-recommends install \
    ca-certificates tzdata apt-utils wget curl bzip2 unzip make gcc g++ pkg-config xsltproc \
    zlib1g-dev libxml2-dev \
    gfortran libopenblas-dev liblapack-dev \
    default-jre coreutils procps && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# --- Create a virtual environment for Python libraries ---
ENV APP_VENV_PATH=/opt/app_venv
# Use the Python from the base image
RUN python${PYTHON_VERSION_FOR_APP} -m venv ${APP_VENV_PATH}

# Install Python packages (libraries and tools) into this venv using the global uv
ENV PYTHON_PACKAGES="scipy polars cutseq snakemake==9.9.0"
RUN echo "Installing Python packages: ${PYTHON_PACKAGES} into ${APP_VENV_PATH}" && \
    uv pip install --python ${APP_VENV_PATH}/bin/python --no-cache ${PYTHON_PACKAGES}

# --- Build samtools (binary will be copied directly) ---
WORKDIR /build/samtools_src
RUN wget -qO- https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 | tar -xjvf - --strip-components 1 && \
    ./configure --without-curses --disable-bz2 --disable-lzma && \
    make -j$(nproc) samtools

# --- Build hisat2 (from hisat-3n branch) ---
WORKDIR /build/hisat2_build
RUN wget -qO- https://github.com/DaehwanKimLab/hisat2/archive/refs/heads/hisat-3n.tar.gz | tar -xzf - --strip-components=1 && \
    make -j$(nproc) # This builds hisat2, hisat2-build, hisat2-build-s, hisat2-build-l, etc.

# --- Prepare UMICollapse ---
WORKDIR /build/umicollapse_build
RUN wget -q -P ./ https://github.com/Daniel-Liu-c0deb0t/UMICollapse/raw/refs/heads/master/umicollapse.jar && \
    mkdir lib && \
    wget -q -P ./lib https://repo1.maven.org/maven2/com/github/samtools/htsjdk/2.19.0/htsjdk-2.19.0.jar && \
    wget -q -P ./lib https://repo1.maven.org/maven2/org/xerial/snappy/snappy-java/1.1.7.3/snappy-java-1.1.7.3.jar

# --- Prepare picard ---
WORKDIR /build/picard_build
RUN wget -qO picard.jar https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar && \
    echo "downloaded picard JAR..."

# -- Download and prepare pbr --
WORKDIR /build/pbr_build
RUN wget https://github.com/y9c/pbr/releases/download/latest-prerelease/pbr-x86_64-unknown-linux-musl.tar.gz && \
    tar zxvf pbr-x86_64-unknown-linux-musl.tar.gz && \
    chmod +x pbr && \
    rm pbr-x86_64-unknown-linux-musl.tar.gz

# ----------- Final Stage -----------
FROM ghcr.io/astral-sh/uv:${UV_BASE_IMAGE_TAG} AS final
# This image already has Python ${PYTHON_VERSION_FOR_APP} and uv.

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC
ENV PIPELINE_HOME=/pipeline
# APP_VENV_PATH is to be copied from builder
ENV APP_VENV_PATH=/opt/app_venv

# PATH order: venv python, then other tool locations, then standard PATH from base image
ENV PATH="${APP_VENV_PATH}/bin:${PIPELINE_HOME}/bin:${PIPELINE_HOME}/samtools:${PIPELINE_HOME}/hisat2-hisat-3n:/usr/local/bin:$PATH"

ENV UMICollapse_JAR_PATH="${PIPELINE_HOME}/umicollapse/umicollapse.jar"
ENV PICARD_JAR_PATH="${PIPELINE_HOME}/picard/picard.jar"

# Install essential runtime dependencies. Python is from the base image.
RUN apt-get update && \
    apt-get -y --no-install-recommends install \
    ca-certificates tzdata default-jre zlib1g libxml2 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Create directory structure
RUN mkdir -p ${PIPELINE_HOME}/bin \
             ${PIPELINE_HOME}/external/trichromat \
             ${PIPELINE_HOME}/umicollapse \
             ${PIPELINE_HOME}/picard \
             ${PIPELINE_HOME}/samtools \
             ${PIPELINE_HOME}/hisat2-hisat-3n \
             ${APP_VENV_PATH} \
             /usr/local/bin \
             /workspace

# Copy the application's Python virtual environment (with scipy, polars, cutseq, snakemake)
COPY --from=builder ${APP_VENV_PATH} ${APP_VENV_PATH}

# Copy samtools binary
COPY --from=builder /build/samtools_src/samtools ${PIPELINE_HOME}/samtools/samtools

# Copy ALL relevant HISAT2 executables (wrappers and compiled binaries)
COPY --from=builder /build/hisat2_build/hisat-3n ${PIPELINE_HOME}/hisat2-hisat-3n/hisat-3n
COPY --from=builder /build/hisat2_build/hisat2-align-s ${PIPELINE_HOME}/hisat2-hisat-3n/hisat2-align-s
COPY --from=builder /build/hisat2_build/hisat2-align-l ${PIPELINE_HOME}/hisat2-hisat-3n/hisat2-align-l
COPY --from=builder /build/hisat2_build/hisat-3n-build ${PIPELINE_HOME}/hisat2-hisat-3n/hisat-3n-build
COPY --from=builder /build/hisat2_build/hisat2-build-s ${PIPELINE_HOME}/hisat2-hisat-3n/hisat2-build-s
COPY --from=builder /build/hisat2_build/hisat2-build-l ${PIPELINE_HOME}/hisat2-hisat-3n/hisat2-build-l
COPY --from=builder /build/hisat2_build/hisat-3n-table ${PIPELINE_HOME}/hisat2-hisat-3n/hisat-3n-table

# Copy UMICollapse and picard
COPY --from=builder /build/umicollapse_build/ ${PIPELINE_HOME}/UMICollapse/
COPY --from=builder /build/picard_build/picard.jar ${PIPELINE_HOME}/picard/picard.jar

# Copy pbr (downloaded binary)
COPY --from=builder /build/pbr_build/pbr ${PIPELINE_HOME}/pbr/pbr

# Copy application-specific files
# COPY ./bin ${PIPELINE_HOME}/bin/
COPY ./VERSION ./Snakefile ./default.yaml ./entrypoint ${PIPELINE_HOME}/
COPY ./external/trichromat/workflow/ ${PIPELINE_HOME}/external/trichromat/workflow/
COPY ./external/trichromat/Snakefile ./external/trichromat/workflow_utils.py ./external/trichromat/default.yaml ./external/trichromat/config.schema.yaml ${PIPELINE_HOME}/external/trichromat/
COPY ./external/trichromat/bin/ ${PIPELINE_HOME}/script/

WORKDIR /workspace

# Executables for cutseq and snakemake will be in ${APP_VENV_PATH}/bin, which is on PATH
RUN chmod +x ${PIPELINE_HOME}/entrypoint && \
    find ${PIPELINE_HOME}/bin/ -type f -exec chmod +x {} \; && \
    chmod +x ${PIPELINE_HOME}/samtools/samtools && \
    chmod +x ${PIPELINE_HOME}/pbr/pbr && \
    find ${PIPELINE_HOME}/hisat2-hisat-3n/ -maxdepth 1 -type f -exec chmod +x {} \; -o -type l -exec chmod +x {} \;

ENTRYPOINT ["/pipeline/entrypoint"]
