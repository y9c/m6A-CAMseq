# m6A-CAMseq final image
# Builds on top of the trichromat base image

# Pin the base image to the immutable commit-specific tag (NOT the mutable
# `latest` tag). The trichromat release pipeline pushes both `latest` and a
# commit-SHA tag. Pinning to the commit SHA makes builds deterministic and
# prevents a race where a concurrent base-image rebuild makes `latest` resolve
# to a stale image (which caused the old fork-based hisat-3n binaries to leak
# into the SIF).
ARG BASE_IMAGE_TAG=a1ca1748c6822508529e602fcf3368391c40c2c4
FROM ghcr.io/y9c/trichromat:${BASE_IMAGE_TAG}

# The base image already contains:
# - A /pipeline directory with samtools, hisat-3n, etc.
# - A /opt/app_venv with snakemake, cutseq, etc.
# - A /workspace working directory

# Copy CAMseq-specific files over
COPY ./VERSION ./Snakefile ./default.yaml ./entrypoint ${PIPELINE_HOME}/
# The base image does not contain the trichromat source, so we copy it.
# This allows CAMseq to use a specific version of trichromat.
COPY ./external/trichromat/ ${PIPELINE_HOME}/external/trichromat/
COPY ./external/trichromat/bin/ ${PIPELINE_HOME}/script/

# Ensure entrypoint is executable
RUN chmod +x ${PIPELINE_HOME}/entrypoint

# The entrypoint from the base image will be used: ENTRYPOINT ["/pipeline/entrypoint"]
# The working directory is already set to /workspace
