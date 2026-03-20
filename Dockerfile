# m6A-CAMseq final image
# Builds on top of the trichromat base image

# You can specify the base image tag to use
ARG BASE_IMAGE_TAG=latest
FROM y9c/trichromat:${BASE_IMAGE_TAG}

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
