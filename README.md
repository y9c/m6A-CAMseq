# m<sup>6</sup>A-CAMseq

![diagram](./docs/diagram.svg)

## Qucik start

### step 1. Prepare configuration file

minimum configuration example:

```yaml
workdir: workspace

# can be "ECLIP10", "ECLIP6", "TAKARAV3", "SACSEQV3", "STRANDED"
libtype: ECLIP10

# TODO: if the genome index is not available, the pipeline will build it automatically
genome_index: /data/reference/genome/Arabidopsis_thaliana/hisat2_tx_3n/TAIR10.release57

reference:
  contamination:
    # Multiple contaminant reference files are allowed
    - ref/Agrobacterium.fa.gz
  genes:
    # Multiple gene reference files are allowed
    - ref/spikein.fa
    - ref/ERCC92.fa
    - ref/cress_rRNA.fa
  genome:
    # Only one genome is allowed
    - /data/reference/genome/Arabidopsis_thaliana/TAIR10.fa

# note that different from v1, a 3rd level of "data:" tag is added
# this is a preserved behavior for future expansion
samples:
  cress1:
    data:
      - R1: data/cress1_R1.fq.gz
        R2: data/cress1_R2.fq.gz
  cress2:
    data:
      - R1: data/cress2_R1.fq.gz
        R2: data/cress2_R2.fq.gz
```

advanced configuration: please refer to [docs/configuration.md](docs/configuration.md)

### step 2: Install apptainer and run

```bash
apptainer run docker://y9ch/camseq -c data.yaml -j 48
```

NOTES:

- apptainer cannot detect mounted partitions by default, you need to add `-B /partition_name` flag right after the `run` command.

- If you work on a remote server with internet access, you can compile the pipeline first

  ```bash
  apptainer build camseq.sif docker://y9ch/camseq:v2

  # then run the pipeline using the sif file
  ptainer run camseq.sif -c data.yaml -j 48
  ```

- For using the previous version of pipeline (v1), try:

  ```bash
  apptainer run docker://y9ch/camseq:v1 -c data.yaml -j 48
  ```

## Customization

System Requirements

This package has been tested on Linux operating systems. It requires the following software dependencies:

- [Python](https://www.python.org/downloads/) 3.7 or higher
- [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) 8.0.0 or higher
- [hisat2-3n](https://github.com/DaehwanKimLab/hisat2/tree/hisat-3n)
- cutseq

## Documentation

The documentation is available at [https://y9c.github.io/m6A-CAMseq/](https://y9c.github.io/m6A-CAMseq/)

## Changelog

- update in May 11, 2025
  - data processing steps is now based on [trichromat](https://github.com/y9c/trichromat)
  - m6A site detection and filtering is now automatically and have been combined with the pipeline

&nbsp;

<p align="center">
<img
  src="https://raw.githubusercontent.com/y9c/y9c/master/resource/footer_line.svg?sanitize=true"
/>
</p>
<p align="center">
Copyright &copy; 2023-present
<a href="https://github.com/y9c" target="_blank">Chang Y</a>
</p>
<p align="center">
