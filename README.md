# m<sup>6</sup>A-CAMseq

![diagram](./docs/diagram.svg)

## Qucik start

- Prepare configuration file

```yaml
reference:
  genes:
    - ../test/spike-in.fa
    - ../test/ERCC92.fa
    - ~/reference/rRNA/Arabidopsis_thaliana.rRNA.fa
  genome:
    - /data/reference/genome/Arabidopsis_thaliana/TAIR10.fa

# TODO: automatically generate index (?)
genome_index: /data/reference/genome/Arabidopsis_thaliana/hisat2_tx_3n/TAIR10.release57

samples:
  test1:
    - R1: ../test/test_R1.fq.gz
      R2: ../test/test_R2.fq.gz
```

- Install apptainer and run

```bash
apptainer run -B /data docker://y9ch/camseq -c config.yaml -j 48
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
