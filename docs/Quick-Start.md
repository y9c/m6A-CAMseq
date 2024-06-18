---
title: Quick Start
nav_exclude: false
nav_order: 2
---

<!-- prettier-ignore-start -->
# Quick Start
{: .fs-9 }
<!-- prettier-ignore-end -->

1. **Specific the path of references (_.fasta_) and samples (_.fastq_) in a configure file (_.YAML_).**

   For example, write down and save the following block into a text file and named it as `data.yaml`.

   ```yaml
   reference:
     genes:
       fa: ./ref/genes.fa
     genome:
       fa: /data/reference/genome/Mus_musculus/GRCm39.fa
       star: /data/reference/genome/Mus_musculus/star/GRCm39.release108

   samples:
     mESCWT-rep1-input:
       data:
         - R1: ./test/IP16.fastq.gz
       group: mESCWT
       treated: false
     mESCWT-rep1-treated:
       data:
         - R1: ./test/IP4.fastq.gz
       group: mESCWT
       treated: true
     mESCWT-rep2-treated:
       data:
         - R1: ./test/IP5.fastq.gz
       group: mESCWT
       treated: true
   ```

2. **Run all the analysis by one command**:

   ```bash
   snakemake -s Snakefile -c config.yaml
   ```
