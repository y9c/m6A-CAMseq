from snakemake.utils import min_version
from collections import defaultdict


min_version("8.0")


BATCH = config["batch"]


workdir: f"workspace_{BATCH}"


config.update(config.get(f"meta_{BATCH}", {}))

REF = config.get("ref_" + config.get("ref"), {})


TEMPDIR = Path(
    os.path.relpath(
        config.get("tempdir", os.path.join(workflow.basedir, ".tmp")), workflow.basedir
    )
)
PATH = config["path"]
INTERNALDIR = Path("internal_files")
MARKDUP = config.get("markdup", True)
SPLICE_GENOME = config.get("splice_genome", True)
WITH_UMI = config.get("with_umi", False)

SAMPLE2DATA = defaultdict(lambda: defaultdict(dict))
GROUP2SAMPLE = defaultdict(lambda: defaultdict(list))
SAMPLE2LIB = defaultdict(dict)
for s, v in config[f"samples_{BATCH}"].items():
    SAMPLE2LIB[s] = v.get("libtype", config.get("libtype", ""))
    if "group" in v and "lib" in v:
        GROUP2SAMPLE[v["group"]][v["lib"]].append(s)
    for i, v2 in enumerate(v["data"], 1):
        r = f"run{i}"
        SAMPLE2DATA[str(s)][r] = {
            k: os.path.expanduser(v3) for k, v3 in dict(v2).items()
        }


rule all:
    input:
        "report_reads/report_cutadapt.html",
        "report_reads/report_trimmed.html",
        expand("report_reads/count/{sample}.tsv", sample=SAMPLE2DATA.keys()),
        expand(
            INTERNALDIR / "pileup_sites/{sample}.{reftype}.tsv.gz",
            sample=SAMPLE2DATA.keys(),
            reftype=["genes", "genome"],
        ),
        expand(
            "report_reads/isize/{sample}.{reftype}.tsv",
            sample=SAMPLE2DATA.keys(),
            reftype=["genes", "genome"],
        ),
        "report_reads/ratio/probe.tsv",
        expand(
            INTERNALDIR / "genes_stat/{sample}.{reftype}.tsv",
            sample=SAMPLE2DATA.keys(),
            reftype=["genes", "genome"],
        ),
        expand(
            INTERNALDIR / "motif_stat/{sample}.{reftype}.tsv",
            sample=SAMPLE2DATA.keys(),
            reftype=["genes", "genome"],
        ),
        "merged_sites/merged.tsv.gz",


ruleorder: cutadapt_PE > cutadapt_SE


rule cutadapt_SE:
    input:
        lambda wildcards: SAMPLE2DATA[wildcards.sample][wildcards.rn].get("R1", "/"),
    output:
        c=temp(TEMPDIR / "trimmed_reads/SE/{sample}_{rn}_R1.fq.gz"),
        s=INTERNALDIR / "discarded_reads/{sample}_{rn}_tooshort_R1.fq.gz",
        report="report_reads/trimming/{sample}_{rn}.json",
    params:
        library=lambda wildcards: SAMPLE2LIB[wildcards.sample],
    threads: 16
    shell:
        """
        cutseq -t {threads} -A {params.library:q} -m 20 --auto-rc -o {output.c} -s {output.s} --json-file {output.report} {input} 
        """


rule cutadapt_PE:
    input:
        lambda wildcards: SAMPLE2DATA[wildcards.sample][wildcards.rn].get("R1", "/"),
        lambda wildcards: SAMPLE2DATA[wildcards.sample][wildcards.rn].get("R2", "/"),
    output:
        c=[
            temp(TEMPDIR / "trimmed_reads/PE/{sample}_{rn}_R1.fq.gz"),
            temp(TEMPDIR / "trimmed_reads/PE/{sample}_{rn}_R2.fq.gz"),
        ],
        s=[
            INTERNALDIR / "discarded_reads/{sample}_{rn}_tooshort_R1.fq.gz",
            INTERNALDIR / "discarded_reads/{sample}_{rn}_tooshort_R2.fq.gz",
        ],
        report="report_reads/trimming/{sample}_{rn}.json",
    params:
        library=lambda wildcards: SAMPLE2LIB[wildcards.sample],
    threads: 16
    shell:
        """
        cutseq -t {threads} -A {params.library:q} -m 20 --auto-rc -o {output.c} -s {output.s} --json-file {output.report} {input} 
        """


rule report_cutadapt:
    input:
        [
            f"report_reads/trimming/{sample}_{rn}.json"
            for sample, d1 in SAMPLE2DATA.items()
            for rn, d2 in d1.items()
        ],
    output:
        "report_reads/report_cutadapt.html",
    shell:
        "multiqc -f -m cutadapt --no-data-dir -n {output} {input}"


rule qc_trimmed:
    input:
        lambda wildcards: (
            TEMPDIR / "trimmed_reads/PE/{sample}_{rn}_{rd}.fq.gz"
            if len(SAMPLE2DATA[wildcards.sample][wildcards.rn]) == 2
            else TEMPDIR / "trimmed_reads/SE/{sample}_{rn}_{rd}.fq.gz"
        ),
    output:
        html="report_reads/trimmed/{sample}_{rn}_{rd}/fastqc_report.html",
        text="report_reads/trimmed/{sample}_{rn}_{rd}/fastqc_data.txt",
        summary="report_reads/trimmed/{sample}_{rn}_{rd}/summary.txt",
    params:
        "report_reads/trimmed/{sample}_{rn}_{rd}",
    shell:
        "falco -o {params} {input}"


rule report_qc_trimmed:
    input:
        [
            f"report_reads/trimmed/{s}_{r}_{i}/fastqc_data.txt"
            for s, v in SAMPLE2DATA.items()
            for r, v2 in v.items()
            for i in v2.keys()
        ],
    output:
        "report_reads/report_trimmed.html",
    shell:
        "multiqc -f -m fastqc -n {output} {input}"


rule prepare_genes_index:
    input:
        GENES_FASTA,
    output:
        fa="genes_index/genes.fa",
        index="genes_index/genes.3n.CT.1.ht2",
    params:
        prefix="genes_index/genes",
    threads: 4
    shell:
        """
        cat {input} | sed 's/(//g' | sed 's/)//g' > {output.fa}
        {PATH[samtools]} faidx {output.fa}
        rm -f genes_index/genes.3n.*.ht2
        {PATH[hisat3nbuild]} -p {threads} --base-change {BASE_CHANGE} {output.fa} {params.prefix}
        """


rule hisat2_3n_mapping_genes_PE:
    input:
        fq1=TEMPDIR / "trimmed_reads_PE/{sample}_{rn}_R1.fq.gz",
        fq2=TEMPDIR / "trimmed_reads_PE/{sample}_{rn}_R2.fq.gz",
        index="genes_index/genes.3n.CT.1.ht2",
    output:
        bam=temp(TEMPDIR / "run_mapping_PE/{sample}_{rn}.genes.bam"),
        unmap=temp(TEMPDIR / "run_mapping_PE/{sample}_{rn}.genes.unmap.bam"),
        summary="report_reads/mapping/{sample}_{rn}.genes.summary",
    params:
        index="genes_index/genes",
        basechange=BASE_CHANGE,
        directional=(
            "--directional-mapping-reverse"
            if STRANDNESS == "R"
            else "--directional-mapping" if STRANDNESS == "F" else ""
        ),
        flt="flag.proper_pair && !flag.unmap && !flag.munmap"
        + (
            " && !flag.read1 == !flag.reverse"
            if (GENE_NORC and STRANDNESS == "R")
            else (
                " && !flag.read1 != !flag.reverse"
                if (GENE_NORC and STRANDNESS == "F")
                else ""
            )
        ),
    threads: 24
    shell:
        """
        {PATH[hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -1 {input.fq1} -2 {input.fq2} --base-change {params.basechange} \
            {params.directional} --avoid-pseudogene --no-softclip --no-spliced-alignment --bowtie2-dp 2 --np 0 --rdg 5,3 --rfg 5,3 --mp 3,1 --score-min L,-3,-0.6 |\
            {PATH[samtools]} view -@ {threads} -O BAM -e {params.flt:q} -U {output.unmap} -o {output.bam}
        """


rule hisat2_3n_mapping_genome_PE:
    input:
        fq1=TEMPDIR / "genes_unmapped/PE/{sample}_{rn}_unmap_R1.fq.gz",
        fq2=TEMPDIR / "genes_unmapped/PE/{sample}_{rn}_unmap_R2.fq.gz",
    output:
        bam=temp(TEMPDIR / "rerun_mapping/PE/{sample}_{rn}.genome.bam"),
        summary="report_reads/remap/{sample}_{rn}_PE.summary",
    params:
        index=REF["genome"]["hisat3n"],
        basechange=config.get("base_change", "A,G"),
        directional="--directional-mapping",
        splice_args=(
            "--pen-noncansplice 20 --min-intronlen 20 --max-intronlen 20"
            if SPLICE_GENOME
            else "--no-spliced-alignment"
        ),
    threads: 24
    shell:
        # --no-discordant --no-mixed
        """
        export TMPDIR="/scratch/midway3/yec"
        {PATH[hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -1 {input.fq1} -2 {input.fq2} --base-change {params.basechange} {params.directional} {params.splice_args} \
            --avoid-pseudogene --bowtie2-dp 1 --np 0 --rdg 5,3 --rfg 5,3 --sp 9,3 --mp 3,1 --score-min L,-3,-0.5 |\
            {PATH[samtools]} view -@ {threads} -O BAM -o {output.bam}
        """


rule hisat2_3n_mapping_genome_SE:
    input:
        fq=TEMPDIR / "genes_unmapped/SE/{sample}_{rn}_unmap_R1.fq.gz",
    output:
        bam=temp(TEMPDIR / "rerun_mapping/SE/{sample}_{rn}.genome.bam"),
        summary="report_reads/remap/{sample}_{rn}_SE.summary",
    params:
        index=REF["genome"]["hisat3n"],
        basechange=config.get("base_change", "A,G"),
        directional="--directional-mapping",
        splice_args=(
            "--pen-noncansplice 20 --min-intronlen 20 --max-intronlen 20"
            if SPLICE_GENOME
            else "--no-spliced-alignment"
        ),
    threads: 24
    shell:
        """
        {PATH[hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -U {input.fq} --base-change {params.basechange} {params.directional} {params.splice_args} \
            --avoid-pseudogene --bowtie2-dp 1 --np 0 --rdg 5,3 --rfg 5,3 --sp 9,3 --mp 3,1 --score-min L,-3,-0.5 |\
            {PATH[samtools]} view -@ {threads} -o {output.bam}
        """


rule add_tag_remap:
    input:
        lambda wildcards: (
            TEMPDIR / "rerun_mapping/PE/{sample}_{rn}.genome.bam"
            if len(SAMPLE2DATA[wildcards.sample][wildcards.rn]) == 2
            else TEMPDIR / "rerun_mapping/SE/{sample}_{rn}.genome.bam"
        ),
    output:
        temp(TEMPDIR / "rerun_bam/{sample}_{rn}.bam"),
    params:
        is_reverse="--forward-lib",
    threads: 20
    shell:
        "python ../src/add_tag.py {input} {output} --threads {threads} {params.is_reverse}"


rule filter_and_sort_bam_remap:
    input:
        TEMPDIR / "rerun_bam/{sample}_{rn}.bam",
    output:
        unmap=temp(TEMPDIR / "rerun_unmap/{sample}_{rn}.bam"),
        mapped=INTERNALDIR / "genome_bam/{sample}_{rn}.bam",
    threads: 16
    shell:
        """
        samtools view -e 'exists([AP]) && [AP]<=5 && !flag.secondary' -@ {threads} -U {output.unmap} -h {input} |\
            samtools sort -m 3G -O BAM -o {output.mapped}
        """


ruleorder: extract_unmapped_reads_PE > extract_unmapped_reads_SE


rule extract_unmapped_reads_PE:
    input:
        un=TEMPDIR / "rerun_unmap/{sample}_{rn}.bam",
    output:
        r1=INTERNALDIR / "discarded_reads/{sample}_{rn}_unmap_R1.fq.gz",
        r2=INTERNALDIR / "discarded_reads/{sample}_{rn}_unmap_R2.fq.gz",
    shell:
        """
        {PATH[samtools]} fastq -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n {input}
        """


rule extract_unmapped_reads_SE:
    input:
        un=TEMPDIR / "rerun_unmap/{sample}_{rn}.bam",
    output:
        r1=INTERNALDIR / "discarded_reads/{sample}_{rn}_unmap_R1.fq.gz",
    shell:
        """
        {PATH[samtools]} fastq -0 {output.r1} -n {input}
        """


rule combine_runs:
    input:
        lambda wildcards: [
            (
                INTERNALDIR / f"genes_bam/{wildcards.sample}_{r}.bam"
                if wildcards.reftype == "genes"
                else INTERNALDIR / f"genome_bam/{wildcards.sample}_{r}.bam"
            )
            for r in SAMPLE2DATA[wildcards.sample]
        ],
    output:
        bam=temp(TEMPDIR / "combined_runs/{sample}.{reftype}.bam"),
        bai=temp(TEMPDIR / "combined_runs/{sample}.{reftype}.bam.bai"),
    threads: 16
    shell:
        """
        samtools merge -@ {threads} -f --write-index -o {output.bam}##idx##{output.bai} {input}
        """


rule stat_combined:
    input:
        bam=TEMPDIR / "combined_runs/{sample}.{reftype}.bam",
    output:
        stat="report_reads/combined/{sample}.{reftype}.txt",
        n="report_reads/combined/{sample}.{reftype}.count",
    threads: 2
    shell:
        """
        samtools flagstat -@ {threads} -O TSV {input} > {output.stat}
        samtools view -@ {threads} -c -F 384 {input} > {output.n}
        """


rule drop_duplicates:
    input:
        bam=TEMPDIR / "combined_runs/{sample}.{reftype}.bam",
    output:
        bam=INTERNALDIR / "aligned_bam/{sample}.{reftype}.bam",
        txt="report_reads/dedup/{sample}.{reftype}.log",
    params:
        path_umicollapse=config["path"]["umicollapse"],
    threads: 16
    run:
        if SAMPLE2LIB[wildcards.sample] in ["ECLIP10", "ECLIP6", "TAKARAV3"]:
            shell(
                """
            export TMPDIR="/scratch/midway3/yec"
            /software/java-15.0.2-el8-x86_64/bin/java -server -Xms8G -Xmx40G -Xss100M -Djava.io.tmpdir=$TMPDIR -jar {params.path_umicollapse} bam \
                -t 2 -T {threads} --data naive --merge avgqual --two-pass -i {input.bam} -o {output.bam} >{output.txt}
            """
            )
        elif MARKDUP:
            shell(
                """
            export TMPDIR="/scratch/midway3/yec"
                ~/tools/jdk8u322-b06-jre/bin/java -Xmx36G -jar ~/tools/gatk-4.2.5.0/gatk-package-4.2.5.0-local.jar MarkDuplicates \
                    -I {input} -O {output.bam} -M {output.txt} \
                    --DUPLICATE_SCORING_STRATEGY SUM_OF_BASE_QUALITIES --REMOVE_DUPLICATES true --VALIDATION_STRINGENCY SILENT --TMP_DIR $TMPDIR
            """
            )
        else:
            shell(
                """
                cp {input.bam} {output.bam}
                touch {output.txt}
            """
            )


rule dedup_index:
    input:
        bam=INTERNALDIR / "aligned_bam/{sample}.{reftype}.bam",
    output:
        bai=INTERNALDIR / "aligned_bam/{sample}.{reftype}.bam.bai",
    threads: 6
    shell:
        """
        samtools index -@ {threads} {input}
        """


rule stat_dedup:
    input:
        bam=INTERNALDIR / "aligned_bam/{sample}.{reftype}.bam",
    output:
        stat="report_reads/dedup/{sample}.{reftype}.txt",
        n="report_reads/dedup/{sample}.{reftype}.count",
    threads: 2
    shell:
        """
        samtools flagstat -@ {threads} -O TSV {input} > {output.stat}
        samtools view -@ {threads} -c -F 384 {input} > {output.n}
        """


rule pileup_base:
    input:
        bam=INTERNALDIR / "aligned_bam/{sample}.{reftype}.bam",
        bai=INTERNALDIR / "aligned_bam/{sample}.{reftype}.bam.bai",
        ref=lambda wildcards: (
            INTERNALDIR / "reference_file/genes.fa"
            if wildcards.reftype == "genes"
            else os.path.expanduser(REF["genome"]["fa"])
        ),
    output:
        INTERNALDIR / "pileup_sites/{sample}.{reftype}.tsv.gz",
    threads: 8
    shell:
        "python ../src/count_pileup.py {input.bam} {input.ref} | bgzip -@ {threads} -c > {output}"
