# calling
TEMPDIR = Path(
    os.path.relpath(
        config.get("tempdir", os.path.join(workflow.basedir, ".tmp")), workflow.basedir
    )
)
INTERNALDIR = Path("internal_files")

PATH = config["path"]

LIBTYPE = config.get("libtype", "")

MARKDUP = config.get("markdup", True)
STRANDNESS = config.get("strandness", True)
GENE_NORC = config.get("gene_norc", True)
BASE_CHANGE = config.get("base_change", "A,G")
SPLICE_GENOME = config.get("splice_genome", True)
SPLICE_CONTAM = config.get("splice_contamination", False)
REF = config.get("_REF", {})
READS = config.get("_READS", {})
SAMPLE_LIB = config.get("_LIB", {})
SAMPLE_ADP = config.get("_ADP", {})
SAMPLE_UMI = config.get("_UMI", {})
SAMPLE_STD = config.get("_STD", {})
SAMPLE_DUP = config.get("_DUP", {})
SAMPLE_NORC = config.get("_NORC", {})

MOTIF_FLANKING = config.get("motif_flanking", 2)
DROP_CLUSTERED = config.get("drop_clustered", True)


rule build_internal_fa_index:
    input:
        INTERNALDIR / "reference_file/{reftype}.fa",
    output:
        INTERNALDIR / "reference_file/{reftype}.fa.fai",
    threads: 8
    shell:
        """
        {PATH[samtools]} faidx -@ {threads} {input} --fai-idx {output}
        """


rule build_genome_fa_index:
    input:
        REF["genome"][0],
    output:
        REF["genome"][0] + ".fai",
    threads: 8
    shell:
        """
        {PATH[samtools]} faidx -@ {threads} {input} --fai-idx {output}
        """


# hisat:
# ref     pos     strand  convertedBaseCount      unconvertedBaseCount
# pbr:
# #chrom  pos0    ref_base        depth   a       c       g       t       n
# fwd: $1, $2 + 1, +, $7, $5
# rev: $1, $2 + 1, -, $6, $8

# {PATH[samtools]} view -@ {threads} -e "rlen < 100000 && [XM] * 20 <= (qlen-sclen) && [Zf] <= 3 && 3 * [Zf] <= [Zf] + [Yf]" -h {input} | \
#     {PATH[hisat3ntable]} -p {threads} --alignments - --ref {params.fa} --output-name /dev/stdout --base-change {params.basechange} | \
#     cut -f 1,2,3,5,7 | \
#     gzip > {output}


rule pileup_by_fwd_strand:
    input:
        bam=INTERNALDIR / "aligned_bam/{sample}.{reftype}.bam",
        bai=INTERNALDIR / "aligned_bam/{sample}.{reftype}.bam.bai",
        fa=lambda wildcards: (
            REF["genome"][0]
            if wildcards.reftype == "genome"
            else INTERNALDIR / f"reference_file/{wildcards.reftype}.fa"
        ),
        fai=lambda wildcards: (
            REF["genome"][0] + ".fai"
            if wildcards.reftype == "genome"
            else INTERNALDIR / f"reference_file/{wildcards.reftype}.fa.fai"
        ),
    output:
        temp(TEMPDIR / "pileup_by_strand/{sample}.{reftype}.fwd.tsv"),
    params:
        basechange=BASE_CHANGE,
        motif_flanking=MOTIF_FLANKING,
        drop_clustered=(
            "&& [Zf] <= 3 && 2 * [Zf] <= [Yf]"
            if DROP_CLUSTERED
            else ""
        ),
    threads: 16
    shell:
        """
        {PATH[samtools]} view -@ {threads} -e "rlen < 100000 && [XM] * 20 <= (qlen-sclen) {params.drop_clustered} && mapq >= 20" -h {input.bam} | \
            {PATH[hisat3ntable]} -p {threads} --alignments - --ref {input.fa} --output-name /dev/stdout --base-change {params.basechange} | \
            awk '$3 == "+"' | cut -f 1,2,3,5,7 | \
            {PATH[annotateMotif]} -r {input.fa} -k {params.motif_flanking} > {output}
        """


rule pileup_by_rev_strand:
    input:
        bam=INTERNALDIR / "aligned_bam/{sample}.{reftype}.bam",
        bai=INTERNALDIR / "aligned_bam/{sample}.{reftype}.bam.bai",
        fa=lambda wildcards: (
            INTERNALDIR / "reference_file/genes.fa"
            if wildcards.reftype == "genes"
            else REF["genome"][0]
        ),
        fai=lambda wildcards: (
            INTERNALDIR / "reference_file/genes.fa.fai"
            if wildcards.reftype == "genes"
            else REF["genome"][0] + ".fai"
        ),
    output:
        temp(TEMPDIR / "pileup_by_strand/{sample}.{reftype}.rev.tsv"),
    params:
        basechange=BASE_CHANGE,
        motif_flanking=MOTIF_FLANKING,
        drop_clustered=(
            "&& [Zf] <= 3 && 2 * [Zf] <= [Yf]"
            if DROP_CLUSTERED
            else ""
        ),
    threads: 16
    shell:
        """
        {PATH[samtools]} view -@ {threads} -e "rlen < 100000 && [XM] * 20 <= (qlen-sclen) {params.drop_clustered} && mapq >= 20" -h {input.bam} | \
            {PATH[hisat3ntable]} -p {threads} --alignments - --ref {input.fa} --output-name /dev/stdout --base-change {params.basechange} | \
            awk '$3 == "-"' | cut -f 1,2,3,5,7 | \
            {PATH[annotateMotif]} -r {input.fa} -k {params.motif_flanking} > {output}
        """


rule pileup_to_arrow:
    input:
        TEMPDIR / "pileup_by_strand/{sample}.{reftype}.{strand}.tsv",
    output:
        INTERNALDIR / "pileup_arrow/{sample}.{reftype}.{strand}.arrow",
    threads: 2
    shell:
        """
        {PATH[convertPileup]} {input} {output}
        """


rule join_pileup_table:
    input:
        expand(
            INTERNALDIR / "pileup_arrow/{sample}.{{reftype}}.{strand}.arrow",
            sample=READS.keys(),
            strand=["fwd", "rev"],
        ),
    output:
        "report_sites/joined/{reftype}.arrow",
    params:
        names=[name for name in READS.keys() for _ in range(2)],
        strand=["+", "-"] * len(READS.keys()),
    threads: 8
    shell:
        """
        {PATH[joinPileup]} -f {input} -n {params.names} -s {params.strand} -o {output}
        """


# usage: filter_sites.py [-h] -f FILE -o OUTPUT [-u MIN_UNCON] [-d MIN_DEPTH] [-r MIN_RATIO] [-p MIN_PVAL]
rule filter_sites:
    input:
        "report_sites/joined/{reftype}.arrow",
    output:
        "report_sites/filtered/{reftype}.tsv.gz",
    params:
        min_uncon=config.get("cutoff", {}).get("min_uncon", 1),
        min_depth=config.get("cutoff", {}).get("min_depth", 10),
        min_ratio=config.get("cutoff", {}).get("min_ratio", 0.05),
        min_pval=config.get("cutoff", {}).get("min_pval", 1),
    threads: 8
    shell:
        """
        {PATH[filterSites]} -f {input} -o {output} -u {params.min_uncon} -d {params.min_depth} -r {params.min_ratio} -p {params.min_pval}
        """
