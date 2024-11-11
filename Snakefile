from collections import defaultdict


WORKDIR = os.path.relpath(
    config.get("workdir", "workspace"), os.path.dirname(workflow.configfiles[-1])
)
TEMPDIR = Path(
    os.path.relpath(
        config.get("tempdir", os.path.join(workflow.basedir, ".tmp")), workflow.basedir
    )
)
INTERNALDIR = Path("internal_files")


workdir: WORKDIR


PATH = config["path"]


LIBTYPE = config["libtype"]
WITH_UMI = LIBTYPE in ["ECLIP10", "ECLIP6", "TAKARAV3", "SACSEQV3"]
MARKDUP = config.get("markdup", True)
STRANDNESS = config.get("strandness", True)
GENE_NORC = config.get("gene_norc", True)
BASE_CHANGE = config.get("base_change", "A,G")
SPLICE_GENOME = config.get("splice_genome", True)


REF = config["reference"]
for k, v in REF.items():
    v3 = []
    for v2 in v:
        v2 = os.path.expanduser(v2)
        v3.append(v2 if os.path.isabs(v2) else os.path.relpath(v2, WORKDIR))
    REF[k] = v3

GENES_FASTA = REF.get("genes", [])
GENOME_FASTA = REF.get("genome", [])

SAMPLE2DATA = defaultdict(lambda: defaultdict(dict))
for s, v in config[f"samples"].items():
    for i, v2 in enumerate(v, 1):
        r = f"run{i}"
        SAMPLE2DATA[str(s)][r] = {
            k: (
                os.path.expanduser(v3)
                if os.path.isabs(os.path.expanduser(v3))
                else os.path.relpath(os.path.expanduser(v3), WORKDIR)
            )
            for k, v3 in dict(v2).items()
        }


rule all:
    input:
        expand(
            "report_reads/dedup/{sample}.{reftype}.count",
            sample=SAMPLE2DATA.keys(),
            reftype=["genes", "genome"],
        ),
        expand(
            "count_sites/pileup/{sample}.{reftype}.tsv.gz",
            sample=SAMPLE2DATA.keys(),
            reftype=["genes", "genome"],
        ),
        [
            INTERNALDIR / f"discarded_reads/{sample}_{rn}_unmap_{rd}.fq.gz"
            for sample, v in SAMPLE2DATA.items()
            for rn, v2 in v.items()
            for rd in (["R1", "R2"] if len(v2) == 2 else ["R1"])
        ],


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
        library=LIBTYPE,
    threads: 16
    shell:
        """
        {PATH[cutseq]} -t {threads} -A {params.library:q} -m 20 --auto-rc -o {output.c} -s {output.s} --json-file {output.report} {input} 
        """


rule prepare_genes_index:
    input:
        [os.path.expanduser(f) for f in GENES_FASTA],
    output:
        fa=INTERNALDIR / "genes_index/genes.fa",
        index=INTERNALDIR / "genes_index/genes.3n.CT.1.ht2",
    params:
        prefix=INTERNALDIR / "genes_index/genes",
    threads: 4
    shell:
        """
        cat {input} | sed 's/(//g' | sed 's/)//g' > {output.fa}
        {PATH[samtools]} faidx {output.fa}
        rm -f {params.prefix}.3n.*.ht2
        {PATH[hisat3nbuild]} -p {threads} --base-change {BASE_CHANGE} {output.fa} {params.prefix}
        """


rule hisat2_3n_mapping_genes_PE:
    input:
        fq1=TEMPDIR / "trimmed_reads/PE/{sample}_{rn}_R1.fq.gz",
        fq2=TEMPDIR / "trimmed_reads/PE/{sample}_{rn}_R2.fq.gz",
        index=INTERNALDIR / "genes_index/genes.3n.CT.1.ht2",
    output:
        bam=temp(TEMPDIR / "genes_mapping/PE/{sample}_{rn}.bam"),
        summary="report_reads/mapping/{sample}_{rn}.genes.summary",
    params:
        index=INTERNALDIR / "genes_index/genes",
        basechange=BASE_CHANGE,
        directional="--directional-mapping" if STRANDNESS else "",
    threads: 24
    shell:
        """
        {PATH[hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -1 {input.fq1} -2 {input.fq2} --base-change {params.basechange} {params.directional} \
            --avoid-pseudogene --no-softclip --no-spliced-alignment --np 0 --rdg 5,3 --rfg 5,3 --mp 3,1 --score-min L,-3,-0.75 |\
            {PATH[samtools]} view -@ {threads} -O BAM -o {output.bam}
        """


rule filter_and_sort_bam_genes:
    input:
        TEMPDIR / "genes_mapping/PE/{sample}_{rn}.bam",
    output:
        unmap=temp(TEMPDIR / "genes_unmapped/PE/{sample}_{rn}.bam"),
        mapped=temp(TEMPDIR / "genes_bam/{sample}_{rn}.bam"),
    params:
        flt="flag.proper_pair && !flag.unmap && !flag.munmap"
        + (" && !flag.read1 != !flag.reverse" if (GENE_NORC and STRANDNESS) else ""),
    threads: 16
    shell:
        """
        {PATH[samtools]} view -e '{params.flt}' -@ {threads} -U {output.unmap} -h {input} |\
            {PATH[samtools]} sort -m 3G -O BAM -o {output.mapped}
        """


rule extract_unmapped_reads_PE_from_genes:
    input:
        un=TEMPDIR / "genes_unmapped/PE/{sample}_{rn}.bam",
    output:
        r1=temp(TEMPDIR / "genes_unmapped/PE/{sample}_{rn}_R1.fq.gz"),
        r2=temp(TEMPDIR / "genes_unmapped/PE/{sample}_{rn}_R2.fq.gz"),
    shell:
        """
        {PATH[samtools]} fastq -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n {input}
        """


rule hisat2_3n_mapping_genome_PE:
    input:
        fq1=TEMPDIR / "genes_unmapped/PE/{sample}_{rn}_R1.fq.gz",
        fq2=TEMPDIR / "genes_unmapped/PE/{sample}_{rn}_R2.fq.gz",
    output:
        bam=temp(TEMPDIR / "genome_mapping/PE/{sample}_{rn}.bam"),
        summary="report_reads/mapping/{sample}_{rn}_genome.summary",
    params:
        index=config.get("genome_index"),
        basechange=BASE_CHANGE,
        directional="--directional-mapping",
        splice_args=(
            "--pen-noncansplice 20 --min-intronlen 20 --max-intronlen 20"
            if SPLICE_GENOME
            else "--no-spliced-alignment"
        ),
    threads: 24
    shell:
        """
        {PATH[hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -1 {input.fq1} -2 {input.fq2} --base-change {params.basechange} {params.directional} {params.splice_args} \
            --avoid-pseudogene --np 0 --rdg 5,3 --rfg 5,3 --sp 9,3 --mp 3,1 --score-min L,-3,-0.75 |\
            {PATH[samtools]} view -@ {threads} -O BAM -o {output.bam}
        """


rule filter_and_sort_bam_genome:
    input:
        TEMPDIR / "genome_mapping/PE/{sample}_{rn}.bam",
    output:
        unmap=temp(TEMPDIR / "genome_unmapped/{sample}_{rn}.bam"),
        mapped=temp(TEMPDIR / "genome_bam/{sample}_{rn}.bam"),
    params:
        flt="flag.proper_pair && !flag.unmap && !flag.munmap",
    threads: 16
    shell:
        """
        {PATH[samtools]} view -e '{params.flt}' -@ {threads} -U {output.unmap} -h {input} |\
            {PATH[samtools]} sort -m 3G -O BAM -o {output.mapped}
        """


rule extract_unmapped_reads_PE_from_genome:
    input:
        un=TEMPDIR / "genome_unmapped/{sample}_{rn}.bam",
    output:
        r1=INTERNALDIR / "discarded_reads/{sample}_{rn}_unmap_R1.fq.gz",
        r2=INTERNALDIR / "discarded_reads/{sample}_{rn}_unmap_R2.fq.gz",
    shell:
        """
        {PATH[samtools]} fastq -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n {input}
        """


rule combine_runs:
    input:
        lambda wildcards: [
            (
                TEMPDIR / f"genes_bam/{wildcards.sample}_{r}.bam"
                if wildcards.reftype == "genes"
                else TEMPDIR / f"genome_bam/{wildcards.sample}_{r}.bam"
            )
            for r in SAMPLE2DATA[wildcards.sample]
        ],
    output:
        bam=temp(TEMPDIR / "combined_runs/{sample}.{reftype}.bam"),
        bai=temp(TEMPDIR / "combined_runs/{sample}.{reftype}.bam.bai"),
    threads: 16
    shell:
        """
        {PATH[samtools]} merge -@ {threads} -f --write-index -o {output.bam}##idx##{output.bai} {input}
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
        {PATH[samtools]} flagstat -@ {threads} -O TSV {input} > {output.stat}
        {PATH[samtools]} view -@ {threads} -c -F 384 {input} > {output.n}
        """


rule drop_duplicates:
    input:
        bam=TEMPDIR / "combined_runs/{sample}.{reftype}.bam",
    output:
        bam=INTERNALDIR / "aligned_bam/{sample}.{reftype}.bam",
        txt="report_reads/dedup/{sample}.{reftype}.log",
    params:
        path_umicollapse=config["path"]["umicollapse"],
        tempdir=TEMPDIR,
    threads: 16
    run:
        if WITH_UMI:
            shell(
                """
            java -server -Xms8G -Xmx36G -Xss100M -Djava.io.tmpdir={params.tempdir} -jar {PATH[umicollapse]} bam \
                -t 2 -T {threads} --data naive --merge avgqual --two-pass -i {input.bam} -o {output.bam} >{output.txt}
            """
            )
        elif MARKDUP:
            shell(
                """
            java -Xms8G -Xmx36G -Xss100M -jar {PATH[gatk]} MarkDuplicates \
                -I {input} -O {output.bam} -M {output.txt} \
                --DUPLICATE_SCORING_STRATEGY SUM_OF_BASE_QUALITIES --REMOVE_DUPLICATES true --VALIDATION_STRINGENCY SILENT --TMP_DIR {params.tempdir}
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
        {PATH[samtools]} index -@ {threads} {input}
        """


rule stat_dedup:
    input:
        bam=INTERNALDIR / "aligned_bam/{sample}.{reftype}.bam",
        bai=INTERNALDIR / "aligned_bam/{sample}.{reftype}.bam.bai",
    output:
        stat="report_reads/dedup/{sample}.{reftype}.txt",
        n="report_reads/dedup/{sample}.{reftype}.count",
    threads: 2
    shell:
        """
        {PATH[samtools]} flagstat -@ {threads} -O TSV {input.bam} > {output.stat}
        {PATH[samtools]} view -@ {threads} -c -F 384 {input.bam} > {output.n}
        """


rule hisat2_3n_calling_filtered:
    input:
        INTERNALDIR / "aligned_bam/{sample}.{reftype}.bam",
    output:
        "count_sites/pileup/{sample}.{reftype}.tsv.gz",
    params:
        fa=lambda wildcards: (
            INTERNALDIR / "genes_index/genes.fa"
            if wildcards.reftype == "genes"
            else GENOME_FASTA[0]
        ),
        basechange=BASE_CHANGE,
    threads: 16
    shell:
        """
        {PATH[samtools]} view -@ {threads} -e "rlen < 100000 && [XM] * 20 <= (qlen-sclen) && [Zf] <= 3 && 3 * [Zf] <= [Zf] + [Yf]" -h {input} | \
            {PATH[hisat3ntable]} -p {threads} --alignments - --ref {params.fa} --output-name /dev/stdout --base-change {params.basechange} | \
            cut -f 1,2,3,5,7 | \
            bgzip -@ {threads} -c > {output}
        """
