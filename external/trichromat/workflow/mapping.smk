# mapping workflow

import shlex

TEMPDIR = Path(
    os.path.relpath(
        config.get("tempdir", os.path.join(workflow.basedir, ".tmp")), workflow.basedir
    )
)
INTERNALDIR = Path("internal_files")

PATH = config["path"]

LIBTYPE = config.get("libtype", "")

MARKDUP = config.get("markdup", True)
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


# rule all:
#     input:
#         expand("report_reads/unmap/{sample}.count", sample=READS.keys()),
#         expand(
#             "report_reads/combined/{sample}.{reftype}.count",
#             sample=READS.keys(),
#             reftype=REF.keys(),
#         ),
#         expand(
#             "report_reads/dedup/{sample}.{reftype}.count",
#             sample=READS.keys(),
#             reftype=REF.keys(),
#         ),


rule cutadapt_SE:
    input:
        lambda wildcards: READS[wildcards.sample][wildcards.rn].get("R1", "/"),
    output:
        c=temp(TEMPDIR / "trimmed_reads/SE/{sample}_{rn}_R1.fq.gz"),
        s=INTERNALDIR / "discarded_reads/{sample}_{rn}_tooshort_R1.fq.gz",
        report=temp(TEMPDIR / "trimmed_reads/SE/{sample}_{rn}.json"),
    params:
        cut=lambda wildcards: (
            f"-A {shlex.quote(SAMPLE_LIB[wildcards.sample])}"
            if SAMPLE_LIB[wildcards.sample]
            else f"-a {shlex.quote(SAMPLE_ADP[wildcards.sample])}"
        ),
    threads: 16
    shell:
        """
        {config[path][cutseq]} -t {threads} {params.cut} -m 20 --auto-rc -o {output.c} -s {output.s} --json-file {output.report} {input} 
        """


rule cutadapt_PE:
    input:
        lambda wildcards: READS[wildcards.sample][wildcards.rn].get("R1", "/"),
        lambda wildcards: READS[wildcards.sample][wildcards.rn].get("R2", "/"),
    output:
        c=[
            temp(TEMPDIR / "trimmed_reads/PE/{sample}_{rn}_R1.fq.gz"),
            temp(TEMPDIR / "trimmed_reads/PE/{sample}_{rn}_R2.fq.gz"),
        ],
        s=[
            INTERNALDIR / "discarded_reads/{sample}_{rn}_tooshort_R1.fq.gz",
            INTERNALDIR / "discarded_reads/{sample}_{rn}_tooshort_R2.fq.gz",
        ],
        report=temp(TEMPDIR / "trimmed_reads/PE/{sample}_{rn}.json"),
    params:
        cut=lambda wildcards: (
            f"-A {shlex.quote(SAMPLE_LIB[wildcards.sample])}"
            if SAMPLE_LIB[wildcards.sample]
            else f"-a {shlex.quote(SAMPLE_ADP[wildcards.sample])}"
        ),
    threads: 16
    shell:
        """
        {config[path][cutseq]} -t {threads} {params.cut} -m 20 --auto-rc -o {output.c} -s {output.s} --json-file {output.report} {input} 
        """


rule move_cutadapt_report:
    input:
        lambda wildcards: (
            TEMPDIR / "trimmed_reads/PE/{sample}_{rn}.json"
            if len(READS[wildcards.sample][wildcards.rn]) == 2
            else TEMPDIR / "trimmed_reads/SE/{sample}_{rn}.json"
        ),
    output:
        "report_reads/trimming/{sample}_{rn}.json",
    shell:
        """
        mv {input} {output}
        """


# map to contamination (drop reads)


rule combine_contamination_fa:
    input:
        REF.get("contamination", []),
    output:
        fa=INTERNALDIR / "reference_file/contamination.fa",
        fai=INTERNALDIR / "reference_file/contamination.fa.fai",
    shell:
        """
        zcat -f {input} > {output.fa}
        {config[path][samtools]} faidx {output.fa} --fai-idx {output.fai}
        """


rule build_contamination_hisat3n_index:
    input:
        INTERNALDIR / "reference_file/contamination.fa",
    output:
        INTERNALDIR / "reference_file/contamination.3n.GA.1.ht2",
    params:
        prefix=str(INTERNALDIR / "reference_file/contamination"),
        base_change=BASE_CHANGE,
    threads: 4
    shell:
        """
        rm -f {params.prefix}*.ht2
        {config[path][hisat3nbuild]} -p {threads} --base-change {params.base_change} {input} {params.prefix}
        """


rule hisat2_3n_mapping_contamination_SE:
    input:
        fq=TEMPDIR / "trimmed_reads/SE/{sample}_{rn}_R1.fq.gz",
        idx=INTERNALDIR / "reference_file/contamination.3n.GA.1.ht2",
    output:
        mapped=temp(
            TEMPDIR / "contamination_mapping/SE/{sample}_{rn}.contamination.bam"
        ),
        unmapped=temp(TEMPDIR / "contamination_mapping/SE/{sample}_{rn}.unmap.bam"),
        summary="report_reads/premap/{sample}_{rn}_SE.summary",
    params:
        index=str(INTERNALDIR / "reference_file/contamination"),
        directional=lambda wildcards: (
            "--directional-mapping" if SAMPLE_STD[wildcards.sample] else ""
        ),
        splice_args=(
            "--pen-noncansplice 20 --min-intronlen 20 --max-intronlen 20"
            if SPLICE_CONTAM
            else "--no-spliced-alignment"
        ),
        base_change=BASE_CHANGE,
    threads: 24
    shell:
        """
        {config[path][hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -U {input.fq} --base-change {params.base_change} {params.directional} {params.splice_args} \
            --np 0 --rdg 5,3 --rfg 5,3 --sp 9,3 --mp 3,1 --score-min L,-2,-0.8 |\
            {config[path][samtools]} view -@ {threads} -e '!flag.unmap && qlen-sclen >= 30 && [XM] * 15 < qlen-sclen' -O BAM -U {output.unmapped} -o {output.mapped}
        """


rule hisat2_3n_mapping_contamination_PE:
    input:
        fq1=TEMPDIR / "trimmed_reads/PE/{sample}_{rn}_R1.fq.gz",
        fq2=TEMPDIR / "trimmed_reads/PE/{sample}_{rn}_R2.fq.gz",
        idx=INTERNALDIR / "reference_file/contamination.3n.GA.1.ht2",
    output:
        mapped=temp(
            TEMPDIR / "contamination_mapping/PE/{sample}_{rn}.contamination.bam"
        ),
        unmapped=temp(TEMPDIR / "contamination_mapping/PE/{sample}_{rn}.unmap.bam"),
        summary="report_reads/premap/{sample}_{rn}_PE.summary",
    params:
        index=str(INTERNALDIR / "reference_file/contamination"),
        directional=lambda wildcards: (
            "--directional-mapping" if SAMPLE_STD[wildcards.sample] else ""
        ),
        splice_args=(
            "--pen-noncansplice 20 --min-intronlen 20 --max-intronlen 20"
            if SPLICE_CONTAM
            else "--no-spliced-alignment"
        ),
        base_change=BASE_CHANGE,
    threads: 24
    shell:
        """
        {config[path][hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -1 {input.fq1} -2 {input.fq2} --base-change {params.base_change} {params.directional} {params.splice_args} \
            --no-discordant --no-mixed \
            --np 0 --rdg 5,3 --rfg 5,3 --sp 9,3 --mp 3,1 --score-min L,-2,-0.8 |\
            {config[path][samtools]} view -@ {threads} -e 'flag.proper_pair && !flag.unmap && !flag.munmap && qlen-sclen >= 30 && [XM] * 15 < (qlen-sclen)' -O BAM -U {output.unmapped} -o {output.mapped}
        """


rule fix_hisat_tlen_bug_premap_SE:
    input:
        TEMPDIR / "contamination_mapping/SE/{sample}_{rn}.contamination.bam",
    output:
        temp(TEMPDIR / "contamination_fixmate/SE/{sample}_{rn}.bam"),
    threads: 4
    shell:
        "cp {input} {output}"


rule fix_hisat_tlen_bug_premap_PE:
    input:
        TEMPDIR / "contamination_mapping/PE/{sample}_{rn}.contamination.bam",
    output:
        temp(TEMPDIR / "contamination_fixmate/PE/{sample}_{rn}.bam"),
    threads: 4
    shell:
        "{config[path][samtools]} fixmate -@ {threads} -m -O BAM {input} {output}"


rule sort_bam_contamination:
    input:
        TEMPDIR / "contamination_fixmate/{mode}/{sample}_{rn}.bam",
    output:
        temp(TEMPDIR / "contamination_bam/{mode}/{sample}_{rn}.bam"),
    threads: 24
    shell:
        """
        {config[path][samtools]} sort -@ {threads} -m 3G -O BAM -o {output} {input}
        """


rule extract_unmapped_from_premap_SE:
    input:
        un=TEMPDIR / "contamination_mapping/SE/{sample}_{rn}.unmap.bam",
    output:
        r1=temp(TEMPDIR / "prefilter_reads/SE/{sample}_{rn}_R1.fq.gz"),
    shell:
        """
        {config[path][samtools]} fastq -0 {output.r1} -n {input}
        """


rule extract_unmapped_from_premap_PE:
    input:
        un=TEMPDIR / "contamination_mapping/PE/{sample}_{rn}.unmap.bam",
    output:
        r1=temp(TEMPDIR / "prefilter_reads/PE/{sample}_{rn}_R1.fq.gz"),
        r2=temp(TEMPDIR / "prefilter_reads/PE/{sample}_{rn}_R2.fq.gz"),
    shell:
        """
        {config[path][samtools]} fastq -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n {input}
        """


# map to genes


rule prepare_genes_index:
    input:
        REF.get("genes", []),
    output:
        fa=INTERNALDIR / "reference_file/genes.fa",
        index=INTERNALDIR / "reference_file/genes.3n.CT.1.ht2",
    params:
        prefix=INTERNALDIR / "reference_file/genes",
        base_change=BASE_CHANGE,
    threads: 4
    shell:
        """
        zcat -f {input} > {output.fa}
        {config[path][samtools]} faidx {output.fa}
        rm -f {params.prefix}.3n.*.ht2
        {config[path][hisat3nbuild]} -p {threads} --base-change {params.base_change} {output.fa} {params.prefix}
        """


rule hisat2_3n_mapping_genes_SE:
    input:
        fq1=(
            TEMPDIR / "prefilter_reads/SE/{sample}_{rn}_R1.fq.gz"
            if "contamination" in REF
            else TEMPDIR / "trimmed_reads/SE/{sample}_{rn}_R1.fq.gz"
        ),
        index=INTERNALDIR / "reference_file/genes.3n.CT.1.ht2",
    output:
        bam=temp(TEMPDIR / "genes_mapping/SE/{sample}_{rn}.bam"),
        summary="report_reads/mapping/{sample}_{rn}.genes.summary",
    params:
        index=INTERNALDIR / "reference_file/genes",
        directional=lambda wildcards: (
            "--directional-mapping" if SAMPLE_STD[wildcards.sample] else ""
        ),
        base_change=BASE_CHANGE,
    threads: 24
    shell:
        """
        {config[path][hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -U {input.fq1} --base-change {params.base_change} {params.directional} \
            --avoid-pseudogene --no-softclip --no-spliced-alignment --np 0 --rdg 5,3 --rfg 5,3 --mp 3,1 --score-min L,-3,-0.75 |\
            {config[path][samtools]} view -@ {threads} -O BAM -o {output.bam}
        """


rule hisat2_3n_mapping_genes_PE:
    input:
        fq1=(
            TEMPDIR / "prefilter_reads/PE/{sample}_{rn}_R1.fq.gz"
            if "contamination" in REF
            else TEMPDIR / "trimmed_reads/PE/{sample}_{rn}_R1.fq.gz"
        ),
        fq2=(
            TEMPDIR / "prefilter_reads/PE/{sample}_{rn}_R2.fq.gz"
            if "contamination" in REF
            else TEMPDIR / "trimmed_reads/PE/{sample}_{rn}_R2.fq.gz"
        ),
        index=INTERNALDIR / "reference_file/genes.3n.CT.1.ht2",
    output:
        bam=temp(TEMPDIR / "genes_mapping/PE/{sample}_{rn}.bam"),
        summary="report_reads/mapping/{sample}_{rn}.genes.summary",
    params:
        index=INTERNALDIR / "reference_file/genes",
        directional=lambda wildcards: (
            "--directional-mapping" if SAMPLE_STD[wildcards.sample] else ""
        ),
        base_change=BASE_CHANGE,
    threads: 24
    shell:
        """
        {config[path][hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -1 {input.fq1} -2 {input.fq2} --base-change {params.base_change} {params.directional} \
            --avoid-pseudogene --no-softclip --no-spliced-alignment --np 0 --rdg 5,3 --rfg 5,3 --mp 3,1 --score-min L,-3,-0.75 |\
            {config[path][samtools]} view -@ {threads} -O BAM -o {output.bam}
        """


rule filter_and_sort_bam_genes:
    input:
        TEMPDIR / "genes_mapping/{mode}/{sample}_{rn}.bam",
    output:
        unmap=temp(TEMPDIR / "genes_unmapped/{mode}/{sample}_{rn}.bam"),
        mapped=temp(TEMPDIR / "genes_bam/{mode}/{sample}_{rn}.bam"),
    params:
        flt=lambda wildcards: (
            "flag.proper_pair && !flag.unmap && !flag.munmap"
            + (
                " && !flag.read1 != !flag.reverse"
                if (SAMPLE_NORC[wildcards.sample] and SAMPLE_STD[wildcards.sample])
                else ""
            )
            if wildcards.mode == "PE"
            else "!flag.unmap"
            + (
                " && !flag.reverse"
                if (SAMPLE_NORC[wildcards.sample] and SAMPLE_STD[wildcards.sample])
                else ""
            )
        ),
    threads: 16
    shell:
        """
        {config[path][samtools]} view -e '{params.flt}' -@ {threads} -u -U {output.unmap} {input} |\
            {config[path][samtools]} fixmate -@ {threads} -m -u - - |\
            {config[path][samtools]} sort -m 3G -O BAM -o {output.mapped}
        """


rule extract_unmapped_reads_SE_from_genes:
    input:
        TEMPDIR / "genes_unmapped/SE/{sample}_{rn}.bam",
    output:
        temp(TEMPDIR / "genes_unmapped/SE/{sample}_{rn}_R1.fq.gz"),
    threads: 4
    shell:
        """
        {config[path][samtools]} fastq -@ {threads} -0 {output} {input}
        """


rule extract_unmapped_reads_PE_from_genes:
    input:
        TEMPDIR / "genes_unmapped/PE/{sample}_{rn}.bam",
    output:
        r1=temp(TEMPDIR / "genes_unmapped/PE/{sample}_{rn}_R1.fq.gz"),
        r2=temp(TEMPDIR / "genes_unmapped/PE/{sample}_{rn}_R2.fq.gz"),
    threads: 4
    shell:
        """
        {config[path][samtools]} fastq -@ {threads} -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n {input}
        """


# map to genome


rule hisat2_3n_mapping_genome_SE:
    input:
        fq1=TEMPDIR / "genes_unmapped/SE/{sample}_{rn}_R1.fq.gz",
    output:
        bam=temp(TEMPDIR / "genome_mapping/SE/{sample}_{rn}.bam"),
        summary="report_reads/mapping/{sample}_{rn}.genome.summary",
    params:
        index=config.get("genome_index"),
        directional=lambda wildcards: (
            "--directional-mapping" if SAMPLE_STD[wildcards.sample] else ""
        ),
        splice_args=(
            "--pen-noncansplice 20 --min-intronlen 20 --max-intronlen 20"
            if SPLICE_GENOME
            else "--no-spliced-alignment"
        ),
        base_change=BASE_CHANGE,
    threads: 24
    shell:
        """
        {config[path][hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -U {input.fq1} --base-change {params.base_change} {params.directional} {params.splice_args} \
            --avoid-pseudogene --np 0 --rdg 5,3 --rfg 5,3 --sp 9,3 --mp 3,1 --score-min L,-3,-0.75 |\
            {config[path][samtools]} view -@ {threads} -O BAM -o {output.bam}
        """


rule hisat2_3n_mapping_genome_PE:
    input:
        fq1=TEMPDIR / "genes_unmapped/PE/{sample}_{rn}_R1.fq.gz",
        fq2=TEMPDIR / "genes_unmapped/PE/{sample}_{rn}_R2.fq.gz",
    output:
        bam=temp(TEMPDIR / "genome_mapping/PE/{sample}_{rn}.bam"),
        summary="report_reads/mapping/{sample}_{rn}.genome.summary",
    params:
        index=config.get("genome_index"),
        directional=lambda wildcards: (
            "--directional-mapping" if SAMPLE_STD[wildcards.sample] else ""
        ),
        splice_args=(
            "--pen-noncansplice 20 --min-intronlen 20 --max-intronlen 20"
            if SPLICE_GENOME
            else "--no-spliced-alignment"
        ),
        base_change=BASE_CHANGE,
    threads: 24
    shell:
        """
        {config[path][hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -1 {input.fq1} -2 {input.fq2} --base-change {params.base_change} {params.directional} {params.splice_args} \
            --avoid-pseudogene --np 0 --rdg 5,3 --rfg 5,3 --sp 9,3 --mp 3,1 --score-min L,-3,-0.75 |\
            {config[path][samtools]} view -@ {threads} -O BAM -o {output.bam}
        """


rule filter_and_sort_bam_genome:
    input:
        TEMPDIR / "genome_mapping/{mode}/{sample}_{rn}.bam",
    output:
        unmap=temp(TEMPDIR / "genome_unmapped/{mode}/{sample}_{rn}.bam"),
        mapped=temp(TEMPDIR / "genome_bam/{mode}/{sample}_{rn}.bam"),
    params:
        flt=lambda wildcards: (
            "flag.proper_pair && !flag.unmap && !flag.munmap"
            if wildcards.mode == "PE"
            else "!flag.unmap"
        ),
    threads: 16
    shell:
        # samtools collate -@ 4 -O -u example.bam | 
        """
        {config[path][samtools]} view -e '{params.flt}' -@ {threads} -u -U {output.unmap} {input} |\
            {config[path][samtools]} fixmate -@ {threads} -m -u - - |\
            {config[path][samtools]} sort -m 3G -O BAM -o {output.mapped}
        """


ruleorder: extract_unmapped_reads_PE_from_genome > extract_unmapped_reads_SE_from_genome


rule extract_unmapped_reads_SE_from_genome:
    input:
        TEMPDIR / "genome_unmapped/SE/{sample}_{rn}.bam",
    output:
        INTERNALDIR / "discarded_reads/{sample}_{rn}_unmap_R1.fq.gz",
    threads: 4
    shell:
        """
        {config[path][samtools]} fastq -@ {threads} -0 {output} {input}
        """


rule extract_unmapped_reads_PE_from_genome:
    input:
        TEMPDIR / "genome_unmapped/PE/{sample}_{rn}.bam",
    output:
        r1=INTERNALDIR / "discarded_reads/{sample}_{rn}_unmap_R1.fq.gz",
        r2=INTERNALDIR / "discarded_reads/{sample}_{rn}_unmap_R2.fq.gz",
    threads: 4
    shell:
        """
        {config[path][samtools]} fastq -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n {input}
        """


rule stat_unmapped:
    input:
        lambda wildcards: [
            INTERNALDIR / f"discarded_reads/{wildcards.sample}_{rn}_unmap_R1.fq.gz"
            for rn in READS[wildcards.sample].keys()
        ],
    output:
        "report_reads/unmap/{sample}.count",
    threads: 2
    shell:
        """
        zcat {input} | {config[path][samtools]} view -@ {threads} -c > {output}
        """


rule combine_runs:
    input:
        lambda wildcards: [
            (
                TEMPDIR / f"{wildcards.reftype}_bam/PE/{wildcards.sample}_{rn}.bam"
                if len(rs) == 2
                else TEMPDIR
                / f"{wildcards.reftype}_bam/SE/{wildcards.sample}_{rn}.bam"
            )
            for rn, rs in READS[wildcards.sample].items()
        ],
    output:
        bam=temp(TEMPDIR / "combined_runs/{sample}.{reftype}.bam"),
        bai=temp(TEMPDIR / "combined_runs/{sample}.{reftype}.bam.bai"),
    threads: 16
    shell:
        """
        {config[path][samtools]} merge -@ {threads} -f --write-index -o {output.bam}##idx##{output.bai} {input}
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
        {config[path][samtools]} flagstat -@ {threads} -O TSV {input} > {output.stat}
        {config[path][samtools]} view -@ {threads} -c -F 384 {input} > {output.n}
        """


# rule drop_duplicates_not_working:
#     input:
#         bam=TEMPDIR / "combined_runs/{sample}.{reftype}.bam",
#     output:
#         bam=INTERNALDIR / "aligned_bam/{sample}.{reftype}.bam",
#         txt="report_reads/dedup/{sample}.{reftype}.log",
#     params:
#         umi="--barcode-rgx '.*_([!-?A-~]+)'",
#         markdup=MARKDUP,
#     threads: 16
#      shell:
#          """
#          if [[ "{params.markdup}" == "True" ]]; then
#              {config[path][samtools]} markdup -@ {threads} -c -S -l 600 --mode t -t --include-fails --duplicate-count --write-index {params.umi} -f {output.txt} {input} {output.bam}
#          else
#              cp {input.bam} {output.bam} && touch {output.txt}
#          fi
#          """


rule drop_duplicates:
    input:
        bam=TEMPDIR / "combined_runs/{sample}.{reftype}.bam",
        bai=TEMPDIR / "combined_runs/{sample}.{reftype}.bam.bai",
    output:
        bam=INTERNALDIR / "aligned_bam/{sample}.{reftype}.bam",
        txt="report_reads/dedup/{sample}.{reftype}.log",
    params:
        with_umi=lambda wildcards: SAMPLE_UMI[wildcards.sample],
        markdup=lambda wildcards: SAMPLE_DUP[wildcards.sample],
        collapse_paired=lambda wildcards: (
            "--paired --remove-unpaired"
            if max(len(rd) for rn, rd in READS[wildcards.sample].items()) == 2
            else ""
        ),
        flowmode_single=lambda wildcards: (
            "--FLOW_MODE true --FLOW_USE_END_IN_UNPAIRED_READS true --FLOW_USE_UNPAIRED_CLIPPED_END false"
            if max(len(rd) for rn, rd in READS[wildcards.sample].items()) == 1
            else ""
        ),
        java_tmp=TEMPDIR,
    threads: 16
    run:
        if params.with_umi:
            shell(
                """
                mkdir -p {params.java_tmp}
                {config[path][umicollapse]} -Djava.io.tmpdir={params.java_tmp} \
                    -T {threads} --remove-chimeric --data bktree --merge avgqual \
                    {params.collapse_paired} \
                    -i {input.bam} -o {output.bam} >{output.txt}
                """
            )
        elif params.markdup:
            shell(
                """
                mkdir -p {params.java_tmp}
                INTER_BAM={params.java_tmp}/{wildcards.sample}.{wildcards.reftype}.bam 
                {config[path][picard]} -Djava.io.tmpdir={params.java_tmp} AddOrReplaceReadGroups \
                    --RGID {wildcards.sample} --RGLB LIB --RGPL ILLUMINA --RGPU LANE --RGSM {wildcards.sample} --VALIDATION_STRINGENCY SILENT \
                    --INPUT {input.bam} --OUTPUT ${{INTER_BAM}}
                {config[path][picard]} -Djava.io.tmpdir={params.java_tmp} MarkDuplicates \
                    --DUPLICATE_SCORING_STRATEGY SUM_OF_BASE_QUALITIES --REMOVE_DUPLICATES true --VALIDATION_STRINGENCY SILENT --TAG_DUPLICATE_SET_MEMBERS true --ADD_PG_TAG_TO_READS false {params.flowmode_single} \
                    --TMP_DIR {params.java_tmp} \
                    --INPUT ${{INTER_BAM}} --OUTPUT {output.bam} --METRICS_FILE {output.txt}
                rm ${{INTER_BAM}}
                """
            )
        else:
            shell(
                """
                cp {input.bam} {output.bam}
                touch {output.txt}
                """
            )


rule write_bam_index:
    input:
        INTERNALDIR / "aligned_bam/{sample}.{reftype}.bam",
    output:
        INTERNALDIR / "aligned_bam/{sample}.{reftype}.bam.bai",
    threads: 4
    shell:
        """
        {config[path][samtools]} index -@ {threads} -b {input} {output}
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
        {config[path][samtools]} flagstat -@ {threads} -O TSV {input.bam} > {output.stat}
        {config[path][samtools]} view -@ {threads} -c -F 384 {input.bam} > {output.n}
        """


rule collect_read_counts:
    input:
        trimming=[
            f"report_reads/trimming/{sample}_{rn}.json"
            for sample, rs in READS.items()
            for rn in rs
        ],
        unmap=expand("report_reads/unmap/{sample}.count", sample=READS.keys()),
        combined=expand(
            "report_reads/combined/{sample}.{reftype}.count",
            sample=READS.keys(),
            reftype=REF.keys(),
        ),
        dedup=expand(
            "report_reads/dedup/{sample}.{reftype}.count",
            sample=READS.keys(),
            reftype=REF.keys(),
        ),
    output:
        "report_reads/read_counts_summary.tsv",
    shell:
        """
        {config[path][collectReadCounts]} --trimming {input.trimming} --unmap {input.unmap} --combined {input.combined} --dedup {input.dedup} --output-name {output}
        """
