from collections import defaultdict
from external.trichromat.workflow_utils import preprocess_config


configfile: "default.yaml"


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
BASE_CHANGE = config.get("base_change", "A,G")

# LIBTYPE = config["libtype"]
# WITH_UMI = LIBTYPE in ["ECLIP10", "ECLIP6", "TAKARAV3", "SACSEQV3"]
# STRANDNESS = config.get("strandness", True)
# MARKDUP = config.get("markdup", True)
# GENE_NORC = config.get("gene_norc", True)
# SPLICE_GENOME = config.get("splice_genome", True)

preprocess_config(config, WORKDIR)
REF = config["_REF"]
READS = config["_READS"]


rule all:
    input:
        expand(
            "report_reads/unmap/{sample}.count",
            sample=READS.keys(),
        ),
        expand(
            "report_reads/dedup/{sample}.{reftype}.count",
            sample=READS.keys(),
            reftype=REF.keys(),
        ),
        expand(
            "report_sites/filtered/{reftype}.tsv.gz",
            reftype=["genes", "genome"],
        ),


module trichromat_workflow:
    snakefile:
        # You can also link to the trichromat workflow directly from github
        # github("y9c/trichromat", path="Snakefile", tag="main")
        # debugging only
        # "../trichromat/Snakefile"
        "external/trichromat/Snakefile"
    config:
        config


use rule * from trichromat_workflow as trichromat_*


rule hisat2_3n_calling_filtered:
    input:
        INTERNALDIR / "aligned_bam/{sample}.{reftype}.bam",
    output:
        "report_sites/pileup/{sample}.{reftype}.tsv.gz",
    params:
        fa=lambda wildcards: (
            INTERNALDIR / "genes_index/genes.fa"
            if wildcards.reftype == "genes"
            else REF["genome"][0]
        ),
        basechange=BASE_CHANGE,
    threads: 16
    shell:
        """
        {PATH[samtools]} view -@ {threads} -e "rlen < 100000 && [XM] * 20 <= (qlen-sclen) && [Zf] <= 3 && 3 * [Zf] <= [Zf] + [Yf]" -h {input} | \
            {PATH[hisat3ntable]} -p {threads} --alignments - --ref {params.fa} --output-name /dev/stdout --base-change {params.basechange} | \
            cut -f 1,2,3,5,7 | \
            gzip -c > {output}
        """


rule join_pileup_table:
    input:
        expand(
            "report_sites/pileup/{sample}.{{reftype}}.tsv.gz",
            sample=READS.keys(),
        ),
    output:
        "report_sites/joined/{reftype}.arrow",
    params:
        names=list(READS.keys()),
    threads: 8
    shell:
        """
        {PATH[joinPileup]} -f {input} -n {params.names} -o {output}
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
