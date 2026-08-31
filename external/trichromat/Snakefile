from snakemake.utils import validate
from workflow_utils import preprocess_config


configfile: "default.yaml"


validate(config, "config.schema.yaml")


if workflow._workdir_handler is not None:
    WORKDIR = workflow._workdir_handler.workdir
else:
    WORKDIR = os.path.relpath(
        config.get("workdir", "workspace"), os.path.dirname(workflow.configfiles[-1])
    )

    workdir: WORKDIR


preprocess_config(config, WORKDIR)


configfile: "default.yaml"


include: "workflow/mapping.smk"
include: "workflow/calling.smk"


rule all:
    input:
        "report_reads/read_counts_summary.tsv",
        expand("report_sites/filtered/{reftype}.tsv.gz", reftype=["genes", "genome"]),
