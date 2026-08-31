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
        "report_reads/read_counts_summary.tsv",
        expand(
            "report_sites/filtered/{reftype}.tsv.gz",
            reftype=["genes", "genome"],
        ),


module trichromat_workflow:
    snakefile:
        # trichromat is vendored (copied) into external/trichromat so this repo
        # is self-contained and needs no git submodules. To refresh it, copy the
        # latest code from https://github.com/y9c/trichromat into external/trichromat/.
        # You can also wire the upstream workflow directly instead:
        # github("y9c/trichromat", path="Snakefile", tag="main")
        "external/trichromat/Snakefile"
    config:
        config


use rule * from trichromat_workflow as trichromat_*
