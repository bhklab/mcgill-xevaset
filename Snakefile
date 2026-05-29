from pathlib import Path


configfile: "config/pipeline.yaml"


rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])


include: "workflow/rules/curation.smk"


rule all:
    input:
        xeva=results / "Xeva_McGill.rds",
        xeva_csv_archive=results / "Xeva_McGill_csv.tar.gz",
        unmatched_treatments=results / "unmatched_treatments.csv",
        unmapped_genes=results / "unmapped_genes.csv",
    localrule: True
