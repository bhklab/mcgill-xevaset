# Outputs

Final pipeline outputs are written under `data/results/`.

| Output | Description |
| --- | --- |
| `McGill_MultiAssayExperiment.rds` | Integrated molecular assay object. |
| `Xeva_McGill.rds` | Curated McGill `XevaSet`. |
| `Xeva_McGill_csv/` | CSV export of curated metadata, response metrics, tumor measurements, molecular profile matrices, molecular feature metadata, and the mapping between models and molecular samples. |
| `Xeva_McGill_csv.tar.gz` | Snakemake-generated archive of the CSV export directory. |
| `unmatched_treatments.csv` | Source treatment records not matched during treatment annotation. |
| `unmapped_genes.csv` | Molecular features or fusion components not mapped to current GENCODE. |
| `qc/mcgill_xeva_qc.html` | Rendered QC walkthrough. |
