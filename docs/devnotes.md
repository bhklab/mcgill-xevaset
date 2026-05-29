# Developer Notes

## Scope Decisions

- The configured workbook in `config/pipeline.yaml` is the treatment-response source of truth.
- The curation includes RNAseq, mutation, CNV, logCNV, fusion, and segment molecular profiles. Segment-level copy number data is stored as a matrix-style `SummarizedExperiment` so it can move through the current `XevaSet` and CSV export path.

## Directory Decisions

- Raw files remain under `data/rawdata/` and are not modified by the workflow.
- Pipeline intermediates are written under `data/procdata/`.
- Final outputs are written under `data/results/`.
- The unmatched treatment report is written directly under `data/results/`.

## Identifier Decisions

- PDX workbook sheets are treated as patient/model group identifiers, using the `GCRC####` sheet name as `patient.id`.
- Mouse-level model IDs follow the legacy McGill convention based on mouse ID, workbook row, GCRC ID, and a two-letter drug code.
- Controls are normalized to `untreated` in curated tables while retaining source labels in provenance columns.
- Molecular profiles map to models through patient-level `GCRC` IDs from `SA_GCRC_ID.xlsx` and clinical metadata.
- WGS-derived molecular profiles use WGS sample columns from `SA_GCRC_ID.xlsx`, with RNAseq sample columns as a fallback for samples present in cbioportal files but absent from the WGS columns.
- Gene-level molecular feature metadata is rebuilt from the current GENCODE release tracked by AnnotationGx instead of carrying older source gene annotations forward. RNAseq keeps the original Ensembl feature ID as provenance, while GENCODE fields supply the current gene symbol, biotype, coordinates, and release metadata. AnnotationGx direct gene helpers are available, but currently make one request per feature; the pipeline uses batched current-GENCODE lookups for runtime practicality.

## Metadata Decisions

- Tissue metadata is set to breast cancer / BRCA for all workbook models.
- Treatment annotation uses AnnotationDB through configurable paths in `config/pipeline.yaml`.
- Drug metadata intentionally avoids workbook-level aggregate provenance columns for `source.dose`, `source.schedule`, and `source.patients`. Per-model source dose and schedule are retained in model metadata because they vary by model record.
- Drug annotation keeps compact compound fields from AnnotationDB when available: `pubchem.cid`, `molecular.formula`, `smiles`, and `inchikey`.
- PDX model and sample metadata are curated from raw workbook, clinical, and crosswalk inputs only. The pipeline does not query or store cell-line metadata for this dataset.
- Annotation misses are retained in the curated dataset and listed in CSV reports with the source identifiers and available raw metadata used for matching.
- Reports intentionally avoid request-level columns such as URLs, endpoints, status codes, or call traces.

## Code Style Decisions

- Scripts are written as readable linear curation steps.
- Small helper functions are used only where repeated parsing or table writing would otherwise obscure the script.
- The workflow generally lets Snakemake and R report failures directly instead of wrapping every possible malformed-input condition.
- CSV outputs write missing values as empty fields rather than literal `NA` strings so readers that do not treat `NA` specially still import blanks.
- Empty source strings are normalized to missing values before final CSV export so missing fields are emitted as bare empty CSV fields, not quoted empty strings.
- Final export CSVs quote-wrap character-like values while leaving numeric values unquoted to make automatic CSV parser type inference less ambiguous.
- CSV exports exclude empty Xeva object scaffolding and session metadata tables that do not help non-R users consume the curated dataset.
- The former batch-definition/model-list export is omitted because the same model membership is already represented by `model.csv` and the sensitivity/model tables.

## Dependency Decisions

- The workflow declares `bioconductor-xeva` through Pixi/Bioconda.
- As of May 8, 2026, `bioconductor-xeva` 1.26.0 depends on `bioconductor-pharmacogx` 3.14.0, which Bioconda publishes for `linux-64` and `osx-64` but not native `osx-arm64`.
- Pixi target platforms are therefore set to `linux-64` and `osx-64`. Native Apple Silicon execution should use a Linux container/runner or an x86_64/Rosetta Pixi environment.
