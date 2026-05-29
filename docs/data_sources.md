# Data Sources

## Raw Inputs

The default raw inputs are configured in `config/pipeline.yaml`.

- `GCRC_PDX_1x1x1_Trial_BHK_April_2024.xlsx`: treatment-response workbook.
- `Morag_Park_clinical_info_April_2024.xlsx`: clinical metadata.
- `SA_GCRC_ID.xlsx`: sample, PDX, RNAseq, and WGS crosswalk.
- `pdxMorag_rnaseq.rda`: RNAseq `ExpressionSet`.
- `data_mutations_extended.txt`: cBioPortal-style mutation data.
- `data_CNA.txt`: discrete CNA matrix.
- `data_log2CNA.txt`: log2 CNA matrix.
- `data_fusions.txt`: cBioPortal-style fusion calls.
- `data_segments.txt`: cBioPortal-style segment-level copy number data.

## AnnotationDB

Treatment AnnotationDB lookup paths are configured through:

- `annotationdb.base_url`
- `annotationdb.compound_many_path`

The pipeline stores the final unmatched treatment report at `data/results/unmatched_treatments.csv`.

## GENCODE Feature Annotation

Gene-level molecular feature metadata is rebuilt from the current GENCODE release used by the pipeline. Species and lookup batch settings are configured under `annotationgx.gencode` in `config/pipeline.yaml`.

The pipeline uses batched current-GENCODE lookups for the feature tables so RNAseq, CNV, logCNV, mutation, and fusion annotations can be regenerated in a practical runtime while storing the release used in feature metadata.

The pipeline stores unmapped molecular features at `data/results/unmapped_genes.csv`.
