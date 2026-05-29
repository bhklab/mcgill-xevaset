# McGill XevaSet Pipeline

This repository curates the McGill Park lab breast cancer PDX data into a `XevaSet`.

The workflow uses Snakemake for execution and pixi for package dependencies. Raw inputs are expected under `data/rawdata/`; the pipeline writes intermediates under `data/procdata/` and final outputs under `data/results/`.

Gene-level molecular feature metadata is rebuilt against the current GENCODE release tracked by AnnotationGx, and unmapped features are reported in `data/results/unmapped_genes.csv`.

The Xeva dependency stack currently resolves on `linux-64` and `osx-64`. On Apple Silicon, Pixi falls back to the `osx-64` environment under Rosetta.

Run from the repository root:

```bash
pixi install
pixi run snakemake --cores <n>
```

Knit the QC report with:

```bash
pixi run qc
```

Pipeline behavior and paths are configured in `config/pipeline.yaml`.
