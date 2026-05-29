# McGill XevaSet Pipeline

This pipeline curates the McGill Park lab breast cancer PDX data into a `XevaSet`, plus CSV exports for non-R users. Molecular feature metadata for gene-level profiles is rebuilt against the current GENCODE release used by the pipeline.

Pixi targets `linux-64` and `osx-64` because the current Bioconda `bioconductor-xeva` dependency stack is not published for native `osx-arm64`.

Run from the repository root:

```bash
pixi run setup
pixi run snakemake --cores <n>
```

The setup task installs AnnotationGx from the configured remote GitHub ref.

Edit `config/pipeline.yaml` to change source paths, output directories, response settings, and the treatment AnnotationDB base URL.
