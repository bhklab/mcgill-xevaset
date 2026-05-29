snakemake@source("helpers.R")
snk <- parse_snakemake()

suppressPackageStartupMessages({
  library(dplyr)
  library(fs)
  library(MultiAssayExperiment)
  library(S4Vectors)
  library(SummarizedExperiment)
})

sample_df <- read_csv(snk$input$sample)
rownames(sample_df) <- sample_df$sampleid

assay_paths <- unname(unlist(snk$input$assays))
assays <- lapply(assay_paths, readRDS)
assay_names <- vapply(
  assays,
  function(se) S4Vectors::metadata(se)$datatype,
  character(1)
)
names(assays) <- assay_names

sample_map <- lapply(names(assays), function(assay_name) {
  se <- assays[[assay_name]]
  col_df <- as.data.frame(
    SummarizedExperiment::colData(se),
    stringsAsFactors = FALSE
  )
  primary <- normalize_chr(col_df$patient.id)
  primary[is.na(primary)] <- normalize_chr(col_df$sampleid)[is.na(primary)]

  data.frame(
    assay = assay_name,
    primary = primary,
    colname = colnames(se),
    stringsAsFactors = FALSE
  )
}) |>
  bind_rows() |>
  filter(!is.na(.data$primary))
sample_map$assay <- factor(sample_map$assay, levels = names(assays))

missing_primary <- setdiff(unique(sample_map$primary), rownames(sample_df))
if (length(missing_primary) > 0) {
  missing_df <- data.frame(
    sampleid = missing_primary,
    patient.id = missing_primary,
    tissue.id = snk$config$metadata$tissue_id,
    tissue = snk$config$metadata$tissue,
    tissue.name = snk$config$metadata$tissue_name,
    stringsAsFactors = FALSE,
    row.names = missing_primary
  )
  for (col_name in setdiff(names(sample_df), names(missing_df))) {
    missing_df[[col_name]] <- NA
  }
  missing_df <- missing_df[, names(sample_df), drop = FALSE]
  sample_df <- rbind(sample_df, missing_df)
}

mae <- MultiAssayExperiment::MultiAssayExperiment(
  experiments = assays,
  colData = S4Vectors::DataFrame(sample_df),
  sampleMap = sample_map,
  metadata = list(dataset = snk$config$dataset$name)
)

dir_create(path_dir(snk$output$mae))
saveRDS(mae, snk$output$mae)
