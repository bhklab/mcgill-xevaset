snakemake@source("helpers.R")
snk <- parse_snakemake()

suppressPackageStartupMessages({
  library(dplyr)
  library(S4Vectors)
  library(SummarizedExperiment)
})

model_df <- read_csv(snk$input$model) |>
  distinct(.data$model.id, .data$patient.id)

extract_mapping <- function(se_path) {
  se <- readRDS(se_path)
  datatype <- S4Vectors::metadata(se)$datatype
  col_df <- as.data.frame(
    SummarizedExperiment::colData(se),
    stringsAsFactors = FALSE
  )
  col_df$colname <- colnames(se)

  col_df |>
    transmute(
      patient.id = normalize_chr(.data$patient.id),
      biobase.id = normalize_chr(coalesce(
        .data$biobase.id,
        .data$sampleid,
        .data$colname
      )),
      source = normalize_chr(.data$Source),
      mDataType = datatype
    ) |>
    filter(!is.na(.data$patient.id), !is.na(.data$biobase.id))
}

mapping_df <- bind_rows(
  extract_mapping(snk$input$rnaseq),
  extract_mapping(snk$input$mutation),
  extract_mapping(snk$input$cnv),
  extract_mapping(snk$input$fusion),
  extract_mapping(snk$input$segments)
) |>
  inner_join(model_df, by = "patient.id", relationship = "many-to-many") |>
  select(
    "model.id",
    "biobase.id",
    "source",
    "mDataType"
  ) |>
  distinct() |>
  arrange(.data$mDataType, .data$model.id, .data$biobase.id)

write_csv(mapping_df, snk$output$modToBiobaseMap)
