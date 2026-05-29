snakemake@source("helpers.R")
snk <- parse_snakemake()

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(S4Vectors)
  library(SummarizedExperiment)
})

mut <- read_tsv(snk$input$mutation) |>
  filter(
    !is.na(.data$Hugo_Symbol),
    nzchar(.data$Hugo_Symbol),
    !is.na(.data$Tumor_Sample_Barcode),
    nzchar(.data$Tumor_Sample_Barcode)
  )

mutation_calls <- mut |>
  distinct(
    .data$Hugo_Symbol,
    .data$Tumor_Sample_Barcode,
    .data$Variant_Classification
  ) |>
  group_by(.data$Hugo_Symbol, .data$Tumor_Sample_Barcode) |>
  summarise(
    variant = paste(sort(unique(.data$Variant_Classification)), collapse = ","),
    .groups = "drop"
  )

genes <- sort(unique(mut$Hugo_Symbol))
barcodes <- sort(unique(mut$Tumor_Sample_Barcode))
sample_ids <- sample_id_from_profile(barcodes)

mutation_mat <- matrix(
  "0",
  nrow = length(genes),
  ncol = length(barcodes),
  dimnames = list(genes, sample_ids)
)
mutation_mat[cbind(
  match(mutation_calls$Hugo_Symbol, genes),
  match(mutation_calls$Tumor_Sample_Barcode, barcodes)
)] <- mutation_calls$variant

crosswalk <- readxl::read_excel(snk$input$crosswalk) |>
  as.data.frame(stringsAsFactors = FALSE) |>
  sample_ids_from_crosswalk()
wgs_map <- molecular_sample_map_from_crosswalk(crosswalk)

col_df <- data.frame(
  sampleid = sample_ids,
  biobase.id = sample_ids,
  Barcode = barcodes,
  Source.from.barcode = ifelse(grepl("_Br_X_", barcodes), "PDX", "Patient"),
  stringsAsFactors = FALSE
) |>
  left_join(wgs_map, by = "sampleid") |>
  mutate(Source = coalesce(.data$Source, .data$Source.from.barcode))
rownames(col_df) <- sample_ids

row_df <- data.frame(gene = genes, stringsAsFactors = FALSE, row.names = genes)

se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(exprs = mutation_mat),
  rowData = S4Vectors::DataFrame(row_df),
  colData = S4Vectors::DataFrame(col_df)
)
S4Vectors::metadata(se) <- list(datatype = "mutation", annotation = "mutation")

fs::dir_create(fs::path_dir(snk$output$se))
saveRDS(se, snk$output$se)
write_matrix_csv(
  matrix = mutation_mat,
  row_ids = row_df$gene,
  row_id_name = "gene",
  output_path = snk$output$matrix
)
