snakemake@source("helpers.R")
snk <- parse_snakemake()

suppressPackageStartupMessages({
  library(Biobase)
  library(dplyr)
  library(readxl)
  library(S4Vectors)
  library(SummarizedExperiment)
})

rna <- readRDS(snk$input$rnaseq)
assay_mat <- Biobase::exprs(rna)

clinical <- readxl::read_excel(snk$input$clinical) |>
  as.data.frame(stringsAsFactors = FALSE)
names(clinical) <- clean_names_simple(names(clinical))
clinical_map <- clinical |>
  transmute(
    sampleid = normalize_chr(.data$sample_id),
    patient.id.clinical = normalize_chr(.data$gcrc_id2),
    Source.clinical = normalize_chr(.data$source)
  ) |>
  filter(!is.na(.data$sampleid))

crosswalk <- readxl::read_excel(snk$input$crosswalk) |>
  as.data.frame(stringsAsFactors = FALSE) |>
  sample_ids_from_crosswalk()
rna_map <- bind_rows(
  crosswalk |>
    transmute(
      sampleid = .data$RNAseq.PDX,
      patient.id.crosswalk = .data$patient.id,
      Source.crosswalk = "PDX"
    ),
  crosswalk |>
    transmute(
      sampleid = .data$RNAseq.Patient,
      patient.id.crosswalk = .data$patient.id,
      Source.crosswalk = "Patient"
    )
) |>
  filter(!is.na(.data$sampleid))

col_df <- as.data.frame(Biobase::pData(rna), stringsAsFactors = FALSE)
col_df$sampleid <- colnames(assay_mat)
col_df$biobase.id <- col_df$sampleid
col_df <- col_df |>
  left_join(clinical_map, by = "sampleid") |>
  left_join(rna_map, by = "sampleid") |>
  mutate(
    patient.id = coalesce(
      .data$patient.id.clinical,
      .data$patient.id.crosswalk
    ),
    Source = coalesce(.data$Source.clinical, .data$Source.crosswalk)
  )
rownames(col_df) <- colnames(assay_mat)

row_df <- as.data.frame(Biobase::fData(rna), stringsAsFactors = FALSE)
if (!nrow(row_df)) {
  row_df <- data.frame(gene_id = rownames(assay_mat), stringsAsFactors = FALSE)
}
rownames(row_df) <- rownames(assay_mat)

se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(exprs = as.matrix(assay_mat)),
  rowData = S4Vectors::DataFrame(row_df),
  colData = S4Vectors::DataFrame(col_df)
)
S4Vectors::metadata(se) <- list(datatype = "RNAseq", annotation = "RNAseq")

fs::dir_create(fs::path_dir(snk$output$se))
saveRDS(se, snk$output$se)
write_matrix_csv(
  matrix = assay_mat,
  row_ids = rownames(assay_mat),
  row_id_name = "feature_id",
  output_path = snk$output$matrix
)
