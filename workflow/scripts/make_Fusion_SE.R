snakemake@source("helpers.R")
snk <- parse_snakemake()

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(S4Vectors)
  library(SummarizedExperiment)
})

fusion_df <- read_tsv(snk$input$fusion)
names(fusion_df) <- clean_names_simple(names(fusion_df))

fusion_df <- fusion_df |>
  transmute(
    fusion = normalize_chr(.data$fusion),
    sampleid = sample_id_from_profile(.data$tumor_sample_barcode),
    Barcode = normalize_chr(.data$tumor_sample_barcode),
    hugo_symbol = normalize_chr(.data$hugo_symbol),
    entrez_gene_id = normalize_chr(.data$entrez_gene_id),
    center = normalize_chr(.data$center),
    dna_support = normalize_chr(.data$dna_support),
    rna_support = normalize_chr(.data$rna_support),
    method = normalize_chr(.data$method),
    frame = normalize_chr(.data$frame),
    fusion_status = normalize_chr(.data$fusion_status)
  ) |>
  filter(!is.na(.data$fusion), !is.na(.data$sampleid))

fusion_df$dna_support_value <- as.integer(
  tolower(fusion_df$dna_support) == "yes"
)
fusion_df$rna_support_value <- as.integer(
  tolower(fusion_df$rna_support) == "yes"
)
fusion_df$dna_support_value[is.na(fusion_df$dna_support_value)] <- 0L
fusion_df$rna_support_value[is.na(fusion_df$rna_support_value)] <- 0L

feature_df <- fusion_df |>
  group_by(.data$fusion) |>
  summarise(
    fusion_name = dplyr::first(.data$fusion),
    left_gene = sub("-.*$", "", dplyr::first(.data$fusion)),
    right_gene = sub("^[^-]+-", "", dplyr::first(.data$fusion)),
    hugo_symbol = collapse_unique(.data$hugo_symbol),
    entrez_gene_id = collapse_unique(.data$entrez_gene_id),
    center = collapse_unique(.data$center),
    method = collapse_unique(.data$method),
    frame = collapse_unique(.data$frame),
    fusion_status = collapse_unique(.data$fusion_status),
    .groups = "drop"
  ) |>
  arrange(.data$fusion)
feature_df <- as.data.frame(feature_df, stringsAsFactors = FALSE)
rownames(feature_df) <- feature_df$fusion

sample_ids <- sort(unique(fusion_df$sampleid))
presence_mat <- matrix(
  0L,
  nrow = nrow(feature_df),
  ncol = length(sample_ids),
  dimnames = list(feature_df$fusion, sample_ids)
)
dna_mat <- presence_mat
rna_mat <- presence_mat

call_df <- fusion_df |>
  group_by(.data$fusion, .data$sampleid) |>
  summarise(
    presence = 1L,
    dna_support = max(.data$dna_support_value),
    rna_support = max(.data$rna_support_value),
    .groups = "drop"
  )

call_idx <- cbind(
  match(call_df$fusion, rownames(presence_mat)),
  match(call_df$sampleid, colnames(presence_mat))
)
presence_mat[call_idx] <- call_df$presence
dna_mat[call_idx] <- call_df$dna_support
rna_mat[call_idx] <- call_df$rna_support

crosswalk <- readxl::read_excel(snk$input$crosswalk) |>
  as.data.frame(stringsAsFactors = FALSE) |>
  sample_ids_from_crosswalk()
sample_map <- molecular_sample_map_from_crosswalk(crosswalk)

barcode_df <- fusion_df |>
  distinct(.data$sampleid, .data$Barcode) |>
  mutate(
    Source.from.barcode = case_when(
      grepl("_Br_X_", .data$Barcode) ~ "PDX",
      grepl("_Br_M_", .data$Barcode) ~ "Metastasis",
      TRUE ~ "Patient"
    )
  )

col_df <- data.frame(
  sampleid = sample_ids,
  biobase.id = sample_ids,
  stringsAsFactors = FALSE
) |>
  left_join(barcode_df, by = "sampleid") |>
  left_join(sample_map, by = "sampleid") |>
  mutate(Source = coalesce(.data$Source, .data$Source.from.barcode))
col_df <- as.data.frame(col_df, stringsAsFactors = FALSE)
rownames(col_df) <- col_df$sampleid

se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(
    exprs = presence_mat,
    dna_support = dna_mat,
    rna_support = rna_mat
  ),
  rowData = S4Vectors::DataFrame(feature_df),
  colData = S4Vectors::DataFrame(col_df)
)
S4Vectors::metadata(se) <- list(datatype = "fusion", annotation = "fusion")

fs::dir_create(fs::path_dir(snk$output$se))
saveRDS(se, snk$output$se)
write_matrix_csv(
  matrix = presence_mat,
  row_ids = feature_df$fusion,
  row_id_name = "fusion",
  output_path = snk$output$matrix
)
