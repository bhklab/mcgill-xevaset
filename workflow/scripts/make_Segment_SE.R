snakemake@source("helpers.R")
snk <- parse_snakemake()

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(S4Vectors)
  library(SummarizedExperiment)
})

segment_df <- read_tsv(snk$input$segments)
names(segment_df) <- clean_names_simple(names(segment_df))

segment_df <- segment_df |>
  transmute(
    Barcode = normalize_chr(.data$id),
    sampleid = sample_id_from_profile(.data$id),
    chrom = sub("^chr", "", normalize_chr(.data$chrom)),
    start = as.integer(.data$loc_start),
    end = as.integer(.data$loc_end),
    num_mark = as.numeric(.data$num_mark),
    seg_mean = as.numeric(.data$seg_mean)
  ) |>
  filter(
    !is.na(.data$sampleid),
    !is.na(.data$chrom),
    !is.na(.data$start),
    !is.na(.data$end),
    !is.na(.data$seg_mean)
  )
segment_df$feature_id <- paste0(
  "chr",
  segment_df$chrom,
  ":",
  segment_df$start,
  "-",
  segment_df$end
)

feature_df <- segment_df |>
  distinct(.data$feature_id, .data$chrom, .data$start, .data$end) |>
  arrange(.data$chrom, .data$start, .data$end)
feature_df <- as.data.frame(feature_df, stringsAsFactors = FALSE)
rownames(feature_df) <- feature_df$feature_id

sample_ids <- sort(unique(segment_df$sampleid))
segment_mat <- matrix(
  NA_real_,
  nrow = nrow(feature_df),
  ncol = length(sample_ids),
  dimnames = list(feature_df$feature_id, sample_ids)
)
segment_idx <- cbind(
  match(segment_df$feature_id, rownames(segment_mat)),
  match(segment_df$sampleid, colnames(segment_mat))
)
segment_mat[segment_idx] <- segment_df$seg_mean

crosswalk <- readxl::read_excel(snk$input$crosswalk) |>
  as.data.frame(stringsAsFactors = FALSE) |>
  sample_ids_from_crosswalk()
sample_map <- molecular_sample_map_from_crosswalk(crosswalk)

barcode_df <- segment_df |>
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
  assays = list(exprs = segment_mat),
  rowData = S4Vectors::DataFrame(feature_df),
  colData = S4Vectors::DataFrame(col_df)
)
S4Vectors::metadata(se) <- list(datatype = "segments", annotation = "segments")

fs::dir_create(fs::path_dir(snk$output$se))
saveRDS(se, snk$output$se)
write_matrix_csv(
  matrix = segment_mat,
  row_ids = feature_df$feature_id,
  row_id_name = "segment",
  output_path = snk$output$matrix
)
