snakemake@source("helpers.R")
snk <- parse_snakemake()

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(S4Vectors)
  library(SummarizedExperiment)
})

datatype <- snk$params$datatype
cnv <- read_tsv(snk$input$cnv)
genes <- normalize_chr(cnv$Hugo_Symbol)
keep <- !is.na(genes)
cnv <- cnv[keep, , drop = FALSE]
genes <- genes[keep]

barcodes <- setdiff(names(cnv), "Hugo_Symbol")
sample_ids <- sample_id_from_profile(barcodes)
assay_df <- cnv[, barcodes, drop = FALSE]
assay_df[] <- lapply(assay_df, function(col) as.numeric(as.character(col)))
assay_mat <- as.matrix(assay_df)
storage.mode(assay_mat) <- "double"
rownames(assay_mat) <- make.unique(genes, sep = "__")
colnames(assay_mat) <- sample_ids

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

row_df <- data.frame(
  gene = genes,
  stringsAsFactors = FALSE,
  row.names = rownames(assay_mat)
)

se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(exprs = assay_mat),
  rowData = S4Vectors::DataFrame(row_df),
  colData = S4Vectors::DataFrame(col_df)
)
S4Vectors::metadata(se) <- list(datatype = datatype, annotation = datatype)

fs::dir_create(fs::path_dir(snk$output$se))
saveRDS(se, snk$output$se)
write_matrix_csv(
  matrix = assay_mat,
  row_ids = row_df$gene,
  row_id_name = "gene",
  output_path = snk$output$matrix
)
