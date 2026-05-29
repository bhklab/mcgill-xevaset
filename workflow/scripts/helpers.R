suppressPackageStartupMessages({
  library(fs)
})

parse_snakemake <- function() {
  if (!exists("snakemake")) {
    stop("This script is intended to be executed by Snakemake.", call. = FALSE)
  }

  if (length(snakemake@log) > 0) {
    log_file <- snakemake@log[[1]]
    dir_create(path_dir(log_file))
    sink(log_file, append = FALSE, split = TRUE)
    message_connection <- file(log_file, open = "at")
    sink(message_connection, append = TRUE, type = "message")
  }

  list(
    input = snakemake@input,
    output = snakemake@output,
    params = snakemake@params,
    threads = snakemake@threads,
    config = snakemake@config
  )
}


clean_names_simple <- function(x) {
  x <- gsub("%", "", x, fixed = TRUE)
  x <- gsub("([a-z0-9])([A-Z])", "\\1_\\2", x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("(^_+|_+$)", "", x)
  tolower(x)
}


normalize_chr <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "NULL", "NaN")] <- NA_character_
  x
}


first_non_missing <- function(x) {
  if (is.character(x) || is.factor(x)) {
    x <- normalize_chr(x)
  }
  idx <- which(!is.na(x))
  if (!length(idx)) {
    return(NA)
  }
  x[[idx[[1]]]]
}


collapse_unique <- function(x, sep = ";") {
  x <- normalize_chr(x)
  x <- sort(unique(x[!is.na(x)]))
  if (!length(x)) {
    return(NA_character_)
  }
  paste(x, collapse = sep)
}


excel_date <- function(x) {
  if (inherits(x, "Date")) {
    return(x)
  }
  if (inherits(x, "POSIXt")) {
    return(as.Date(x))
  }
  numeric_value <- suppressWarnings(as.numeric(x))
  if (!is.na(numeric_value)) {
    date_value <- suppressWarnings(as.Date(
      numeric_value,
      origin = "1899-12-30"
    ))
    if (!is.na(date_value) && date_value < as.Date("1902-01-01")) {
      return(as.Date(NA))
    }
    return(date_value)
  }
  x <- normalize_chr(x)
  if (is.na(x) || !grepl("^\\d{4}-\\d{2}-\\d{2}", x)) {
    return(as.Date(NA))
  }
  suppressWarnings(as.Date(x))
}


drug_code <- function(x) {
  x <- gsub("[^A-Za-z0-9]", "", as.character(x))
  x <- ifelse(nchar(x) >= 2, substr(x, 1, 2), x)
  x <- paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, 2)))
  x[!nzchar(x)] <- "NA"
  x
}


read_tsv <- function(path) {
  data.table::fread(
    path,
    sep = "\t",
    data.table = FALSE,
    check.names = FALSE
  )
}


read_csv <- function(path) {
  data.table::fread(
    path,
    sep = ",",
    data.table = FALSE,
    check.names = FALSE
  )
}


write_csv <- function(df, path, quote = "auto") {
  dir_create(path_dir(path))
  data.table::fwrite(df, path, sep = ",", na = "", quote = quote)
}


quote_csv_value <- function(x) {
  paste0('"', gsub('"', '""', x, fixed = TRUE), '"')
}


write_matrix_csv <- function(matrix, row_ids, row_id_name, output_path) {
  out_df <- data.frame(
    row_ids = row_ids,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  names(out_df)[[1]] <- row_id_name
  out_df <- cbind(out_df, as.data.frame(matrix, check.names = FALSE))
  write_csv(out_df, output_path)
}


sample_ids_from_crosswalk <- function(crosswalk_df) {
  names(crosswalk_df) <- clean_names_simple(names(crosswalk_df))
  get_col <- function(name) {
    if (name %in% names(crosswalk_df)) {
      return(normalize_chr(crosswalk_df[[name]]))
    }
    rep(NA_character_, nrow(crosswalk_df))
  }

  crosswalk_df |>
    dplyr::transmute(
      patient.id = get_col("pdx_id"),
      SA.ID = get_col("sa_id"),
      WGS.PDX = get_col("wgs_pdx"),
      WGS.Patient = get_col("wgs_patient"),
      WGS.Germline = get_col("wgs_germline"),
      RNAseq.PDX = get_col("rnaseq_pdx"),
      RNAseq.Patient = get_col("rnaseq_patient"),
      Site = get_col("site"),
      Histology = get_col("histology"),
      ER = get_col("er"),
      PR = get_col("pr"),
      HER2 = get_col("her2"),
      Grade = get_col("grade"),
      BRCA1.2 = get_col("brca1_2"),
      AIMS = get_col("aims")
    ) |>
    dplyr::filter(!is.na(.data$patient.id))
}


molecular_sample_map_from_crosswalk <- function(crosswalk_df) {
  dplyr::bind_rows(
    crosswalk_df |>
      dplyr::transmute(
        sampleid = .data$WGS.PDX,
        patient.id = .data$patient.id,
        Source = "PDX"
      ),
    crosswalk_df |>
      dplyr::transmute(
        sampleid = .data$WGS.Patient,
        patient.id = .data$patient.id,
        Source = "Patient"
      ),
    crosswalk_df |>
      dplyr::transmute(
        sampleid = .data$WGS.Germline,
        patient.id = .data$patient.id,
        Source = "Germline"
      ),
    crosswalk_df |>
      dplyr::transmute(
        sampleid = .data$RNAseq.PDX,
        patient.id = .data$patient.id,
        Source = "PDX"
      ),
    crosswalk_df |>
      dplyr::transmute(
        sampleid = .data$RNAseq.Patient,
        patient.id = .data$patient.id,
        Source = "Patient"
      )
  ) |>
    tidyr::separate_rows(sampleid, sep = "\\s*,\\s*") |>
    dplyr::mutate(sampleid = normalize_chr(.data$sampleid)) |>
    dplyr::filter(
      !is.na(.data$sampleid),
      !tolower(.data$sampleid) %in% c("no sample", "none")
    ) |>
    dplyr::distinct(.data$sampleid, .keep_all = TRUE)
}


sample_id_from_profile <- function(x) {
  x <- basename(as.character(x))
  x <- sub("^WG_", "", x)
  x <- sub("\\..*$", "", x)
  x <- sub("^.*_(SA[[:alnum:]_]+)$", "\\1", x)
  normalize_chr(x)
}


expression_set_to_se <- function(eset, datatype) {
  assay_mat <- Biobase::exprs(eset)

  row_df <- as.data.frame(Biobase::fData(eset), stringsAsFactors = FALSE)
  if (!nrow(row_df)) {
    row_df <- data.frame(feature_id = rownames(assay_mat))
  }
  rownames(row_df) <- rownames(assay_mat)

  col_df <- as.data.frame(Biobase::pData(eset), stringsAsFactors = FALSE)
  if (!nrow(col_df)) {
    col_df <- data.frame(sampleid = colnames(assay_mat))
  }
  rownames(col_df) <- colnames(assay_mat)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(exprs = as.matrix(assay_mat)),
    rowData = S4Vectors::DataFrame(row_df),
    colData = S4Vectors::DataFrame(col_df)
  )
  S4Vectors::metadata(se) <- list(datatype = datatype, annotation = datatype)
  se
}


se_to_expression_set <- function(se) {
  assay_mat <- SummarizedExperiment::assay(se)
  row_df <- as.data.frame(
    SummarizedExperiment::rowData(se),
    stringsAsFactors = FALSE
  )
  col_df <- as.data.frame(
    SummarizedExperiment::colData(se),
    stringsAsFactors = FALSE
  )

  if (!nrow(row_df)) {
    row_df <- data.frame(feature_id = rownames(assay_mat))
  }
  if (!nrow(col_df)) {
    col_df <- data.frame(sampleid = colnames(assay_mat))
  }
  rownames(row_df) <- rownames(assay_mat)
  rownames(col_df) <- colnames(assay_mat)

  Biobase::ExpressionSet(
    assayData = as.matrix(assay_mat),
    phenoData = Biobase::AnnotatedDataFrame(col_df),
    featureData = Biobase::AnnotatedDataFrame(row_df)
  )
}


serialize_csv_value <- function(x) {
  if (is.null(x) || !length(x)) {
    return(NA_character_)
  }
  if (inherits(x, "sessionInfo")) {
    x <- paste(capture.output(print(x)), collapse = " | ")
  }
  if (inherits(x, "POSIXt")) {
    return(format(x, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
  }
  if (inherits(x, "Date")) {
    return(as.character(x))
  }
  if (is.factor(x)) {
    x <- as.character(x)
  }
  if (is.data.frame(x)) {
    x <- lapply(x, function(col) if (is.factor(col)) as.character(col) else col)
  }

  value <- if (is.list(x) || length(x) > 1) {
    jsonlite::toJSON(
      x,
      auto_unbox = TRUE,
      null = "null",
      dataframe = "rows",
      na = "null"
    )
  } else {
    as.character(x[[1]])
  }

  gsub("[\r\n]+", " ", value)
}


normalize_for_csv <- function(df) {
  df[] <- lapply(
    df,
    function(col) {
      if (is.factor(col)) {
        return(normalize_chr(as.character(col)))
      }
      if (inherits(col, "POSIXt")) {
        return(format(col, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
      }
      if (inherits(col, "Date")) {
        return(as.character(col))
      }
      if (is.logical(col)) {
        return(as.character(col))
      }
      if (is.list(col)) {
        return(normalize_chr(vapply(col, serialize_csv_value, character(1))))
      }
      if (is.character(col)) {
        return(normalize_chr(col))
      }
      col
    }
  )
  df
}


write_csv_export <- function(df, output_path, rownames_col = NULL) {
  out_df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)

  if (!is.null(rownames_col)) {
    out_df <- data.frame(
      stats::setNames(
        data.frame(rownames = rownames(out_df), stringsAsFactors = FALSE),
        rownames_col
      ),
      out_df,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  out_df <- normalize_for_csv(out_df)
  encoded_df <- as.data.frame(
    lapply(out_df, function(col) {
      if (is.numeric(col) || is.integer(col)) {
        return(ifelse(is.na(col), "", as.character(col)))
      }

      col <- as.character(col)
      ifelse(is.na(col), "", quote_csv_value(col))
    }),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  header <- paste(quote_csv_value(names(encoded_df)), collapse = ",")
  dir_create(path_dir(output_path))
  writeLines(header, output_path, useBytes = TRUE)
  data.table::fwrite(
    encoded_df,
    output_path,
    sep = ",",
    quote = FALSE,
    na = "",
    col.names = FALSE,
    append = TRUE
  )
}
