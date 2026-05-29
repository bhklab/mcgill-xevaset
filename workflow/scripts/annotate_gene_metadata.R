snakemake@source("helpers.R")
snk <- parse_snakemake()

suppressPackageStartupMessages({
  library(data.table)
  library(httr)
  library(jsonlite)
  library(S4Vectors)
  library(SummarizedExperiment)
})

gencode_cfg <- snk$config$annotationgx$gencode
species <- gencode_cfg$species
id_batch_size <- as.integer(gencode_cfg$id_batch_size)
request_timeout_seconds <- as.integer(gencode_cfg$request_timeout_seconds)
symbol_fallback <- isTRUE(gencode_cfg$symbol_fallback)

get_current_release_info <- function(species) {
  response <- httr::GET(
    "https://rest.ensembl.org/info/data",
    query = list("content-type" = "application/json"),
    httr::accept_json(),
    httr::timeout(request_timeout_seconds)
  )
  httr::stop_for_status(response)
  payload <- jsonlite::fromJSON(
    httr::content(response, as = "text", encoding = "UTF-8")
  )

  data.table(
    species = species,
    ensembl_release = as.integer(payload$releases[[1]]),
    queried_at = as.character(Sys.time())
  )
}

release_info <- get_current_release_info(species)
ensembl_release <- release_info$ensembl_release[[1]]

empty_gene_table <- function() {
  data.table(
    query_type = character(),
    query_value = character(),
    gene_id = character(),
    gene_version = numeric(),
    symbol = character(),
    description = character(),
    biotype = character(),
    seq_region_name = character(),
    start = integer(),
    end = integer(),
    strand = integer(),
    strand_char = character(),
    assembly_name = character(),
    canonical_transcript = character(),
    object_type = character(),
    species = character(),
    release = integer()
  )
}

strand_to_char <- function(x) {
  dplyr::case_when(
    as.integer(x) == 1L ~ "+",
    as.integer(x) == -1L ~ "-",
    TRUE ~ "*"
  )
}

field <- function(record, name, default = NA) {
  value <- record[[name]]
  if (is.null(value) || !length(value)) {
    return(default)
  }
  value[[1]]
}

record_to_row <- function(record, query_value, query_type) {
  if (is.null(record) || !length(record)) {
    return(NULL)
  }

  strand <- field(record, "strand", NA_integer_)
  data.table(
    query_type = query_type,
    query_value = query_value,
    gene_id = field(record, "id", NA_character_),
    gene_version = field(record, "version", NA_real_),
    symbol = field(record, "display_name", NA_character_),
    description = field(record, "description", NA_character_),
    biotype = field(record, "biotype", NA_character_),
    seq_region_name = field(record, "seq_region_name", NA_character_),
    start = field(record, "start", NA_integer_),
    end = field(record, "end", NA_integer_),
    strand = strand,
    strand_char = strand_to_char(strand),
    assembly_name = field(record, "assembly_name", NA_character_),
    canonical_transcript = field(record, "canonical_transcript", NA_character_),
    object_type = field(record, "object_type", NA_character_),
    species = field(record, "species", NA_character_),
    release = ensembl_release
  )
}

as_gene_table <- function(result) {
  if (is.null(result) || !length(result$genes) || !nrow(result$genes)) {
    return(empty_gene_table())
  }

  out <- as.data.table(result$genes)
  missing_cols <- setdiff(names(empty_gene_table()), names(out))
  for (col_name in missing_cols) {
    out[[col_name]] <- NA
  }
  out[, names(empty_gene_table()), with = FALSE]
}

lookup_ids <- function(ids) {
  ids <- sort(unique(na.omit(ids)))
  if (!length(ids)) {
    return(empty_gene_table())
  }

  id_batches <- split(ids, ceiling(seq_along(ids) / id_batch_size))
  rows <- lapply(seq_along(id_batches), function(batch_idx) {
    batch_ids <- id_batches[[batch_idx]]
    message(
      "GENCODE ID batch ",
      batch_idx,
      "/",
      length(id_batches),
      " (",
      length(batch_ids),
      " IDs)"
    )
    response <- httr::POST(
      "https://rest.ensembl.org/lookup/id",
      query = list("content-type" = "application/json"),
      body = list(ids = batch_ids, expand = 0),
      encode = "json",
      httr::accept_json(),
      httr::timeout(request_timeout_seconds)
    )
    httr::stop_for_status(response)
    payload <- jsonlite::fromJSON(
      httr::content(response, as = "text", encoding = "UTF-8"),
      simplifyVector = FALSE
    )

    data.table::rbindlist(
      lapply(names(payload), function(id) {
        record_to_row(payload[[id]], id, "id")
      }),
      fill = TRUE
    )
  })

  data.table::rbindlist(rows, fill = TRUE)
}

lookup_symbols <- function(symbols) {
  symbols <- sort(unique(na.omit(symbols)))
  if (!length(symbols)) {
    return(empty_gene_table())
  }

  message(
    "Looking up ",
    length(symbols),
    " remaining gene symbols with current GENCODE"
  )
  symbol_batches <- split(symbols, ceiling(seq_along(symbols) / id_batch_size))
  rows <- lapply(seq_along(symbol_batches), function(batch_idx) {
    batch_symbols <- symbol_batches[[batch_idx]]
    message(
      "GENCODE symbol batch ",
      batch_idx,
      "/",
      length(symbol_batches),
      " (",
      length(batch_symbols),
      " symbols)"
    )
    response <- httr::POST(
      paste0("https://rest.ensembl.org/lookup/symbol/", species),
      query = list("content-type" = "application/json"),
      body = list(symbols = batch_symbols, expand = 0),
      encode = "json",
      httr::accept_json(),
      httr::timeout(request_timeout_seconds)
    )
    httr::stop_for_status(response)
    payload <- jsonlite::fromJSON(
      httr::content(response, as = "text", encoding = "UTF-8"),
      simplifyVector = FALSE
    )

    data.table::rbindlist(
      lapply(names(payload), function(symbol) {
        record_to_row(payload[[symbol]], symbol, "symbol")
      }),
      fill = TRUE
    )
  })

  data.table::rbindlist(rows, fill = TRUE)
}

id_from_feature <- function(x) {
  x <- normalize_chr(x)
  x <- sub("\\.[0-9]+$", "", x)
  x[!grepl("^ENS[A-Z]*G[0-9]+$", x)] <- NA_character_
  x
}

gene_annotation_cols <- c(
  "gencode.query.type" = "query_type",
  "gencode.query.value" = "query_value",
  "gencode.mapped" = "gencode.mapped",
  "gencode.gene.id" = "gene_id",
  "gencode.gene.version" = "gene_version",
  "gencode.symbol" = "symbol",
  "gencode.description" = "description",
  "gencode.biotype" = "biotype",
  "gencode.seq.region" = "seq_region_name",
  "gencode.start" = "start",
  "gencode.end" = "end",
  "gencode.strand" = "strand_char",
  "gencode.assembly" = "assembly_name",
  "gencode.canonical.transcript" = "canonical_transcript",
  "gencode.object.type" = "object_type",
  "gencode.species" = "species",
  "gencode.release" = "release"
)

add_gencode_columns <- function(row_df, annotation_df, prefix = "gencode") {
  col_names <- names(gene_annotation_cols)
  if (!identical(prefix, "gencode")) {
    col_names <- sub("^gencode", paste0("gencode.", prefix), col_names)
  }

  for (idx in seq_along(gene_annotation_cols)) {
    row_df[[col_names[[idx]]]] <- annotation_df[[gene_annotation_cols[[idx]]]]
  }
  row_df
}

annotate_features <- function(feature_df, id_hits, symbol_hits) {
  query_df <- as.data.table(feature_df)

  id_part <- id_hits[, names(empty_gene_table()), with = FALSE]
  setnames(id_part, names(id_part), paste0("id.", names(id_part)))
  annotated <- merge(
    query_df,
    id_part,
    by.x = "id_query",
    by.y = "id.query_value",
    all.x = TRUE,
    sort = FALSE
  )

  symbol_part <- symbol_hits[, names(empty_gene_table()), with = FALSE]
  setnames(
    symbol_part,
    names(symbol_part),
    paste0("symbol.", names(symbol_part))
  )
  annotated <- merge(
    annotated,
    symbol_part,
    by.x = "symbol_query",
    by.y = "symbol.query_value",
    all.x = TRUE,
    sort = FALSE
  )

  setorder(annotated, row_index)
  id_mapped <- !is.na(annotated$id.gene_id) & nzchar(annotated$id.gene_id)
  symbol_mapped <- !is.na(annotated$symbol.gene_id) &
    nzchar(annotated$symbol.gene_id)

  out <- data.table(
    row_index = annotated$row_index,
    feature_rowname = annotated$feature_rowname,
    feature_column = annotated$feature_column,
    feature_value = annotated$feature_value,
    query_type = ifelse(
      id_mapped,
      "id",
      ifelse(symbol_mapped, "symbol", annotated$primary_query_type)
    ),
    query_value = ifelse(
      id_mapped,
      annotated$id_query,
      ifelse(
        symbol_mapped,
        annotated$symbol_query,
        annotated$primary_query_value
      )
    ),
    gencode.mapped = id_mapped | symbol_mapped
  )

  for (col_name in setdiff(
    names(empty_gene_table()),
    c("query_type", "query_value")
  )) {
    out[[col_name]] <- ifelse(
      id_mapped,
      annotated[[paste0("id.", col_name)]],
      annotated[[paste0("symbol.", col_name)]]
    )
  }

  out
}

build_feature_queries <- function(
  row_df,
  feature_col,
  profile_name,
  id_col = NULL,
  symbol_col = NULL
) {
  feature_value <- normalize_chr(row_df[[feature_col]])
  id_query <- if (!is.null(id_col)) {
    id_from_feature(row_df[[id_col]])
  } else {
    id_from_feature(feature_value)
  }
  symbol_query <- if (!is.null(symbol_col)) {
    normalize_chr(row_df[[symbol_col]])
  } else {
    feature_value
  }
  primary_type <- ifelse(!is.na(id_query), "id", "symbol")
  primary_value <- ifelse(!is.na(id_query), id_query, symbol_query)

  data.table(
    profile = profile_name,
    row_index = seq_len(nrow(row_df)),
    feature_rowname = rownames(row_df),
    feature_column = feature_col,
    feature_value = feature_value,
    id_query = id_query,
    symbol_query = symbol_query,
    primary_query_type = primary_type,
    primary_query_value = primary_value
  )
}

ses <- list(
  RNAseq = readRDS(snk$input$rnaseq),
  mutation = readRDS(snk$input$mutation),
  CNV = readRDS(snk$input$cnv),
  fusion = readRDS(snk$input$fusion),
  segments = readRDS(snk$input$segments)
)

rna_row_df <- as.data.frame(rowData(ses$RNAseq), stringsAsFactors = FALSE)
rna_row_df$feature_id <- rownames(rna_row_df)
rna_queries <- build_feature_queries(
  rna_row_df,
  feature_col = "feature_id",
  profile_name = "RNAseq",
  id_col = "feature_id",
  symbol_col = "Symbol"
)

gene_profiles <- list(
  mutation = list(se = ses$mutation, feature_col = "gene"),
  CNV = list(se = ses$CNV, feature_col = "gene")
)
gene_queries <- lapply(names(gene_profiles), function(profile_name) {
  row_df <- as.data.frame(
    rowData(gene_profiles[[profile_name]]$se),
    stringsAsFactors = FALSE
  )
  build_feature_queries(
    row_df,
    feature_col = gene_profiles[[profile_name]]$feature_col,
    profile_name = profile_name
  )
})

fusion_row_df <- as.data.frame(rowData(ses$fusion), stringsAsFactors = FALSE)
fusion_left_queries <- build_feature_queries(
  fusion_row_df,
  feature_col = "left_gene",
  profile_name = "fusion"
)
fusion_right_queries <- build_feature_queries(
  fusion_row_df,
  feature_col = "right_gene",
  profile_name = "fusion"
)

all_query_df <- rbindlist(
  c(
    list(rna_queries),
    gene_queries,
    list(fusion_left_queries, fusion_right_queries)
  ),
  fill = TRUE
)
id_queries <- sort(unique(na.omit(all_query_df$id_query)))

message(
  "Looking up ",
  length(id_queries),
  " Ensembl gene IDs with current GENCODE"
)
id_hits <- lookup_ids(id_queries)

symbol_hits_from_ids <- id_hits[
  !is.na(symbol) & nzchar(symbol),
]
symbol_hits_from_ids[, query_type := "symbol"]
symbol_hits_from_ids[, query_value := symbol]

symbol_queries <- sort(unique(na.omit(all_query_df$symbol_query)))
missing_symbols <- setdiff(symbol_queries, symbol_hits_from_ids$query_value)
symbol_hits_fallback <- if (symbol_fallback) {
  lookup_symbols(missing_symbols)
} else {
  empty_gene_table()
}
symbol_hits <- rbindlist(
  list(symbol_hits_fallback, symbol_hits_from_ids),
  fill = TRUE
)
symbol_hits <- symbol_hits[!duplicated(paste(query_type, query_value))]

unmapped_rows <- list()
save_release_metadata <- function(se) {
  se_metadata <- metadata(se)
  se_metadata$gencode <- as.list(release_info[1, ])
  metadata(se) <- se_metadata
  se
}

annotate_regular_profile <- function(profile_name, se, query_df, output_path) {
  source_row_df <- as.data.frame(rowData(se), stringsAsFactors = FALSE)
  annotation_df <- annotate_features(query_df, id_hits, symbol_hits)

  if (identical(profile_name, "RNAseq")) {
    row_df <- data.frame(
      feature_id = rownames(source_row_df),
      stringsAsFactors = FALSE,
      row.names = rownames(source_row_df)
    )
  } else {
    row_df <- source_row_df
  }

  row_df <- add_gencode_columns(row_df, annotation_df)
  rowData(se) <- S4Vectors::DataFrame(row_df)
  se <- save_release_metadata(se)

  unmapped <- annotation_df[gencode.mapped == FALSE]
  if (nrow(unmapped)) {
    source_cols <- as.data.table(source_row_df[
      unmapped$row_index,
      ,
      drop = FALSE
    ])
    setnames(
      source_cols,
      names(source_cols),
      paste0("source.", names(source_cols))
    )
    unmapped_rows[[profile_name]] <<- cbind(
      data.table(profile = profile_name),
      unmapped[,
        .(
          feature_rowname,
          feature_column,
          feature_value,
          query_type,
          query_value
        )
      ],
      source_cols
    )
  }

  fs::dir_create(fs::path_dir(output_path))
  saveRDS(se, output_path)
}

annotate_regular_profile("RNAseq", ses$RNAseq, rna_queries, snk$output$rnaseq)
for (profile_name in names(gene_profiles)) {
  annotate_regular_profile(
    profile_name,
    gene_profiles[[profile_name]]$se,
    gene_queries[[which(names(gene_profiles) == profile_name)]],
    snk$output[[tolower(profile_name)]]
  )
}

left_annotation <- annotate_features(fusion_left_queries, id_hits, symbol_hits)
right_annotation <- annotate_features(
  fusion_right_queries,
  id_hits,
  symbol_hits
)
fusion_out_df <- fusion_row_df
fusion_out_df <- add_gencode_columns(
  fusion_out_df,
  left_annotation,
  prefix = "left"
)
fusion_out_df <- add_gencode_columns(
  fusion_out_df,
  right_annotation,
  prefix = "right"
)
rowData(ses$fusion) <- S4Vectors::DataFrame(fusion_out_df)
ses$fusion <- save_release_metadata(ses$fusion)
fs::dir_create(fs::path_dir(snk$output$fusion))
saveRDS(ses$fusion, snk$output$fusion)

fusion_unmapped <- rbindlist(
  list(
    left_annotation[gencode.mapped == FALSE],
    right_annotation[gencode.mapped == FALSE]
  ),
  fill = TRUE
)
if (nrow(fusion_unmapped)) {
  source_cols <- as.data.table(
    fusion_row_df[fusion_unmapped$row_index, , drop = FALSE]
  )
  setnames(
    source_cols,
    names(source_cols),
    paste0("source.", names(source_cols))
  )
  unmapped_rows$fusion <- cbind(
    data.table(profile = "fusion"),
    fusion_unmapped[,
      .(
        feature_rowname,
        feature_column,
        feature_value,
        query_type,
        query_value
      )
    ],
    source_cols
  )
}

ses$segments <- save_release_metadata(ses$segments)
fs::dir_create(fs::path_dir(snk$output$segments))
saveRDS(ses$segments, snk$output$segments)

unmapped_genes <- if (length(unmapped_rows)) {
  rbindlist(unmapped_rows, fill = TRUE)
} else {
  data.table(
    profile = character(),
    feature_rowname = character(),
    feature_column = character(),
    feature_value = character(),
    query_type = character(),
    query_value = character()
  )
}
write_csv_export(unmapped_genes, snk$output$unmapped)
