snakemake@source("helpers.R")
snk <- parse_snakemake()

suppressPackageStartupMessages({
  library(dplyr)
  library(jsonlite)
})

model_df <- read_csv(snk$input$model)
sample_df <- read_csv(snk$input$sample)
treatment_df <- read_csv(snk$input$treatment)

base_url <- sub("/+$", "", snk$config$annotationdb$base_url)
compound_path <- snk$config$annotationdb$compound_many_path

fetch_json <- function(path, params) {
  query <- paste(
    paste0(names(params), "=", utils::URLencode(params, reserved = TRUE)),
    collapse = "&"
  )
  tryCatch(
    jsonlite::fromJSON(
      paste0(base_url, path, "?", query, "&format=json"),
      simplifyVector = FALSE
    ),
    error = function(e) list()
  )
}

first_record <- function(records) {
  if (is.data.frame(records)) {
    if (!nrow(records)) {
      return(list())
    }
    return(as.list(records[1, , drop = FALSE]))
  }
  if (!length(records)) {
    return(list())
  }
  records[[1]]
}

field <- function(record, names) {
  for (name in names) {
    if (!is.null(record[[name]]) && length(record[[name]]) > 0) {
      return(serialize_csv_value(record[[name]]))
    }
  }
  NA_character_
}

treatment_annotations <- lapply(seq_len(nrow(treatment_df)), function(idx) {
  treatment <- treatment_df[idx, , drop = FALSE]

  if (treatment$drug.id == snk$config$workbook$control_drug) {
    return(data.frame(
      drug.id = treatment$drug.id,
      annotationdb.match = FALSE,
      submitted.identifier = treatment$source_treatment_name,
      identifier.type = "source_treatment_name",
      compound.name = NA_character_,
      pubchem.cid = NA_character_,
      molecular.formula = NA_character_,
      smiles = NA_character_,
      inchikey = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  records <- fetch_json(
    compound_path,
    c(compound = treatment$source_treatment_name)
  )
  record <- first_record(records)

  data.frame(
    drug.id = treatment$drug.id,
    annotationdb.match = length(record) > 0,
    submitted.identifier = treatment$source_treatment_name,
    identifier.type = "source_treatment_name",
    compound.name = field(record, c("title", "name", "compound_name")),
    pubchem.cid = field(
      record,
      c("cid", "pubchem.cid", "pubchem_cid", "compound_id")
    ),
    molecular.formula = field(
      record,
      c(
        "molecular_formula",
        "molecular.formula",
        "molecularFormula",
        "MolecularFormula"
      )
    ),
    smiles = field(record, c("canonical_smiles", "isomeric_smiles", "smiles")),
    inchikey = field(record, c("inchikey", "inchi_key", "InChIKey")),
    stringsAsFactors = FALSE
  )
}) |>
  bind_rows()

treatment_annotated <- treatment_df |>
  left_join(treatment_annotations, by = "drug.id")

unmatched_treatments <- treatment_annotated |>
  filter(
    .data$drug.id != snk$config$workbook$control_drug,
    !.data$annotationdb.match
  ) |>
  transmute(
    drug.id = .data$drug.id,
    query_treatment_name = .data$source_treatment_name,
    standard.name = .data$standard.name,
    source.patient.count = .data$source.patient.count,
    source.model.count = .data$source.model.count
  )

write_csv(model_df, snk$output$model)
write_csv(sample_df, snk$output$sample)
write_csv(treatment_annotated, snk$output$treatment)
write_csv_export(unmatched_treatments, snk$output$unmatched_treatments)
