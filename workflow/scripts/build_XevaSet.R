snakemake@source("helpers.R")
snk <- parse_snakemake()

suppressPackageStartupMessages({
  library(Biobase)
  library(dplyr)
  library(fs)
  library(MultiAssayExperiment)
  library(SummarizedExperiment)
  library(Xeva)
})

model_df <- read_csv(snk$input$model) |>
  as.data.frame(stringsAsFactors = FALSE)
drug_df <- read_csv(snk$input$treatment) |>
  as.data.frame(stringsAsFactors = FALSE)
experiment_df <- read_csv(snk$input$experiment) |>
  as.data.frame(stringsAsFactors = FALSE)
mod_map <- read_csv(snk$input$modToBiobaseMap) |>
  as.data.frame(stringsAsFactors = FALSE)
exp_design <- readRDS(snk$input$expDesign)
mae <- readRDS(snk$input$multiAssayExperiment)

experiment_df <- experiment_df |>
  select(
    "model.id",
    "drug",
    "time",
    "length",
    "width",
    "volume",
    "volume.normal",
    "patient.id",
    "date",
    everything()
  )

parse_dose_value <- function(x) {
  x <- normalize_chr(x)
  dose <- rep(NA_real_, length(x))
  has_number <- !is.na(x) & grepl("[0-9]", x)
  dose[has_number] <- as.numeric(sub(
    "^.*?([0-9]+(?:\\.[0-9]+)?).*$",
    "\\1",
    x[has_number]
  ))
  dose
}

xeva_experiment_df <- experiment_df |>
  mutate(dose = parse_dose_value(.data$dose)) |>
  select(
    "model.id",
    "drug",
    "time",
    "volume",
    "width",
    "length",
    "date",
    "dose",
    "patient.id"
  )

set_response_quietly <- function(xeva, ...) {
  withCallingHandlers(
    Xeva::setResponse(xeva, ...),
    warning = function(warning_condition) {
      warning_text <- conditionMessage(warning_condition)
      expected_response_warning <- grepl(
        "insufficient data after time",
        warning_text
      ) ||
        grepl(
          "start volume zero, adding 1 to compute volume.normal",
          warning_text,
          fixed = TRUE
        ) ||
        grepl(
          "collapsing to unique 'x' values",
          warning_text,
          fixed = TRUE
        )

      if (expected_response_warning) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

rownames(model_df) <- model_df$model.id
rownames(drug_df) <- drug_df$drug.id
molecular_profiles <- lapply(
  MultiAssayExperiment::experiments(mae),
  se_to_expression_set
)

mcgill <- Xeva::createXevaSet(
  name = snk$config$dataset$name,
  model = model_df,
  drug = drug_df,
  experiment = xeva_experiment_df,
  expDesign = exp_design,
  molecularProfiles = molecular_profiles,
  modToBiobaseMap = mod_map
)

as_integer_or_na <- function(x) {
  suppressWarnings(as.integer(normalize_chr(x)))
}

as_logical_or_na <- function(x) {
  x <- tolower(normalize_chr(x))
  dplyr::case_when(
    x == "true" ~ TRUE,
    x == "false" ~ FALSE,
    TRUE ~ NA
  )
}

mcgill@drug$source.patient.count <- as_integer_or_na(
  mcgill@drug$source.patient.count
)
mcgill@drug$source.model.count <- as_integer_or_na(
  mcgill@drug$source.model.count
)
mcgill@drug$pubchem.cid <- as_integer_or_na(mcgill@drug$pubchem.cid)
mcgill@drug$annotationdb.match <- as_logical_or_na(
  mcgill@drug$annotationdb.match
)

for (metric in snk$config$response$metrics) {
  message("Computing response metric: ", metric)
  mcgill <- set_response_quietly(
    mcgill,
    res.measure = metric,
    min.time = snk$config$response$min_time,
    max.time = snk$config$response$max_time,
    treatment.only = FALSE,
    verbose = TRUE
  )
}

build_experiment_measurements <- function(experiment_list) {
  lapply(experiment_list, function(model_obj) {
    measurement_df <- as.data.frame(model_obj@data, stringsAsFactors = FALSE)
    measurement_df$model.id <- model_obj@model.id
    measurement_df$measurement_index <- seq_len(nrow(measurement_df))
    measurement_df |>
      select(
        "model.id",
        "measurement_index",
        everything()
      ) |>
      normalize_for_csv()
  }) |>
    data.table::rbindlist(fill = TRUE) |>
    as.data.frame(stringsAsFactors = FALSE)
}

export_component <- function(
  output_dir,
  relative_path,
  df,
  rownames_col = NULL
) {
  output_path <- fs::path(output_dir, relative_path)
  write_csv_export(df, output_path, rownames_col = rownames_col)
  data.frame(
    relative_path = relative_path,
    rows = nrow(as.data.frame(df)),
    columns = ncol(as.data.frame(df)),
    stringsAsFactors = FALSE
  )
}

export_xeva_csv <- function(xeva, output_dir) {
  if (dir_exists(output_dir)) {
    dir_delete(output_dir)
  }
  dir_create(output_dir)

  manifest <- list(
    model = export_component(output_dir, "model.csv", xeva@model),
    drug = export_component(output_dir, "drug.csv", xeva@drug),
    sensitivity_model = export_component(
      output_dir,
      "sensitivity_model.csv",
      xeva@sensitivity$model
    ),
    sensitivity_batch = export_component(
      output_dir,
      "sensitivity_batch.csv",
      xeva@sensitivity$batch
    ),
    mod_to_biobase_map = export_component(
      output_dir,
      "mod_to_biobase_map.csv",
      xeva@modToBiobaseMap
    )
  )

  manifest$experiment_measurements <- export_component(
    output_dir,
    "experiment_measurements.csv",
    build_experiment_measurements(xeva@experiment)
  )

  for (profile_name in names(xeva@molecularProfiles)) {
    profile <- xeva@molecularProfiles[[profile_name]]
    profile_manifest <- export_component(
      output_dir,
      fs::path("molecular_profiles", paste0(profile_name, ".csv")),
      as.data.frame(Biobase::exprs(profile), check.names = FALSE),
      rownames_col = "feature_rowname"
    )
    manifest[[paste0("molecular_profile_", profile_name)]] <- profile_manifest

    feature_manifest <- export_component(
      output_dir,
      fs::path("molecular_profile_features", paste0(profile_name, ".csv")),
      as.data.frame(Biobase::fData(profile), check.names = FALSE),
      rownames_col = "feature_rowname"
    )
    manifest[[paste0(
      "molecular_profile_features_",
      profile_name
    )]] <- feature_manifest
  }

  manifest_df <- bind_rows(manifest, .id = "component")
  write_csv_export(manifest_df, fs::path(output_dir, "file_manifest.csv"))
}

dir_create(path_dir(snk$output$xeva))
saveRDS(mcgill, snk$output$xeva)
export_xeva_csv(mcgill, snk$output$csvDir)
