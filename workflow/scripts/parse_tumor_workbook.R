snakemake@source("helpers.R")
snk <- parse_snakemake()

suppressPackageStartupMessages({
  library(dplyr)
  library(fs)
  library(readxl)
})

workbook <- snk$input$workbook
ignored_sheets <- snk$config$workbook$ignored_sheets
control_label <- snk$config$workbook$source_control_label
control_drug <- snk$config$workbook$control_drug

sheets <- setdiff(readxl::excel_sheets(workbook), ignored_sheets)

cell_chr <- function(df, row, col) {
  if (row > nrow(df) || col > ncol(df)) {
    return(NA_character_)
  }
  normalize_chr(df[[col]][[row]])
}

cell_num <- function(df, row, col) {
  suppressWarnings(as.numeric(cell_chr(df, row, col)))
}

parse_block <- function(df, sheet, start_row, end_row) {
  block <- df[start_row:end_row, , drop = FALSE]
  col2 <- normalize_chr(block[[2]])
  mouse_header_rows <- which(tolower(col2) == "mouse #")
  if (!length(mouse_header_rows)) {
    return(data.frame())
  }
  mouse_header_row <- mouse_header_rows[[1]]
  data_row <- mouse_header_row + 1L

  arm <- cell_chr(block, mouse_header_row, 1)
  source_treatment_name <- if (identical(arm, control_label)) {
    control_label
  } else {
    cell_chr(block, data_row, 1)
  }
  drug <- if (identical(arm, control_label)) {
    control_drug
  } else {
    source_treatment_name
  }

  dose <- cell_chr(block, data_row + 1L, 2)
  schedule <- cell_chr(block, data_row + 2L, 2)
  mouse_id <- cell_chr(block, data_row, 2)
  source_data_row <- start_row + data_row - 1L
  if (is.na(mouse_id)) {
    mouse_id <- paste0("missing_mouse_", source_data_row)
  }

  header_values <- tolower(normalize_chr(as.character(block[
    mouse_header_row,
  ])))
  volume_cols <- which(header_values == "volume")
  if (!length(volume_cols)) {
    return(data.frame())
  }
  length_cols <- volume_cols - 2L
  width_cols <- volume_cols - 1L

  patient_num <- gsub("^GCRC", "", sheet)
  model_id <- sprintf(
    "m%s.%02d.%s%s",
    sub("\\.0$", "", mouse_id),
    source_data_row,
    patient_num,
    drug_code(source_treatment_name)
  )

  out <- lapply(seq_along(volume_cols), function(idx) {
    volume_col <- volume_cols[[idx]]
    length_col <- length_cols[[idx]]
    width_col <- width_cols[[idx]]

    length_value <- cell_num(block, data_row, length_col)
    width_value <- cell_num(block, data_row, width_col)
    source_volume <- cell_num(block, data_row, volume_col)
    if (!any(is.na(c(width_value, length_value)))) {
      dims <- sort(c(width_value, length_value))
      width_value <- dims[[1]]
      length_value <- dims[[2]]
    }
    calculated_volume <- if (!any(is.na(c(width_value, length_value)))) {
      length_value * width_value^2 / 2
    } else {
      source_volume
    }

    data.frame(
      model.id = model_id,
      patient.id = sheet,
      sheet = sheet,
      source_block_start_row = start_row,
      source_data_row = source_data_row,
      mouse.id = mouse_id,
      source_arm = arm,
      source_treatment_name = source_treatment_name,
      drug = drug,
      dose = dose,
      schedule = schedule,
      time = cell_num(block, 2, volume_col),
      date = excel_date(cell_chr(block, 1, volume_col)),
      length = length_value,
      width = width_value,
      source_volume = source_volume,
      volume = calculated_volume,
      stringsAsFactors = FALSE
    )
  }) |>
    bind_rows() |>
    filter(!is.na(.data$time), .data$time >= 0, !is.na(.data$volume))

  if (!nrow(out)) {
    return(out)
  }

  out$volume.normal <- (out$volume - out$volume[[1]]) / out$volume[[1]]
  out$measurement_index <- seq_len(nrow(out))
  out
}

measurement_tables <- list()
for (sheet in sheets) {
  message("Parsing sheet: ", sheet)
  df <- readxl::read_excel(
    workbook,
    sheet = sheet,
    col_names = FALSE,
    col_types = "text",
    .name_repair = "minimal",
    na = c("", "NA")
  )
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)

  start_rows <- which(tolower(normalize_chr(df[[2]])) == "implantation date")
  end_rows <- c(start_rows[-1] - 1L, nrow(df))

  sheet_tables <- lapply(
    seq_along(start_rows),
    function(idx) parse_block(df, sheet, start_rows[[idx]], end_rows[[idx]])
  )
  measurement_tables[[sheet]] <- bind_rows(sheet_tables)
}

measurements <- bind_rows(measurement_tables) |>
  arrange(.data$patient.id, .data$source_data_row, .data$time)

block_summary <- measurements |>
  group_by(
    .data$patient.id,
    .data$model.id,
    .data$source_block_start_row,
    .data$source_data_row,
    .data$mouse.id,
    .data$source_arm,
    .data$source_treatment_name,
    .data$drug,
    .data$dose,
    .data$schedule
  ) |>
  summarise(
    n_measurements = n(),
    min_time = min(.data$time, na.rm = TRUE),
    max_time = max(.data$time, na.rm = TRUE),
    baseline_volume = first(.data$volume),
    final_volume = last(.data$volume),
    .groups = "drop"
  ) |>
  arrange(.data$patient.id, .data$source_data_row)

write_csv(measurements, snk$output$measurements)
write_csv(block_summary, snk$output$block_summary)
