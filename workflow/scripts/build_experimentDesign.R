snakemake@source("helpers.R")
snk <- parse_snakemake()

suppressPackageStartupMessages({
  library(dplyr)
  library(fs)
})

model_df <- read_csv(snk$input$model)
control_drug <- snk$config$workbook$control_drug
drugs <- setdiff(sort(unique(model_df$drug)), control_drug)

exp_design <- list()
for (patient_id in sort(unique(model_df$patient.id))) {
  for (drug_id in drugs) {
    treatment <- model_df$model.id[
      model_df$patient.id == patient_id & model_df$drug == drug_id
    ]
    control <- model_df$model.id[
      model_df$patient.id == patient_id & model_df$drug == control_drug
    ]

    if (length(treatment) > 0 && length(control) > 0) {
      batch_name <- paste(patient_id, drug_id, sep = ".")
      exp_design[[batch_name]] <- list(
        batch.name = batch_name,
        treatment = treatment,
        control = control
      )
    }
  }
}

dir_create(path_dir(snk$output$expDesign))
saveRDS(exp_design, snk$output$expDesign)
