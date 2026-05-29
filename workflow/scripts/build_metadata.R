snakemake@source("helpers.R")
snk <- parse_snakemake()

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
})

measurements <- read_csv(snk$input$measurements)

clinical <- readxl::read_excel(snk$input$clinical) |>
  as.data.frame(stringsAsFactors = FALSE)
names(clinical) <- clean_names_simple(names(clinical))

clinical_get <- function(name) {
  if (name %in% names(clinical)) {
    return(normalize_chr(clinical[[name]]))
  }
  rep(NA_character_, nrow(clinical))
}

clinical <- clinical |>
  mutate(
    patient.id = clinical_get("gcrc_id2"),
    Sample.ID = clinical_get("sample_id"),
    Source = clinical_get("source"),
    SA.ID = clinical_get("sa_id"),
    Site.Breast = clinical_get("site_breast"),
    Metastatic.Sample = clinical_get("metastatic_sample"),
    Histology.clinical = clinical_get("histology"),
    ER.clinical = clinical_get("er"),
    PR.clinical = clinical_get("pr"),
    HER2.clinical = clinical_get("her2"),
    AIMS.clinical = clinical_get("aims"),
    BRCA1.2_Mutation = clinical_get("brca1_2_mutation"),
    Systemic.Treatment = clinical_get("systemic_treatment"),
    Patient.Metastasis = clinical_get("patient_metastasis"),
    Vital.Status = clinical_get("vital_status")
  ) |>
  filter(!is.na(.data$patient.id)) |>
  arrange(.data$patient.id, desc(.data$Source == "PDX"))

clinical_one <- clinical |>
  group_by(.data$patient.id) |>
  slice(1L) |>
  ungroup() |>
  select(
    "patient.id",
    "Sample.ID",
    "Source",
    "SA.ID",
    "Site.Breast",
    "Metastatic.Sample",
    "Histology.clinical",
    "ER.clinical",
    "PR.clinical",
    "HER2.clinical",
    "AIMS.clinical",
    "BRCA1.2_Mutation",
    "Systemic.Treatment",
    "Patient.Metastasis",
    "Vital.Status"
  )

crosswalk <- readxl::read_excel(snk$input$crosswalk) |>
  as.data.frame(stringsAsFactors = FALSE) |>
  sample_ids_from_crosswalk()

patient_ids <- measurements |>
  distinct(.data$patient.id) |>
  arrange(.data$patient.id)

sample_df <- patient_ids |>
  left_join(crosswalk, by = "patient.id") |>
  left_join(
    clinical_one,
    by = "patient.id",
    suffix = c(".crosswalk", ".clinical")
  ) |>
  mutate(
    sampleid = .data$patient.id,
    tissue.id = snk$config$metadata$tissue_id,
    tissue = snk$config$metadata$tissue,
    tissue.name = snk$config$metadata$tissue_name,
    oncotree.code = snk$config$metadata$oncotree_code,
    oncotree.name = snk$config$metadata$oncotree_name,
    oncotree.mainType = snk$config$metadata$oncotree_main_type,
    SA.ID = coalesce(.data$SA.ID.crosswalk, .data$SA.ID.clinical),
    Site = coalesce(.data$Site, .data$Site.Breast),
    Histology = coalesce(.data$Histology, .data$Histology.clinical),
    ER = coalesce(.data$ER, .data$ER.clinical),
    PR = coalesce(.data$PR, .data$PR.clinical),
    HER2 = coalesce(.data$HER2, .data$HER2.clinical),
    AIMS = coalesce(.data$AIMS, .data$AIMS.clinical)
  ) |>
  select(
    "sampleid",
    "patient.id",
    "tissue.id",
    "tissue",
    "tissue.name",
    "oncotree.code",
    "oncotree.name",
    "oncotree.mainType",
    everything()
  )

treatment_df <- measurements |>
  distinct(
    .data$drug,
    .data$source_treatment_name,
    .data$patient.id,
    .data$model.id
  ) |>
  group_by(.data$drug, .data$source_treatment_name) |>
  summarise(
    standard.name = first(.data$drug),
    treatment.type = if_else(
      first(.data$drug) == snk$config$workbook$control_drug,
      "control",
      "drug"
    ),
    source.patient.count = n_distinct(.data$patient.id),
    source.model.count = n_distinct(.data$model.id),
    .groups = "drop"
  ) |>
  transmute(
    drug.id = .data$drug,
    standard.name = if_else(
      .data$drug == snk$config$workbook$control_drug,
      snk$config$workbook$control_drug,
      .data$standard.name
    ),
    source_treatment_name = .data$source_treatment_name,
    treatment.type = .data$treatment.type,
    source.patient.count = .data$source.patient.count,
    source.model.count = .data$source.model.count
  ) |>
  arrange(.data$treatment.type, .data$drug.id)

model_df <- measurements |>
  distinct(
    .data$model.id,
    .data$patient.id,
    .data$drug,
    .data$source_treatment_name,
    .data$dose,
    .data$schedule,
    .data$source_data_row,
    .data$mouse.id
  ) |>
  left_join(sample_df, by = "patient.id") |>
  transmute(
    model.id = .data$model.id,
    tissue.id = .data$tissue.id,
    tissue = .data$tissue,
    tissue.name = .data$tissue.name,
    oncotree.code = .data$oncotree.code,
    oncotree.name = .data$oncotree.name,
    oncotree.mainType = .data$oncotree.mainType,
    patient.id = .data$patient.id,
    drug = .data$drug,
    source_treatment_name = .data$source_treatment_name,
    source.dose = .data$dose,
    source.schedule = .data$schedule,
    source_data_row = .data$source_data_row,
    mouse.id = .data$mouse.id,
    sampleid = .data$sampleid,
    SA.ID = .data$SA.ID,
    RNAseq.PDX = .data$RNAseq.PDX,
    WGS.PDX = .data$WGS.PDX,
    Site = .data$Site,
    Histology = .data$Histology,
    ER = .data$ER,
    PR = .data$PR,
    HER2 = .data$HER2,
    AIMS = .data$AIMS
  ) |>
  arrange(.data$patient.id, .data$source_data_row)

write_csv(model_df, snk$output$model)
write_csv(sample_df, snk$output$sample)
write_csv(treatment_df, snk$output$treatment)
