# ==============================================================================
# 02_check_bulk_tcr_depth_and_metadata.R
# ==============================================================================
# Goal:
#   Check sample metadata and sequencing depth for the Curie003 bulk TCRβ
#   RepSeqData object before downstream repertoire analyses.
#
# Input:
#   data/processed/repseq_objects/RepSeqData_Curie003.RData
#
# Outputs:
#   results/qc/bulk_tcr_depth_by_sample.csv
#   results/qc/bulk_tcr_depth_by_group.csv
#   results/qc/bulk_tcr_paired_mouse_check.csv
#
# Notes:
#   This script does not modify the RepSeqData object.
#   Internal path columns are not exported.
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Load packages
# ------------------------------------------------------------------------------

library(AnalyzAIRR)
library(data.table)
library(dplyr)
library(readr)


# ------------------------------------------------------------------------------
# 1. Define input and output paths
# ------------------------------------------------------------------------------

repseq_path <- file.path(
  "data",
  "processed",
  "repseq_objects",
  "RepSeqData_Curie003.RData"
)

qc_dir <- file.path("results", "qc")

dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# 2. Load RepSeqData object
# ------------------------------------------------------------------------------

load(repseq_path)

if (!exists("RepSeqData_Curie003")) {
  stop("Object RepSeqData_Curie003 was not found after loading.")
}


# ------------------------------------------------------------------------------
# 3. Extract assayData and metaData
# ------------------------------------------------------------------------------

assay_data <- as.data.table(RepSeqData_Curie003@assayData)
metadata <- as.data.frame(RepSeqData_Curie003@metaData)


# ------------------------------------------------------------------------------
# 4. Check required columns
# ------------------------------------------------------------------------------

required_assay_cols <- c(
  "sample_id",
  "aaCDR3",
  "ntCDR3",
  "aaClone",
  "ntClone",
  "count"
)

required_meta_cols <- c(
  "sample_id",
  "mouse",
  "cell_subset",
  "injection",
  "nSequences"
)

missing_assay_cols <- setdiff(required_assay_cols, colnames(assay_data))
missing_meta_cols <- setdiff(required_meta_cols, colnames(metadata))

if (length(missing_assay_cols) > 0) {
  stop(
    "Missing assayData columns: ",
    paste(missing_assay_cols, collapse = ", ")
  )
}

if (length(missing_meta_cols) > 0) {
  stop(
    "Missing metaData columns: ",
    paste(missing_meta_cols, collapse = ", ")
  )
}


# ------------------------------------------------------------------------------
# 5. Clean key metadata columns
# ------------------------------------------------------------------------------

assay_data <- assay_data %>%
  mutate(sample_id = as.character(sample_id))

metadata <- metadata %>%
  mutate(
    sample_id = as.character(sample_id),
    mouse = as.character(mouse),
    cell_subset = as.character(cell_subset),
    injection = as.character(injection)
  )


# ------------------------------------------------------------------------------
# 6. Check sample ID matching between assayData and metaData
# ------------------------------------------------------------------------------

assay_samples <- sort(unique(assay_data$sample_id))
metadata_samples <- sort(unique(metadata$sample_id))

samples_only_in_assay <- setdiff(assay_samples, metadata_samples)
samples_only_in_metadata <- setdiff(metadata_samples, assay_samples)

if (length(samples_only_in_assay) > 0) {
  warning(
    "Samples present in assayData but missing from metaData: ",
    paste(samples_only_in_assay, collapse = ", ")
  )
}

if (length(samples_only_in_metadata) > 0) {
  warning(
    "Samples present in metaData but missing from assayData: ",
    paste(samples_only_in_metadata, collapse = ", ")
  )
}


# ------------------------------------------------------------------------------
# 7. Calculate sequencing depth and clonotype counts from assayData
# ------------------------------------------------------------------------------

depth_from_assay <- assay_data %>%
  group_by(sample_id) %>%
  summarise(
    total_count_from_assay = sum(count, na.rm = TRUE),
    n_clonotype_rows = n(),
    n_unique_aaCDR3_from_assay = n_distinct(aaCDR3),
    n_unique_ntCDR3_from_assay = n_distinct(ntCDR3),
    n_unique_aaClone_from_assay = n_distinct(aaClone),
    n_unique_ntClone_from_assay = n_distinct(ntClone),
    .groups = "drop"
  )


# ------------------------------------------------------------------------------
# 8. Create sample-level QC table
# ------------------------------------------------------------------------------

metadata_export_cols <- c(
  "sample_id",
  "mouse",
  "cell_subset",
  "organ",
  "injection",
  "protein",
  "dose",
  "perc_Tregs_spleen",
  "nSequences",
  "aaCDR3",
  "ntCDR3",
  "aaClone",
  "ntClone",
  "chao1",
  "iChao"
)

metadata_export_cols <- intersect(metadata_export_cols, colnames(metadata))

depth_by_sample <- metadata %>%
  select(all_of(metadata_export_cols)) %>%
  left_join(depth_from_assay, by = "sample_id") %>%
  mutate(
    depth_difference = total_count_from_assay - nSequences,
    depth_matches_metadata = depth_difference == 0
  ) %>%
  arrange(cell_subset, injection, mouse)


# ------------------------------------------------------------------------------
# 9. Create group-level QC table
# ------------------------------------------------------------------------------

depth_by_group <- depth_by_sample %>%
  group_by(cell_subset, injection) %>%
  summarise(
    n_samples = n(),
    n_mice = n_distinct(mouse),
    total_sequences = sum(nSequences, na.rm = TRUE),
    min_sequences = min(nSequences, na.rm = TRUE),
    median_sequences = median(nSequences, na.rm = TRUE),
    mean_sequences = mean(nSequences, na.rm = TRUE),
    max_sequences = max(nSequences, na.rm = TRUE),
    min_aaCDR3 = min(aaCDR3, na.rm = TRUE),
    median_aaCDR3 = median(aaCDR3, na.rm = TRUE),
    mean_aaCDR3 = mean(aaCDR3, na.rm = TRUE),
    max_aaCDR3 = max(aaCDR3, na.rm = TRUE),
    .groups = "drop"
  )


# ------------------------------------------------------------------------------
# 10. Check paired Treg and Tconv samples by mouse
# ------------------------------------------------------------------------------

paired_mouse_check <- metadata %>%
  group_by(mouse, cell_subset) %>%
  summarise(n_samples = n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = cell_subset,
    values_from = n_samples,
    values_fill = 0
  )

if (!"Treg" %in% colnames(paired_mouse_check)) {
  paired_mouse_check$Treg <- 0
}

if (!"Tconv" %in% colnames(paired_mouse_check)) {
  paired_mouse_check$Tconv <- 0
}

paired_mouse_check <- paired_mouse_check %>%
  mutate(
    has_paired_treg_tconv = Treg > 0 & Tconv > 0
  ) %>%
  arrange(mouse)


# ------------------------------------------------------------------------------
# 11. Save QC tables
# ------------------------------------------------------------------------------

write_csv(
  depth_by_sample,
  file.path(qc_dir, "bulk_tcr_depth_by_sample.csv")
)

write_csv(
  depth_by_group,
  file.path(qc_dir, "bulk_tcr_depth_by_group.csv")
)

write_csv(
  paired_mouse_check,
  file.path(qc_dir, "bulk_tcr_paired_mouse_check.csv")
)