# ==============================================================================
# 04_calculate_bulk_treg_clonotype_frequencies.R
# ==============================================================================
# Goal:
#   Calculate within-sample clonotype frequencies and frequency-bin summaries
#   for Curie003 bulk Treg TCRβ repertoires.
#
# Input:
#   data/processed/repseq_objects/RepSeqData_Curie003.RData
#
# Outputs:
#   results/intermediate/bulk_treg_clonotype_frequencies.csv
#   results/intermediate/bulk_treg_frequency_bin_summary.csv
#
# Notes:
#   This script is focused on Treg samples, as these cells are the 
#   focus of the project
#
#   Frequencies are calculated as:
#     clonotype count / total TCRβ sequences in the corresponding sample
#
#   Frequency bins:
#     ]0.000001, 0.00001]
#     ]0.00001, 0.0001]
#     ]0.0001, 0.001]
#     ]0.001, 0.01]
#     ]0.01, 1]
#
#   The notation ]a, b] means greater than a and less than or equal to b.
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

output_dir <- file.path("results", "intermediate")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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
# 4. Check columns needed for this script
# ------------------------------------------------------------------------------

required_assay_cols <- c(
  "sample_id",
  "V",
  "J",
  "ntCDR3",
  "aaCDR3",
  "aaClone",
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
    "Missing assayData columns needed for this script: ",
    paste(missing_assay_cols, collapse = ", ")
  )
}

if (length(missing_meta_cols) > 0) {
  stop(
    "Missing metaData columns needed for this script: ",
    paste(missing_meta_cols, collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# 5. Clean key columns
# ------------------------------------------------------------------------------

assay_data <- assay_data %>%
  mutate(
    sample_id = as.character(sample_id),
    V = as.character(V),
    J = as.character(J),
    ntCDR3 = as.character(ntCDR3),
    aaCDR3 = as.character(aaCDR3),
    aaClone = as.character(aaClone),
    count = as.numeric(count)
  )

metadata <- metadata %>%
  mutate(
    sample_id = as.character(sample_id),
    mouse = as.character(mouse),
    cell_subset = as.character(cell_subset),
    injection = as.character(injection),
    nSequences = as.numeric(nSequences)
  )

# ------------------------------------------------------------------------------
# 6. Keep Treg metadata only
# ------------------------------------------------------------------------------

metadata_treg <- metadata %>%
  filter(cell_subset == "Treg") %>%
  select(
    sample_id,
    mouse,
    cell_subset,
    injection,
    nSequences
  )

if (nrow(metadata_treg) == 0) {
  stop("No Treg samples were found in metaData.")
}

# ------------------------------------------------------------------------------
# 7. Keep Treg clonotypes only
# ------------------------------------------------------------------------------

treg_assay_data <- assay_data %>%
  filter(sample_id %in% metadata_treg$sample_id)

if (nrow(treg_assay_data) == 0) {
  stop("No Treg clonotypes were found in assayData.")
}

# ------------------------------------------------------------------------------
# 8. Calculate total sequence count per Treg sample
# ------------------------------------------------------------------------------

treg_sample_depth <- treg_assay_data %>%
  group_by(sample_id) %>%
  summarise(
    total_count_from_assay = sum(count, na.rm = TRUE),
    n_clonotype_rows = n(),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# 9. Calculate clonotype frequencies
# ------------------------------------------------------------------------------

treg_clonotype_frequencies <- treg_assay_data %>%
  left_join(treg_sample_depth, by = "sample_id") %>%
  left_join(metadata_treg, by = "sample_id") %>%
  mutate(
    frequency = count / total_count_from_assay,
    frequency_percent = frequency * 100
  )

# ------------------------------------------------------------------------------
# 10. Check frequency sums per sample
# ------------------------------------------------------------------------------

frequency_check <- treg_clonotype_frequencies %>%
  group_by(sample_id) %>%
  summarise(
    total_frequency = sum(frequency, na.rm = TRUE),
    total_frequency_percent = sum(frequency_percent, na.rm = TRUE),
    .groups = "drop"
  )

bad_frequency_sums <- frequency_check %>%
  filter(abs(total_frequency - 1) > 1e-6)

if (nrow(bad_frequency_sums) > 0) {
  warning(
    "Some Treg samples have frequency sums different from 1. Check frequency_check."
  )
}

# ------------------------------------------------------------------------------
# 11. Add frequency bins
# ------------------------------------------------------------------------------

treg_clonotype_frequencies <- treg_clonotype_frequencies %>%
  mutate(
    frequency_bin = case_when(
      frequency > 0.01 & frequency <= 1 ~ "]0.01,1]",
      frequency > 0.001 & frequency <= 0.01 ~ "]0.001,0.01]",
      frequency > 0.0001 & frequency <= 0.001 ~ "]0.0001,0.001]",
      frequency > 0.00001 & frequency <= 0.0001 ~ "]0.00001,0.0001]",
      frequency > 0.000001 & frequency <= 0.00001 ~ "]0.000001,0.00001]",
      TRUE ~ NA_character_
    ),
    frequency_bin = factor(
      frequency_bin,
      levels = c(
        "]0.000001,0.00001]",
        "]0.00001,0.0001]",
        "]0.0001,0.001]",
        "]0.001,0.01]",
        "]0.01,1]"
      )
    )
  )

# ------------------------------------------------------------------------------
# 12. Create frequency-bin summary table
# ------------------------------------------------------------------------------

treg_frequency_bin_summary <- treg_clonotype_frequencies %>%
  filter(!is.na(frequency_bin)) %>%
  group_by(
    sample_id,
    mouse,
    cell_subset,
    injection,
    frequency_bin
  ) %>%
  summarise(
    n_clonotypes = n(),
    total_count = sum(count, na.rm = TRUE),
    total_frequency = sum(frequency, na.rm = TRUE),
    total_frequency_percent = total_frequency * 100,
    .groups = "drop"
  ) %>%
  arrange(injection, mouse, frequency_bin)

# ------------------------------------------------------------------------------
# 13. Save output tables
# ------------------------------------------------------------------------------

write_csv(
  treg_clonotype_frequencies,
  file.path(output_dir, "bulk_treg_clonotype_frequencies.csv")
)

write_csv(
  treg_frequency_bin_summary,
  file.path(output_dir, "bulk_treg_frequency_bin_summary.csv")
)

# ------------------------------------------------------------------------------
# 14. Print quick checks
# ------------------------------------------------------------------------------

message("Bulk Treg clonotype frequency calculation complete.")
message("Output directory: ", output_dir)

message("Number of Treg samples: ", n_distinct(treg_clonotype_frequencies$sample_id))
message("Number of Treg clonotype rows: ", nrow(treg_clonotype_frequencies))

message("Frequency sum check:")
print(summary(frequency_check$total_frequency))

message("Frequency bins:")
print(table(treg_clonotype_frequencies$frequency_bin, useNA = "ifany"))