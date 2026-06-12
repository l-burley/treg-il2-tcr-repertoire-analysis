# ==============================================================================
# 01_create_repseq_objects_from_mixcr.R
# ==============================================================================
# Goal:
#   Create AnalyzAIRR RepSeqData objects from MiXCR TRB clonotype tables.
#
# Input:
#   data/raw/metadata.csv
#
# Required metadata columns:
#   sample_id
#   cell_subset
#   injection
#   virtual_filepath
#
# Outputs:
#   data/processed/repseq_objects/RepSeqData_Curie003_treg.RData
#   data/processed/repseq_objects/RepSeqData_Curie003_tconv.RData
#   data/processed/repseq_objects/RepSeqData_Curie003_all.RData
#
# Notes:
#   Raw MiXCR files and full metadata are not included in the GitHub repository.
#   The metadata file should contain paths to the MiXCR clonotype tables on the
#   machine where the analysis is run.
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Load packages
# ------------------------------------------------------------------------------

library(AnalyzAIRR)
library(dplyr)


# ------------------------------------------------------------------------------
# 1. Define relative input and output paths
# ------------------------------------------------------------------------------

metadata_path <- file.path("data", "raw", "metadata.csv")

output_dir <- file.path("data", "processed", "repseq_objects")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# 2. Load metadata
# ------------------------------------------------------------------------------

metadata <- read.csv(metadata_path, header = TRUE, stringsAsFactors = FALSE)


# ------------------------------------------------------------------------------
# 3. Check required metadata columns
# ------------------------------------------------------------------------------

required_cols <- c(
  "sample_id",
  "cell_subset",
  "injection",
  "virtual_filepath"
)

missing_cols <- setdiff(required_cols, colnames(metadata))

if (length(missing_cols) > 0) {
  stop(
    "The metadata file is missing the following required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}


# ------------------------------------------------------------------------------
# 4. Check file paths listed in metadata
# ------------------------------------------------------------------------------

missing_files <- metadata$virtual_filepath[!file.exists(metadata$virtual_filepath)]

if (length(missing_files) > 0) {
  stop(
    "Some MiXCR files listed in metadata$virtual_filepath were not found:\n",
    paste(missing_files, collapse = "\n")
  )
}


# ------------------------------------------------------------------------------
# 5. Prepare metadata columns
# ------------------------------------------------------------------------------

metadata <- metadata %>%
  mutate(
    sample_id = as.character(sample_id),
    cell_subset = factor(cell_subset),
    injection = factor(injection)
  )


# ------------------------------------------------------------------------------
# 6. Split metadata into Treg and Tconv samples
# ------------------------------------------------------------------------------

metadata_treg <- metadata %>%
  filter(cell_subset == "Treg")

metadata_tconv <- metadata %>%
  filter(cell_subset == "Tconv")


# ------------------------------------------------------------------------------
# 7. Check sample numbers
# ------------------------------------------------------------------------------

if (nrow(metadata_treg) == 0) {
  stop("No Treg samples were found in metadata$cell_subset.")
}

if (nrow(metadata_tconv) == 0) {
  stop("No Tconv samples were found in metadata$cell_subset.")
}

message("Number of Treg samples: ", nrow(metadata_treg))
message("Number of Tconv samples: ", nrow(metadata_tconv))


# ------------------------------------------------------------------------------
# 8. Prepare sample metadata for AnalyzAIRR
# ------------------------------------------------------------------------------

rownames(metadata_treg) <- metadata_treg$sample_id
rownames(metadata_tconv) <- metadata_tconv$sample_id

metadata_treg_airr <- metadata_treg %>%
  select(-sample_id)

metadata_tconv_airr <- metadata_tconv %>%
  select(-sample_id)


# ------------------------------------------------------------------------------
# 9. Create RepSeqData object for Treg samples
# ------------------------------------------------------------------------------

RepSeqData_Curie003_treg <- readAIRRSet(
  fileList = metadata_treg_airr$virtual_filepath,
  fileFormat = "MiXCR",
  chain = "TRB",
  sampleinfo = metadata_treg_airr,
  keep.ambiguous = FALSE,
  keep.unproductive = FALSE,
  filter.singletons = FALSE,
  aa.th = 8,
  outFiltered = FALSE,
  cores = 1L
)

save(
  RepSeqData_Curie003_treg,
  file = file.path(output_dir, "RepSeqData_Curie003_treg.RData")
)


# ------------------------------------------------------------------------------
# 10. Create RepSeqData object for Tconv samples
# ------------------------------------------------------------------------------

RepSeqData_Curie003_tconv <- readAIRRSet(
  fileList = metadata_tconv_airr$virtual_filepath,
  fileFormat = "MiXCR",
  chain = "TRB",
  sampleinfo = metadata_tconv_airr,
  keep.ambiguous = FALSE,
  keep.unproductive = FALSE,
  filter.singletons = FALSE,
  aa.th = 8,
  outFiltered = FALSE,
  cores = 1L
)

save(
  RepSeqData_Curie003_tconv,
  file = file.path(output_dir, "RepSeqData_Curie003_tconv.RData")
)


# ------------------------------------------------------------------------------
# 11. Merge Treg and Tconv RepSeqData objects
# ------------------------------------------------------------------------------

RepSeqData_Curie003 <- mergeRepSeq(
  a = RepSeqData_Curie003_treg,
  b = RepSeqData_Curie003_tconv
)

save(
  RepSeqData_Curie003,
  file = file.path(output_dir, "RepSeqData_Curie003.RData")
)


# ------------------------------------------------------------------------------
# 12. Save cleaned metadata used for object creation
# ------------------------------------------------------------------------------

write.csv(
  metadata,
  file = file.path(output_dir, "metadata_Curie003_used_for_repseq_objects.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 13. Print final summary
# ------------------------------------------------------------------------------

message("RepSeqData objects created successfully.")
message("Output directory: ", output_dir)

message("Saved files:")
message("  - RepSeqData_Curie003_treg.RData")
message("  - RepSeqData_Curie003_tconv.RData")
message("  - RepSeqData_Curie003_all.RData")
message("  - metadata_Curie003_used_for_repseq_objects.csv")