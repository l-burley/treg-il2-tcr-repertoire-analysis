# ==============================================================================
# 12_prepare_and_run_gliph2_treg_tripod_reference.R
# ==============================================================================
# Goal:
#   Prepare GLIPH2 input from Curie003 Treg CDR3β sequences and run GLIPH2 using
#   a TRiPoD untreated B6 Treg reference repertoire.
#
#   The GLIPH2 input is restricted to Treg CDR3β clonotypes detected at
#   frequency >= 0.01% within each sample.
#
# Inputs:
#   data/processed/repseq_objects/RepSeqData_Curie003.RData
#   data/external/gliph2_reference/tripod_treg/ref_tripod_Treg_beta.txt
#   data/external/gliph2_reference/tripod_treg/ref_L_tripod_Treg.txt
#   data/external/gliph2_reference/tripod_treg/ref_V_tripod_Treg.txt
#
# Outputs:
#   results/gliph2/input/gliph2_input.tsv
#   results/gliph2/input/treg_filtered_0.01pct_with_freq.tsv
#   results/gliph2/input/tripod_reference_used.tsv
#   results/gliph2/input/tripod_length_freq_used.tsv
#   results/gliph2/input/tripod_vgene_freq_used.tsv
#   results/gliph2/input/run_info_tripod_treg_ref.tsv
#   results/gliph2/raw_output/
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Load packages
# ------------------------------------------------------------------------------

library(AnalyzAIRR)
library(turboGliph)
library(dplyr)
library(data.table)


# ------------------------------------------------------------------------------
# 1. Define input and output paths
# ------------------------------------------------------------------------------

repseq_path <- file.path(
  "data",
  "processed",
  "repseq_objects",
  "RepSeqData_Curie003.RData"
)

reference_dir <- file.path(
  "data",
  "external",
  "gliph2_reference",
  "tripod_treg"
)

ref_beta_path <- file.path(reference_dir, "ref_tripod_Treg_beta.txt")
ref_length_path <- file.path(reference_dir, "ref_L_tripod_Treg.txt")
ref_v_path <- file.path(reference_dir, "ref_V_tripod_Treg.txt")

gliph2_root <- file.path("results", "gliph2")
gliph2_input_dir <- file.path(gliph2_root, "input")
gliph2_raw_dir <- file.path(gliph2_root, "raw_output")
gliph2_processed_dir <- file.path(gliph2_root, "processed")
gliph2_results_dir <- file.path(gliph2_root, "results")
gliph2_figure_dir <- file.path(gliph2_root, "figures")

dir.create(gliph2_root, recursive = TRUE, showWarnings = FALSE)
dir.create(gliph2_input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gliph2_raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gliph2_processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gliph2_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gliph2_figure_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# 2. Set analysis parameters
# ------------------------------------------------------------------------------

frequency_threshold <- 0.0001
frequency_threshold_label <- "0.01%"

test_run_sim_depth <- 50
full_run_sim_depth <- 1000
n_cores <- 4

run_test_sample <- TRUE


# ------------------------------------------------------------------------------
# 3. Check that required input files exist
# ------------------------------------------------------------------------------

required_files <- c(
  repseq_path,
  ref_beta_path,
  ref_length_path,
  ref_v_path
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "The following required files are missing:\n",
    paste(missing_files, collapse = "\n")
  )
}


# ------------------------------------------------------------------------------
# 4. Load TRiPoD Treg reference files
# ------------------------------------------------------------------------------

ref_beta_raw <- fread(
  ref_beta_path,
  header = FALSE,
  col.names = c("CDR3b", "TRBV", "TRBJ")
)

ref_beta <- ref_beta_raw %>%
  select(
    CDR3b,
    TRBV
  ) %>%
  mutate(
    CDR3b = as.character(CDR3b),
    TRBV = as.character(TRBV)
  ) %>%
  filter(
    !is.na(CDR3b),
    CDR3b != "",
    !is.na(TRBV),
    TRBV != ""
  )

ref_length_freq <- fread(
  ref_length_path,
  header = FALSE,
  col.names = c("cdr3_length", "frequency")
) %>%
  mutate(
    cdr3_length = as.numeric(cdr3_length),
    frequency = as.numeric(frequency)
  )

ref_v_freq <- fread(
  ref_v_path,
  header = FALSE,
  col.names = c("TRBV", "frequency")
) %>%
  mutate(
    TRBV = as.character(TRBV),
    frequency = as.numeric(frequency)
  )

if (nrow(ref_beta) == 0) {
  stop("The TRiPoD beta-chain reference file contains no valid rows.")
}

if (nrow(ref_length_freq) == 0) {
  stop("The TRiPoD CDR3 length frequency file contains no valid rows.")
}

if (nrow(ref_v_freq) == 0) {
  stop("The TRiPoD V-gene frequency file contains no valid rows.")
}

fwrite(
  ref_beta,
  file.path(gliph2_input_dir, "tripod_reference_used.tsv"),
  sep = "\t"
)

fwrite(
  ref_length_freq,
  file.path(gliph2_input_dir, "tripod_length_freq_used.tsv"),
  sep = "\t"
)

fwrite(
  ref_v_freq,
  file.path(gliph2_input_dir, "tripod_vgene_freq_used.tsv"),
  sep = "\t"
)


# ------------------------------------------------------------------------------
# 5. Load Curie003 RepSeqData object
# ------------------------------------------------------------------------------

load(repseq_path)

if (!exists("RepSeqData_Curie003")) {
  stop("Object RepSeqData_Curie003 was not found after loading.")
}

repseq_object <- RepSeqData_Curie003

assay_data <- repseq_object@assayData %>%
  as.data.frame()

metadata <- repseq_object@metaData %>%
  as.data.frame()


# ------------------------------------------------------------------------------
# 6. Check required RepSeqData columns
# ------------------------------------------------------------------------------

required_assay_cols <- c(
  "sample_id",
  "aaCDR3",
  "V",
  "count"
)

required_meta_cols <- c(
  "sample_id",
  "mouse",
  "cell_subset",
  "injection",
  "protein"
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
# 7. Filter to Treg samples
# ------------------------------------------------------------------------------

treg_metadata <- metadata %>%
  mutate(
    sample_id = as.character(sample_id),
    mouse = as.character(mouse),
    cell_subset = as.character(cell_subset),
    injection = as.character(injection),
    protein = as.character(protein)
  ) %>%
  filter(cell_subset == "Treg") %>%
  distinct(sample_id, .keep_all = TRUE)

if (nrow(treg_metadata) == 0) {
  stop("No Treg samples were found in the metadata.")
}

treg_sample_ids <- treg_metadata$sample_id

treg_data <- assay_data %>%
  mutate(
    sample_id = as.character(sample_id),
    aaCDR3 = as.character(aaCDR3),
    V = as.character(V),
    count = as.numeric(count)
  ) %>%
  filter(sample_id %in% treg_sample_ids)

if (nrow(treg_data) == 0) {
  stop("No assayData rows were found for Treg samples.")
}


# ------------------------------------------------------------------------------
# 8. Apply 0.01% within-sample frequency threshold
# ------------------------------------------------------------------------------

treg_filtered <- treg_data %>%
  group_by(sample_id) %>%
  mutate(
    sample_total_reads = sum(count, na.rm = TRUE),
    freq = count / sample_total_reads
  ) %>%
  filter(freq >= frequency_threshold) %>%
  ungroup()

if (nrow(treg_filtered) == 0) {
  stop("No Treg clonotypes remained after applying the frequency threshold.")
}

per_sample_filter_summary <- treg_filtered %>%
  group_by(sample_id) %>%
  summarise(
    n_clones_kept = n(),
    pct_reads_kept = round(sum(freq, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  ) %>%
  left_join(
    treg_metadata %>%
      select(
        sample_id,
        mouse,
        injection,
        protein
      ),
    by = "sample_id"
  ) %>%
  arrange(injection, mouse)


# ------------------------------------------------------------------------------
# 9. Format input for GLIPH2
# ------------------------------------------------------------------------------

# The TRiPoD mouse reference V genes use names such as mTRBV13-1.
# If the Curie003 V genes are stored as TRBV13-1, an "m" prefix is added.

gliph2_input <- treg_filtered %>%
  mutate(
    V_mouse = ifelse(
      grepl("^mTRBV", V),
      V,
      paste0("m", V)
    )
  ) %>%
  select(
    CDR3b = aaCDR3,
    TRBV = V_mouse,
    counts = count,
    patient = sample_id
  ) %>%
  mutate(
    CDR3b = as.character(CDR3b),
    TRBV = as.character(TRBV),
    counts = as.numeric(counts),
    patient = as.character(patient)
  )

if (any(is.na(gliph2_input$CDR3b)) || any(gliph2_input$CDR3b == "")) {
  stop("NA or empty values found in GLIPH2 input column CDR3b.")
}

if (any(is.na(gliph2_input$TRBV)) || any(gliph2_input$TRBV == "")) {
  stop("NA or empty values found in GLIPH2 input column TRBV.")
}

if (any(is.na(gliph2_input$counts))) {
  stop("NA values found in GLIPH2 input column counts.")
}

if (any(is.na(gliph2_input$patient)) || any(gliph2_input$patient == "")) {
  stop("NA or empty values found in GLIPH2 input column patient.")
}


# ------------------------------------------------------------------------------
# 10. Check V-gene overlap between input and reference
# ------------------------------------------------------------------------------

input_v_genes <- sort(unique(gliph2_input$TRBV))
reference_v_genes <- sort(unique(ref_v_freq$TRBV))

v_genes_missing_from_reference <- setdiff(
  input_v_genes,
  reference_v_genes
)

v_gene_check <- data.frame(
  input_v_gene = input_v_genes,
  present_in_reference = input_v_genes %in% reference_v_genes
)

fwrite(
  v_gene_check,
  file.path(gliph2_input_dir, "vgene_reference_overlap_check.tsv"),
  sep = "\t"
)

if (length(v_genes_missing_from_reference) > 0) {
  warning(
    "Some input V genes were not found in the TRiPoD reference V-gene file: ",
    paste(v_genes_missing_from_reference, collapse = ", ")
  )
}


# ------------------------------------------------------------------------------
# 11. Save GLIPH2 input and filtering summary
# ------------------------------------------------------------------------------

fwrite(
  gliph2_input,
  file.path(gliph2_input_dir, "gliph2_input.tsv"),
  sep = "\t"
)

fwrite(
  treg_filtered,
  file.path(gliph2_input_dir, "treg_filtered_0.01pct_with_freq.tsv"),
  sep = "\t"
)

fwrite(
  per_sample_filter_summary,
  file.path(gliph2_input_dir, "treg_filtering_summary_0.01pct.tsv"),
  sep = "\t"
)


# ------------------------------------------------------------------------------
# 12. Optional test run on one sample
# ------------------------------------------------------------------------------

if (run_test_sample) {
  
  test_sample_id <- gliph2_input$patient[1]
  
  test_input <- gliph2_input %>%
    filter(patient == test_sample_id)
  
  message("Starting GLIPH2 test run on sample: ", test_sample_id)
  message("Sequences in test run: ", nrow(test_input))
  
  gliph2_test_result <- gliph2(
    cdr3_sequences = test_input,
    refdb_beta = ref_beta,
    v_usage_freq = ref_v_freq,
    cdr3_length_freq = ref_length_freq,
    sim_depth = test_run_sim_depth,
    n_cores = 1,
    all_aa_interchangeable = TRUE
  )
  
  message("GLIPH2 test run complete.")
  
  if (!is.null(gliph2_test_result$cluster_properties)) {
    message("Test run cluster properties:")
    print(head(gliph2_test_result$cluster_properties))
  }
}


# ------------------------------------------------------------------------------
# 13. Full GLIPH2 run with TRiPoD Treg reference
# ------------------------------------------------------------------------------

message("Starting full GLIPH2 run.")
message("Input rows: ", nrow(gliph2_input))
message("Input samples: ", n_distinct(gliph2_input$patient))
message("Reference rows: ", nrow(ref_beta))
message("Output folder: ", gliph2_raw_dir)

gliph2_full_result <- gliph2(
  cdr3_sequences = gliph2_input,
  refdb_beta = ref_beta,
  v_usage_freq = ref_v_freq,
  cdr3_length_freq = ref_length_freq,
  sim_depth = full_run_sim_depth,
  n_cores = n_cores,
  result_folder = gliph2_raw_dir,
  all_aa_interchangeable = TRUE
)

message("Full GLIPH2 run complete.")

if (!is.null(gliph2_full_result$cluster_properties)) {
  message("Total clusters found: ", nrow(gliph2_full_result$cluster_properties))
  message("First rows of cluster properties:")
  print(head(gliph2_full_result$cluster_properties))
}


# ------------------------------------------------------------------------------
# 14. Save run information
# ------------------------------------------------------------------------------

run_info <- data.frame(
  item = c(
    "repseq_path",
    "reference_beta_file",
    "reference_length_file",
    "reference_v_file",
    "gliph2_input_file",
    "raw_output_folder",
    "frequency_threshold_label",
    "frequency_threshold_fraction",
    "cell_subset",
    "test_run_sim_depth",
    "full_run_sim_depth",
    "n_cores",
    "all_aa_interchangeable",
    "n_treg_samples",
    "n_input_rows"
  ),
  value = c(
    repseq_path,
    ref_beta_path,
    ref_length_path,
    ref_v_path,
    file.path(gliph2_input_dir, "gliph2_input.tsv"),
    gliph2_raw_dir,
    frequency_threshold_label,
    as.character(frequency_threshold),
    "Treg",
    as.character(test_run_sim_depth),
    as.character(full_run_sim_depth),
    as.character(n_cores),
    "TRUE",
    as.character(n_distinct(gliph2_input$patient)),
    as.character(nrow(gliph2_input))
  )
)

fwrite(
  run_info,
  file.path(gliph2_input_dir, "run_info_tripod_treg_ref.tsv"),
  sep = "\t"
)


# ------------------------------------------------------------------------------
# 15. Print final summary
# ------------------------------------------------------------------------------

message("GLIPH2 input preparation and run complete.")
message("Input files saved to: ", gliph2_input_dir)
message("Raw GLIPH2 output saved to: ", gliph2_raw_dir)

message("Filtering summary:")
print(per_sample_filter_summary)

message("Run information:")
print(run_info)