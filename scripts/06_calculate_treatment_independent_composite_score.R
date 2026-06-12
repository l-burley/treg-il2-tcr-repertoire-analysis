# ==============================================================================
# 06_calculate_treatment_independent_composite_score.R
# ==============================================================================
# Goal:
#   Calculate treatment-independent composite scores for shared expanded Treg
#   CDR3β sequences.
#
# Input:
#   data/processed/repseq_objects/RepSeqData_Curie003.RData
#
# Outputs:
#   results/tables/cdr3_scores_full_treatment_independent.csv
#   results/tables/cdr3_scores_full_treatment_independent.tsv
#   results/tables/cdr3_scores_top20_treatment_independent.csv
#   results/tables/cdr3_scores_top100_treatment_independent.csv
#   results/tables/cdr3_scores_top100_treatment_independent.tsv
#   results/tables/cdr3_score_components_raw_treatment_independent.csv
#   results/tables/cdr3_score_components_raw_treatment_independent.tsv
#
# Filtering criteria:
#   1. Treg samples only
#   2. CDR3β frequency >= 0.01% within a Treg sample
#   3. Exact amino-acid CDR3β detected in at least 2 mice
#
# Composite score criteria:
#   1. Number of mice containing the exact CDR3β
#   2. Mean frequency across detected mice
#   3. LD=1 neighbour support
#   4. Frequency consistency across detected mice
#   5. Number of mice with frequency >= 0.1%
#
# Notes:
#   This score is treatment-independent because it does not include IL-2/PBS
#   enrichment criteria. It is used for comparison with the IL-2-prioritised
#   composite score in Figure 3.
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Load packages
# ------------------------------------------------------------------------------

library(AnalyzAIRR)
library(data.table)
library(dplyr)
library(tidyr)
library(stringdist)
library(stringr)


# ------------------------------------------------------------------------------
# 1. Define input and output paths
# ------------------------------------------------------------------------------

repseq_path <- file.path(
  "data",
  "processed",
  "repseq_objects",
  "RepSeqData_Curie003.RData"
)

table_dir <- file.path("results", "tables")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# 2. Load RepSeqData object
# ------------------------------------------------------------------------------

load(repseq_path)

if (!exists("RepSeqData_Curie003")) {
  stop("Object RepSeqData_Curie003 was not found after loading.")
}

repseq_object <- RepSeqData_Curie003


# ------------------------------------------------------------------------------
# 3. Build Treg metadata table
# ------------------------------------------------------------------------------

sc_meta <- repseq_object@metaData %>%
  as.data.frame() %>%
  filter(cell_subset == "Treg") %>%
  mutate(
    sample_id = as.character(sample_id),
    mouse = as.character(mouse),
    injection = as.character(injection),
    protein = as.character(protein),
    protein_clean = case_when(
      str_detect(protein, "MCB") ~ "mutant_IL2",
      str_detect(protein, "CB") ~ "IL2",
      str_detect(protein, "PBS") ~ "PBS",
      TRUE ~ protein
    )
  ) %>%
  distinct(sample_id, .keep_all = TRUE) %>%
  select(
    sample_id,
    mouse,
    injection,
    protein_clean
  )

sc_treg_ids <- sc_meta$sample_id

if (length(sc_treg_ids) == 0) {
  stop("No Treg samples were found in the metadata.")
}


# ------------------------------------------------------------------------------
# 4. Compute per-sample frequency for all Treg samples
# ------------------------------------------------------------------------------

sc_treg_raw <- repseq_object@assayData %>%
  as.data.frame() %>%
  mutate(
    sample_id = as.character(sample_id),
    aaCDR3 = as.character(aaCDR3),
    count = as.numeric(count)
  ) %>%
  filter(sample_id %in% sc_treg_ids) %>%
  group_by(sample_id) %>%
  mutate(
    sample_total = sum(count, na.rm = TRUE),
    freq = count / sample_total
  ) %>%
  ungroup() %>%
  filter(freq >= 0.0001) %>%
  left_join(sc_meta, by = "sample_id") %>%
  select(
    sample_id,
    mouse,
    injection,
    protein_clean,
    cdr3_aa = aaCDR3,
    count,
    freq,
    sample_total
  )

if (nrow(sc_treg_raw) == 0) {
  stop("No Treg CDR3β rows remained after applying the 0.01% frequency filter.")
}


# ------------------------------------------------------------------------------
# 5. Filter to CDR3β sequences present in at least 2 mice
# ------------------------------------------------------------------------------

sc_mice_per_cdr3 <- sc_treg_raw %>%
  group_by(cdr3_aa) %>%
  summarise(
    n_mice_total = n_distinct(mouse),
    .groups = "drop"
  )

sc_shared_cdr3s <- sc_mice_per_cdr3 %>%
  filter(n_mice_total >= 2) %>%
  pull(cdr3_aa)

sc_treg_filt <- sc_treg_raw %>%
  filter(cdr3_aa %in% sc_shared_cdr3s)

if (length(sc_shared_cdr3s) == 0) {
  stop("No CDR3β sequences were detected in at least 2 mice.")
}


# ------------------------------------------------------------------------------
# 6. Criterion 1: number of mice with exact CDR3β
# ------------------------------------------------------------------------------

sc_c1 <- sc_mice_per_cdr3 %>%
  filter(cdr3_aa %in% sc_shared_cdr3s) %>%
  rename(c1_n_mice_exact = n_mice_total)


# ------------------------------------------------------------------------------
# 7. Criterion 2: mean frequency across detected mice
# ------------------------------------------------------------------------------

sc_c2 <- sc_treg_filt %>%
  group_by(cdr3_aa, mouse) %>%
  summarise(
    freq_in_mouse = sum(freq, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(cdr3_aa) %>%
  summarise(
    c2_mean_freq = mean(freq_in_mouse, na.rm = TRUE),
    c2_median_freq = median(freq_in_mouse, na.rm = TRUE),
    c2_max_freq = max(freq_in_mouse, na.rm = TRUE),
    c2_sd_freq = sd(freq_in_mouse, na.rm = TRUE),
    .groups = "drop"
  )


# ------------------------------------------------------------------------------
# 8. Criterion 3: LD=1 neighbour support
# ------------------------------------------------------------------------------

message("Computing pairwise Levenshtein distances for ",
        length(sc_shared_cdr3s), " CDR3β sequences...")

sc_ld_matrix <- stringdistmatrix(
  sc_shared_cdr3s,
  sc_shared_cdr3s,
  method = "lv"
)

rownames(sc_ld_matrix) <- sc_shared_cdr3s
colnames(sc_ld_matrix) <- sc_shared_cdr3s

sc_ld1_neighbors <- lapply(sc_shared_cdr3s, function(cdr3) {
  dists <- sc_ld_matrix[cdr3, ]
  names(dists[dists == 1])
})

names(sc_ld1_neighbors) <- sc_shared_cdr3s

sc_cdr3_to_mice <- sc_treg_filt %>%
  distinct(cdr3_aa, mouse) %>%
  group_by(cdr3_aa) %>%
  summarise(
    mice_set = list(unique(mouse)),
    .groups = "drop"
  )

sc_cdr3_mice_lookup <- setNames(
  sc_cdr3_to_mice$mice_set,
  sc_cdr3_to_mice$cdr3_aa
)

sc_c3 <- data.frame(
  cdr3_aa = sc_shared_cdr3s,
  c3_n_neighbors = sapply(sc_ld1_neighbors, length),
  c3_ld1_neighbor_mice = sapply(sc_shared_cdr3s, function(cdr3) {
    neighbors <- sc_ld1_neighbors[[cdr3]]
    
    if (length(neighbors) == 0) {
      return(0L)
    }
    
    length(unique(unlist(sc_cdr3_mice_lookup[neighbors])))
  }),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------------------------
# 9. Criterion 4: frequency consistency across detected mice
# ------------------------------------------------------------------------------

sc_c4 <- sc_treg_filt %>%
  group_by(cdr3_aa, mouse) %>%
  summarise(
    freq_in_mouse = sum(freq, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(cdr3_aa) %>%
  summarise(
    c4_cv = ifelse(
      mean(freq_in_mouse) == 0,
      0,
      sd(freq_in_mouse) / mean(freq_in_mouse)
    ),
    .groups = "drop"
  ) %>%
  mutate(
    c4_cv = pmin(c4_cv, 2)
  )


# ------------------------------------------------------------------------------
# 10. Criterion 5: number of mice above 0.1% frequency
# ------------------------------------------------------------------------------

sc_c5 <- sc_treg_filt %>%
  group_by(cdr3_aa, mouse) %>%
  summarise(
    freq_in_mouse = sum(freq, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(cdr3_aa) %>%
  summarise(
    c5_n_mice_above_0.1pct = sum(freq_in_mouse >= 0.001, na.rm = TRUE),
    .groups = "drop"
  )


# ------------------------------------------------------------------------------
# 11. Combine all raw score components
# ------------------------------------------------------------------------------

sc_raw_scores <- sc_c1 %>%
  left_join(sc_c2, by = "cdr3_aa") %>%
  left_join(sc_c3, by = "cdr3_aa") %>%
  left_join(sc_c4, by = "cdr3_aa") %>%
  left_join(sc_c5, by = "cdr3_aa") %>%
  mutate(
    c3_ld1_neighbor_mice = replace_na(c3_ld1_neighbor_mice, 0),
    c3_n_neighbors = replace_na(c3_n_neighbors, 0),
    c5_n_mice_above_0.1pct = replace_na(c5_n_mice_above_0.1pct, 0)
  )


# ------------------------------------------------------------------------------
# 12. Normalize criteria and calculate treatment-independent score
# ------------------------------------------------------------------------------

normalize_01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  
  if (rng[1] == rng[2]) {
    return(rep(0.5, length(x)))
  }
  
  (x - rng[1]) / (rng[2] - rng[1])
}

sc_scored <- sc_raw_scores %>%
  mutate(
    norm_c1 = normalize_01(c1_n_mice_exact),
    norm_c2 = normalize_01(log10(c2_mean_freq + 1e-9)),
    norm_c3 = normalize_01(c3_ld1_neighbor_mice),
    norm_c4 = 1 - normalize_01(c4_cv),
    norm_c5 = normalize_01(c5_n_mice_above_0.1pct),
    
    treatment_independent_score =
      0.35 * norm_c1 +
      0.25 * norm_c2 +
      0.20 * norm_c3 +
      0.15 * norm_c4 +
      0.05 * norm_c5
  ) %>%
  arrange(desc(treatment_independent_score))


# ------------------------------------------------------------------------------
# 13. Add annotation columns
# ------------------------------------------------------------------------------

sc_mouse_list <- sc_treg_filt %>%
  distinct(cdr3_aa, mouse) %>%
  group_by(cdr3_aa) %>%
  summarise(
    mouse_list = paste(sort(unique(mouse)), collapse = ", "),
    .groups = "drop"
  )

sc_neighbor_list <- data.frame(
  cdr3_aa = names(sc_ld1_neighbors),
  ld1_neighbors = sapply(sc_ld1_neighbors, function(neighbors) {
    if (length(neighbors) == 0) {
      return("none")
    }
    
    paste(neighbors, collapse = ", ")
  }),
  stringsAsFactors = FALSE
)

assign_bin <- function(freq) {
  case_when(
    freq > 0.01 ~ "bin1_>1pct",
    freq > 0.001 ~ "bin2_0.1-1pct",
    freq > 0.0001 ~ "bin3_0.01-0.1pct",
    TRUE ~ "bin4_<0.01pct"
  )
}

sc_bin_summary <- sc_treg_filt %>%
  group_by(cdr3_aa, mouse) %>%
  summarise(
    freq_in_mouse = sum(freq, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    freq_bin = assign_bin(freq_in_mouse)
  ) %>%
  group_by(cdr3_aa) %>%
  summarise(
    bin_summary = paste(paste0(mouse, ":", freq_bin), collapse = "; "),
    .groups = "drop"
  )

sc_final <- sc_scored %>%
  left_join(sc_mouse_list, by = "cdr3_aa") %>%
  left_join(sc_neighbor_list, by = "cdr3_aa") %>%
  left_join(sc_bin_summary, by = "cdr3_aa") %>%
  mutate(
    treatment_independent_rank = row_number()
  ) %>%
  select(
    treatment_independent_rank,
    cdr3_aa,
    treatment_independent_score,
    c1_n_mice_exact,
    c2_mean_freq,
    c2_max_freq,
    c3_n_neighbors,
    c3_ld1_neighbor_mice,
    c4_cv,
    c5_n_mice_above_0.1pct,
    norm_c1,
    norm_c2,
    norm_c3,
    norm_c4,
    norm_c5,
    mouse_list,
    ld1_neighbors,
    bin_summary
  )


# ------------------------------------------------------------------------------
# 14. Save outputs
# ------------------------------------------------------------------------------

fwrite(
  sc_final,
  file.path(table_dir, "cdr3_scores_full_treatment_independent.tsv"),
  sep = "\t"
)

fwrite(
  sc_final,
  file.path(table_dir, "cdr3_scores_full_treatment_independent.csv")
)

fwrite(
  sc_final %>% slice_head(n = 20),
  file.path(table_dir, "cdr3_scores_top20_treatment_independent.csv")
)

fwrite(
  sc_final %>% slice_head(n = 100),
  file.path(table_dir, "cdr3_scores_top100_treatment_independent.tsv"),
  sep = "\t"
)

fwrite(
  sc_final %>% slice_head(n = 100),
  file.path(table_dir, "cdr3_scores_top100_treatment_independent.csv")
)

fwrite(
  sc_raw_scores,
  file.path(table_dir, "cdr3_score_components_raw_treatment_independent.tsv"),
  sep = "\t"
)

fwrite(
  sc_raw_scores,
  file.path(table_dir, "cdr3_score_components_raw_treatment_independent.csv")
)


# ------------------------------------------------------------------------------
# 15. Print quick checks
# ------------------------------------------------------------------------------

message("Treatment-independent composite scoring complete.")
message("Tables saved to: ", table_dir)

message("Treg samples: ", length(sc_treg_ids))
message("CDR3β sequences scored: ", nrow(sc_final))

message("Score range: ",
        round(min(sc_final$treatment_independent_score, na.rm = TRUE), 4),
        " to ",
        round(max(sc_final$treatment_independent_score, na.rm = TRUE), 4))

message("Top 20 treatment-independent CDR3β sequences:")
print(
  sc_final %>%
    select(
      treatment_independent_rank,
      cdr3_aa,
      treatment_independent_score,
      c1_n_mice_exact,
      c2_mean_freq,
      c3_ld1_neighbor_mice,
      c4_cv,
      c5_n_mice_above_0.1pct
    ) %>%
    slice_head(n = 20),
  n = 20
)