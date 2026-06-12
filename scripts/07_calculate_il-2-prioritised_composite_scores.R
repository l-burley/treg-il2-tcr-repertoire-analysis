# ==============================================================================
# 07_calculate_il-2-prioritised_composite_scores.R
# ==============================================================================
# Goal:
#   Calculate IL-2-prioritised composite scores for shared expanded Treg CDR3β
#   sequences.
#
# Input:
#   data/processed/repseq_objects/RepSeqData_Curie003.RData
#
# Outputs:
#   results/tables/cdr3_scores_full_il2.csv
#   results/tables/cdr3_scores_full_il2.tsv
#   results/tables/cdr3_scores_top20_il2.csv
#   results/tables/cdr3_scores_top100_il2.csv
#   results/tables/cdr3_scores_top100_il2.tsv
#   results/tables/cdr3_score_components_raw_il2.csv
#   results/tables/cdr3_score_components_raw_il2.tsv
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
#   6. IL-2/PBS fold-change
#   7. IL-2/PBS Cliff's delta
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
# 6. Build complete CDR3β × sample table for IL-2/PBS criteria
# ------------------------------------------------------------------------------

# This follows the original thesis calculation:
# the table is completed directly from sc_treg_filt without collapsing duplicate
# aaCDR3-sample rows.

sc_sample_totals <- sc_treg_raw %>%
  distinct(sample_id, sample_total)

sc_complete <- sc_treg_filt %>%
  select(cdr3_aa, sample_id, count, freq) %>%
  complete(
    cdr3_aa = sc_shared_cdr3s,
    sample_id = sc_treg_ids,
    fill = list(count = 0, freq = 0)
  ) %>%
  left_join(sc_meta, by = "sample_id") %>%
  left_join(sc_sample_totals, by = "sample_id")


# ------------------------------------------------------------------------------
# 7. Criterion 1: number of mice with exact CDR3β
# ------------------------------------------------------------------------------

sc_c1 <- sc_mice_per_cdr3 %>%
  filter(cdr3_aa %in% sc_shared_cdr3s) %>%
  rename(c1_n_mice_exact = n_mice_total)


# ------------------------------------------------------------------------------
# 8. Criterion 2: mean frequency across detected mice
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
# 9. Criterion 3: LD=1 neighbour support
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
# 10. Criterion 4: frequency consistency across detected mice
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
# 11. Criterion 5: number of mice above 0.1% frequency
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
# 12. Criterion 6: IL-2/PBS mean frequency fold-change
# ------------------------------------------------------------------------------

PSEUDOCOUNT <- 1e-6

sc_c6 <- sc_complete %>%
  mutate(
    treatment = ifelse(protein_clean == "PBS", "PBS", "IL2")
  ) %>%
  group_by(cdr3_aa, treatment) %>%
  summarise(
    mean_freq = mean(freq + PSEUDOCOUNT, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = treatment,
    values_from = mean_freq,
    names_prefix = "mean_freq_"
  ) %>%
  mutate(
    mean_freq_IL2 = replace_na(mean_freq_IL2, PSEUDOCOUNT),
    mean_freq_PBS = replace_na(mean_freq_PBS, PSEUDOCOUNT),
    c6_fold_change = mean_freq_IL2 / mean_freq_PBS,
    c6_log2_fold_change = log2(c6_fold_change)
  ) %>%
  select(
    cdr3_aa,
    c6_fold_change,
    c6_log2_fold_change,
    c6_mean_freq_il2 = mean_freq_IL2,
    c6_mean_freq_pbs = mean_freq_PBS
  )


# ------------------------------------------------------------------------------
# 13. Criterion 7: IL-2/PBS Cliff's delta
# ------------------------------------------------------------------------------

sc_complete_il2 <- sc_complete %>%
  filter(protein_clean != "PBS") %>%
  select(cdr3_aa, sample_id, freq)

sc_complete_pbs <- sc_complete %>%
  filter(protein_clean == "PBS") %>%
  select(cdr3_aa, sample_id, freq)

compute_cliffs_delta <- function(cdr3) {
  il2_freqs <- sc_complete_il2$freq[sc_complete_il2$cdr3_aa == cdr3]
  pbs_freqs <- sc_complete_pbs$freq[sc_complete_pbs$cdr3_aa == cdr3]
  
  n_pairs <- length(il2_freqs) * length(pbs_freqs)
  
  if (n_pairs == 0) {
    return(NA_real_)
  }
  
  pair_scores <- outer(
    il2_freqs,
    pbs_freqs,
    FUN = function(a, b) sign(a - b)
  )
  
  sum(pair_scores) / n_pairs
}

sc_c7 <- data.frame(
  cdr3_aa = sc_shared_cdr3s,
  c7_cliffs_delta = sapply(sc_shared_cdr3s, compute_cliffs_delta),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------------------------
# 14. Combine all raw score components
# ------------------------------------------------------------------------------

sc_raw_scores <- sc_c1 %>%
  left_join(sc_c2, by = "cdr3_aa") %>%
  left_join(sc_c3, by = "cdr3_aa") %>%
  left_join(sc_c4, by = "cdr3_aa") %>%
  left_join(sc_c5, by = "cdr3_aa") %>%
  left_join(sc_c6, by = "cdr3_aa") %>%
  left_join(sc_c7, by = "cdr3_aa") %>%
  mutate(
    c3_ld1_neighbor_mice = replace_na(c3_ld1_neighbor_mice, 0),
    c3_n_neighbors = replace_na(c3_n_neighbors, 0),
    c5_n_mice_above_0.1pct = replace_na(c5_n_mice_above_0.1pct, 0),
    c7_cliffs_delta = replace_na(c7_cliffs_delta, 0)
  )


# ------------------------------------------------------------------------------
# 15. Normalize criteria and calculate composite score
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
    norm_c6 = normalize_01(c6_log2_fold_change),
    norm_c7 = normalize_01(c7_cliffs_delta),
    
    composite_score =
      0.25 * norm_c1 +
      0.20 * norm_c2 +
      0.15 * norm_c3 +
      0.10 * norm_c4 +
      0.05 * norm_c5 +
      0.15 * norm_c6 +
      0.10 * norm_c7
  ) %>%
  arrange(desc(composite_score))


# ------------------------------------------------------------------------------
# 16. Add annotation columns
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
    il2_prioritised_rank = row_number()
  ) %>%
  select(
    il2_prioritised_rank,
    cdr3_aa,
    composite_score,
    c1_n_mice_exact,
    c2_mean_freq,
    c2_max_freq,
    c3_n_neighbors,
    c3_ld1_neighbor_mice,
    c4_cv,
    c5_n_mice_above_0.1pct,
    c6_fold_change,
    c6_log2_fold_change,
    c6_mean_freq_il2,
    c6_mean_freq_pbs,
    c7_cliffs_delta,
    norm_c1,
    norm_c2,
    norm_c3,
    norm_c4,
    norm_c5,
    norm_c6,
    norm_c7,
    mouse_list,
    ld1_neighbors,
    bin_summary
  )


# ------------------------------------------------------------------------------
# 17. Save outputs
# ------------------------------------------------------------------------------

fwrite(
  sc_final,
  file.path(table_dir, "cdr3_scores_full_il2.tsv"),
  sep = "\t"
)

fwrite(
  sc_final,
  file.path(table_dir, "cdr3_scores_full_il2.csv")
)

fwrite(
  sc_final %>% slice_head(n = 20),
  file.path(table_dir, "cdr3_scores_top20_il2.csv")
)

fwrite(
  sc_final %>% slice_head(n = 100),
  file.path(table_dir, "cdr3_scores_top100_il2.tsv"),
  sep = "\t"
)

fwrite(
  sc_final %>% slice_head(n = 100),
  file.path(table_dir, "cdr3_scores_top100_il2.csv")
)

fwrite(
  sc_raw_scores,
  file.path(table_dir, "cdr3_score_components_raw_il2.tsv"),
  sep = "\t"
)

fwrite(
  sc_raw_scores,
  file.path(table_dir, "cdr3_score_components_raw_il2.csv")
)


# ------------------------------------------------------------------------------
# 18. Print quick checks
# ------------------------------------------------------------------------------

message("Composite scoring complete.")
message("Tables saved to: ", table_dir)

message("Treg samples: ", length(sc_treg_ids))
message("IL-2 samples: ", sum(sc_meta$protein_clean %in% c("IL2", "mutant_IL2")))
message("PBS samples: ", sum(sc_meta$protein_clean == "PBS"))

message("CDR3β sequences scored: ", nrow(sc_final))

message("Score range: ",
        round(min(sc_final$composite_score, na.rm = TRUE), 4),
        " to ",
        round(max(sc_final$composite_score, na.rm = TRUE), 4))

message("Top 20 IL-2-prioritised CDR3β sequences:")
print(
  sc_final %>%
    select(
      il2_prioritised_rank,
      cdr3_aa,
      composite_score,
      c1_n_mice_exact,
      c2_mean_freq,
      c3_ld1_neighbor_mice,
      c5_n_mice_above_0.1pct,
      c6_fold_change,
      c7_cliffs_delta
    ) %>%
    slice_head(n = 20),
  n = 20
)