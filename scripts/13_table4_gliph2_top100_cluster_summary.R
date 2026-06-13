# ==============================================================================
# 13_table4_gliph2_top100_cluster_summary.R
# ==============================================================================
# Goal:
#   Create GLIPH2 cluster summary table for clusters containing top-100
#   IL-2-prioritised Treg CDR3β sequences.
#
#   Clusters are retained if they contain:
#     1. at least 3 unique sample-derived CDR3β sequences
#     2. at least 2 CDR3β sequences from the top-100 IL-2-prioritised list
#
# Inputs:
#   results/gliph2/raw_output/convergence_groups.txt
#   results/gliph2/raw_output/cluster_member_details.txt
#   results/tables/cdr3_scores_top100_il2.csv
#   data/processed/repseq_objects/RepSeqData_Curie003.RData
#
# Outputs:
#   results/tables/table4_gliph2_top100_cluster_summary.csv
#   results/tables/table4_gliph2_top100_cluster_summary.tsv
#   results/gliph2/results/top100_gliph2_wide_summary.csv
#   results/gliph2/results/top100_gliph2_wide_summary.tsv
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Load packages
# ------------------------------------------------------------------------------

library(AnalyzAIRR)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)


# ------------------------------------------------------------------------------
# 1. Define input and output paths
# ------------------------------------------------------------------------------

gliph2_raw_dir <- file.path(
  "results",
  "gliph2",
  "raw_output"
)

gliph2_results_dir <- file.path(
  "results",
  "gliph2",
  "results"
)

table_dir <- file.path(
  "results",
  "tables"
)

repseq_path <- file.path(
  "data",
  "processed",
  "repseq_objects",
  "RepSeqData_Curie003.RData"
)

gliph2_groups_path <- file.path(
  gliph2_raw_dir,
  "convergence_groups.txt"
)

gliph2_members_path <- file.path(
  gliph2_raw_dir,
  "cluster_member_details.txt"
)

top100_path <- file.path(
  table_dir,
  "cdr3_scores_top100_il2.csv"
)

dir.create(gliph2_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# 2. Check input files
# ------------------------------------------------------------------------------

required_files <- c(
  gliph2_groups_path,
  gliph2_members_path,
  top100_path,
  repseq_path
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "The following required files are missing:\n",
    paste(missing_files, collapse = "\n")
  )
}


# ------------------------------------------------------------------------------
# 3. Load GLIPH2 raw output files
# ------------------------------------------------------------------------------

g2_raw_groups <- fread(gliph2_groups_path)
g2_raw_members <- fread(gliph2_members_path)


# ------------------------------------------------------------------------------
# 4. Check required GLIPH2 columns
# ------------------------------------------------------------------------------

required_group_cols <- c(
  "tag",
  "type",
  "cluster_size",
  "unique_cdr3_sample",
  "OvE",
  "fisher.score",
  "clonal.expansion.score",
  "cdr3.length.score",
  "vgene.score",
  "total.score"
)

required_member_cols <- c(
  "tag",
  "CDR3b",
  "TRBV",
  "patient",
  "counts"
)

missing_group_cols <- setdiff(required_group_cols, colnames(g2_raw_groups))
missing_member_cols <- setdiff(required_member_cols, colnames(g2_raw_members))

if (length(missing_group_cols) > 0) {
  stop(
    "Missing columns in convergence_groups.txt: ",
    paste(missing_group_cols, collapse = ", ")
  )
}

if (length(missing_member_cols) > 0) {
  stop(
    "Missing columns in cluster_member_details.txt: ",
    paste(missing_member_cols, collapse = ", ")
  )
}


# ------------------------------------------------------------------------------
# 5. Standardise GLIPH2 cluster and member tables
# ------------------------------------------------------------------------------

cluster_stats <- g2_raw_groups %>%
  rename(
    cluster_id = tag
  ) %>%
  mutate(
    cluster_id = as.character(cluster_id),
    type = as.character(type),
    cluster_size = as.numeric(cluster_size),
    unique_cdr3_sample = as.numeric(unique_cdr3_sample),
    OvE = as.numeric(OvE),
    fisher.score = as.numeric(fisher.score),
    clonal.expansion.score = as.numeric(clonal.expansion.score),
    cdr3.length.score = as.numeric(cdr3.length.score),
    vgene.score = as.numeric(vgene.score),
    total.score = as.numeric(total.score)
  ) %>%
  select(
    cluster_id,
    type,
    cluster_size,
    unique_cdr3_sample,
    OvE,
    fisher.score,
    clonal.expansion.score,
    cdr3.length.score,
    vgene.score,
    total.score
  )

g2_members <- g2_raw_members %>%
  rename(
    cluster_id = tag,
    cdr3_aa = CDR3b,
    v_gene = TRBV,
    sample_id = patient,
    count = counts
  ) %>%
  mutate(
    cluster_id = as.character(cluster_id),
    cdr3_aa = as.character(cdr3_aa),
    v_gene = as.character(v_gene),
    sample_id = as.character(sample_id),
    count = as.numeric(count)
  )


# ------------------------------------------------------------------------------
# 6. Load Treg metadata
# ------------------------------------------------------------------------------

load(repseq_path)

if (!exists("RepSeqData_Curie003")) {
  stop("Object RepSeqData_Curie003 was not found after loading.")
}

metadata <- RepSeqData_Curie003@metaData %>%
  as.data.frame()

required_meta_cols <- c(
  "sample_id",
  "mouse",
  "cell_subset",
  "injection",
  "protein"
)

missing_meta_cols <- setdiff(required_meta_cols, colnames(metadata))

if (length(missing_meta_cols) > 0) {
  stop(
    "Missing columns in RepSeqData metadata: ",
    paste(missing_meta_cols, collapse = ", ")
  )
}

treg_metadata <- metadata %>%
  mutate(
    sample_id = as.character(sample_id),
    mouse = as.character(mouse),
    cell_subset = as.character(cell_subset),
    injection = as.character(injection),
    protein = as.character(protein),
    protein_clean = case_when(
      str_detect(protein, "MCB") ~ "mutant_IL2",
      str_detect(protein, "CB") ~ "IL2",
      str_detect(protein, "PBS") ~ "PBS",
      TRUE ~ protein
    )
  ) %>%
  filter(cell_subset == "Treg") %>%
  distinct(sample_id, .keep_all = TRUE) %>%
  select(
    sample_id,
    mouse,
    injection,
    protein_clean
  )


# ------------------------------------------------------------------------------
# 7. Annotate GLIPH2 members with metadata
# ------------------------------------------------------------------------------

g2_annot <- g2_members %>%
  left_join(
    treg_metadata,
    by = "sample_id"
  )

if (sum(is.na(g2_annot$mouse)) > 0) {
  warning(
    "Some GLIPH2 member rows did not match Treg metadata. Unmatched sample IDs: ",
    paste(unique(g2_annot$sample_id[is.na(g2_annot$mouse)]), collapse = ", ")
  )
}


# ------------------------------------------------------------------------------
# 8. Load top-100 IL-2-prioritised CDR3β sequences
# ------------------------------------------------------------------------------

top100_scores <- fread(top100_path)

required_top100_cols <- c(
  "il2_prioritised_rank",
  "cdr3_aa",
  "composite_score"
)

missing_top100_cols <- setdiff(required_top100_cols, colnames(top100_scores))

if (length(missing_top100_cols) > 0) {
  stop(
    "Missing columns in top-100 IL-2 score table: ",
    paste(missing_top100_cols, collapse = ", ")
  )
}

top100_ranked <- top100_scores %>%
  mutate(
    il2_prioritised_rank = as.numeric(il2_prioritised_rank),
    cdr3_aa = as.character(cdr3_aa),
    composite_score = as.numeric(composite_score)
  ) %>%
  arrange(il2_prioritised_rank) %>%
  slice_head(n = 100) %>%
  transmute(
    top100_rank = il2_prioritised_rank,
    cdr3_aa,
    composite_score,
    c1_n_mice_exact,
    c6_log2_fold_change,
    c7_cliffs_delta,
    mouse_list
  )

if (nrow(top100_ranked) == 0) {
  stop("No top-100 IL-2-prioritised CDR3β sequences were found.")
}


# ------------------------------------------------------------------------------
# 9. Create one-row-per-CDR3 wide GLIPH2 summary
# ------------------------------------------------------------------------------

top100_long <- top100_ranked %>%
  left_join(
    g2_annot %>%
      distinct(
        cdr3_aa,
        cluster_id
      ),
    by = "cdr3_aa"
  ) %>%
  left_join(
    cluster_stats,
    by = "cluster_id"
  )

cluster_counts_per_cdr3 <- top100_long %>%
  filter(!is.na(cluster_id)) %>%
  group_by(cdr3_aa) %>%
  summarise(
    n_clusters = n_distinct(cluster_id),
    .groups = "drop"
  )

top100_wide <- top100_long %>%
  left_join(
    cluster_counts_per_cdr3,
    by = "cdr3_aa"
  ) %>%
  mutate(
    n_clusters = replace_na(n_clusters, 0)
  ) %>%
  group_by(
    top100_rank,
    cdr3_aa,
    composite_score,
    c1_n_mice_exact,
    c6_log2_fold_change,
    c7_cliffs_delta,
    mouse_list,
    n_clusters
  ) %>%
  summarise(
    in_gliph2 = any(!is.na(cluster_id)),
    cluster_ids = ifelse(
      all(is.na(cluster_id)),
      "none",
      paste(na.omit(cluster_id), collapse = " | ")
    ),
    cluster_types = ifelse(
      all(is.na(cluster_id)),
      "none",
      paste(na.omit(type), collapse = " | ")
    ),
    cluster_sizes = ifelse(
      all(is.na(cluster_id)),
      NA_character_,
      paste(na.omit(cluster_size), collapse = " | ")
    ),
    cluster_unique_cdr3 = ifelse(
      all(is.na(cluster_id)),
      NA_character_,
      paste(na.omit(unique_cdr3_sample), collapse = " | ")
    ),
    cluster_OvE = ifelse(
      all(is.na(cluster_id)),
      NA_character_,
      paste(round(na.omit(OvE), 2), collapse = " | ")
    ),
    cluster_fisher_pval = ifelse(
      all(is.na(cluster_id)),
      NA_character_,
      paste(signif(na.omit(fisher.score), 3), collapse = " | ")
    ),
    cluster_fisher_sig = ifelse(
      all(is.na(cluster_id)),
      FALSE,
      any(fisher.score <= 0.05, na.rm = TRUE)
    ),
    cluster_total_score = ifelse(
      all(is.na(cluster_id)),
      NA_character_,
      paste(signif(na.omit(total.score), 3), collapse = " | ")
    ),
    .groups = "drop"
  ) %>%
  arrange(top100_rank)


# ------------------------------------------------------------------------------
# 10. Create cluster-level summary for Table 4
# ------------------------------------------------------------------------------

top100_cluster_members <- g2_annot %>%
  inner_join(
    top100_ranked %>%
      select(
        top100_rank,
        cdr3_aa,
        composite_score
      ),
    by = "cdr3_aa"
  ) %>%
  distinct(
    cluster_id,
    top100_rank,
    cdr3_aa,
    composite_score
  )

cluster_member_summary <- g2_annot %>%
  group_by(cluster_id) %>%
  summarise(
    n_unique_sample_cdr3 = n_distinct(cdr3_aa),
    n_samples = n_distinct(sample_id),
    n_mice = n_distinct(mouse[!is.na(mouse)]),
    mouse_list = paste(sort(unique(mouse[!is.na(mouse)])), collapse = ", "),
    all_cluster_cdr3 = paste(sort(unique(cdr3_aa)), collapse = " "),
    .groups = "drop"
  )

top100_cluster_summary <- top100_cluster_members %>%
  group_by(cluster_id) %>%
  summarise(
    n_top100_cdr3 = n_distinct(cdr3_aa),
    top100_ranks = paste(sort(unique(top100_rank)), collapse = ", "),
    top100_cdr3 = paste(
      cdr3_aa[order(top100_rank)],
      collapse = " "
    ),
    top100_cdr3_with_rank = paste(
      paste0(top100_rank[order(top100_rank)], ": ", cdr3_aa[order(top100_rank)]),
      collapse = " | "
    ),
    .groups = "drop"
  )

table4_full <- cluster_stats %>%
  left_join(
    cluster_member_summary,
    by = "cluster_id"
  ) %>%
  left_join(
    top100_cluster_summary,
    by = "cluster_id"
  ) %>%
  mutate(
    n_top100_cdr3 = replace_na(n_top100_cdr3, 0),
    passes_table4_filter = n_unique_sample_cdr3 >= 3 & n_top100_cdr3 >= 2
  ) %>%
  filter(passes_table4_filter) %>%
  arrange(
    fisher.score,
    desc(n_top100_cdr3),
    desc(n_unique_sample_cdr3)
  )

if (nrow(table4_full) == 0) {
  warning("No GLIPH2 clusters passed the Table 4 filtering criteria.")
}


# ------------------------------------------------------------------------------
# 11. Create selected-column Table 4 output
# ------------------------------------------------------------------------------

table4_selected <- table4_full %>%
  transmute(
    `GLIPH2 cluster ID` = cluster_id,
    `Cluster type` = type,
    `Cluster size` = cluster_size,
    `Unique CDR3β` = unique_cdr3_sample,
    `Fisher p-value` = fisher.score,
    `Total score` = total.score,
    `Top-100 ranks` = top100_ranks,
    `Top-100 CDR3β sequences` = top100_cdr3
  )


# ------------------------------------------------------------------------------
# 12. Save output files
# ------------------------------------------------------------------------------

fwrite(
  top100_wide,
  file.path(gliph2_results_dir, "top100_gliph2_wide_summary.tsv"),
  sep = "\t"
)

fwrite(
  top100_wide,
  file.path(gliph2_results_dir, "top100_gliph2_wide_summary.csv")
)

fwrite(
  table4_full,
  file.path(table_dir, "table4_gliph2_top100_cluster_summary_full.csv")
)

fwrite(
  table4_full,
  file.path(table_dir, "table4_gliph2_top100_cluster_summary_full.tsv"),
  sep = "\t"
)

fwrite(
  table4_selected,
  file.path(table_dir, "table4_gliph2_top100_cluster_summary.csv")
)

fwrite(
  table4_selected,
  file.path(table_dir, "table4_gliph2_top100_cluster_summary.tsv"),
  sep = "\t"
)


# ------------------------------------------------------------------------------
# 13. Print quick summary
# ------------------------------------------------------------------------------

message("Table 4 GLIPH2 summary complete.")
message("Tables saved to: ", table_dir)
message("Additional top-100 wide summary saved to: ", gliph2_results_dir)

message("Number of top-100 CDR3β sequences in at least one GLIPH2 cluster: ",
        sum(top100_wide$in_gliph2))

message("Number of GLIPH2 clusters passing Table 4 filters: ",
        nrow(table4_selected))

message("Clusters passing Table 4 filters:")
print(table4_selected)