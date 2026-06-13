# ==============================================================================
# 14_figure6_gliph2_network_and_logos.R
# ==============================================================================
# Goal:
#   Create Figure 6-related outputs from retained GLIPH2 clusters containing
#   top-100 IL-2-prioritised CDR3β sequences.
#
#   Part A:
#     Export Cytoscape-ready node and edge tables for retained GLIPH2 clusters.
#     Each retained GLIPH2 cluster is represented as one connected component.
#     Nodes represent CDR3β sequences and edges connect sequences belonging to
#     the same GLIPH2 cluster.
#
#   Part B:
#     Generate sequence logo plots for the selected GLIPH2 clusters shown in
#     Figure 6B. The selected clusters are manually listed to make the output
#     reproducible and easy to understand.
#
# Inputs:
#   results/tables/table4_gliph2_top100_cluster_summary_full.csv
#   data/processed/repseq_objects/RepSeqData_Curie003.RData
#
# Outputs:
#   results/gliph2/figures/figure6A_gliph2_cytoscape_nodes.csv
#   results/gliph2/figures/figure6A_gliph2_cytoscape_edges.csv
#   results/gliph2/figures/figure6A_gliph2_cluster_summary.csv
#   results/gliph2/figures/figure6B_selected_gliph2_logo_grid.pdf
#   results/gliph2/figures/figure6B_selected_gliph2_logo_grid.png
#   results/gliph2/figures/figure6B_selected_gliph2_logo_length_summary.csv
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Load packages
# ------------------------------------------------------------------------------

library(AnalyzAIRR)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(ggseqlogo)
library(patchwork)


# ------------------------------------------------------------------------------
# 1. Define input and output paths
# ------------------------------------------------------------------------------

table4_full_path <- file.path(
  "results",
  "tables",
  "table4_gliph2_top100_cluster_summary_full.csv"
)

repseq_path <- file.path(
  "data",
  "processed",
  "repseq_objects",
  "RepSeqData_Curie003.RData"
)

figure_dir <- file.path(
  "results",
  "gliph2",
  "figures"
)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# 2. Analysis settings
# ------------------------------------------------------------------------------

subset_to_use <- "Treg"

presence_freq_threshold <- 0

make_complete_cluster_edges <- TRUE

# GLIPH2 clusters selected for the Figure 6B sequence logo panel.
# These are manually listed so the logo grid exactly matches the final figure.

selected_logo_cluster_ids <- c(
  "SLAGG%YE_ADEFGHLRSTV",
  "S%SYE_ADEFLRSVWY",
  "S%GGQNT_ADEFGILPRS",
  "S%TANSD_DGILQRSV",
  "SL%GQNT_ADGPQRST",
  "S%GSQNT_AGILPQRSW",
  "S%DWGYE_PQRSW",
  "S%DWGNYAE_ALPQRS"
)


# ------------------------------------------------------------------------------
# 3. Check input files
# ------------------------------------------------------------------------------

required_files <- c(
  table4_full_path,
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
# 4. Load retained GLIPH2 cluster table
# ------------------------------------------------------------------------------

table4_full <- fread(table4_full_path)

required_table_cols <- c(
  "cluster_id",
  "type",
  "cluster_size",
  "unique_cdr3_sample",
  "fisher.score",
  "total.score",
  "n_top100_cdr3",
  "top100_ranks",
  "top100_cdr3",
  "all_cluster_cdr3"
)

missing_table_cols <- setdiff(required_table_cols, colnames(table4_full))

if (length(missing_table_cols) > 0) {
  stop(
    "Missing columns in table4_gliph2_top100_cluster_summary_full.csv: ",
    paste(missing_table_cols, collapse = ", ")
  )
}

gliph2_cluster_summary <- table4_full %>%
  mutate(
    cluster_id = as.character(cluster_id),
    cluster_type = as.character(type),
    cluster_size_display = as.numeric(cluster_size),
    unique_cdr3_ref = as.numeric(unique_cdr3_sample),
    cluster_fisher_pval = as.numeric(fisher.score),
    cluster_total_score = as.numeric(total.score),
    n_top100_in_cluster = as.numeric(n_top100_cdr3),
    top100_ranks_in_cluster = as.character(top100_ranks),
    top100_cdr3s_in_cluster = as.character(top100_cdr3),
    all_sample_cdr3s_in_cluster = as.character(all_cluster_cdr3),
    fisher_display = formatC(cluster_fisher_pval, format = "e", digits = 2),
    total_score_display = signif(cluster_total_score, 3)
  ) %>%
  arrange(
    cluster_fisher_pval,
    desc(n_top100_in_cluster),
    desc(cluster_size_display)
  )

if (nrow(gliph2_cluster_summary) == 0) {
  stop("No retained GLIPH2 clusters were found in the Table 4 full summary.")
}


# ------------------------------------------------------------------------------
# 5. Load RepSeqData object for node attributes
# ------------------------------------------------------------------------------

load(repseq_path)

if (!exists("RepSeqData_Curie003")) {
  stop("Object RepSeqData_Curie003 was not found after loading.")
}

rep_assay <- RepSeqData_Curie003@assayData %>%
  as.data.frame()

rep_meta <- RepSeqData_Curie003@metaData %>%
  as.data.frame()

required_assay_cols <- c(
  "sample_id",
  "aaCDR3",
  "count"
)

required_meta_cols <- c(
  "sample_id",
  "mouse",
  "cell_subset",
  "injection",
  "protein"
)

missing_assay_cols <- setdiff(required_assay_cols, colnames(rep_assay))
missing_meta_cols <- setdiff(required_meta_cols, colnames(rep_meta))

if (length(missing_assay_cols) > 0) {
  stop(
    "Missing RepSeqData assayData columns: ",
    paste(missing_assay_cols, collapse = ", ")
  )
}

if (length(missing_meta_cols) > 0) {
  stop(
    "Missing RepSeqData metaData columns: ",
    paste(missing_meta_cols, collapse = ", ")
  )
}

rep_meta <- rep_meta %>%
  mutate(
    sample_id = as.character(sample_id),
    mouse = as.character(mouse),
    cell_subset = as.character(cell_subset),
    injection = as.character(injection),
    protein = str_trim(as.character(protein))
  )

if (!is.null(subset_to_use)) {
  keep_samples <- rep_meta %>%
    filter(cell_subset == subset_to_use) %>%
    pull(sample_id)
  
  rep_assay <- rep_assay %>%
    filter(sample_id %in% keep_samples)
  
  rep_meta <- rep_meta %>%
    filter(sample_id %in% keep_samples)
}

sample_totals <- rep_assay %>%
  mutate(
    sample_id = as.character(sample_id),
    count = as.numeric(count)
  ) %>%
  group_by(sample_id) %>%
  summarise(
    sample_total_count = sum(count, na.rm = TRUE),
    .groups = "drop"
  )

rep_long <- rep_assay %>%
  mutate(
    sample_id = as.character(sample_id),
    aaCDR3 = as.character(aaCDR3),
    count = as.numeric(count)
  ) %>%
  left_join(
    sample_totals,
    by = "sample_id"
  ) %>%
  left_join(
    rep_meta %>%
      select(
        sample_id,
        mouse,
        cell_subset,
        injection,
        protein
      ),
    by = "sample_id"
  ) %>%
  mutate(
    freq = count / sample_total_count
  )


# ------------------------------------------------------------------------------
# 6. Helper functions
# ------------------------------------------------------------------------------

split_cdr3_list <- function(x) {
  if (is.na(x) || x == "") {
    character(0)
  } else if (str_detect(x, "\\|")) {
    str_split(x, "\\s*\\|\\s*")[[1]] %>%
      str_trim() %>%
      discard(~ .x == "")
  } else {
    str_split(x, "\\s+")[[1]] %>%
      str_trim() %>%
      discard(~ .x == "")
  }
}

split_rank_list <- function(x) {
  if (is.na(x) || x == "") {
    character(0)
  } else if (str_detect(x, "\\|")) {
    str_split(x, "\\s*\\|\\s*")[[1]] %>%
      str_trim() %>%
      discard(~ .x == "")
  } else if (str_detect(x, ",")) {
    str_split(x, "\\s*,\\s*")[[1]] %>%
      str_trim() %>%
      discard(~ .x == "")
  } else {
    str_split(x, "\\s+")[[1]] %>%
      str_trim() %>%
      discard(~ .x == "")
  }
}

prepare_logo_input <- function(cluster_table) {
  cluster_table %>%
    mutate(
      cluster_fisher_pval = as.numeric(cluster_fisher_pval),
      cluster_size_display = as.numeric(cluster_size_display),
      n_top100_in_cluster = as.numeric(n_top100_in_cluster),
      all_cdr3_vector = map(all_sample_cdr3s_in_cluster, split_cdr3_list),
      n_total_cdr3 = map_int(all_cdr3_vector, length),
      lengths_present = map_chr(
        all_cdr3_vector,
        ~ paste(sort(unique(nchar(.x))), collapse = ", ")
      ),
      most_common_length = map_int(
        all_cdr3_vector,
        ~ as.integer(names(sort(table(nchar(.x)), decreasing = TRUE)[1]))
      ),
      cdr3_for_logo = map2(
        all_cdr3_vector,
        most_common_length,
        ~ .x[nchar(.x) == .y]
      ),
      n_logo_cdr3 = map_int(cdr3_for_logo, length),
      fisher_display = formatC(
        cluster_fisher_pval,
        format = "e",
        digits = 2
      )
    )
}

make_logo_plot <- function(seqs,
                           cluster_id,
                           ranks,
                           n_logo,
                           n_total,
                           seq_length) {
  
  ggseqlogo(seqs, method = "bits") +
    labs(
      title = as.character(cluster_id),
      subtitle = paste0(
        "Ranks: ", ranks,
        " | ", n_logo, "/", n_total,
        " seqs, length ", seq_length
      ),
      x = NULL,
      y = NULL
    ) +
    theme_classic(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", size = 7),
      plot.subtitle = element_text(size = 6),
      axis.text.x = element_text(size = 5),
      axis.text.y = element_text(size = 5),
      axis.title = element_blank(),
      legend.position = "none",
      plot.margin = margin(3, 3, 3, 3)
    )
}


# ------------------------------------------------------------------------------
# 7. Expand retained GLIPH2 clusters to one row per cluster-CDR3β sequence
# ------------------------------------------------------------------------------

gliph_nodes_long <- gliph2_cluster_summary %>%
  mutate(
    cluster_row_id = row_number(),
    all_cdr3_list = map(all_sample_cdr3s_in_cluster, split_cdr3_list),
    top100_cdr3_list = map(top100_cdr3s_in_cluster, split_cdr3_list),
    top100_rank_list = map(top100_ranks_in_cluster, split_rank_list)
  ) %>%
  select(
    cluster_row_id,
    cluster_id,
    cluster_type,
    cluster_size_display,
    unique_cdr3_ref,
    cluster_fisher_pval,
    fisher_display,
    cluster_total_score,
    total_score_display,
    n_top100_in_cluster,
    all_cdr3_list,
    top100_cdr3_list,
    top100_rank_list
  ) %>%
  unnest_longer(
    all_cdr3_list,
    values_to = "cdr3"
  ) %>%
  mutate(
    cdr3 = str_trim(cdr3)
  ) %>%
  filter(
    !is.na(cdr3),
    cdr3 != ""
  )


# ------------------------------------------------------------------------------
# 8. Add top-100 annotation to nodes
# ------------------------------------------------------------------------------

top100_lookup <- gliph2_cluster_summary %>%
  mutate(
    cluster_row_id = row_number(),
    top100_cdr3_list = map(top100_cdr3s_in_cluster, split_cdr3_list),
    top100_rank_list = map(top100_ranks_in_cluster, split_rank_list)
  ) %>%
  select(
    cluster_row_id,
    top100_cdr3_list,
    top100_rank_list
  ) %>%
  mutate(
    top100_lookup_tbl = map2(
      top100_cdr3_list,
      top100_rank_list,
      ~ tibble(
        cdr3 = .x,
        top100_rank_in_cluster = .y
      )
    )
  ) %>%
  select(
    cluster_row_id,
    top100_lookup_tbl
  ) %>%
  unnest(top100_lookup_tbl)

gliph_nodes_long <- gliph_nodes_long %>%
  left_join(
    top100_lookup,
    by = c("cluster_row_id", "cdr3")
  ) %>%
  mutate(
    is_top100 = !is.na(top100_rank_in_cluster),
    top100_rank_in_cluster = ifelse(
      is.na(top100_rank_in_cluster),
      "",
      top100_rank_in_cluster
    )
  )


# ------------------------------------------------------------------------------
# 9. Summarise RepSeq contribution per CDR3β sequence
# ------------------------------------------------------------------------------

cdr3_contribution <- rep_long %>%
  filter(aaCDR3 %in% unique(gliph_nodes_long$cdr3)) %>%
  group_by(
    aaCDR3,
    sample_id,
    mouse,
    cell_subset,
    injection,
    protein
  ) %>%
  summarise(
    count = sum(count, na.rm = TRUE),
    freq = sum(freq, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(freq > presence_freq_threshold) %>%
  group_by(aaCDR3) %>%
  summarise(
    contributing_mice = paste(sort(unique(mouse)), collapse = " | "),
    contributing_samples = paste(sort(unique(sample_id)), collapse = " | "),
    contributing_injections = paste(sort(unique(injection)), collapse = " | "),
    contributing_proteins = paste(sort(unique(protein)), collapse = " | "),
    n_contributing_mice = n_distinct(mouse),
    n_contributing_samples = n_distinct(sample_id),
    total_count_across_samples = sum(count, na.rm = TRUE),
    mean_freq_detected_samples = mean(freq, na.rm = TRUE),
    max_freq_detected_samples = max(freq, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(
    cdr3 = aaCDR3
  )


# ------------------------------------------------------------------------------
# 10. Build Cytoscape node table for Figure 6A
# ------------------------------------------------------------------------------

nodes <- gliph_nodes_long %>%
  left_join(
    cdr3_contribution,
    by = "cdr3"
  ) %>%
  mutate(
    node_id = paste(cluster_id, cdr3, sep = "__"),
    node_label = cdr3,
    contributing_mice = ifelse(is.na(contributing_mice), "", contributing_mice),
    contributing_samples = ifelse(is.na(contributing_samples), "", contributing_samples),
    contributing_injections = ifelse(is.na(contributing_injections), "", contributing_injections),
    contributing_proteins = ifelse(is.na(contributing_proteins), "", contributing_proteins),
    n_contributing_mice = ifelse(is.na(n_contributing_mice), 0, n_contributing_mice),
    n_contributing_samples = ifelse(is.na(n_contributing_samples), 0, n_contributing_samples),
    total_count_across_samples = ifelse(is.na(total_count_across_samples), 0, total_count_across_samples),
    mean_freq_detected_samples = ifelse(is.na(mean_freq_detected_samples), 0, mean_freq_detected_samples),
    max_freq_detected_samples = ifelse(is.na(max_freq_detected_samples), 0, max_freq_detected_samples)
  ) %>%
  select(
    node_id,
    node_label,
    cdr3,
    cluster_id,
    cluster_row_id,
    cluster_type,
    cluster_size_display,
    unique_cdr3_ref,
    cluster_fisher_pval,
    fisher_display,
    cluster_total_score,
    total_score_display,
    n_top100_in_cluster,
    is_top100,
    top100_rank_in_cluster,
    contributing_mice,
    n_contributing_mice,
    contributing_samples,
    n_contributing_samples,
    contributing_injections,
    contributing_proteins,
    total_count_across_samples,
    mean_freq_detected_samples,
    max_freq_detected_samples
  ) %>%
  distinct(
    node_id,
    .keep_all = TRUE
  )


# ------------------------------------------------------------------------------
# 11. Build Cytoscape edge table for Figure 6A
# ------------------------------------------------------------------------------

make_edges_for_cluster <- function(df) {
  node_ids <- df$node_id %>%
    as.character() %>%
    str_replace_all('"', "")
  
  if (length(node_ids) < 2) {
    return(tibble())
  }
  
  pairs <- combn(node_ids, 2, simplify = FALSE)
  
  tibble(
    source = map_chr(pairs, 1),
    target = map_chr(pairs, 2),
    interaction = "same_gliph2_cluster",
    cluster_id = unique(df$cluster_id),
    cluster_row_id = unique(df$cluster_row_id),
    edge_weight = 1
  )
}

if (make_complete_cluster_edges) {
  edges <- nodes %>%
    group_by(
      cluster_row_id,
      cluster_id
    ) %>%
    group_split() %>%
    map_dfr(make_edges_for_cluster)
} else {
  stop("Only complete within-cluster edges are implemented.")
}


# ------------------------------------------------------------------------------
# 12. Create network cluster summary
# ------------------------------------------------------------------------------

cluster_summary_for_network <- nodes %>%
  group_by(
    cluster_id,
    cluster_row_id,
    cluster_type,
    cluster_size_display,
    unique_cdr3_ref,
    cluster_fisher_pval,
    fisher_display,
    cluster_total_score,
    total_score_display
  ) %>%
  summarise(
    n_nodes = n(),
    n_top100_nodes = sum(is_top100),
    top100_ranks = paste(
      top100_rank_in_cluster[top100_rank_in_cluster != ""],
      collapse = " | "
    ),
    top100_cdr3s = paste(cdr3[is_top100], collapse = " | "),
    all_cdr3s = paste(cdr3, collapse = " | "),
    contributing_mice_all = paste(
      sort(unique(unlist(str_split(
        contributing_mice[contributing_mice != ""],
        "\\s*\\|\\s*"
      )))),
      collapse = " | "
    ),
    n_contributing_mice_all = n_distinct(
      unlist(str_split(
        contributing_mice[contributing_mice != ""],
        "\\s*\\|\\s*"
      ))
    ),
    .groups = "drop"
  )


# ------------------------------------------------------------------------------
# 13. Save Cytoscape files for Figure 6A
# ------------------------------------------------------------------------------

fwrite(
  nodes,
  file.path(figure_dir, "figure6A_gliph2_cytoscape_nodes.csv")
)

fwrite(
  edges,
  file.path(figure_dir, "figure6A_gliph2_cytoscape_edges.csv")
)

fwrite(
  cluster_summary_for_network,
  file.path(figure_dir, "figure6A_gliph2_cluster_summary.csv")
)


# ------------------------------------------------------------------------------
# 14. Select clusters for Figure 6B sequence logos
# ------------------------------------------------------------------------------

selected_logo_clusters <- gliph2_cluster_summary %>%
  filter(cluster_id %in% selected_logo_cluster_ids) %>%
  mutate(
    cluster_id = factor(cluster_id, levels = selected_logo_cluster_ids)
  ) %>%
  arrange(cluster_id) %>%
  mutate(
    cluster_id = as.character(cluster_id)
  )

missing_selected_clusters <- setdiff(
  selected_logo_cluster_ids,
  selected_logo_clusters$cluster_id
)

if (length(missing_selected_clusters) > 0) {
  warning(
    "The following selected logo clusters were not found in the GLIPH2 summary: ",
    paste(missing_selected_clusters, collapse = ", ")
  )
}

selected_logo_input <- selected_logo_clusters %>%
  prepare_logo_input()

if (nrow(selected_logo_input) == 0) {
  stop("No selected GLIPH2 clusters were available for logo plotting.")
}


# ------------------------------------------------------------------------------
# 15. Save sequence-length summary for Figure 6B
# ------------------------------------------------------------------------------

selected_logo_length_summary <- selected_logo_input %>%
  select(
    cluster_id,
    cluster_size_display,
    n_top100_in_cluster,
    top100_ranks_in_cluster,
    n_total_cdr3,
    lengths_present,
    most_common_length,
    n_logo_cdr3,
    fisher_display
  )

fwrite(
  selected_logo_length_summary,
  file.path(figure_dir, "figure6B_selected_gliph2_logo_length_summary.csv")
)


# ------------------------------------------------------------------------------
# 16. Generate Figure 6B logo grid
# ------------------------------------------------------------------------------

selected_logo_plots <- pmap(
  list(
    selected_logo_input$cdr3_for_logo,
    selected_logo_input$cluster_id,
    selected_logo_input$top100_ranks_in_cluster,
    selected_logo_input$n_logo_cdr3,
    selected_logo_input$n_total_cdr3,
    selected_logo_input$most_common_length
  ),
  make_logo_plot
)

selected_logo_grid <- wrap_plots(
  selected_logo_plots,
  ncol = 2
) +
  plot_annotation(
    title = "Sequence logos for selected GLIPH2 clusters",
    subtitle = "Logos were generated using the most common CDR3β length subset within each retained cluster."
  ) &
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10)
  )


# ------------------------------------------------------------------------------
# 17. Save Figure 6B logo grid
# ------------------------------------------------------------------------------

ggsave(
  filename = file.path(figure_dir, "figure6B_selected_gliph2_logo_grid.pdf"),
  plot = selected_logo_grid,
  width = 8,
  height = 12
)

ggsave(
  filename = file.path(figure_dir, "figure6B_selected_gliph2_logo_grid.png"),
  plot = selected_logo_grid,
  width = 8,
  height = 12,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 18. Print quick checks
# ------------------------------------------------------------------------------

message("Figure 6 GLIPH2 network and logo outputs complete.")
message("Outputs saved to: ", figure_dir)

message("Network summary:")
message("Retained clusters: ", n_distinct(nodes$cluster_id))
message("Cluster-specific nodes: ", nrow(nodes))
message("Edges: ", nrow(edges))
message("Top-100 nodes: ", sum(nodes$is_top100))
message("Original unique CDR3β sequences: ", n_distinct(nodes$cdr3))

overlapping_cdr3s <- nodes %>%
  distinct(
    cdr3,
    cluster_id
  ) %>%
  count(
    cdr3,
    name = "n_clusters"
  ) %>%
  filter(n_clusters > 1) %>%
  arrange(
    desc(n_clusters),
    cdr3
  )

message("CDR3β sequences present in more than one retained GLIPH2 cluster:")
print(overlapping_cdr3s)

message("Selected logo clusters:")
print(
  selected_logo_length_summary %>%
    select(
      cluster_id,
      cluster_size_display,
      n_top100_in_cluster,
      top100_ranks_in_cluster,
      n_total_cdr3,
      lengths_present,
      most_common_length,
      n_logo_cdr3
    )
)