# ==============================================================================
# 18_figure9_umap_treg_states_top_cdr3.R
# ==============================================================================
# Goal:
#   Create Figure 9 showing annotated single-cell Treg states and the UMAP
#   distribution of cells carrying top-ranked IL-2-prioritised CDR3β sequences.
#
#   Panel A:
#     UMAP coloured by annotated Treg cell state.
#
#   Panel B:
#     UMAP highlighting cells carrying:
#       - top 20 IL-2-prioritised CDR3β sequences
#       - top 21-100 IL-2-prioritised CDR3β sequences
#       - other detected CDR3β sequences
#
# Inputs:
#   data/processed/single_cell/seurat_objects/soht_annotated_tcr_filtered.rds
#   results/tables/cdr3_scores_top100_il2.csv
#
# Outputs:
#   results/figures/figure9_umap_treg_states_top_cdr3.png
#   results/figures/figure9_umap_treg_states_top_cdr3.pdf
#   results/tables/figure9_top_cdr3_status_counts.csv
#   results/tables/figure9_top_cdr3_status_by_cell_type.csv
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Load packages
# ------------------------------------------------------------------------------

library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tibble)


# ------------------------------------------------------------------------------
# 1. Define input and output paths
# ------------------------------------------------------------------------------

seurat_path <- file.path(
  "data",
  "processed",
  "single_cell",
  "seurat_objects",
  "soht_annotated_tcr_filtered.rds"
)

top100_path <- file.path(
  "results",
  "tables",
  "cdr3_scores_top100_il2.csv"
)

figure_dir <- file.path("results", "figures")
table_dir <- file.path("results", "tables")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# 2. Load input files
# ------------------------------------------------------------------------------

if (!file.exists(seurat_path)) {
  stop("TCR-filtered Seurat object not found: ", seurat_path)
}

if (!file.exists(top100_path)) {
  stop("Top-100 IL-2 score table not found: ", top100_path)
}

soht_tcr <- readRDS(seurat_path)
top100_scores <- fread(top100_path)

DefaultAssay(soht_tcr) <- "RNA"


# ------------------------------------------------------------------------------
# 3. Check required columns
# ------------------------------------------------------------------------------

required_score_cols <- c(
  "il2_prioritised_rank",
  "cdr3_aa",
  "composite_score"
)

missing_score_cols <- setdiff(required_score_cols, colnames(top100_scores))

if (length(missing_score_cols) > 0) {
  stop(
    "Missing columns in top-100 IL-2 score table: ",
    paste(missing_score_cols, collapse = ", ")
  )
}

required_metadata_cols <- c(
  "cell_type",
  "beta_cdr3"
)

missing_metadata_cols <- setdiff(
  required_metadata_cols,
  colnames(soht_tcr@meta.data)
)

if (length(missing_metadata_cols) > 0) {
  stop(
    "Missing metadata columns in TCR-filtered Seurat object: ",
    paste(missing_metadata_cols, collapse = ", ")
  )
}

if (!"umap" %in% names(soht_tcr@reductions)) {
  stop("The TCR-filtered Seurat object does not contain a UMAP reduction.")
}


# ------------------------------------------------------------------------------
# 4. Prepare top 20 and top 100 IL-2-prioritised CDR3β sequences
# ------------------------------------------------------------------------------

top100_ranked <- top100_scores %>%
  mutate(
    il2_prioritised_rank = as.numeric(il2_prioritised_rank),
    cdr3_aa = as.character(cdr3_aa),
    composite_score = as.numeric(composite_score)
  ) %>%
  arrange(il2_prioritised_rank) %>%
  slice_head(n = 100)

top20_cdr3 <- top100_ranked %>%
  slice_head(n = 20) %>%
  pull(cdr3_aa)

top100_cdr3 <- top100_ranked %>%
  pull(cdr3_aa)

if (length(top20_cdr3) == 0) {
  stop("No top-20 IL-2-prioritised CDR3β sequences were found.")
}

if (length(top100_cdr3) == 0) {
  stop("No top-100 IL-2-prioritised CDR3β sequences were found.")
}


# ------------------------------------------------------------------------------
# 5. Add top-ranked CDR3β status to metadata
# ------------------------------------------------------------------------------

soht_tcr$beta_cdr3 <- as.character(soht_tcr$beta_cdr3)

soht_tcr$top_cdr3_status <- case_when(
  soht_tcr$beta_cdr3 %in% top20_cdr3 ~ "Top 20",
  soht_tcr$beta_cdr3 %in% top100_cdr3 ~ "Top 21-100",
  !is.na(soht_tcr$beta_cdr3) & soht_tcr$beta_cdr3 != "" ~ "Other CDR3β",
  TRUE ~ "No beta CDR3"
)

soht_tcr$top_cdr3_status <- factor(
  soht_tcr$top_cdr3_status,
  levels = c(
    "Top 20",
    "Top 21-100",
    "Other CDR3β",
    "No beta CDR3"
  )
)


# ------------------------------------------------------------------------------
# 6. Set plotting colours
# ------------------------------------------------------------------------------

celltype_cols <- c(
  "Resting memory-like Tregs" = "#60DBB5",
  "Activated memory-like Tregs" = "#FF7442",
  "Activated checkpoint-high Tregs" = "#6D8FD6",
  "Cytotoxic-like Tregs" = "#E742B7",
  "CD25-high suppressive Tregs" = "#89D51A",
  "Proliferating Tregs" = "#FFE04F",
  "B/plasma contamination" = "#A366EF",
  "Myeloid/NK contamination" = "#52C9F0"
)

celltype_cols <- celltype_cols[
  names(celltype_cols) %in% unique(as.character(soht_tcr$cell_type))
]

top_cdr3_cols <- c(
  "Top 20" = "#B13A35",
  "Top 21-100" = "#F4A582",
  "Other CDR3β" = "#D9D9D9",
  "No beta CDR3" = "#F0F0F0"
)


# ------------------------------------------------------------------------------
# 7. Create Panel A: annotated Treg states
# ------------------------------------------------------------------------------

p_celltype_umap <- DimPlot(
  soht_tcr,
  reduction = "umap",
  group.by = "cell_type",
  cols = celltype_cols,
  pt.size = 0.3
) +
  labs(
    title = NULL
  ) +
  theme(
    legend.position = "right"
  )


# ------------------------------------------------------------------------------
# 8. Create Panel B: top-ranked CDR3β status
# ------------------------------------------------------------------------------

# Plot other CDR3β cells as a grey background and overlay top-100 cells so that
# highlighted cells are not hidden behind the background.

umap_df <- FetchData(
  soht_tcr,
  vars = c(
    "UMAP_1",
    "UMAP_2",
    "top_cdr3_status"
  )
) %>%
  rownames_to_column("cell_barcode") %>%
  mutate(
    top_cdr3_status = factor(
      top_cdr3_status,
      levels = c(
        "Other CDR3β",
        "No beta CDR3",
        "Top 21-100",
        "Top 20"
      )
    )
  )

background_df <- umap_df %>%
  filter(top_cdr3_status %in% c("Other CDR3β", "No beta CDR3"))

highlight_df <- umap_df %>%
  filter(top_cdr3_status %in% c("Top 21-100", "Top 20"))

p_top_cdr3_umap <- ggplot() +
  geom_point(
    data = background_df,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      color = top_cdr3_status
    ),
    size = 0.3,
    alpha = 0.55
  ) +
  geom_point(
    data = highlight_df,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      color = top_cdr3_status
    ),
    size = 0.45,
    alpha = 0.95
  ) +
  scale_color_manual(
    values = top_cdr3_cols,
    breaks = c(
      "Top 20",
      "Top 21-100",
      "Other CDR3β",
      "No beta CDR3"
    ),
    name = NULL
  ) +
  coord_equal() +
  theme_classic() +
  labs(
    x = "UMAP_1",
    y = "UMAP_2",
    title = NULL
  ) +
  theme(
    legend.position = "right"
  )


# ------------------------------------------------------------------------------
# 9. Combine panels
# ------------------------------------------------------------------------------

p_figure9 <- p_celltype_umap + p_top_cdr3_umap +
  plot_layout(
    ncol = 2,
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 16
    ),
    legend.position = "right"
  )


# ------------------------------------------------------------------------------
# 10. Export summary tables
# ------------------------------------------------------------------------------

top_cdr3_status_counts <- soht_tcr@meta.data %>%
  count(top_cdr3_status, name = "n_cells") %>%
  mutate(
    percent_cells = n_cells / sum(n_cells) * 100
  )

top_cdr3_status_by_cell_type <- soht_tcr@meta.data %>%
  count(cell_type, top_cdr3_status, name = "n_cells") %>%
  group_by(cell_type) %>%
  mutate(
    cell_type_total = sum(n_cells),
    percent_within_cell_type = n_cells / cell_type_total * 100
  ) %>%
  ungroup() %>%
  arrange(cell_type, top_cdr3_status)

fwrite(
  top_cdr3_status_counts,
  file.path(table_dir, "figure9_top_cdr3_status_counts.csv")
)

fwrite(
  top_cdr3_status_by_cell_type,
  file.path(table_dir, "figure9_top_cdr3_status_by_cell_type.csv")
)


# ------------------------------------------------------------------------------
# 11. Save Figure 9
# ------------------------------------------------------------------------------

ggsave(
  filename = file.path(figure_dir, "figure9_umap_treg_states_top_cdr3.pdf"),
  plot = p_figure9,
  width = 12,
  height = 5.5
)

ggsave(
  filename = file.path(figure_dir, "figure9_umap_treg_states_top_cdr3.png"),
  plot = p_figure9,
  width = 12,
  height = 5.5,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 12. Print quick checks
# ------------------------------------------------------------------------------

message("Figure 9 complete.")
message("Figure saved to: ", figure_dir)
message("Tables saved to: ", table_dir)

message("Top CDR3β status counts:")
print(top_cdr3_status_counts)

message("Top CDR3β status by cell type:")
print(top_cdr3_status_by_cell_type)