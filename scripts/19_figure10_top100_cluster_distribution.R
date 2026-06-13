# ==============================================================================
# 19_figure10_top100_cluster_distribution.R
# ==============================================================================
# Goal:
#   Create Figure 10 showing the percentage of cells within each annotated
#   Treg state that carry one of the top-100 IL-2-prioritised CDR3β sequences.
#
# Inputs:
#   data/processed/single_cell/seurat_objects/soht_annotated_tcr_filtered.rds
#   results/tables/cdr3_scores_top100_il2.csv
#
# Outputs:
#   results/figures/figure10_top100_cluster_distribution.png
#   results/figures/figure10_top100_cluster_distribution.pdf
#   results/tables/figure10_pct_cells_top100_per_cell_type.csv
#   results/tables/figure10_top100_detection_per_cell_type.csv
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Load packages
# ------------------------------------------------------------------------------

library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
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


# ------------------------------------------------------------------------------
# 4. Prepare top-100 IL-2-prioritised CDR3β sequences
# ------------------------------------------------------------------------------

top100_ranked <- top100_scores %>%
  mutate(
    il2_prioritised_rank = as.numeric(il2_prioritised_rank),
    cdr3_aa = as.character(cdr3_aa),
    composite_score = as.numeric(composite_score)
  ) %>%
  arrange(il2_prioritised_rank) %>%
  slice_head(n = 100) %>%
  select(
    rank = il2_prioritised_rank,
    cdr3_aa,
    composite_score
  )

top100_cdr3 <- top100_ranked$cdr3_aa

if (length(top100_cdr3) == 0) {
  stop("No top-100 IL-2-prioritised CDR3β sequences were found.")
}


# ------------------------------------------------------------------------------
# 5. Prepare single-cell metadata
# ------------------------------------------------------------------------------

sc_meta <- soht_tcr@meta.data %>%
  rownames_to_column("cell_barcode") %>%
  mutate(
    cell_type = as.character(cell_type),
    beta_cdr3 = as.character(beta_cdr3)
  ) %>%
  filter(
    !is.na(beta_cdr3),
    beta_cdr3 != "",
    beta_cdr3 != "NA"
  ) %>%
  mutate(
    in_top100 = beta_cdr3 %in% top100_cdr3
  )

if (nrow(sc_meta) == 0) {
  stop("No cells with beta CDR3 metadata were found.")
}


# ------------------------------------------------------------------------------
# 6. Calculate percentage of cells carrying a top-100 CDR3β by cell type
# ------------------------------------------------------------------------------

pct_cells_top100 <- sc_meta %>%
  group_by(cell_type) %>%
  summarise(
    total_cells = n(),
    cells_with_top100 = sum(in_top100, na.rm = TRUE),
    pct_cells_with_top100 = round(cells_with_top100 / total_cells * 100, 2),
    .groups = "drop"
  ) %>%
  filter(
    !cell_type %in% c(
      "B/plasma contamination",
      "Myeloid/NK contamination"
    )
  ) %>%
  arrange(desc(pct_cells_with_top100))

top100_detection_per_cell_type <- sc_meta %>%
  filter(in_top100) %>%
  group_by(cell_type) %>%
  summarise(
    n_unique_top100_detected = n_distinct(beta_cdr3),
    n_cells_with_top100 = n(),
    .groups = "drop"
  ) %>%
  right_join(
    pct_cells_top100 %>%
      select(cell_type, total_cells),
    by = "cell_type"
  ) %>%
  mutate(
    n_unique_top100_detected = ifelse(
      is.na(n_unique_top100_detected),
      0,
      n_unique_top100_detected
    ),
    n_cells_with_top100 = ifelse(
      is.na(n_cells_with_top100),
      0,
      n_cells_with_top100
    ),
    pct_top100_sequences_detected = round(n_unique_top100_detected / 100 * 100, 1)
  ) %>%
  arrange(desc(n_unique_top100_detected))


# ------------------------------------------------------------------------------
# 7. Save summary tables
# ------------------------------------------------------------------------------

fwrite(
  pct_cells_top100,
  file.path(table_dir, "figure10_pct_cells_top100_per_cell_type.csv")
)

fwrite(
  top100_detection_per_cell_type,
  file.path(table_dir, "figure10_top100_detection_per_cell_type.csv")
)


# ------------------------------------------------------------------------------
# 8. Create Figure 10
# ------------------------------------------------------------------------------

p_figure10 <- ggplot(
  pct_cells_top100,
  aes(
    x = reorder(cell_type, pct_cells_with_top100),
    y = pct_cells_with_top100
  )
) +
  geom_col(
    fill = "#24325F",
    color = "black",
    width = 0.7
  ) +
  geom_text(
    aes(label = paste0(pct_cells_with_top100, "%")),
    hjust = -0.2,
    size = 3.5
  ) +
  coord_flip() +
  theme_classic(base_size = 11) +
  labs(
    x = "",
    y = "% of cells carrying a top-100 CDR3β sequence"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15))
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_blank(),
    panel.grid = element_blank()
  )


# ------------------------------------------------------------------------------
# 9. Save Figure 10
# ------------------------------------------------------------------------------

ggsave(
  filename = file.path(figure_dir, "figure10_top100_cluster_distribution.pdf"),
  plot = p_figure10,
  width = 8,
  height = 4.5
)

ggsave(
  filename = file.path(figure_dir, "figure10_top100_cluster_distribution.png"),
  plot = p_figure10,
  width = 8,
  height = 4.5,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 10. Print quick checks
# ------------------------------------------------------------------------------

message("Figure 10 complete.")
message("Figure saved to: ", figure_dir)
message("Tables saved to: ", table_dir)

message("Percentage of cells carrying top-100 CDR3β by cell type:")
print(pct_cells_top100)

message("Top-100 CDR3β detection per cell type:")
print(top100_detection_per_cell_type)