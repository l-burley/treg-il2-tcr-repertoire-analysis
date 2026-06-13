# ==============================================================================
# 16_figure7_chain_pairing_top20.R
# ==============================================================================
# Goal:
#   Create Figure 7 showing TCR alpha/beta chain-pairing configurations among
#   single-cell barcodes matching the top 20 IL-2-prioritised CDR3β sequences.
#
#   Each bar represents one top-20 CDR3β sequence. Bar segments show the
#   percentage of matched single-cell barcodes with each alpha/beta chain
#   configuration. The number above each bar indicates the total number of
#   matched barcodes for that CDR3β sequence.
#
# Inputs:
#   data/processed/single_cell/seurat_objects/soht_annotated_all_cells.rds
#   results/tables/cdr3_scores_top100_il2.csv
#
# Outputs:
#   results/figures/figure7_chain_pairing_top20.png
#   results/figures/figure7_chain_pairing_top20.pdf
#   results/tables/figure7_chain_pairing_top20_summary.csv
#   results/tables/figure7_chain_pairing_top20_long.csv
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Load packages
# ------------------------------------------------------------------------------

library(Seurat)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
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
  "soht_annotated_all_cells.rds"
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
  stop("Annotated Seurat object not found: ", seurat_path)
}

if (!file.exists(top100_path)) {
  stop("Top-100 IL-2 score table not found: ", top100_path)
}

soht <- readRDS(seurat_path)

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
  "tra_cdr3",
  "trb_cdr3"
)

missing_metadata_cols <- setdiff(
  required_metadata_cols,
  colnames(soht@meta.data)
)

if (length(missing_metadata_cols) > 0) {
  stop(
    "Missing metadata columns in Seurat object: ",
    paste(missing_metadata_cols, collapse = ", ")
  )
}


# ------------------------------------------------------------------------------
# 4. Prepare top 20 IL-2-prioritised CDR3β sequences
# ------------------------------------------------------------------------------

top20_rank_table <- top100_scores %>%
  mutate(
    il2_prioritised_rank = as.numeric(il2_prioritised_rank),
    cdr3_aa = as.character(cdr3_aa),
    composite_score = as.numeric(composite_score)
  ) %>%
  arrange(il2_prioritised_rank) %>%
  slice_head(n = 20) %>%
  transmute(
    rank = il2_prioritised_rank,
    cdr3_aa,
    composite_score,
    cdr3_rank_label = paste0(rank, ". ", cdr3_aa)
  )

top20_cdr3 <- top20_rank_table$cdr3_aa

if (length(top20_cdr3) == 0) {
  stop("No top-20 CDR3β sequences were found.")
}


# ------------------------------------------------------------------------------
# 5. Helper functions for TCR parsing
# ------------------------------------------------------------------------------

split_cdr3 <- function(x) {
  x <- as.character(x)
  
  if (is.na(x) || x == "" || x == "NA") {
    return(character(0))
  }
  
  parts <- unlist(strsplit(x, split = ";|,|\\|"))
  parts <- trimws(parts)
  parts <- parts[parts != "" & parts != "NA"]
  
  unique(parts)
}

count_cdr3 <- function(x) {
  sapply(x, function(y) length(split_cdr3(y)))
}

contains_cdr3 <- function(cdr3_string, query_cdr3) {
  query_cdr3 %in% split_cdr3(cdr3_string)
}


# ------------------------------------------------------------------------------
# 6. Prepare single-cell metadata
# ------------------------------------------------------------------------------

sc_meta <- soht@meta.data %>%
  rownames_to_column("cell_barcode") %>%
  mutate(
    tra_cdr3 = as.character(tra_cdr3),
    trb_cdr3 = as.character(trb_cdr3)
  )

# Recalculate n_alpha and n_beta if they were not already stored in the object.
if (!"n_alpha" %in% colnames(sc_meta)) {
  sc_meta$n_alpha <- unname(count_cdr3(sc_meta$tra_cdr3))
}

if (!"n_beta" %in% colnames(sc_meta)) {
  sc_meta$n_beta <- unname(count_cdr3(sc_meta$trb_cdr3))
}

sc_meta <- sc_meta %>%
  mutate(
    n_alpha = as.numeric(n_alpha),
    n_beta = as.numeric(n_beta),
    chain_combo = paste0(
      "a",
      pmin(n_alpha, 2),
      "b",
      pmin(n_beta, 2)
    )
  )


# ------------------------------------------------------------------------------
# 7. Match top-20 CDR3β sequences to single-cell barcodes
# ------------------------------------------------------------------------------

matched_cells <- lapply(top20_cdr3, function(query_cdr3) {
  
  sc_meta %>%
    filter(
      vapply(
        trb_cdr3,
        contains_cdr3,
        logical(1),
        query_cdr3 = query_cdr3
      )
    ) %>%
    mutate(
      cdr3_aa = query_cdr3
    )
  
}) %>%
  bind_rows() %>%
  left_join(
    top20_rank_table,
    by = "cdr3_aa"
  )

if (nrow(matched_cells) == 0) {
  stop("No single-cell barcodes matched the top-20 CDR3β sequences.")
}


# ------------------------------------------------------------------------------
# 8. Summarise chain-pairing configurations
# ------------------------------------------------------------------------------

combo_levels <- c(
  "a0b1",
  "a0b2",
  "a1b1",
  "a1b2",
  "a2b1",
  "a2b2"
)

combo_labels <- c(
  "a0b1" = "0α 1β",
  "a0b2" = "0α 2β",
  "a1b1" = "1α 1β",
  "a1b2" = "1α 2β",
  "a2b1" = "2α 1β",
  "a2b2" = "2α 2β"
)

combo_colours <- c(
  "a0b1" = "#E07B54",
  "a0b2" = "#C0392B",
  "a1b1" = "#2E86AB",
  "a1b2" = "#5BA4CF",
  "a2b1" = "#27AE60",
  "a2b2" = "#8E44AD"
)

cdr3_order <- top20_rank_table %>%
  arrange(rank) %>%
  pull(cdr3_rank_label)

chain_pairing_long <- matched_cells %>%
  filter(chain_combo %in% combo_levels) %>%
  count(
    rank,
    cdr3_aa,
    cdr3_rank_label,
    composite_score,
    chain_combo,
    name = "n_barcodes_combo"
  ) %>%
  complete(
    rank,
    cdr3_aa,
    cdr3_rank_label,
    composite_score,
    chain_combo = combo_levels,
    fill = list(n_barcodes_combo = 0)
  ) %>%
  group_by(
    rank,
    cdr3_aa,
    cdr3_rank_label,
    composite_score
  ) %>%
  mutate(
    n_barcodes = sum(n_barcodes_combo),
    pct = ifelse(
      n_barcodes > 0,
      n_barcodes_combo / n_barcodes * 100,
      0
    )
  ) %>%
  ungroup() %>%
  mutate(
    chain_combo = factor(chain_combo, levels = combo_levels),
    cdr3_rank_label = factor(cdr3_rank_label, levels = cdr3_order)
  ) %>%
  arrange(rank, chain_combo)

chain_pairing_summary <- chain_pairing_long %>%
  select(
    rank,
    cdr3_aa,
    cdr3_rank_label,
    composite_score,
    n_barcodes,
    chain_combo,
    n_barcodes_combo,
    pct
  )

chain_pairing_wide <- chain_pairing_long %>%
  select(
    rank,
    cdr3_aa,
    cdr3_rank_label,
    composite_score,
    n_barcodes,
    chain_combo,
    pct
  ) %>%
  mutate(
    chain_combo = paste0(chain_combo, "_pct")
  ) %>%
  pivot_wider(
    names_from = chain_combo,
    values_from = pct,
    values_fill = 0
  ) %>%
  arrange(rank)


# ------------------------------------------------------------------------------
# 9. Create labels for total barcode counts
# ------------------------------------------------------------------------------

label_df <- chain_pairing_long %>%
  distinct(
    rank,
    cdr3_rank_label,
    n_barcodes
  ) %>%
  mutate(
    cdr3_rank_label = factor(cdr3_rank_label, levels = cdr3_order),
    label = as.character(n_barcodes)
  )


# ------------------------------------------------------------------------------
# 10. Plot Figure 7
# ------------------------------------------------------------------------------

p_figure7 <- ggplot(
  chain_pairing_long,
  aes(
    x = cdr3_rank_label,
    y = pct,
    fill = chain_combo
  )
) +
  geom_col(
    width = 0.8,
    colour = "white",
    linewidth = 0.25
  ) +
  geom_text(
    data = label_df,
    aes(
      x = cdr3_rank_label,
      y = 104,
      label = label
    ),
    inherit.aes = FALSE,
    angle = 45,
    hjust = 0,
    vjust = 0,
    size = 2.8
  ) +
  scale_fill_manual(
    values = combo_colours,
    labels = combo_labels,
    name = "Chain\ncombination"
  ) +
  scale_y_continuous(
    limits = c(0, 118),
    breaks = c(0, 25, 50, 75, 100),
    labels = function(x) paste0(x, "%"),
    expand = c(0, 0)
  ) +
  labs(
    x = "Ranked CDR3β sequence",
    y = "Percentage of matched barcodes"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 7
    ),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 8.5),
    plot.margin = margin(12, 12, 30, 12)
  )


# ------------------------------------------------------------------------------
# 11. Save figure and output tables
# ------------------------------------------------------------------------------

ggsave(
  filename = file.path(figure_dir, "figure7_chain_pairing_top20.pdf"),
  plot = p_figure7,
  width = 9.5,
  height = 5.5
)

ggsave(
  filename = file.path(figure_dir, "figure7_chain_pairing_top20.png"),
  plot = p_figure7,
  width = 9.5,
  height = 5.5,
  dpi = 300
)

fwrite(
  chain_pairing_summary,
  file.path(table_dir, "figure7_chain_pairing_top20_summary.csv")
)

fwrite(
  chain_pairing_wide,
  file.path(table_dir, "figure7_chain_pairing_top20_wide.csv")
)

fwrite(
  matched_cells %>%
    select(
      cell_barcode,
      cdr3_aa,
      rank,
      composite_score,
      tra_cdr3,
      trb_cdr3,
      n_alpha,
      n_beta,
      chain_combo,
      everything()
    ),
  file.path(table_dir, "figure7_matched_top20_single_cell_barcodes.csv")
)


# ------------------------------------------------------------------------------
# 12. Print quick checks
# ------------------------------------------------------------------------------

message("Figure 7 complete.")
message("Figure saved to: ", figure_dir)
message("Tables saved to: ", table_dir)

message("Number of top-20 CDR3β sequences with matched barcodes: ",
        n_distinct(matched_cells$cdr3_aa))

message("Matched barcode counts per top-20 CDR3β:")
print(
  label_df %>%
    arrange(as.numeric(rank))
)

message("Chain-pairing summary:")
print(chain_pairing_wide)