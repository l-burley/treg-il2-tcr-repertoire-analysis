# ==============================================================================
# 09_figure4_ld1_network_top100.R
# ==============================================================================
# Goal:
#   Create Cytoscape-compatible node and edge tables for the top 100
#   IL-2-prioritised Treg CDR3β sequences.
#
#   Nodes are top 100 IL-2-prioritised CDR3β sequences.
#   Edges connect CDR3β sequences with Levenshtein distance = 1.
#
# Input:
#   results/tables/cdr3_scores_top100_il2.csv
#
# Outputs:
#   results/tables/top100_il2_cdr3_ld1_nodes_cytoscape.csv
#   results/tables/top100_il2_cdr3_ld1_edges_cytoscape.csv
#
# Notes:
#   These files can be imported into Cytoscape:
#     - edge table: source, target, interaction
#     - node table: node_id and node attributes
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Load packages
# ------------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(stringdist)


# ------------------------------------------------------------------------------
# 1. Define input and output paths
# ------------------------------------------------------------------------------

input_path <- file.path(
  "results",
  "tables",
  "cdr3_scores_top100_il2.csv"
)

table_dir <- file.path("results", "tables")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# 2. Load top 100 IL-2-prioritised CDR3β table
# ------------------------------------------------------------------------------

top100_scores <- fread(input_path)


# ------------------------------------------------------------------------------
# 3. Check required columns
# ------------------------------------------------------------------------------

required_cols <- c(
  "cdr3_aa",
  "composite_score",
  "il2_prioritised_rank"
)

missing_cols <- setdiff(required_cols, colnames(top100_scores))

if (length(missing_cols) > 0) {
  stop(
    "Missing columns needed for this script: ",
    paste(missing_cols, collapse = ", ")
  )
}


# ------------------------------------------------------------------------------
# 4. Prepare node table
# ------------------------------------------------------------------------------

top100_nodes <- top100_scores %>%
  mutate(
    cdr3_aa = as.character(cdr3_aa),
    composite_score = as.numeric(composite_score),
    il2_prioritised_rank = as.numeric(il2_prioritised_rank)
  ) %>%
  distinct(cdr3_aa, .keep_all = TRUE) %>%
  arrange(il2_prioritised_rank) %>%
  mutate(
    node_id = cdr3_aa,
    label = cdr3_aa
  )

if (nrow(top100_nodes) == 0) {
  stop("No top 100 CDR3β sequences were found.")
}


# ------------------------------------------------------------------------------
# 5. Calculate all pairwise Levenshtein distances
# ------------------------------------------------------------------------------

cdr3s <- top100_nodes$cdr3_aa

pair_grid <- expand.grid(
  source = cdr3s,
  target = cdr3s,
  stringsAsFactors = FALSE
)

# Keep each pair only once and remove self-comparisons.
pair_grid <- pair_grid %>%
  filter(source < target)

ld1_edges <- pair_grid %>%
  mutate(
    levenshtein_distance = stringdist(source, target, method = "lv")
  ) %>%
  filter(levenshtein_distance == 1) %>%
  mutate(
    interaction = "LD1"
  ) %>%
  select(
    source,
    target,
    interaction,
    levenshtein_distance
  )


# ------------------------------------------------------------------------------
# 6. Add node attributes
# ------------------------------------------------------------------------------

connected_nodes <- unique(c(ld1_edges$source, ld1_edges$target))

node_table <- top100_nodes %>%
  mutate(
    has_ld1_connection = node_id %in% connected_nodes,
    degree_ld1 = sapply(node_id, function(x) {
      sum(ld1_edges$source == x | ld1_edges$target == x)
    })
  ) %>%
  select(
    node_id,
    label,
    cdr3_aa,
    il2_prioritised_rank,
    composite_score,
    has_ld1_connection,
    degree_ld1,
    everything()
  )


# ------------------------------------------------------------------------------
# 7. Save Cytoscape input files
# ------------------------------------------------------------------------------

fwrite(
  node_table,
  file.path(table_dir, "top100_il2_cdr3_ld1_nodes_cytoscape.csv")
)

fwrite(
  ld1_edges,
  file.path(table_dir, "top100_il2_cdr3_ld1_edges_cytoscape.csv")
)


# ------------------------------------------------------------------------------
# 8. Print quick checks
# ------------------------------------------------------------------------------

message("Cytoscape LD=1 network files created.")
message("Tables saved to: ", table_dir)

message("Number of nodes: ", nrow(node_table))
message("Number of LD=1 edges: ", nrow(ld1_edges))
message("Number of connected nodes: ", length(connected_nodes))
message("Number of unconnected nodes: ", sum(!node_table$has_ld1_connection))

message("Top connected nodes by LD=1 degree:")
print(
  node_table %>%
    arrange(desc(degree_ld1), il2_prioritised_rank) %>%
    select(
      il2_prioritised_rank,
      cdr3_aa,
      composite_score,
      degree_ld1,
      has_ld1_connection
    ) %>%
    slice_head(n = 10)
)