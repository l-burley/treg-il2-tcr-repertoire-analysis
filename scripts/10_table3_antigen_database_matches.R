# ==============================================================================
# 10_table3_antigen_database_matches.R
# ==============================================================================
# Goal:
#   Match the top 100 IL-2-prioritised Treg CDR3β sequences to an external
#   antigen-associated TCR database.
#
#   Exact CDR3β amino acid matches are searched against a combined antigen
#   database. The script exports:
#     1. A full match table with database annotations
#     2. A compact summary table
#     3. A selected-column table containing:
#        Rank, CDR3B, Database, Antigen, Epitope, Organism, Cell subset
#
# Inputs:
#   results/tables/cdr3_scores_top100_il2.csv
#   data/external/antigen_database/antigen_database_combined.csv
#
# Outputs:
#   results/tables/top100_il2_cdr3_database_matches_full.csv
#   results/tables/top100_il2_cdr3_database_matches_summary.csv
#   results/tables/top100_il2_cdr3_database_matches_selected_columns.csv
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Load packages
# ------------------------------------------------------------------------------

library(data.table)
library(dplyr)


# ------------------------------------------------------------------------------
# 1. Define input and output paths
# ------------------------------------------------------------------------------

top100_path <- file.path(
  "results",
  "tables",
  "cdr3_scores_top100_il2.csv"
)

antigen_db_path <- file.path(
  "data",
  "external",
  "antigen_database",
  "antigen_database_combined.csv"
)

table_dir <- file.path("results", "tables")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# 2. Load input files
# ------------------------------------------------------------------------------

top100_scores <- fread(top100_path)
antigen_db <- fread(antigen_db_path)


# ------------------------------------------------------------------------------
# 3. Check required top 100 score columns
# ------------------------------------------------------------------------------

required_top100_cols <- c(
  "il2_prioritised_rank",
  "cdr3_aa",
  "composite_score"
)

missing_top100_cols <- setdiff(required_top100_cols, colnames(top100_scores))

if (length(missing_top100_cols) > 0) {
  stop(
    "Missing columns in top 100 IL-2 score table: ",
    paste(missing_top100_cols, collapse = ", ")
  )
}


# ------------------------------------------------------------------------------
# 4. Prepare top 100 CDR3β sequences
# ------------------------------------------------------------------------------

top100 <- top100_scores %>%
  mutate(
    il2_prioritised_rank = as.numeric(il2_prioritised_rank),
    cdr3_aa = as.character(cdr3_aa),
    composite_score = as.numeric(composite_score)
  ) %>%
  arrange(il2_prioritised_rank) %>%
  slice_head(n = 100) %>%
  mutate(
    rank = il2_prioritised_rank
  ) %>%
  select(
    rank,
    cdr3_aa,
    composite_score,
    everything()
  )

if (nrow(top100) == 0) {
  stop("No top 100 IL-2-prioritised CDR3β sequences were found.")
}


# ------------------------------------------------------------------------------
# 5. Check and standardise antigen database columns
# ------------------------------------------------------------------------------

# These are the annotation columns expected from the combined antigen database.
# Columns that are not present are added as NA so the script can still run with
# slightly different database versions.

standard_db_cols <- c(
  "CDR3_beta",
  "CDR3_alpha",
  "V_beta",
  "J_beta",
  "Epitope",
  "Antigen",
  "Antigen_organism",
  "Cell_subset",
  "PubMed_ID",
  "Database",
  "Verified_score",
  "Identification_score",
  "grouped_antigen",
  "antigen_group",
  "antigen_category",
  "label"
)

add_missing_columns <- function(x, required_cols) {
  missing_cols <- setdiff(required_cols, colnames(x))
  
  for (col in missing_cols) {
    x[[col]] <- NA_character_
  }
  
  return(x)
}

antigen_db_clean <- antigen_db %>%
  add_missing_columns(standard_db_cols) %>%
  mutate(
    CDR3_beta = as.character(CDR3_beta),
    CDR3_alpha = as.character(CDR3_alpha),
    V_beta = as.character(V_beta),
    J_beta = as.character(J_beta),
    Epitope = as.character(Epitope),
    Antigen = as.character(Antigen),
    Antigen_organism = as.character(Antigen_organism),
    Cell_subset = as.character(Cell_subset),
    PubMed_ID = as.character(PubMed_ID),
    Database = as.character(Database),
    Verified_score = as.character(Verified_score),
    Identification_score = as.character(Identification_score),
    grouped_antigen = as.character(grouped_antigen),
    antigen_group = as.character(antigen_group),
    antigen_category = as.character(antigen_category),
    label = as.character(label)
  ) %>%
  select(all_of(standard_db_cols)) %>%
  filter(
    !is.na(CDR3_beta),
    CDR3_beta != ""
  )

if (nrow(antigen_db_clean) == 0) {
  stop("No valid CDR3_beta entries were found in the antigen database.")
}


# ------------------------------------------------------------------------------
# 6. Match top 100 CDR3β sequences to antigen database
# ------------------------------------------------------------------------------

top100_db_matches <- top100 %>%
  left_join(
    antigen_db_clean,
    by = c("cdr3_aa" = "CDR3_beta")
  ) %>%
  filter(
    !is.na(Database) |
      !is.na(Antigen) |
      !is.na(Epitope) |
      !is.na(Antigen_organism)
  ) %>%
  arrange(
    rank,
    Database,
    Antigen,
    Epitope
  )


# ------------------------------------------------------------------------------
# 7. Make full match table
# ------------------------------------------------------------------------------

top100_db_matches_full <- top100_db_matches %>%
  select(
    rank,
    cdr3_aa,
    composite_score,
    Database,
    CDR3_alpha,
    V_beta,
    J_beta,
    Epitope,
    Antigen,
    Antigen_organism,
    Cell_subset,
    PubMed_ID,
    Verified_score,
    Identification_score,
    grouped_antigen,
    antigen_group,
    antigen_category,
    label,
    everything()
  )


# ------------------------------------------------------------------------------
# 8. Make compact summary table
# ------------------------------------------------------------------------------

top100_match_summary <- top100_db_matches %>%
  group_by(
    rank,
    cdr3_aa,
    composite_score
  ) %>%
  summarise(
    n_database_hits = n(),
    source_databases = paste(
      unique(na.omit(Database[Database != ""])),
      collapse = "; "
    ),
    antigens = paste(
      unique(na.omit(Antigen[Antigen != ""])),
      collapse = "; "
    ),
    epitopes = paste(
      unique(na.omit(Epitope[Epitope != ""])),
      collapse = "; "
    ),
    antigen_organisms = paste(
      unique(na.omit(Antigen_organism[Antigen_organism != ""])),
      collapse = "; "
    ),
    cell_subsets = paste(
      unique(na.omit(Cell_subset[Cell_subset != ""])),
      collapse = "; "
    ),
    pubmed_ids = paste(
      unique(na.omit(PubMed_ID[PubMed_ID != ""])),
      collapse = "; "
    ),
    .groups = "drop"
  ) %>%
  arrange(rank)


# ------------------------------------------------------------------------------
# 9. Make selected-column table
# ------------------------------------------------------------------------------

top100_match_selected_columns <- top100_match_summary %>%
  transmute(
    Rank = rank,
    CDR3B = cdr3_aa,
    Database = source_databases,
    Antigen = antigens,
    Epitope = epitopes,
    Organism = antigen_organisms,
    `Cell subset` = cell_subsets
  ) %>%
  arrange(Rank)


# ------------------------------------------------------------------------------
# 10. Save output files
# ------------------------------------------------------------------------------

fwrite(
  top100_db_matches_full,
  file.path(table_dir, "top100_il2_cdr3_database_matches_full.csv")
)

fwrite(
  top100_match_summary,
  file.path(table_dir, "top100_il2_cdr3_database_matches_summary.csv")
)

fwrite(
  top100_match_selected_columns,
  file.path(table_dir, "top100_il2_cdr3_database_matches_selected_columns.csv")
)


# ------------------------------------------------------------------------------
# 11. Print quick summary
# ------------------------------------------------------------------------------

message("Antigen database matching complete.")
message("Tables saved to: ", table_dir)

message(
  "Number of top 100 CDR3β sequences with at least one database match: ",
  n_distinct(top100_db_matches$cdr3_aa)
)

message(
  "Total number of database match rows: ",
  nrow(top100_db_matches)
)

message("Selected-column antigen match table:")
print(top100_match_selected_columns)