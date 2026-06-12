# ==============================================================================
# 08_figure3_compare_composite_scores.R
# ==============================================================================
# Goal:
#   Create Figure 3 comparing the top 100 CDR3β sequences from:
#     1. Treatment-independent composite scoring
#     2. IL-2-prioritised composite scoring
#
#   For each scoring method, count how many of the top 100 CDR3β sequences are
#   detected in each Treg sample at frequency >= 0.01%.
#
# Inputs:
#   data/processed/repseq_objects/RepSeqData_Curie003.RData
#   results/tables/cdr3_scores_full_treatment_independent.csv
#   results/tables/cdr3_scores_full_il2.csv
#
# Outputs:
#   results/tables/figure3_top100_presence_summary.csv
#   results/figures/figure3_treatment_independent_vs_il2_top100_presence_side_by_side.png
#   results/figures/figure3_treatment_independent_vs_il2_top100_presence_side_by_side.pdf
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Load packages
# ------------------------------------------------------------------------------

library(AnalyzAIRR)
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)


# ------------------------------------------------------------------------------
# 1. Define input and output paths
# ------------------------------------------------------------------------------

repseq_path <- file.path(
  "data",
  "processed",
  "repseq_objects",
  "RepSeqData_Curie003.RData"
)

treatment_independent_path <- file.path(
  "results",
  "tables",
  "cdr3_scores_full_treatment_independent.csv"
)

il2_prioritised_path <- file.path(
  "results",
  "tables",
  "cdr3_scores_full_il2.csv"
)

figure_dir <- file.path("results", "figures")
table_dir <- file.path("results", "tables")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# 2. Load RepSeqData object
# ------------------------------------------------------------------------------

load(repseq_path)

if (!exists("RepSeqData_Curie003")) {
  stop("Object RepSeqData_Curie003 was not found after loading.")
}

repseq_object <- RepSeqData_Curie003

dat <- repseq_object@assayData %>%
  as.data.frame()

meta <- repseq_object@metaData %>%
  as.data.frame()


# ------------------------------------------------------------------------------
# 3. Load composite score tables
# ------------------------------------------------------------------------------

treatment_independent_scores <- fread(treatment_independent_path)

il2_prioritised_scores <- fread(il2_prioritised_path)


# ------------------------------------------------------------------------------
# 4. Check required columns
# ------------------------------------------------------------------------------

required_dat_cols <- c(
  "sample_id",
  "aaCDR3",
  "count"
)

required_meta_cols <- c(
  "sample_id",
  "mouse",
  "cell_subset",
  "injection",
  "protein",
  "nSequences"
)

required_treatment_cols <- c(
  "cdr3_aa",
  "treatment_independent_score"
)

required_il2_cols <- c(
  "cdr3_aa",
  "composite_score"
)

missing_dat_cols <- setdiff(required_dat_cols, colnames(dat))
missing_meta_cols <- setdiff(required_meta_cols, colnames(meta))
missing_treatment_cols <- setdiff(
  required_treatment_cols,
  colnames(treatment_independent_scores)
)
missing_il2_cols <- setdiff(
  required_il2_cols,
  colnames(il2_prioritised_scores)
)

if (length(missing_dat_cols) > 0) {
  stop("Missing assayData columns: ", paste(missing_dat_cols, collapse = ", "))
}

if (length(missing_meta_cols) > 0) {
  stop("Missing metaData columns: ", paste(missing_meta_cols, collapse = ", "))
}

if (length(missing_treatment_cols) > 0) {
  stop(
    "Missing treatment-independent score columns: ",
    paste(missing_treatment_cols, collapse = ", ")
  )
}

if (length(missing_il2_cols) > 0) {
  stop(
    "Missing IL-2-prioritised score columns: ",
    paste(missing_il2_cols, collapse = ", ")
  )
}


# ------------------------------------------------------------------------------
# 5. Clean input tables
# ------------------------------------------------------------------------------

dat <- dat %>%
  mutate(
    sample_id = as.character(sample_id),
    aaCDR3 = as.character(aaCDR3),
    count = as.numeric(count)
  )

meta <- meta %>%
  mutate(
    sample_id = as.character(sample_id),
    mouse = as.character(mouse),
    cell_subset = as.character(cell_subset),
    injection = as.character(injection),
    protein = as.character(protein),
    nSequences = as.numeric(nSequences)
  )

treatment_independent_scores <- treatment_independent_scores %>%
  mutate(
    cdr3_aa = as.character(cdr3_aa),
    treatment_independent_score = as.numeric(treatment_independent_score)
  )

il2_prioritised_scores <- il2_prioritised_scores %>%
  mutate(
    cdr3_aa = as.character(cdr3_aa),
    composite_score = as.numeric(composite_score)
  )


# ------------------------------------------------------------------------------
# 6. Select top 100 CDR3β sequences from each scoring method
# ------------------------------------------------------------------------------

treatment_top100_cdr3 <- treatment_independent_scores %>%
  arrange(desc(treatment_independent_score)) %>%
  slice_head(n = 100) %>%
  pull(cdr3_aa) %>%
  unique()

il2_top100_cdr3 <- il2_prioritised_scores %>%
  arrange(desc(composite_score)) %>%
  slice_head(n = 100) %>%
  pull(cdr3_aa) %>%
  unique()

if (length(treatment_top100_cdr3) == 0) {
  stop("No treatment-independent top 100 CDR3β sequences were found.")
}

if (length(il2_top100_cdr3) == 0) {
  stop("No IL-2-prioritised top 100 CDR3β sequences were found.")
}


# ------------------------------------------------------------------------------
# 7. Prepare Treg sample table and clonotype frequencies
# ------------------------------------------------------------------------------

treg_meta <- meta %>%
  filter(cell_subset == "Treg") %>%
  select(
    sample_id,
    mouse,
    injection,
    protein,
    nSequences
  )

if (nrow(treg_meta) == 0) {
  stop("No Treg samples were found in the metadata.")
}

dat_treg <- dat %>%
  inner_join(treg_meta, by = "sample_id") %>%
  mutate(
    frequency = count / nSequences
  )

all_treg_samples <- treg_meta %>%
  distinct(
    sample_id,
    mouse,
    injection
  )


# ------------------------------------------------------------------------------
# 8. Function to count top 100 CDR3β sequences per sample
# ------------------------------------------------------------------------------

count_top100_per_sample <- function(top100_cdr3, count_column_name, percent_column_name) {
  
  top100_present <- dat_treg %>%
    filter(
      aaCDR3 %in% top100_cdr3,
      frequency >= 0.0001
    )
  
  top100_per_sample <- top100_present %>%
    distinct(
      sample_id,
      mouse,
      injection,
      aaCDR3
    ) %>%
    count(
      sample_id,
      mouse,
      injection,
      name = count_column_name
    ) %>%
    mutate(
      "{percent_column_name}" := .data[[count_column_name]] /
        length(top100_cdr3) * 100
    )
  
  top100_per_sample <- all_treg_samples %>%
    left_join(
      top100_per_sample,
      by = c("sample_id", "mouse", "injection")
    ) %>%
    mutate(
      "{count_column_name}" := ifelse(
        is.na(.data[[count_column_name]]),
        0,
        .data[[count_column_name]]
      ),
      "{percent_column_name}" := ifelse(
        is.na(.data[[percent_column_name]]),
        0,
        .data[[percent_column_name]]
      )
    ) %>%
    arrange(injection, mouse)
  
  return(top100_per_sample)
}


# ------------------------------------------------------------------------------
# 9. Count top 100 presence for each scoring method
# ------------------------------------------------------------------------------

treatment_top100_per_sample <- count_top100_per_sample(
  top100_cdr3 = treatment_top100_cdr3,
  count_column_name = "n_treatment_top100_present",
  percent_column_name = "percent_treatment_top100_present"
)

il2_top100_per_sample <- count_top100_per_sample(
  top100_cdr3 = il2_top100_cdr3,
  count_column_name = "n_il2_top100_present",
  percent_column_name = "percent_il2_top100_present"
)


# ------------------------------------------------------------------------------
# 10. Set mouse plotting order
# ------------------------------------------------------------------------------

# Each panel is ordered by treatment group and then by the number of detected
# top 100 CDR3β sequences, matching the original plotting approach.

treatment_top100_per_sample <- treatment_top100_per_sample %>%
  mutate(
    injection = factor(injection, levels = c("CTRL", "AAV")),
    mouse = factor(
      mouse,
      levels = mouse[order(injection, n_treatment_top100_present)]
    )
  )

il2_top100_per_sample <- il2_top100_per_sample %>%
  mutate(
    injection = factor(injection, levels = c("CTRL", "AAV")),
    mouse = factor(
      mouse,
      levels = mouse[order(injection, n_il2_top100_present)]
    )
  )


# ------------------------------------------------------------------------------
# 11. Create treatment-independent top 100 plot
# ------------------------------------------------------------------------------

p_treatment_top100_presence <- ggplot(
  treatment_top100_per_sample,
  aes(
    x = mouse,
    y = n_treatment_top100_present,
    fill = injection
  )
) +
  geom_col(
    color = "black",
    width = 0.75,
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = n_treatment_top100_present),
    vjust = -0.4,
    size = 3.5
  ) +
  scale_fill_manual(
    values = c(
      "CTRL" = "grey75",
      "AAV" = "#24325F"
    ),
    labels = c(
      "CTRL" = "AAV-PBS",
      "AAV" = "AAV-IL-2"
    )
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.title = element_blank()
  ) +
  labs(
    x = "Mouse",
    y = "Number of treatment-independent top 100\nCDR3β sequences detected at ≥0.01%",
    fill = "Treatment"
  ) +
  ylim(
    0,
    max(treatment_top100_per_sample$n_treatment_top100_present) + 5
  )


# ------------------------------------------------------------------------------
# 12. Create IL-2-prioritised top 100 plot
# ------------------------------------------------------------------------------

p_il2_top100_presence <- ggplot(
  il2_top100_per_sample,
  aes(
    x = mouse,
    y = n_il2_top100_present,
    fill = injection
  )
) +
  geom_col(
    color = "black",
    width = 0.75,
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = n_il2_top100_present),
    vjust = -0.4,
    size = 3.5
  ) +
  scale_fill_manual(
    values = c(
      "CTRL" = "grey75",
      "AAV" = "#24325F"
    ),
    labels = c(
      "CTRL" = "AAV-PBS",
      "AAV" = "AAV-IL-2"
    )
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.title = element_blank()
  ) +
  labs(
    x = "Mouse",
    y = "Number of IL-2-prioritised top 100\nCDR3β sequences detected at ≥0.01%",
    fill = "Treatment"
  ) +
  ylim(
    0,
    max(il2_top100_per_sample$n_il2_top100_present) + 5
  )


# ------------------------------------------------------------------------------
# 13. Combine panels side by side
# ------------------------------------------------------------------------------

p_figure3_combined <- p_treatment_top100_presence | p_il2_top100_presence


# ------------------------------------------------------------------------------
# 14. Save output table
# ------------------------------------------------------------------------------

figure3_summary <- treatment_top100_per_sample %>%
  mutate(
    mouse = as.character(mouse),
    injection = as.character(injection)
  ) %>%
  select(
    sample_id,
    mouse,
    injection,
    n_treatment_top100_present,
    percent_treatment_top100_present
  ) %>%
  left_join(
    il2_top100_per_sample %>%
      mutate(
        mouse = as.character(mouse),
        injection = as.character(injection)
      ) %>%
      select(
        sample_id,
        n_il2_top100_present,
        percent_il2_top100_present
      ),
    by = "sample_id"
  )

fwrite(
  figure3_summary,
  file.path(table_dir, "figure3_top100_presence_summary.csv")
)


# ------------------------------------------------------------------------------
# 15. Save combined figure
# ------------------------------------------------------------------------------

ggsave(
  filename = file.path(
    figure_dir,
    "figure3_treatment_independent_vs_il2_top100_presence_side_by_side.png"
  ),
  plot = p_figure3_combined,
  width = 16,
  height = 5,
  dpi = 300
)

ggsave(
  filename = file.path(
    figure_dir,
    "figure3_treatment_independent_vs_il2_top100_presence_side_by_side.pdf"
  ),
  plot = p_figure3_combined,
  width = 16,
  height = 5
)


# ------------------------------------------------------------------------------
# 16. Print quick checks
# ------------------------------------------------------------------------------

message("Figure 3 complete.")
message("Figure saved to: ", figure_dir)
message("Table saved to: ", table_dir)

message("Treatment-independent top 100 length: ", length(treatment_top100_cdr3))
message("IL-2-prioritised top 100 length: ", length(il2_top100_cdr3))
message("Top 100 overlap: ", length(intersect(treatment_top100_cdr3, il2_top100_cdr3)))

message("Treatment-independent top 100 per sample:")
print(treatment_top100_per_sample)

message("IL-2-prioritised top 100 per sample:")
print(il2_top100_per_sample)