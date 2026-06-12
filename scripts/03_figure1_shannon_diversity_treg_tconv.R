# ==============================================================================
# 03_figure1_shannon_diversity_treg_tconv.R
# ==============================================================================
# Goal:
#   Calculate Shannon diversity for each bulk TCRβ sample and compare matched
#   Treg and Tconv repertoires from the same mice.
#
# Input:
#   data/processed/repseq_objects/RepSeqData_Curie003.RData
#
# Outputs:
#   results/tables/figure1_shannon_diversity_by_sample.csv
#   results/tables/figure1_shannon_diversity_paired_wide.csv
#   results/tables/figure1_shannon_wilcoxon_result.csv
#   results/figures/figure1_shannon_diversity_treg_tconv.png
#   results/figures/figure1_shannon_diversity_treg_tconv.pdf
#
# Notes:
#   Shannon diversity is calculated as:
#     H = -sum(p_i * log(p_i))
#
#   where p_i is the frequency of clonotype i within a sample.
#
#   Since each mouse has matched Treg and Tconv samples, a paired Wilcoxon
#   signed-rank test is used.
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Load packages
# ------------------------------------------------------------------------------

library(AnalyzAIRR)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)


# ------------------------------------------------------------------------------
# 1. Define input and output paths
# ------------------------------------------------------------------------------

repseq_path <- file.path(
  "data",
  "processed",
  "repseq_objects",
  "RepSeqData_Curie003.RData"
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


# ------------------------------------------------------------------------------
# 3. Extract assayData and metaData
# ------------------------------------------------------------------------------

assay_data <- as.data.table(RepSeqData_Curie003@assayData)
metadata <- as.data.frame(RepSeqData_Curie003@metaData)


# ------------------------------------------------------------------------------
# 4. Check columns needed for this script
# ------------------------------------------------------------------------------

required_assay_cols <- c(
  "sample_id",
  "count"
)

required_meta_cols <- c(
  "sample_id",
  "mouse",
  "cell_subset",
  "injection"
)

missing_assay_cols <- setdiff(required_assay_cols, colnames(assay_data))
missing_meta_cols <- setdiff(required_meta_cols, colnames(metadata))

if (length(missing_assay_cols) > 0) {
  stop(
    "Missing assayData columns needed for this script: ",
    paste(missing_assay_cols, collapse = ", ")
  )
}

if (length(missing_meta_cols) > 0) {
  stop(
    "Missing metaData columns needed for this script: ",
    paste(missing_meta_cols, collapse = ", ")
  )
}


# ------------------------------------------------------------------------------
# 5. Clean key columns
# ------------------------------------------------------------------------------

assay_data <- assay_data %>%
  mutate(
    sample_id = as.character(sample_id),
    count = as.numeric(count)
  )

metadata <- metadata %>%
  mutate(
    sample_id = as.character(sample_id),
    mouse = as.character(mouse),
    cell_subset = as.character(cell_subset),
    injection = as.character(injection)
  )


# ------------------------------------------------------------------------------
# 6. Calculate Shannon diversity per sample
# ------------------------------------------------------------------------------

shannon_by_sample <- assay_data %>%
  group_by(sample_id) %>%
  mutate(
    total_count = sum(count, na.rm = TRUE),
    frequency = count / total_count
  ) %>%
  summarise(
    shannon_index = -sum(frequency * log(frequency), na.rm = TRUE),
    total_count = first(total_count),
    n_clonotypes = n(),
    .groups = "drop"
  )


# ------------------------------------------------------------------------------
# 7. Add sample metadata
# ------------------------------------------------------------------------------

shannon_with_metadata <- shannon_by_sample %>%
  left_join(
    metadata %>%
      select(
        sample_id,
        mouse,
        cell_subset,
        injection
      ),
    by = "sample_id"
  )


# ------------------------------------------------------------------------------
# 8. Keep matched Treg and Tconv samples
# ------------------------------------------------------------------------------

paired_shannon <- shannon_with_metadata %>%
  filter(cell_subset %in% c("Treg", "Tconv")) %>%
  group_by(mouse) %>%
  filter(all(c("Treg", "Tconv") %in% cell_subset)) %>%
  ungroup()

if (n_distinct(paired_shannon$mouse) == 0) {
  stop("No mice with paired Treg and Tconv samples were found.")
}


# ------------------------------------------------------------------------------
# 9. Convert to wide format for paired Wilcoxon test
# ------------------------------------------------------------------------------

paired_wide <- paired_shannon %>%
  select(mouse, cell_subset, shannon_index) %>%
  pivot_wider(
    names_from = cell_subset,
    values_from = shannon_index
  ) %>%
  mutate(
    difference_Treg_minus_Tconv = Treg - Tconv
  )

if (!all(c("Treg", "Tconv") %in% colnames(paired_wide))) {
  stop("Paired table does not contain both Treg and Tconv columns.")
}


# ------------------------------------------------------------------------------
# 10. Perform paired Wilcoxon signed-rank test
# ------------------------------------------------------------------------------

wilcox_result <- wilcox.test(
  paired_wide$Treg,
  paired_wide$Tconv,
  paired = TRUE,
  exact = FALSE
)

wilcox_table <- tibble(
  test = "paired Wilcoxon signed-rank test",
  comparison = "Treg vs Tconv",
  n_paired_mice = nrow(paired_wide),
  statistic = unname(wilcox_result$statistic),
  p_value = wilcox_result$p.value,
  median_Treg = median(paired_wide$Treg, na.rm = TRUE),
  median_Tconv = median(paired_wide$Tconv, na.rm = TRUE),
  median_difference_Treg_minus_Tconv = median(
    paired_wide$difference_Treg_minus_Tconv,
    na.rm = TRUE
  )
)


# ------------------------------------------------------------------------------
# 11. Create p-value label
# ------------------------------------------------------------------------------

p_value <- wilcox_result$p.value

p_label <- ifelse(
  p_value < 0.001,
  "p < 0.001",
  paste0("p = ", signif(p_value, 3))
)


# ------------------------------------------------------------------------------
# 12. Set plotting order
# ------------------------------------------------------------------------------

paired_shannon <- paired_shannon %>%
  mutate(
    cell_subset = factor(cell_subset, levels = c("Treg", "Tconv")),
    injection = factor(injection)
  )


# ------------------------------------------------------------------------------
# 13. Make paired Shannon diversity plot
# ------------------------------------------------------------------------------

p_shannon_paired <- ggplot(
  paired_shannon,
  aes(
    x = cell_subset,
    y = shannon_index,
    group = mouse
  )
) +
  geom_line(
    alpha = 0.5,
    linewidth = 0.6
  ) +
  geom_point(
    aes(shape = injection),
    size = 3
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(paired_shannon$shannon_index, na.rm = TRUE) * 1.05,
    label = p_label,
    size = 4
  ) +
  coord_cartesian(
    ylim = c(
      min(paired_shannon$shannon_index, na.rm = TRUE),
      max(paired_shannon$shannon_index, na.rm = TRUE) * 1.1
    )
  ) +
  theme_bw() +
  labs(
    title = "Shannon diversity of paired Treg and Tconv repertoires",
    x = "Cell subset",
    y = "Shannon diversity index",
    shape = "Injection"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 13, face = "bold")
  )


# ------------------------------------------------------------------------------
# 14. Save tables
# ------------------------------------------------------------------------------

write_csv(
  paired_shannon,
  file.path(table_dir, "figure1_shannon_diversity_by_sample.csv")
)

write_csv(
  paired_wide,
  file.path(table_dir, "figure1_shannon_diversity_paired_wide.csv")
)

write_csv(
  wilcox_table,
  file.path(table_dir, "figure1_shannon_wilcoxon_result.csv")
)


# ------------------------------------------------------------------------------
# 15. Save figure
# ------------------------------------------------------------------------------

ggsave(
  filename = file.path(figure_dir, "figure1_shannon_diversity_treg_tconv.png"),
  plot = p_shannon_paired,
  width = 6,
  height = 5,
  dpi = 300
)

ggsave(
  filename = file.path(figure_dir, "figure1_shannon_diversity_treg_tconv.pdf"),
  plot = p_shannon_paired,
  width = 6,
  height = 5
)


# ------------------------------------------------------------------------------
# 16. Print quick checks
# ------------------------------------------------------------------------------

message("Figure 1 Shannon diversity analysis complete.")
message("Number of paired mice: ", nrow(paired_wide))
message("Wilcoxon p-value: ", signif(wilcox_result$p.value, 3))
message("Figures saved to: ", figure_dir)
message("Tables saved to: ", table_dir)