# ==============================================================================
# 11_figure5_top20_il2_pbs_tripod.R
# ==============================================================================
# Goal:
#   Create Figure 5 comparing the frequencies of the top 20 IL-2-prioritised
#   Treg CDR3β sequences across AAV-PBS, AAV-IL-2, and an external TRIPOD
#   untreated B6 Treg reference dataset.
#
#   Panel A:
#     Layered dot plot showing the frequency of each top 20 CDR3β sequence
#     across individual AAV-PBS, AAV-IL-2, and TRIPOD Treg repertoires.
#
#   Panel B:
#     Mean-frequency lollipop plot comparing AAV-PBS, AAV-IL-2, and TRIPOD
#     untreated B6 repertoires for each top 20 CDR3β sequence.
#
# Inputs:
#   data/processed/repseq_objects/RepSeqData_Curie003.RData
#   data/external/tripod/RepSeqData_tripod_b6_tregs.RData
#   results/tables/cdr3_scores_top100_il2.csv
#
# Outputs:
#   results/figures/figure5_top20_il2_pbs_tripod_panel.png
#   results/figures/figure5_top20_il2_pbs_tripod_panel.pdf
#   results/tables/figure5_top20_dot_plot_table.csv
#   results/tables/figure5_top20_lollipop_summary.csv
#   results/tables/figure5_top20_lollipop_wide.csv
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Load packages
# ------------------------------------------------------------------------------

library(AnalyzAIRR)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(patchwork)


# ------------------------------------------------------------------------------
# 1. Define input and output paths
# ------------------------------------------------------------------------------

curie_repseq_path <- file.path(
  "data",
  "processed",
  "repseq_objects",
  "RepSeqData_Curie003.RData"
)

tripod_repseq_path <- file.path(
  "data",
  "external",
  "tripod",
  "RepSeqData_tripod_b6_tregs.RData"
)

il2_scores_path <- file.path(
  "results",
  "tables",
  "cdr3_scores_top100_il2.csv"
)

figure_dir <- file.path("results", "figures")
table_dir <- file.path("results", "tables")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# 2. Load Curie RepSeqData object
# ------------------------------------------------------------------------------

load(curie_repseq_path)

if (!exists("RepSeqData_Curie003")) {
  stop("Object RepSeqData_Curie003 was not found after loading Curie data.")
}

curie_meta <- RepSeqData_Curie003@metaData %>%
  as.data.frame()

curie_dat <- RepSeqData_Curie003@assayData %>%
  as.data.frame()


# ------------------------------------------------------------------------------
# 3. Load TRIPOD RepSeqData object
# ------------------------------------------------------------------------------

load(tripod_repseq_path)

if (!exists("RepSeqData_tripod")) {
  stop("Object RepSeqData_tripod was not found after loading TRIPOD data.")
}

tripod_meta <- RepSeqData_tripod@metaData %>%
  as.data.frame()

tripod_dat <- RepSeqData_tripod@assayData %>%
  as.data.frame()


# ------------------------------------------------------------------------------
# 4. Load IL-2-prioritised scoring table and extract top 20 CDR3β sequences
# ------------------------------------------------------------------------------

il2_scores <- fread(il2_scores_path)

required_score_cols <- c(
  "il2_prioritised_rank",
  "cdr3_aa",
  "composite_score"
)

missing_score_cols <- setdiff(required_score_cols, colnames(il2_scores))

if (length(missing_score_cols) > 0) {
  stop(
    "Missing columns in IL-2 score table: ",
    paste(missing_score_cols, collapse = ", ")
  )
}

top20_rank_table <- il2_scores %>%
  mutate(
    il2_prioritised_rank = as.numeric(il2_prioritised_rank),
    cdr3_aa = as.character(cdr3_aa),
    composite_score = as.numeric(composite_score)
  ) %>%
  arrange(il2_prioritised_rank) %>%
  slice_head(n = 20) %>%
  transmute(
    composite_rank = il2_prioritised_rank,
    cdr3_aa,
    composite_score
  )

top20_cdr3 <- top20_rank_table$cdr3_aa

if (length(top20_cdr3) == 0) {
  stop("No top 20 CDR3β sequences were found.")
}


# ------------------------------------------------------------------------------
# 5. Check required Curie columns
# ------------------------------------------------------------------------------

required_curie_dat_cols <- c(
  "sample_id",
  "aaCDR3",
  "count"
)

required_curie_meta_cols <- c(
  "sample_id",
  "mouse",
  "cell_subset",
  "injection",
  "protein"
)

missing_curie_dat_cols <- setdiff(required_curie_dat_cols, colnames(curie_dat))
missing_curie_meta_cols <- setdiff(required_curie_meta_cols, colnames(curie_meta))

if (length(missing_curie_dat_cols) > 0) {
  stop(
    "Missing Curie assayData columns: ",
    paste(missing_curie_dat_cols, collapse = ", ")
  )
}

if (length(missing_curie_meta_cols) > 0) {
  stop(
    "Missing Curie metaData columns: ",
    paste(missing_curie_meta_cols, collapse = ", ")
  )
}


# ------------------------------------------------------------------------------
# 6. Prepare Curie AAV-PBS and AAV-IL-2 Treg data
# ------------------------------------------------------------------------------

curie_meta_clean <- curie_meta %>%
  mutate(
    sample_id = as.character(sample_id),
    mouse = as.character(mouse),
    cell_subset = as.character(cell_subset),
    injection = as.character(injection),
    protein = as.character(protein)
  ) %>%
  filter(cell_subset == "Treg") %>%
  mutate(
    analysis_group = case_when(
      protein == "PBS" | str_detect(protein, "PBS") ~ "PBS",
      protein != "PBS" & !str_detect(protein, "PBS") ~ "AAV_IL2",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(analysis_group)) %>%
  select(
    sample_id,
    mouse,
    cell_subset,
    injection,
    protein,
    analysis_group
  )

curie_dat_meta <- curie_dat %>%
  mutate(
    sample_id = as.character(sample_id),
    aaCDR3 = as.character(aaCDR3),
    count = as.numeric(count)
  ) %>%
  inner_join(
    curie_meta_clean,
    by = "sample_id"
  )

if (nrow(curie_dat_meta) == 0) {
  stop("No Curie Treg data were found after joining assayData and metaData.")
}


# ------------------------------------------------------------------------------
# 7. Calculate Curie top 20 frequencies per sample
# ------------------------------------------------------------------------------

curie_sample_totals <- curie_dat_meta %>%
  group_by(
    sample_id,
    analysis_group
  ) %>%
  summarise(
    total_sample_count = sum(count, na.rm = TRUE),
    .groups = "drop"
  )

curie_top20_grid <- curie_sample_totals %>%
  select(
    sample_id,
    analysis_group,
    total_sample_count
  ) %>%
  crossing(
    aaCDR3 = top20_cdr3
  )

curie_top20_counts <- curie_dat_meta %>%
  filter(aaCDR3 %in% top20_cdr3) %>%
  group_by(
    sample_id,
    aaCDR3
  ) %>%
  summarise(
    cdr3_count = sum(count, na.rm = TRUE),
    .groups = "drop"
  )

curie_top20_freq_long <- curie_top20_grid %>%
  left_join(
    curie_top20_counts,
    by = c("sample_id", "aaCDR3")
  ) %>%
  mutate(
    cdr3_count = replace_na(cdr3_count, 0),
    frequency = cdr3_count / total_sample_count,
    frequency_percent = frequency * 100,
    detected = cdr3_count > 0
  )


# ------------------------------------------------------------------------------
# 8. Prepare TRIPOD metadata matching
# ------------------------------------------------------------------------------

make_sample_key <- function(x) {
  x %>%
    as.character() %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]", "")
}

required_tripod_dat_cols <- c(
  "sample_id",
  "aaCDR3",
  "count"
)

missing_tripod_dat_cols <- setdiff(required_tripod_dat_cols, colnames(tripod_dat))

if (length(missing_tripod_dat_cols) > 0) {
  stop(
    "Missing TRIPOD assayData columns: ",
    paste(missing_tripod_dat_cols, collapse = ", ")
  )
}

if (!"sample_id" %in% colnames(tripod_meta)) {
  stop("Missing sample_id column in TRIPOD metadata.")
}

if (!"sample_id_full" %in% colnames(tripod_meta)) {
  tripod_meta$sample_id_full <- NA_character_
}

if (!"project" %in% colnames(tripod_meta)) {
  tripod_meta$project <- NA_character_
}

if (!"cell_subset" %in% colnames(tripod_meta)) {
  tripod_meta$cell_subset <- NA_character_
}

tripod_dat <- tripod_dat %>%
  mutate(
    sample_id = as.character(sample_id),
    aaCDR3 = as.character(aaCDR3),
    count = as.numeric(count),
    sample_key = make_sample_key(sample_id)
  )

tripod_meta <- tripod_meta %>%
  mutate(
    sample_id = as.character(sample_id),
    sample_id_full = as.character(sample_id_full),
    sample_key_1 = make_sample_key(sample_id),
    sample_key_2 = make_sample_key(sample_id_full)
  )

tripod_meta_keyed <- bind_rows(
  tripod_meta %>% mutate(sample_key = sample_key_1),
  tripod_meta %>% mutate(sample_key = sample_key_2)
) %>%
  filter(
    !is.na(sample_key),
    sample_key != ""
  ) %>%
  distinct(
    sample_key,
    .keep_all = TRUE
  )

tripod_meta_cols_to_keep <- intersect(
  c(
    "sample_key",
    "sample_id",
    "sample_id_full",
    "id",
    "project",
    "disease_acronym",
    "cell_subset",
    "cell_phenotype",
    "anatomic_site",
    "exp_group",
    "batch"
  ),
  colnames(tripod_meta_keyed)
)

tripod_dat_meta <- tripod_dat %>%
  left_join(
    tripod_meta_keyed %>%
      select(all_of(tripod_meta_cols_to_keep)) %>%
      rename(sample_id_meta = sample_id),
    by = "sample_key"
  ) %>%
  mutate(
    analysis_group = case_when(
      project == "TRiPoD" | is.na(project) ~ "TRIPOD",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(analysis_group == "TRIPOD") %>%
  filter(cell_subset %in% c("Treg", "Tregs", "amTregs", "nTregs") | is.na(cell_subset))

if (nrow(tripod_dat_meta) == 0) {
  stop("No TRIPOD Treg data were found after metadata matching.")
}


# ------------------------------------------------------------------------------
# 9. Calculate TRIPOD top 20 frequencies per sample
# ------------------------------------------------------------------------------

tripod_sample_totals <- tripod_dat_meta %>%
  group_by(
    sample_id,
    analysis_group
  ) %>%
  summarise(
    total_sample_count = sum(count, na.rm = TRUE),
    .groups = "drop"
  )

tripod_top20_grid <- tripod_sample_totals %>%
  select(
    sample_id,
    analysis_group,
    total_sample_count
  ) %>%
  crossing(
    aaCDR3 = top20_cdr3
  )

tripod_top20_counts <- tripod_dat_meta %>%
  filter(aaCDR3 %in% top20_cdr3) %>%
  group_by(
    sample_id,
    aaCDR3
  ) %>%
  summarise(
    cdr3_count = sum(count, na.rm = TRUE),
    .groups = "drop"
  )

tripod_top20_freq_long <- tripod_top20_grid %>%
  left_join(
    tripod_top20_counts,
    by = c("sample_id", "aaCDR3")
  ) %>%
  mutate(
    cdr3_count = replace_na(cdr3_count, 0),
    frequency = cdr3_count / total_sample_count,
    frequency_percent = frequency * 100,
    detected = cdr3_count > 0
  )


# ------------------------------------------------------------------------------
# 10. Combine Curie and TRIPOD top 20 frequency tables
# ------------------------------------------------------------------------------

combined_top20_freq_long <- bind_rows(
  curie_top20_freq_long,
  tripod_top20_freq_long
) %>%
  mutate(
    analysis_group = factor(
      analysis_group,
      levels = c("PBS", "AAV_IL2", "TRIPOD")
    )
  )


# ------------------------------------------------------------------------------
# 11. Prepare dot plot table
# ------------------------------------------------------------------------------

top20_dot_all <- combined_top20_freq_long %>%
  filter(analysis_group %in% c("PBS", "AAV_IL2", "TRIPOD")) %>%
  left_join(
    top20_rank_table %>%
      select(
        composite_rank,
        cdr3_aa,
        composite_score
      ),
    by = c("aaCDR3" = "cdr3_aa")
  ) %>%
  mutate(
    analysis_group = factor(
      analysis_group,
      levels = c("PBS", "AAV_IL2", "TRIPOD")
    ),
    aaCDR3_label = paste0(composite_rank, ". ", aaCDR3),
    aaCDR3_label = factor(
      aaCDR3_label,
      levels = top20_rank_table %>%
        mutate(aaCDR3_label = paste0(composite_rank, ". ", cdr3_aa)) %>%
        pull(aaCDR3_label)
    ),
    freq_for_plot = frequency_percent + 1e-6
  )


# ------------------------------------------------------------------------------
# 12. Create panel A: layered dot plot
# ------------------------------------------------------------------------------

p_top20_dot_all_layered <- ggplot() +
  geom_jitter(
    data = top20_dot_all %>%
      filter(analysis_group == "TRIPOD"),
    aes(
      x = aaCDR3_label,
      y = freq_for_plot
    ),
    width = 0.25,
    height = 0,
    alpha = 0.5,
    size = 1.1,
    shape = 21,
    color = "grey40",
    fill = "white",
    stroke = 0.35
  ) +
  geom_jitter(
    data = top20_dot_all %>%
      filter(analysis_group %in% c("PBS", "AAV_IL2")),
    aes(
      x = aaCDR3_label,
      y = freq_for_plot,
      color = analysis_group
    ),
    width = 0.25,
    height = 0,
    alpha = 0.95,
    size = 2.2
  ) +
  scale_y_log10() +
  scale_color_manual(
    values = c(
      "PBS" = "#5E8CC4",
      "AAV_IL2" = "#E1574C"
    ),
    labels = c(
      "PBS" = "AAV-PBS",
      "AAV_IL2" = "AAV-IL-2"
    )
  ) +
  theme_bw() +
  labs(
    x = "Top 20 CDR3β sequences",
    y = "CDR3β frequency in Treg repertoire (%)",
    color = "Group"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 7
    ),
    axis.ticks.x = element_line(),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "bottom"
  )


# ------------------------------------------------------------------------------
# 13. Summarise mean frequency for panel B
# ------------------------------------------------------------------------------

top20_lollipop_summary <- combined_top20_freq_long %>%
  filter(analysis_group %in% c("PBS", "AAV_IL2", "TRIPOD")) %>%
  left_join(
    top20_rank_table %>%
      select(
        composite_rank,
        cdr3_aa,
        composite_score
      ),
    by = c("aaCDR3" = "cdr3_aa")
  ) %>%
  group_by(
    composite_rank,
    aaCDR3,
    composite_score,
    analysis_group
  ) %>%
  summarise(
    n_samples = n_distinct(sample_id),
    n_detected = sum(detected, na.rm = TRUE),
    detection_percent = n_detected / n_samples * 100,
    mean_frequency_percent = mean(frequency_percent, na.rm = TRUE),
    median_frequency_percent = median(frequency_percent, na.rm = TRUE),
    max_frequency_percent = max(frequency_percent, na.rm = TRUE),
    pct_above_0.01 = sum(frequency_percent >= 0.01, na.rm = TRUE) /
      n_samples * 100,
    .groups = "drop"
  )


# ------------------------------------------------------------------------------
# 14. Create wide and long tables for panel B
# ------------------------------------------------------------------------------

top20_lollipop_wide <- top20_lollipop_summary %>%
  select(
    composite_rank,
    aaCDR3,
    composite_score,
    analysis_group,
    mean_frequency_percent,
    detection_percent,
    pct_above_0.01
  ) %>%
  pivot_wider(
    names_from = analysis_group,
    values_from = c(
      mean_frequency_percent,
      detection_percent,
      pct_above_0.01
    ),
    values_fill = 0
  ) %>%
  mutate(
    log2_aav_vs_pbs = log2(
      (mean_frequency_percent_AAV_IL2 + 1e-6) /
        (mean_frequency_percent_PBS + 1e-6)
    ),
    log2_aav_vs_tripod = log2(
      (mean_frequency_percent_AAV_IL2 + 1e-6) /
        (mean_frequency_percent_TRIPOD + 1e-6)
    ),
    cdr3_label = paste0(composite_rank, ". ", aaCDR3)
  ) %>%
  arrange(composite_rank)

top20_lollipop_long <- top20_lollipop_wide %>%
  select(
    composite_rank,
    aaCDR3,
    cdr3_label,
    composite_score,
    mean_frequency_percent_PBS,
    mean_frequency_percent_AAV_IL2,
    mean_frequency_percent_TRIPOD
  ) %>%
  pivot_longer(
    cols = c(
      mean_frequency_percent_PBS,
      mean_frequency_percent_AAV_IL2,
      mean_frequency_percent_TRIPOD
    ),
    names_to = "group",
    values_to = "mean_frequency_percent"
  ) %>%
  mutate(
    group = recode(
      group,
      mean_frequency_percent_PBS = "AAV-PBS",
      mean_frequency_percent_AAV_IL2 = "AAV-IL-2",
      mean_frequency_percent_TRIPOD = "TRIPOD untreated B6"
    ),
    group = factor(
      group,
      levels = c("TRIPOD untreated B6", "AAV-PBS", "AAV-IL-2")
    ),
    cdr3_label = factor(
      cdr3_label,
      levels = rev(top20_lollipop_wide$cdr3_label)
    )
  )

top20_lollipop_wide <- top20_lollipop_wide %>%
  mutate(
    cdr3_label = factor(
      cdr3_label,
      levels = rev(top20_lollipop_wide$cdr3_label)
    )
  )


# ------------------------------------------------------------------------------
# 15. Create line segments for panel B
# ------------------------------------------------------------------------------

segment_tripod_to_aav <- top20_lollipop_wide %>%
  transmute(
    cdr3_label,
    x_start = mean_frequency_percent_TRIPOD + 1e-6,
    x_end = mean_frequency_percent_AAV_IL2 + 1e-6,
    comparison = "TRIPOD to AAV-IL-2"
  )

segment_pbs_to_aav <- top20_lollipop_wide %>%
  transmute(
    cdr3_label,
    x_start = mean_frequency_percent_PBS + 1e-6,
    x_end = mean_frequency_percent_AAV_IL2 + 1e-6,
    comparison = "AAV-PBS to AAV-IL-2"
  )

segment_table <- bind_rows(
  segment_tripod_to_aav,
  segment_pbs_to_aav
)


# ------------------------------------------------------------------------------
# 16. Create panel B: mean-frequency lollipop plot
# ------------------------------------------------------------------------------

p_top20_lollipop_three_group <- ggplot() +
  geom_segment(
    data = segment_table,
    aes(
      y = cdr3_label,
      yend = cdr3_label,
      x = x_start,
      xend = x_end,
      linetype = comparison
    ),
    linewidth = 0.55,
    alpha = 0.65,
    color = "grey45"
  ) +
  geom_point(
    data = top20_lollipop_long,
    aes(
      x = mean_frequency_percent + 1e-6,
      y = cdr3_label,
      color = group
    ),
    size = 2.8,
    alpha = 0.95
  ) +
  scale_x_log10() +
  scale_color_manual(
    values = c(
      "TRIPOD untreated B6" = "grey45",
      "AAV-PBS" = "#5E8CC4",
      "AAV-IL-2" = "#E1574C"
    ),
    name = ""
  ) +
  scale_linetype_manual(
    values = c(
      "TRIPOD to AAV-IL-2" = "solid",
      "AAV-PBS to AAV-IL-2" = "dashed"
    ),
    name = "Comparison"
  ) +
  theme_bw() +
  labs(
    x = "Mean CDR3β frequency in Treg repertoire (%)",
    y = "Top 20 CDR3β sequences"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "bottom"
  )


# ------------------------------------------------------------------------------
# 17. Combine panels
# ------------------------------------------------------------------------------

p_panel_A <- p_top20_dot_all_layered +
  labs(
    title = NULL,
    subtitle = NULL
  ) +
  theme(
    legend.position = "bottom"
  )

p_panel_B <- p_top20_lollipop_three_group +
  labs(
    title = NULL,
    subtitle = NULL
  ) +
  theme(
    legend.position = "bottom"
  )

p_figure5_panel <- p_panel_A / p_panel_B +
  plot_layout(
    heights = c(1, 1.15),
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
    legend.position = "bottom"
  )


# ------------------------------------------------------------------------------
# 18. Save figure and output tables
# ------------------------------------------------------------------------------

ggsave(
  filename = file.path(
    figure_dir,
    "figure5_top20_il2_pbs_tripod_panel.pdf"
  ),
  plot = p_figure5_panel,
  width = 10,
  height = 10
)

ggsave(
  filename = file.path(
    figure_dir,
    "figure5_top20_il2_pbs_tripod_panel.png"
  ),
  plot = p_figure5_panel,
  width = 10,
  height = 10,
  dpi = 300
)

fwrite(
  top20_dot_all,
  file.path(table_dir, "figure5_top20_dot_plot_table.csv")
)

fwrite(
  top20_lollipop_summary,
  file.path(table_dir, "figure5_top20_lollipop_summary.csv")
)

fwrite(
  top20_lollipop_wide,
  file.path(table_dir, "figure5_top20_lollipop_wide.csv")
)


# ------------------------------------------------------------------------------
# 19. Print quick checks
# ------------------------------------------------------------------------------

message("Figure 5 complete.")
message("Figure saved to: ", figure_dir)
message("Tables saved to: ", table_dir)

message("Top 20 CDR3β sequences:")
print(top20_rank_table)

message("Number of samples per group:")
print(
  combined_top20_freq_long %>%
    group_by(analysis_group) %>%
    summarise(
      n_samples = n_distinct(sample_id),
      n_rows = n(),
      .groups = "drop"
    )
)