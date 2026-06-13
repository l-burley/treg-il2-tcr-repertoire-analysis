# ==============================================================================
# 17_figure8_marker_dotplot.R
# ==============================================================================
# Goal:
#   Create Figure 8 showing marker gene expression across annotated single-cell
#   Treg states.
#
#   Dot size represents the percentage of cells expressing each marker.
#   Dot colour represents scaled average expression.
#
# Input:
#   data/processed/single_cell/seurat_objects/soht_annotated_all_cells.rds
#
# Outputs:
#   results/figures/figure8_marker_dotplot.png
#   results/figures/figure8_marker_dotplot.pdf
#   results/tables/figure8_marker_dotplot_values.csv
#   results/tables/figure8_missing_markers.csv
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Load packages
# ------------------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(data.table)
library(ggplot2)


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

figure_dir <- file.path("results", "figures")
table_dir <- file.path("results", "tables")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# 2. Load annotated Seurat object
# ------------------------------------------------------------------------------

if (!file.exists(seurat_path)) {
  stop("Annotated Seurat object not found: ", seurat_path)
}

soht <- readRDS(seurat_path)

DefaultAssay(soht) <- "RNA"


# ------------------------------------------------------------------------------
# 3. Check required metadata
# ------------------------------------------------------------------------------

if (!"cell_type" %in% colnames(soht@meta.data)) {
  stop("The metadata column 'cell_type' was not found in the Seurat object.")
}

Idents(soht) <- "cell_type"


# ------------------------------------------------------------------------------
# 4. Define marker gene panel
# ------------------------------------------------------------------------------

markers <- c(
  "Cd4", "Cd8a", "Cd8b1", "Foxp3", "Il2ra",
  "Igfbp4", "Tgfbr3", "Myo6", "Cdkn1a",
  "Ccr7", "Sell", "Tcf7", "Il7r", "Cxcr3",
  "Gzmb", "Thy1", "Trat1", "Ccl5",
  "S100a6", "S100a13", "S100a4",
  "Cd7", "Cd81", "Tigit", "Tnfrsf9", "Lag3",
  "Coro2a", "Fos", "Jun", "Zbtb20", "Mctp1", "Lrba", "Cd69",
  "Stmn1", "Pclaf", "Mki67", "Birc5",
  "Igkc", "Jchain", "Cd74", "Iglc2",
  "Ctla4", "Tbx21", "Icos", "Klf2", "Tnfrsf18",
  "Tox", "Ikzf2", "Lyz2"
)

markers_present <- markers[markers %in% rownames(soht)]
markers_missing <- setdiff(markers, rownames(soht))

if (length(markers_present) == 0) {
  stop("None of the requested marker genes were found in the Seurat object.")
}

missing_marker_table <- data.frame(
  marker = markers_missing
)

fwrite(
  missing_marker_table,
  file.path(table_dir, "figure8_missing_markers.csv")
)


# ------------------------------------------------------------------------------
# 5. Set cell type order
# ------------------------------------------------------------------------------

cell_type_order <- c(
  "Activated checkpoint-high Tregs",
  "CD25-high suppressive Tregs",
  "Resting memory-like Tregs",
  "Activated memory-like Tregs",
  "Cytotoxic-like Tregs",
  "B/plasma contamination",
  "Proliferating Tregs",
  "Myeloid/NK contamination"
)

cell_type_order_present <- cell_type_order[
  cell_type_order %in% unique(as.character(soht$cell_type))
]

remaining_cell_types <- setdiff(
  unique(as.character(soht$cell_type)),
  cell_type_order_present
)

final_cell_type_order <- c(
  cell_type_order_present,
  sort(remaining_cell_types)
)

soht$cell_type <- factor(
  as.character(soht$cell_type),
  levels = final_cell_type_order
)

Idents(soht) <- "cell_type"


# ------------------------------------------------------------------------------
# 6. Create marker dot plot
# ------------------------------------------------------------------------------

p_marker_dotplot <- DotPlot(
  soht,
  features = markers_present,
  group.by = "cell_type",
  dot.scale = 6
) +
  scale_color_gradient2(
    low = "#2F6FAE",
    mid = "white",
    high = "#B13A35",
    midpoint = 0
  ) +
  RotatedAxis() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "plain",
      size = 8
    ),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ) +
  labs(
    x = "Features",
    y = "Identity"
  )


# ------------------------------------------------------------------------------
# 7. Export dot plot values
# ------------------------------------------------------------------------------

dotplot_values <- p_marker_dotplot$data %>%
  select(
    cell_type = id,
    gene = features.plot,
    percent_expressed = pct.exp,
    average_expression = avg.exp,
    scaled_average_expression = avg.exp.scaled
  ) %>%
  mutate(
    gene = factor(gene, levels = markers_present),
    cell_type = factor(cell_type, levels = final_cell_type_order),
    percent_expressed = round(percent_expressed, 1),
    average_expression = round(average_expression, 3),
    scaled_average_expression = round(scaled_average_expression, 3)
  ) %>%
  arrange(cell_type, gene)

fwrite(
  dotplot_values,
  file.path(table_dir, "figure8_marker_dotplot_values.csv")
)


# ------------------------------------------------------------------------------
# 8. Save Figure 8
# ------------------------------------------------------------------------------

ggsave(
  filename = file.path(figure_dir, "figure8_marker_dotplot.pdf"),
  plot = p_marker_dotplot,
  width = 13,
  height = 7
)

ggsave(
  filename = file.path(figure_dir, "figure8_marker_dotplot.png"),
  plot = p_marker_dotplot,
  width = 13,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 9. Print quick checks
# ------------------------------------------------------------------------------

message("Figure 8 marker dot plot complete.")
message("Figure saved to: ", figure_dir)
message("Dot plot values saved to: ", table_dir)

message("Markers requested: ", length(markers))
message("Markers present: ", length(markers_present))
message("Markers missing: ", length(markers_missing))

if (length(markers_missing) > 0) {
  message("Missing markers:")
  print(markers_missing)
}

message("Cell type order:")
print(final_cell_type_order)