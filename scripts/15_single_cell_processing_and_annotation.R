# ==============================================================================
# 15_single_cell_processing_and_annotation.R
# ==============================================================================
# Goal:
#   Process the merged Curie003 single-cell Treg Seurat object, perform QC,
#   normalization, Harmony integration, clustering, manual cell-state annotation,
#   and TCR metadata processing.
#
#   This script produces the annotated Seurat objects used for downstream
#   single-cell figures:
#     Figure 7: TCR chain pairing among top CDR3β sequences
#     Figure 8: marker gene dot plot
#     Figure 9: UMAP of Treg states and top CDR3β sequences
#     Figure 10: top-100 CDR3β distribution across Treg states
#
# Input:
#   data/processed/single_cell/seurat_objects/merged_seurat_with_vdj_metadata.rds
#
# Outputs:
#   data/processed/single_cell/seurat_objects/soht_annotated_all_cells.rds
#   data/processed/single_cell/seurat_objects/soht_annotated_tcr_filtered.rds
#   results/tables/single_cell_qc_summary.csv
#   results/tables/single_cell_cell_type_counts.csv
#   results/tables/single_cell_cell_type_distribution_by_sample.csv
#   results/tables/single_cell_chain_pairing_summary.csv
#   results/tables/single_cell_tcr_filtered_metadata.csv
#   results/tables/single_cell_beta_clonotype_counts.csv
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Load packages
# ------------------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(harmony)
library(data.table)
library(tibble)


# ------------------------------------------------------------------------------
# 1. Define input and output paths
# ------------------------------------------------------------------------------

input_rds <- file.path(
  "data",
  "processed",
  "single_cell",
  "seurat_objects",
  "merged_seurat_with_vdj_metadata.rds"
)

seurat_output_dir <- file.path(
  "data",
  "processed",
  "single_cell",
  "seurat_objects"
)

table_dir <- file.path(
  "results",
  "tables"
)

dir.create(seurat_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# 2. Analysis settings
# ------------------------------------------------------------------------------

qc_min_features <- 300
qc_max_features <- 6300
qc_max_percent_mt <- 20

harmony_group_variable <- "sample"
dims_to_use <- 1:14
cluster_resolution <- 0.2
cluster_algorithm <- 4

contaminant_cell_types <- c(
  "B/plasma contamination",
  "Myeloid/NK contamination"
)


# ------------------------------------------------------------------------------
# 3. Load merged Seurat object
# ------------------------------------------------------------------------------

if (!file.exists(input_rds)) {
  stop("Input Seurat object not found: ", input_rds)
}

merged_seurat <- readRDS(input_rds)

DefaultAssay(merged_seurat) <- "RNA"


# ------------------------------------------------------------------------------
# 4. Check required metadata columns
# ------------------------------------------------------------------------------

required_metadata_cols <- c(
  "sample",
  "tra_cdr3",
  "trb_cdr3"
)

missing_metadata_cols <- setdiff(
  required_metadata_cols,
  colnames(merged_seurat@meta.data)
)

if (length(missing_metadata_cols) > 0) {
  stop(
    "Missing required metadata columns in the merged Seurat object: ",
    paste(missing_metadata_cols, collapse = ", ")
  )
}


# ------------------------------------------------------------------------------
# 5. Quality control
# ------------------------------------------------------------------------------

merged_seurat[["percent.mt"]] <- PercentageFeatureSet(
  merged_seurat,
  pattern = "^mt-"
)

qc_before <- data.frame(
  step = "before_qc",
  n_cells = ncol(merged_seurat),
  n_genes = nrow(merged_seurat),
  median_nFeature_RNA = median(merged_seurat$nFeature_RNA, na.rm = TRUE),
  median_nCount_RNA = median(merged_seurat$nCount_RNA, na.rm = TRUE),
  median_percent_mt = median(merged_seurat$percent.mt, na.rm = TRUE)
)

so <- subset(
  merged_seurat,
  subset = nFeature_RNA > qc_min_features &
    nFeature_RNA < qc_max_features &
    percent.mt < qc_max_percent_mt
)

qc_after <- data.frame(
  step = "after_qc",
  n_cells = ncol(so),
  n_genes = nrow(so),
  median_nFeature_RNA = median(so$nFeature_RNA, na.rm = TRUE),
  median_nCount_RNA = median(so$nCount_RNA, na.rm = TRUE),
  median_percent_mt = median(so$percent.mt, na.rm = TRUE)
)

qc_summary <- bind_rows(qc_before, qc_after)

fwrite(
  qc_summary,
  file.path(table_dir, "single_cell_qc_summary.csv")
)


# ------------------------------------------------------------------------------
# 6. Join layers, normalize, integrate, and cluster
# ------------------------------------------------------------------------------

# Seurat v5 merged objects often store one counts layer per sample.
# JoinLayers combines these layers before normalization and downstream analysis.

soht <- JoinLayers(so)

soht <- NormalizeData(
  soht,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

soht <- FindVariableFeatures(soht)

soht <- ScaleData(
  soht,
  features = VariableFeatures(soht)
)

soht <- RunPCA(
  soht,
  features = VariableFeatures(soht)
)

soht <- RunHarmony(
  soht,
  group.by.vars = harmony_group_variable
)

soht <- RunUMAP(
  soht,
  reduction = "harmony",
  dims = dims_to_use
)

soht <- FindNeighbors(
  soht,
  reduction = "harmony",
  dims = dims_to_use
)

soht <- FindClusters(
  soht,
  resolution = cluster_resolution,
  algorithm = cluster_algorithm,
  verbose = FALSE
)


# ------------------------------------------------------------------------------
# 7. Assign cell-state annotations
# ------------------------------------------------------------------------------

# Adjust this map if cluster numbering changes after rerunning clustering.

cluster_annotations <- c(
  "1" = "Resting memory-like Tregs",
  "2" = "Activated memory-like Tregs",
  "3" = "Activated checkpoint-high Tregs",
  "4" = "Cytotoxic-like Tregs",
  "5" = "CD25-high suppressive Tregs",
  "6" = "Proliferating Tregs",
  "7" = "B/plasma contamination",
  "8" = "Myeloid/NK contamination"
)

soht$seurat_cluster_original <- as.character(soht$seurat_clusters)

soht$cell_type <- unname(
  cluster_annotations[as.character(soht$seurat_clusters)]
)

if (any(is.na(soht$cell_type))) {
  unmapped_clusters <- sort(unique(soht$seurat_cluster_original[is.na(soht$cell_type)]))
  
  stop(
    "Some Seurat clusters were not assigned cell_type labels: ",
    paste(unmapped_clusters, collapse = ", "),
    "\nUpdate cluster_annotations in this script."
  )
}

Idents(soht) <- "cell_type"


# ------------------------------------------------------------------------------
# 8. Export cell-type count summaries
# ------------------------------------------------------------------------------

cell_type_counts <- soht@meta.data %>%
  count(cell_type, name = "n_cells") %>%
  arrange(desc(n_cells))

fwrite(
  cell_type_counts,
  file.path(table_dir, "single_cell_cell_type_counts.csv")
)

cell_type_distribution_by_sample <- soht@meta.data %>%
  count(sample, cell_type, name = "n_cells") %>%
  group_by(sample) %>%
  mutate(
    sample_total_cells = sum(n_cells),
    percent_cells = n_cells / sample_total_cells * 100
  ) %>%
  ungroup() %>%
  arrange(sample, cell_type)

fwrite(
  cell_type_distribution_by_sample,
  file.path(table_dir, "single_cell_cell_type_distribution_by_sample.csv")
)


# ------------------------------------------------------------------------------
# 9. Add TCR chain-pairing metadata
# ------------------------------------------------------------------------------

count_cdr3 <- function(x) {
  x <- as.character(x)
  
  sapply(x, function(y) {
    if (is.na(y) || y == "" || y == "NA") {
      return(0)
    }
    
    parts <- unlist(strsplit(y, split = ";|,|\\|"))
    parts <- trimws(parts)
    parts <- parts[parts != "" & parts != "NA"]
    
    length(unique(parts))
  })
}

soht$n_alpha <- unname(count_cdr3(soht$tra_cdr3))
soht$n_beta <- unname(count_cdr3(soht$trb_cdr3))

soht$chain_pairing <- case_when(
  soht$n_alpha == 0 & soht$n_beta == 0 ~ "no TCR",
  soht$n_alpha == 1 & soht$n_beta == 1 ~ "single pair",
  soht$n_alpha == 0 & soht$n_beta == 1 ~ "orphan VDJ",
  soht$n_alpha == 1 & soht$n_beta == 0 ~ "orphan VJ",
  soht$n_alpha > 1 & soht$n_beta == 1 ~ "extra VJ",
  soht$n_alpha == 1 & soht$n_beta > 1 ~ "extra VDJ",
  soht$n_alpha > 1 & soht$n_beta > 1 ~ "two full chains",
  TRUE ~ NA_character_
)

chain_pairing_summary <- soht@meta.data %>%
  count(chain_pairing, name = "n_cells") %>%
  mutate(
    percent_cells = n_cells / sum(n_cells) * 100
  ) %>%
  arrange(desc(n_cells))

fwrite(
  chain_pairing_summary,
  file.path(table_dir, "single_cell_chain_pairing_summary.csv")
)


# ------------------------------------------------------------------------------
# 10. Create TCR-filtered object
# ------------------------------------------------------------------------------

# Downstream beta-chain analysis uses cells with exactly one detected TCRβ CDR3.
# Non-Treg contaminating clusters are removed.

soht_tcr <- subset(
  soht,
  subset = n_beta == 1 &
    !cell_type %in% contaminant_cell_types
)

if (ncol(soht_tcr) == 0) {
  stop("No cells remained after TCR filtering.")
}


# ------------------------------------------------------------------------------
# 11. Add beta-chain clonotype and clone-size metadata
# ------------------------------------------------------------------------------

soht_tcr$beta_cdr3 <- as.character(soht_tcr$trb_cdr3)

soht_tcr$beta_cdr3[
  is.na(soht_tcr$beta_cdr3) |
    soht_tcr$beta_cdr3 == "" |
    soht_tcr$beta_cdr3 == "NA"
] <- NA

beta_clone_counts <- soht_tcr@meta.data %>%
  filter(!is.na(beta_cdr3)) %>%
  count(beta_cdr3, name = "beta_clone_count") %>%
  arrange(desc(beta_clone_count))

soht_tcr@meta.data <- soht_tcr@meta.data %>%
  rownames_to_column("cell_barcode") %>%
  left_join(
    beta_clone_counts,
    by = "beta_cdr3"
  ) %>%
  column_to_rownames("cell_barcode")

soht_tcr$cloneSize_beta <- case_when(
  is.na(soht_tcr$beta_clone_count) ~ "No beta",
  soht_tcr$beta_clone_count == 1 ~ "Single",
  soht_tcr$beta_clone_count > 1 & soht_tcr$beta_clone_count <= 5 ~ "Small",
  soht_tcr$beta_clone_count > 5 & soht_tcr$beta_clone_count <= 20 ~ "Medium",
  soht_tcr$beta_clone_count > 20 & soht_tcr$beta_clone_count <= 100 ~ "Large",
  soht_tcr$beta_clone_count > 100 ~ "Hyperexpanded",
  TRUE ~ NA_character_
)

soht_tcr$cloneSize_beta <- factor(
  soht_tcr$cloneSize_beta,
  levels = c(
    "Single",
    "Small",
    "Medium",
    "Large",
    "Hyperexpanded",
    "No beta"
  )
)

soht_tcr$beta_cdr3_length <- nchar(as.character(soht_tcr$beta_cdr3))


# ------------------------------------------------------------------------------
# 12. Export TCR-filtered summaries
# ------------------------------------------------------------------------------

fwrite(
  beta_clone_counts,
  file.path(table_dir, "single_cell_beta_clonotype_counts.csv")
)

tcr_filtered_metadata <- soht_tcr@meta.data %>%
  rownames_to_column("cell_barcode")

fwrite(
  tcr_filtered_metadata,
  file.path(table_dir, "single_cell_tcr_filtered_metadata.csv")
)

tcr_filtered_cell_type_counts <- soht_tcr@meta.data %>%
  count(cell_type, name = "n_cells") %>%
  arrange(desc(n_cells))

fwrite(
  tcr_filtered_cell_type_counts,
  file.path(table_dir, "single_cell_tcr_filtered_cell_type_counts.csv")
)

tcr_filtered_clone_size_summary <- soht_tcr@meta.data %>%
  count(cloneSize_beta, name = "n_cells") %>%
  mutate(
    percent_cells = n_cells / sum(n_cells) * 100
  )

fwrite(
  tcr_filtered_clone_size_summary,
  file.path(table_dir, "single_cell_tcr_filtered_clone_size_summary.csv")
)


# ------------------------------------------------------------------------------
# 13. Save annotated Seurat objects
# ------------------------------------------------------------------------------

saveRDS(
  soht,
  file.path(seurat_output_dir, "soht_annotated_all_cells.rds")
)

saveRDS(
  soht_tcr,
  file.path(seurat_output_dir, "soht_annotated_tcr_filtered.rds")
)


# ------------------------------------------------------------------------------
# 14. Print quick checks
# ------------------------------------------------------------------------------

message("Single-cell processing and annotation complete.")

message("Annotated all-cell object saved to: ",
        file.path(seurat_output_dir, "soht_annotated_all_cells.rds"))

message("TCR-filtered object saved to: ",
        file.path(seurat_output_dir, "soht_annotated_tcr_filtered.rds"))

message("QC summary:")
print(qc_summary)

message("Cell type counts:")
print(cell_type_counts)

message("Chain pairing summary:")
print(chain_pairing_summary)

message("TCR-filtered clone size summary:")
print(tcr_filtered_clone_size_summary)