# treg-il2-tcr-repertoire-analysis
Analysis scripts for my M2 thesis project investigating regulatory T cell (Treg) TCRβ repertoire changes following long-term AAV-IL-2 treatment in C57BL/6 mice.

This repository contains the R scripts used for bulk TCRβ repertoire analysis, CDR3β prioritisation, GLIPH2 sequence-similarity analysis, antigen database matching, and integration with single-cell RNA-seq/TCR-seq data.

Project overview

The goal of this project was to assess whether sustained IL-2 expression reshapes the Treg TCRβ repertoire and to identify recurrent or enriched CDR3β sequences that may represent selectively expanded Treg clonotypes.

The analysis includes:

* bulk TCRβ repertoire diversity and clonotype frequency analysis
* treatment-independent and IL-2-prioritised composite scoring of shared expanded Treg CDR3β sequences
* comparison with TRiPoD untreated B6 Treg reference repertoires
* GLIPH2 clustering of sequence-similar CDR3β motifs
* antigen database matching of top-ranked CDR3β sequences
* single-cell Treg state annotation and integration with TCR metadata

Repository structure

scripts/
  01_create_repseq_objects_from_mixcr.R
  02_check_bulk_tcr_depth_and_metadata.R
  03_figure1_shannon_diversity_treg_tconv.R
  04_calculate_bulk_treg_clonotype_frequencies.R
  05_figure2_treg_frequency_bins.R
  06_calculate_treatment_independent_composite_score.R
  07_calculate_il2_prioritised_composite_score.R
  08_figure3_compare_composite_scores.R
  09_create_ld1_cytoscape_network_top100_il2.R
  10_match_top100_il2_cdr3_to_antigen_database.R
  11_figure5_top20_il2_pbs_tripod.R
  12_prepare_and_run_gliph2_treg_tripod_reference.R
  13_create_table4_gliph2_top100_cluster_summary.R
  14_figure6_gliph2_network_and_logos.R
  15_single_cell_processing_and_annotation.R
  16_figure7_chain_pairing_top20.R
  17_figure8_marker_dotplot.R
  18_figure9_umap_treg_states_top_cdr3.R
  19_figure10_top100_cluster_distribution.R

Data availability

Raw sequencing data, processed RepSeqData objects, Seurat objects, and external reference databases are not included in this repository because they may be large or subject to data-sharing restrictions.
