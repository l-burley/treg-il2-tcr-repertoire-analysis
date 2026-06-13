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


Data availability

Raw sequencing data, processed RepSeqData objects, Seurat objects, and external reference databases are not included in this repository because they may be large or subject to data-sharing restrictions.
