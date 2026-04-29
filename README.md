# Pancreatic-Cancer-scRNAseq
Pipeline for analyzing pancreatic ductal adenocarcinoma (PDAC) single-cell RNA sequencing data across disease stages, tissue types, and treatment conditions.

## Introduction
This repository contains both an R script and a COMET-based analysis workflow used to process and analyze scRNA-seq data from three publicly available GEO datasets (GSE205013, GSE154778, GSE148073). The R code is written in R and was run in version 4.5.3 using Seurat v5.4.0. In addition, a COMET execution pipeline was performed on processed Seurat-derived outputs to identify marker genes and higher-order gene combinations distinguishing cellular states. COMET analysis was performed across multiple condition comparisons (including metastatic vs primary, treated vs untreated, and early vs advanced), using identical parameter settings unless otherwise specified.

The analysis includes quality control, normalization, batch integration with Harmony, clustering, cell type annotation, differential gene expression, pathway enrichment, and figure generation including UMAP and chord diagram plots. The COMET component specifically evaluates cluster-specific gene enrichment patterns using hypergeometric testing and XL-mHG scoring across single genes and gene combinations.

## RStudio and package requirements
- RStudio (R v4.5.3)
- Seurat (v5.4.0)
- SeuratObject
- ggplot2
- harmony
- clustree
- dplyr
- tidyverse
- data.table
- glmGamPoi
- clusterProfiler
- org.Hs.eg.db
- msigdbr
- circlize
- ggrepel
- reshape2
- hdf5r
- plyr
- BiocParallel

## Content summary
- ah_pancreas_preprocessing_Spring2026.R code used to preprocess raw count matrices, including quality control and normalization, as well as cell type annotation, differential gene expression analysis, pathway enrichment analysis using MSigDB Hallmark gene sets, and COMET input preparation. This also includes code to generate chord plots, UMAPs, and feature plots
- COMET_met_vs_prim_pipeline.ipynb notebook documenting the COMET analysis workflow for metastatic versus primary pancreatic cancer comparison. This includes execution of COMET on preprocessed Seurat-derived inputs (marker matrix, UMAP embeddings, cluster assignments, and gene list), evaluation of single-gene and gene-combination enrichment across clusters, and generation of statistical outputs
- Raw count matrices metadata includes age, sex, stage, treatment status, treatment type, tissue site, and tissue type (primary versus metastatic)
