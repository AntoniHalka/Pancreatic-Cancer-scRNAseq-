# Pancreatic-Cancer-scRNAseq
Pipeline for analyzing pancreatic ductal adenocarcinoma (PDAC) single-cell RNA sequencing data across disease stages, tissue types, and treatment conditions.

## Introduction
This repository contains the R script used to preprocess and analyze scRNA-seq data from three publicly available GEO datasets (GSE205013, GSE154778, GSE148073). The code is written in R and was run in version 4.5.3 using Seurat v5.4.0.

The analysis includes quality control, normalization, batch integration with Harmony, clustering, cell type annotation, differential gene expression, pathway enrichment, and figure generation including UMAP and chord diagram plots.

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
- ah_pancreas_preprocessing_Spring2026.R code used to preprocess raw count matrices, including quality control and normalization, as well as cell type annotation, differential gene expression analysis, pathway enrichment analysis using MSigDB Hallmark gene sets, and COMET input preparation; code to generate chord plots, UMAPs, and feature plots
- Raw count matrices metadata includes age, sex, stage, treatment status, treatment type, tissue site, and tissue type (primary versus metastatic)
