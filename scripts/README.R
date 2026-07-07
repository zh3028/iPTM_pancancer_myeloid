# Scripts for iPTM Module Validation in Pan-Cancer Myeloid Cells

## File Descriptions

| Script | Description |
| :--- | :--- |
| 01_load_data.R | Load GSE154763_MYE raw data and create Seurat object |
| 02_build_seurat.R | QC filtering, SCTransform normalization, PCA, UMAP, clustering |
| 03_iptm_scoring.R | iPTM module scoring using AddModuleScore |
| 04_deg_enrichment.R | Differential expression analysis (M2-High vs M2-Low) |
| 05_cellchat.R | Cell-cell communication analysis using CellChat |
| 06_pathway_analysis.R | Hallmark pathway enrichment analysis |

## Requirements

R version 4.5.3 or higher with the following packages:
- Seurat (v5.0.3)
- tidyverse
- ggplot2
- pheatmap
- clusterProfiler
- org.Hs.eg.db
- CellChat (v1.6.1)
- msigdbr
