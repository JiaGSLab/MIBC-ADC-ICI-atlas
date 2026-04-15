# MIBC-ADC-ICI-atlas scRNA-seq analysis scripts

## Overview

This repository contains the organized analysis scripts used for scRNA-seq, Visium HD spatial transcriptomics, and figure generation.


## Directory structure

```
├── 00_utils/
│   ├── plot_helpers.R                     # Shared plotting helper functions
│   └── setup_scRNAseq_environment.R       # Package loading, shared variables, colors, markers, and utility functions
│
├── 01_preprocessing/
│   ├── 01_merge_qc_integration_major_celltypes.R   # QC, integration, clustering, and major cell type assignment
│   └── 02_reintegrate_annotated_objects.R          # Recombine annotated lineage-level objects for downstream analyses
│
├── 02_subclustering/
│   ├── 01_bcell_subclustering_annotation.R         # B cell and plasma cell subclustering and annotation
│   ├── 02_tcell_subclustering_annotation.R         # T cell subclustering and annotation
│   ├── 03_myeloid_subclustering_annotation.R       # Myeloid cell subclustering and annotation
│   └── 04_stromal_subclustering_annotation.R       # Stromal cell subclustering and annotation
│
├── 03_downstream/
│   ├── cellchat_analysis.R                 # CellChat-based cell–cell communication analysis
│   ├── celltype_proportion_analysis.R      # Cell type proportion analysis across groups
│   ├── copykat_per_sample.R                # Sample-wise CopyKAT analysis for malignant cell inference
│   ├── deg_mp4_high_pre_treatment.R        # Differential expression analysis in pre-treatment MP4-high tumour cells
│   ├── kegg_enrichment_plot.R              # KEGG enrichment plotting
│   ├── mast_analysis_malignant.R           # MAST-based differential expression analysis in malignant cells
│   ├── nmf_analysis_malignant_multi.R      # NMF-based malignant-cell meta-program analysis
│   └── progeny_metabolism.R                # PROGENy and scMetabolism pathway analysis
│
├── 04_figures/
│   ├── main/
│   │   ├── Fig1C_1D.R                      # Main figure 1 plots
│   │   ├── Fig2.R                          # Main figure 2 plots
│   │   ├── Fig3.R                          # Main figure 3 plots
│   │   ├── Fig5.R                          # Main figure 5 plots
│   │   └── Figure4.R                       # Main figure 4 plots
│   │
│   └── supplementary/
│       ├── FigS1.R                         # Supplementary figure S1 plots
│       ├── FigS2.R                         # Supplementary figure S2 plots
│       ├── FigS3.R                         # Supplementary figure S3 plots
│       ├── FigS4.R                         # Supplementary figure S4 plots
│       └── FigS5.R                         # Supplementary figure S5 plots
│
├── 05_visiumHD/
│   ├── 00_visium_hd_figure_entry_points.R          # Connect Visium HD outputs to downstream figure scripts
│   ├── 01_build_rctd_reference_from_scrna.R        # Build RCTD reference from scRNA-seq data
│   ├── 02_run_rctd_slide01.R                       # Run RCTD deconvolution for slide 1
│   ├── 03_run_rctd_slide02.R                       # Run RCTD deconvolution for slide 2
│   ├── 04_run_rctd_slide03.R                       # Run RCTD deconvolution for slide 3
│   ├── 05_attach_rctd_to_visium_objects_and_merge.R # Add RCTD annotations and merge Visium HD objects
│   └── 06_visium_hd_spatial_plots_and_celltype_highlights.R # Generate spatial plots and highlighted cell-type views
│
└──README.md 

```

## Citation

Tumour plasticity and spatial immune architecture define response to HER2-ADC plus anti–PD-1 in muscle-invasive bladder cancer; Wenhui Chen, Yan Mei, Lingmin Kong, Weibiao Wang, Huoming Li, Yichen Li, Yanjie Shen, Yuhong Wang, Zhenhuan Feng, Cheng Ma, Jinxin Guo, Liu Huang, Junjie He,Meiqin Li, Yan Guo, Junxing Chen, Guangshuai Jia, Huanjun Wang; In Submission