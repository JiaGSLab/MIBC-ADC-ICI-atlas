# MIBC HER2-ADC plus anti-PD-1 project

This repository contains analysis and figure-generation scripts for the manuscript:

**Tumor plasticity and spatial immune architecture shape clinical response to HER2-ADC plus anti-PD-1 in muscle-invasive bladder cancer**

Study: multimodal analysis of HER2-positive muscle-invasive bladder cancer treated with neoadjuvant disitamab vedotin plus tislelizumab.

## Directory overview

- `scRNAseq/`: scFFPE-seq and Visium HD analyses, together with scripts for Figures 1–5 and related supplementary figures.
- `fig6-spatial-analysis/`: multiplex-immunohistochemistry image processing, cell phenotyping, cellular-neighbourhood analyses and response modelling for Figure 6. See the subdirectory README and source manifest for details.
- `Fig7_HE_prediction/`: H&E response prediction and regional phenotype analyses using frozen Prov-GigaPath embeddings for Figure 7.

## Environments and data

Environment files are provided in the relevant subdirectories. Raw clinical, sequencing and whole-slide-image data are not included in this repository; controlled-access sequencing data are described in the manuscript. Some scripts retain project-specific input paths and should be configured for the local data location before execution.

For analysis details, sample inclusion criteria and statistical methods, refer to the manuscript and the README/source-manifest files within each analysis directory.

# Citation

Tumor plasticity and spatial immune architecture shape clinical response to HER2-ADC plus anti-PD-1 in muscle-invasive bladder cancer; Wenhui Chen, Yan Mei, Lingmin Kong, Weibiao Wang, Huoming Li, Yichen Li, Yanjie Shen, Yuhong Wang, Zhenhuan Feng, Cheng Ma, Jinxin Guo, Liu Huang, Junjie He,Meiqin Li, Yan Guo, Junxing Chen, Guangshuai Jia, Huanjun Wang; In Submission
