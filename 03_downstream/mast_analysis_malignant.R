# ==================================================
# Script: MAST differential expression for malignant cells
# Author: Wenhui Chen
# Purpose: Performs MAST differential expression testing between specified malignant cell groups.
# Input: Main Seurat or summary objects already generated in the project, plus any file paths referenced in the script.
# Output: Intermediate result files and/or figure files written to output_dir.
# Note: Review file paths before running. Keep original analysis logic unchanged.
# ==================================================
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(muscatWrapper)
library(Seurat)
library(qs2)
library(future)
library(future.apply)
options(future.globals.maxSize = 800 * 1024^3)
set.seed(123)
output_dir = "~/scRNAseq/GPU_IT/QC"
if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
setwd(output_dir)

malignant.obj = qs_read("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/malignant_obj_intergration.qs")
Idents(malignant.obj) = "Group"
DefaultAssay(malignant.obj) = "RNA"
# DEG.Pre_NR_vs_Pre_CR = FindMarkers(malignant.obj,
#                                   ident.1= "Pre_NR",
#                                   ident.2 = "Pre_CR",
#                                   test.use="MAST",
#                                   logfc.threshold=0.05,
#                                   min.pct = 0.02,
#                                   latent.vars=c("sample"),
#                                   random.seed=123)
# qs_save(DEG.Pre_NR_vs_Pre_CR,"DEG.Pre_NR_vs_Pre_CR_MAST.qs")

DEG.Pre_NR_vs_Post_NR = FindMarkers(malignant.obj,
                                  ident.1= "Pre_NR",
                                  ident.2 = "Post_NR",
                                  test.use="MAST",
                                  logfc.threshold=0.05,
                                  min.pct = 0.02,
                                  latent.vars=c("sample"),
                                  random.seed=123)

qs_save(DEG.Pre_NR_vs_Post_NR,"DEG.Pre_NR_vs_Post_NR_MAST.qs")