# ==================================================
# Script: Run RCTD for Visium HD slide02
# Author: Wenhui Chen
# Purpose: Runs spacexr RCTD on Visium HD slide02 using the prebuilt scRNA-seq reference.
# Input: Saved SCreference object and Visium HD slide02 binned output directory.
# Output: RCTD result object for slide02.
# Note: Review file paths before running. Keep original analysis logic unchanged.
# ==================================================
library(arrow)
library(spacexr)
library(ggplot2)
library(Seurat)
library(dplyr)
library(forcats)
options(future.globals.maxSize = 100 * 1024^3)
set.seed(123)

load("~/visiumHD/SCreference_celltype.RData")

VisiumData_tissue <- read.VisiumSpatialRNA("~/visiumHD/FWSP20252181/FWSP20252181-2/binned_outputs/square_008um/")
#non_mito_idx <- !grepl("^MT-", rownames(VisiumData_tissue@counts))
#VisiumData_tissue@counts <- VisiumData_tissue@counts[non_mito_idx, ]
Sys.setenv("OPENBLAS_NUM_THREADS"=10)
myRCTD <- create.RCTD(VisiumData_tissue, SCreference, max_cores = 10,UMI_min = 100)
myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")
save(myRCTD, file = "~/visiumHD/RCTD_tissue2_celltype.RData")
