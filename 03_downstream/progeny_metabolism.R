# ==================================================
# Script: PROGENy and metabolism scoring
# Author: Wenhui Chen
# Purpose: Calculates pathway activity and metabolic activity scores for the integrated Seurat object.
# Input: Main Seurat or summary objects already generated in the project, plus any file paths referenced in the script.
# Output: Intermediate result files and/or figure files written to output_dir.
# Note: Review file paths before running. Keep original analysis logic unchanged.
# ==================================================
library(Seurat)
library(dplyr)
library(progeny)
library(scMetabolism)
library(qs2)
output_dir = "~/scRNAseq/GPU_IT/QC"
setwd(output_dir)

seu.combined.IT = qs_read(file.path(output_dir, 'seu_combined_IT_filtered_obj_intergration.qs'))
DefaultAssay(seu.combined.IT) = 'RNA'
seu.combined.IT.copy=seu.combined.IT
seu.combined.IT.copy[['RNA']] = as(seu.combined.IT.copy[['RNA']], 'Assay')
seu.combined.IT.copy<-sc.metabolism.Seurat(obj = seu.combined.IT.copy, method = "AUCell", imputation = F, ncores = 24, metabolism.type = "KEGG")
# meta.obj = CreateSeuratObject(counts = seu.combined.IT.copy@assays$METABOLISM$score,
#                             meta.data = seu.combined.IT.copy@meta.data, min.cells=0,min.features=0)
qs_save(seu.combined.IT.copy,file.path(output_dir, 'seu_combined_IT_filtered_obj_intergration_metabolism.qs'))

seu.combined.IT.progeny = progeny(seu.combined.IT, scale=FALSE, organism="Human", top=500, perm=1, return_assay = TRUE)
qs_save(seu.combined.IT.progeny,file.path(output_dir, 'seu_combined_IT_filtered_obj_intergration_progeny.qs'))