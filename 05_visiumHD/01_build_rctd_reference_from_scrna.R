# ==================================================
# Script: Build RCTD reference from scRNA-seq
# Author: Wenhui Chen
# Purpose: Builds the spacexr reference object from the annotated scRNA-seq dataset for downstream Visium HD deconvolution.
# Input: Annotated scRNA-seq Seurat object and metadata used to define celltype labels.
# Output: SCreference object saved for later RCTD runs.
# Note: Review file paths before running. Keep original analysis logic unchanged.
# ==================================================
#RCTD
library(arrow)
library(spacexr)
library(ggplot2)
library(Seurat)
library(dplyr)
library(forcats)
library(tidyr)
library(purrr)
library(tibble)
options(future.globals.maxSize = 200 * 1024^3)
set.seed(123)

downsample_by_celltype_sample <- function(
    obj,
    celltype_col = "celltype",   # column name in meta.data representing cell type
    sample_col   = "sample",     # column name in meta.data representing sample ID
    cap_per_celltype = 10000,    # maximum number retained for each cell type
    seed = 123                  # random seed
){
  stopifnot(celltype_col %in% colnames(obj@meta.data))
  stopifnot(sample_col   %in% colnames(obj@meta.data))
  set.seed(seed)
  
  meta <- tibble::rownames_to_column(obj@meta.data, "cell") %>%
    dplyr::select(cell, dplyr::all_of(celltype_col), dplyr::all_of(sample_col)) %>%
    dplyr::rename(
      celltype = !!rlang::sym(celltype_col),
      sample   = !!rlang::sym(sample_col)
    )
  
  # total number of cells in each cell type
  ct_counts <- meta %>%
    dplyr::count(celltype, name = "n_ct")
  
  # stratified counts by cell type × sample
  strata_counts <- meta %>%
    dplyr::count(celltype, sample, name = "n_cs") %>%
    dplyr::left_join(ct_counts, by = "celltype") %>%
    dplyr::group_by(celltype) %>%
    dplyr::mutate(
      cap = pmin(n_ct, cap_per_celltype),
      quota = cap * (n_cs / n_ct),
      target_floor = floor(quota),
      remainder = cap - sum(target_floor),
      frac = quota - target_floor
    ) %>%
    dplyr::arrange(celltype, dplyr::desc(frac)) %>%
    dplyr::mutate(add = as.integer(dplyr::row_number() <= remainder)) %>%
    dplyr::arrange(celltype, sample) %>%
    dplyr::mutate(
      target = dplyr::if_else(n_ct <= cap, n_cs, target_floor + add)
    ) %>%
    dplyr::ungroup()
  
  # sample within each stratum
  pick_ids <- meta %>%
    dplyr::left_join(
      strata_counts %>% dplyr::select(celltype, sample, n_cs, target),
      by = c("celltype", "sample")
    ) %>%
    dplyr::group_by(celltype, sample) %>%
    dplyr::group_modify(~{
      n_take <- min(.x$target[1], nrow(.x))
      dplyr::slice_sample(.x, n = n_take)
    }) %>%
    dplyr::ungroup() %>%
    dplyr::pull(cell)
  
  obj_ds <- subset(obj, cells = pick_ids)
  
  return(obj_ds)
}

load("/home/data/chenwenhui/scRNAseq/QC/seurat_object_subcelltype_copykat_filtered_doublet_and_samples_progeny_step8.RData")
intergration.obj = intergration.obj.filtered.progeny
rm(intergration.obj.filtered.progeny)
meta = intergration.obj@meta.data
meta$celltype = as.character(meta$celltype)
meta$subcelltype = as.character(meta$subcelltype)
if (FALSE){
  meta$celltype[grepl('_B_',meta$subcelltype)] = 'B'
  meta$celltype[grepl('_Pla_',meta$subcelltype)] = 'Plasma'
  meta$celltype[grepl('Epithelial',meta$subcelltype)] = 'Epithelial'
  meta$celltype[grepl('Malignant',meta$subcelltype)] = 'Malignant'
  meta$celltype[grepl('_Macro_',meta$subcelltype)] = 'Macrophage'
  meta$celltype[grepl('_Mono_',meta$subcelltype)] = 'Monocyte'
  meta$celltype[grepl('_Mast_',meta$subcelltype)] = 'Mast'
  meta$celltype[grepl('_Netru_',meta$subcelltype)] = 'Neutrophil'
  meta$celltype[grepl('_cDC1_',meta$subcelltype)] = 'cDC1'
  meta$celltype[grepl('_cDC2_',meta$subcelltype)] = 'cDC2'
  meta$celltype[grepl('_cDC3_',meta$subcelltype)] = 'cDC3'
  meta$celltype[grepl('_pDC_',meta$subcelltype)] = 'pDC'
  meta$celltype[grepl('_Fib_',meta$subcelltype)] = 'Fibroblast'
  meta$celltype[grepl('_Endo_',meta$subcelltype)] = 'Endothelial'
  meta$celltype[grepl('_Per_',meta$subcelltype)] = 'Pericyte'
  meta$celltype[grepl('_SMC_',meta$subcelltype)] = 'SMC'
  meta$celltype[grepl('_CD8_',meta$subcelltype)] = 'CD8 T'
  meta$celltype[grepl('_CD4_',meta$subcelltype)] = 'CD4 T'
  meta$celltype[grepl('_CD4\\+CD8\\+_',meta$subcelltype)] = 'Tprof'
  meta[colnames(malignant.obj),'celltype'] = malignant.obj$MP_label
  meta$celltype[meta$celltype=='Malignant']='Epithelial'
}
meta[meta$celltype=='Epithelial','celltype']=meta[meta$celltype=='Epithelial','subcelltype']
table(meta$celltype)
intergration.obj@meta.data = meta
intergration.obj = downsample_by_celltype_sample(intergration.obj,cap_per_celltype=20000)

mat = intergration.obj@assays$RNA$counts
#mat = mat[,rownames(submeta)]
sc_umis <- Matrix::colSums(mat)
celltypes = as.factor(intergration.obj$celltype)
names(celltypes) = colnames(intergration.obj)
unique(celltypes)
SCreference <- Reference(mat, celltypes, sc_umis,n_max_cells=20000,min_UMI=600)
save(SCreference,file = '~/visiumHD/SCreference_celltype.RData')
