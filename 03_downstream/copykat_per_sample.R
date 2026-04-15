# ==================================================
# Script: CopyKAT analysis by sample
# Author: Wenhui Chen
# Purpose: Runs CopyKAT on each sample separately and saves per-sample and merged copy number prediction results.
# Input: Main Seurat or summary objects already generated in the project, plus any file paths referenced in the script.
# Output: Intermediate result files and/or figure files written to output_dir.
# Note: Review file paths before running. Keep original analysis logic unchanged.
# ==================================================
library(Seurat)
library(patchwork)
library(copykat)
library(dplyr)
library(future)
library(future.apply)
library(BiocParallel)
library(foreach)
library(doParallel)
library(furrr)
library(qs2)
outputdir <- "~/scRNAseq/GPU_IT/QC/copykat_output"
options(future.globals.maxSize = 700 * 1024^3)
set.seed(123)
if (!dir.exists(outputdir)) {
  dir.create(outputdir)
}
setwd(outputdir)

copykat_sample <- function(obj,samp) {
  rawmat <- obj@assays$RNA$counts
  Idents(obj) = 'celltype'
  norm.cells = ""
  if (length(unique(obj$celltype)) > 1) {
    norm.cells <- WhichCells(obj,idents='Epithelial',invert=TRUE)
    if (length(norm.cells) < 50) {
      norm.cells = ""
    }
  }
  copykat.samp <- copykat(rawmat = rawmat, 
                          id.type="S",
                          cell.line = "no",
                          norm.cell.names = norm.cells,
                          win.size = 15,
                          KS.cut = 0.1,
                          LOW.DR = 0.02,
                          UP.DR=0.1, 
                          sam.name = samp,
                          distance = 'spearman',
                          n.cores = 1)
  gc()
  return(copykat.samp)
}

#load data
intergration.obj = qs_read('/home/data/chenwenhui/scRNAseq/GPU_IT/QC/seu_combined_IT_filtered_obj_intergration.qs')
samplelist <- unique(intergration.obj$sample)
intergration.obj = SplitObject(intergration.obj, split.by = "sample")
print(paste0("sample numbers: ",length(samplelist)))

#samplelist = setdiff(samplelist,finished_samples)
nworkers <- 15
future::plan(future::multisession, workers = nworkers)

# Parallel processing using future.apply
copykat_list <- future.apply::future_lapply(
  intergration.obj,
  function(obj) {
    samp = unique(obj$sample)
    print(samp)
    #obj <- subset(intergration.obj, subset = sample == samp)
    copykat_sample(obj, samp = samp)
    gc()
  },
  future.seed = TRUE
)
future::plan(sequential)
names(copykat_list) <- samplelist

# Save results
finished_samples_txt = list.files(outputdir,full.names = F,pattern = '*_copykat_prediction.txt')
finished_samples = gsub("_copykat_prediction.txt","",finished_samples_txt)
df.list = lapply(finished_samples_txt,function(x) {
  read.table(file.path(outputdir,x),header = T)
})
names(df.list) = finished_samples
df = do.call(rbind,df.list)
qs_save(df,file = file.path(outputdir,'copykat_prediction.qs'))
qs_save(copykat_list, file = file.path('/home/data/chenwenhui/scRNAseq/GPU_IT/QC', 'copykat_result_IT.qs'))
rm(list = ls())
gc()
print('All samples finished')
