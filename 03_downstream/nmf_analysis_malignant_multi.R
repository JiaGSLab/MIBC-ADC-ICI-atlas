# ==================================================
# Script: NMF analysis for malignant cells
# Author: Wenhui Chen
# Purpose: Runs multi-rank NMF on malignant cells for each sample to identify recurrent transcriptional programs.
# Input: Main Seurat or summary objects already generated in the project, plus any file paths referenced in the script.
# Output: Intermediate result files and/or figure files written to output_dir.
# Note: Review file paths before running. Keep original analysis logic unchanged.
# ==================================================
## --- 1. Environment setup (minimal loading) ---
suppressPackageStartupMessages({
  library(Seurat)
  library(NMF)
  library(future)
  library(future.apply)
  library(qs2)
})

# Set the global variable size limit (for a 1 TB machine, 500-800 GB is recommended)
options(future.globals.maxSize = 800 * 1024^3)
set.seed(123)

output_dir <- "~/scRNAseq/GPU_IT/QC/NMF_output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

## --- 2. Data loading and preprocessing ---
# Load the large object only once
malignant.obj <- qs_read("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/malignant_obj_intergration.qs")

# Filter out MT and HB genes
genes_to_remove <- rownames(malignant.obj)[grepl('^MT-|^HB[^(PES)]', rownames(malignant.obj))]
malignant.obj <- malignant.obj[!rownames(malignant.obj) %in% genes_to_remove, ]

# Split by sample and filter out samples with fewer than 200 cells
obj.list <- Seurat::SplitObject(malignant.obj, split.by = 'sample')
obj.list <- obj.list[sapply(obj.list, function(x) ncol(x) > 200)]

# Key step: clean memory promptly and keep only the list object
rm(malignant.obj); gc()

## --- 3. Core computation function ---
nmf_worker <- function(obj, ranks = 4:9, method = "snmf/r", seed = 123) {
  sample_id <- unique(obj$sample)
  
  # Retrieve the expression matrix (Seurat v5 recommends layer="data")
  mat_raw <- Seurat::GetAssayData(obj, assay = "RNA", layer = "data")
  
  # Gene filtering: keep genes with expression >2 in more than 2% of cells
  keep <- apply(mat_raw, 1, function(x) sum(x > 2) > (ncol(mat_raw) * 0.02))
  if (length(keep)<2000) {
    keep <- apply(mat_raw, 1, function(x) sum(x > 1.7) > (ncol(mat_raw) * 0.02))
  }
  mat <- mat_raw[keep, ]

  
  # Center the matrix
  mat <- mat - rowMeans(mat)
  mat[mat < 0] <- 0
  mat <- as.matrix(mat[rowSums(mat) > 0, ]) # convert to a standard matrix and remove all-zero rows
  
  # Run NMF for each rank
  sample_res <- lapply(ranks, function(r) {
    # Use tryCatch so that failure at a single rank does not stop the whole job
    fit <- tryCatch({
      NMF::nmf(mat, rank = r, method = method, seed = seed, nrun = 1)
    }, error = function(e) return(NULL))
    
    if (is.null(fit)) return(NULL)
    
    return(list(
      rank = r,
      w = NMF::basis(fit),
      h = t(NMF::coef(fit))
    ))
  })
  
  names(sample_res) <- paste0("rank_", ranks)
  return(sample_res)
}

## --- 4. Run in parallel ---
nworkers <- 10  # with 1 TB memory, using 10–15 workers is generally more stable
future::plan(future::multisession, workers = nworkers)

message("Starting NMF for all samples...")
nmf_results <- future.apply::future_lapply(obj.list, nmf_worker, future.seed = TRUE)
message("Calculation finished.")

## --- 5. Save results ---
qs_save(nmf_results, file.path(output_dir, "nmf_results_all_ranks_multi.qs"))

# Clean the environment
rm(list = ls()); gc()