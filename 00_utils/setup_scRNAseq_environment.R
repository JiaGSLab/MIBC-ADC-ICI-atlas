# ==================================================
# Script: Setup scRNA-seq analysis environment
# Author: Wenhui Chen
# Purpose: Loads packages, defines shared colors and helper functions, and prepares base paths used by the scRNA-seq workflow.
# Input: Package libraries, sample_info file, raw count directories, and output paths referenced in the code.
# Output: Shared functions and initialized variables for downstream scRNA-seq scripts.
# ==================================================
if (TRUE) {
  suppressPackageStartupMessages({library(Seurat)
    library(patchwork)
    library(reticulate)
    library(ggrepel)
    library(SoupX)
    library(Matrix)
    library(dplyr)
    library(future)
    library(scCustomize)
    library(BiocParallel)
    library(SeuratWrappers)
    library(ggplot2)
    library(circlize)
    library(RColorBrewer)
    library(ComplexHeatmap)
    library(DescTools)
    library(ggsci)
    library(COSG)
    library(monocle3)
    #library(ROGUE)
    library(foreach)
    library(doParallel)
    library(furrr)
    library(clustree)
    library(future.apply)
    library(gplots)
    library(biomaRt)
    library(ggalluvial)
    library(msigdbr)
    library(ggsignif)
    library(ggthemr)
    library(ggunchull)
    library(corrplot)
    library(Hmisc)
    library(sscVis)
    library(data.table)
    library(circlize)
    library(grid)
    library(Startrac)
    library(readr)
    library(sscClust)
    library(pheatmap)
    library(scRNAtoolVis)
    library(DESeq2)
    library(miloR)
    library(ggpubr)
    library(enrichR)
    library(qs2)
    library(celda)
    library(sceasy)
    library(progeny)
    library(scMetabolism)
    library(rsvd)
    library(spacexr)
    library(SCNT)
    library(scStarPlot)
    library(scop)
    library(scproportion)
    library(CellChat)
    library(scPCselect)}
  )
  options(future.globals.maxSize = 600 * 1024^3)
  set.seed(123)
  colors_15 <- c("#A52A2A", "#053E7A", "#E46B4F", "#8B8B00", "#B8860B",  
                "#CD661D", "#8B5A2B", "#CDAD00", "#2F4F4F", "#483D8B", 
                "#7A378B", "#4A7088", "#556B2F", "#008B8B", "#CD7054")
  load('~/scRNAseq/colors_90.RData') #this is a color list and can be commented out if not needed
  colors_90 = colors_91
  #load('~/scRNAseq/QC/blacklist.RData')  #TCR 1:247 and BCR 248:683 correspond to TCR and immunoglobulin genes, with a few mitochondrial and red blood cell genes; removing them from HVGs can help clustering, but the effect is not large, so this line can be commented out
  markers=list(Epi=c('KRT7','KRT19','EPCAM','ERBB2','ERBB3','TP63','S100P'),
              Tcell=c('CD2','CD3D','CD3E','CD4','CD8A','NKG7'),
              B=c('CD79A','MS4A1','CD19','MZB1','JCHAIN'),
              Macro=c('LYZ','CD68','CSF1R','VCAN','FCGR3A','CD14'),
              pDC=c('LILRA4','GZMB','IL3RA'),
              cDC=c('CLEC9A','FLT3','IDO1','CD1C','FCER1A','LAMP3','CCR7','FSCN1'),
              Mast=c('MS4A2','CPA3','KIT'),
              Nutr=c('FCGR3B','S100A8','S100A9','CSF3R'),
              Fib=c('COL1A1','COL3A1','DCN','LUM'),
              Endo=c('VWF','PECAM1','CDH5'),
              SMC=c('ACTA2','TAGLN','MYH11'),
              Pericyte=c('RGS5','PDGFRB'),
              Prolif=c('MKI67','TOP2A'))
  non_immune.samples = c('LLZ-1','LS-2','LHF-1','ZZS-2') #when analyzing non-malignant cells, these samples should be filtered out first because the cell numbers are too low
  PR.samples = c('pt22','pt33','pt34','pt40','pt55','pt56','pt25','pt1','pt10') #patients confirmed to have partial response
  output_dir <- "~/scRNAseq/GPU_IT/QC"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  setwd(output_dir)
  source("/home/data/chenwenhui/scRNAseq/GPU_IT/scripts/function_for_plot.R")
}

#
sample_info = read.csv('/home/data/chenwenhui/sample_info.csv', header = TRUE, stringsAsFactors = FALSE,fill = TRUE)
sample_info = sample_info[sample_info$scFFPE.seq==1,]
sample_info$Time_point = ifelse(sample_info$Time_point == 'pre-treatment','Pre','Post')
sample_info$Outcome = ifelse(sample_info$Outcome == 'response','CR','NR')
sample_info$idx = NULL

calculate_simple_gene_set_score <- function(seurat_obj, gene_set, assay = "RNA", slot = "data", score_name = "sig_score") {
  #check whether gene_set is a list; input must be provided as a list
  if (!is.list(gene_set)) {
    stop("gene_set must be a list.")
  }

  #check whether score_name has length 1 or the same length as gene_set; otherwise throw an error
  if (length(score_name) != 1 && length(score_name) != length(gene_set)) {
    stop("score_name must be a single element or have the same length as gene_set.")
  }
  # 1. Check whether the gene set exists in the data
  # Retrieve all expressed genes in the current assay
  expressed_genes <- rownames(Seurat::GetAssayData(object = seurat_obj, assay = assay, slot = slot))

  if (length(score_name) ==1) {
      score_name = paste0(score_name,1:length(gene_set),'_AveExp')
  } else {
      score_name = paste0(score_name,'_AveExp')
  }
  
  for (i in 1:length(gene_set)) {
    # Find the overlapping genes
    valid_genes <- intersect(gene_set[[i]], expressed_genes)

    # Check whether any genes matched
    if (length(valid_genes) == 0) {
      print("No genes from the provided gene set were found in the Seurat object assay. Please check the gene symbols.")
      next
    }

    if (length(valid_genes) < length(gene_set[[i]])) {
      warning(
        paste0(
          "Note: only part of the provided gene set is present in the data (", length(gene_set[[i]]), " genes); only ", 
          length(valid_genes), "  genes are present and will be used for scoring", 
          length(valid_genes), "  genes will be used for calculation."
        )
      )
    }

    # 2. Extract the expression matrix
    expression_matrix <- Seurat::GetAssayData(object = seurat_obj, assay = assay, slot = slot)[valid_genes,]
    print(dim(expression_matrix))
    # 3. Calculate the average expression for each cell
    # matrixStats::rowMeans2 or base::rowMeans can be faster
    # Because the matrix is transposed, compute the mean of each column
    # colMeans calculates the mean of each column, that is, the average for each cell
    gene_set_scores <- colMeans(expression_matrix)

    # 4. Add the score to Seurat object metadata
    seurat_obj@meta.data[score_name[i]] <- gene_set_scores

    message(paste0("Successfully calculated and added score: ", score_name[i], "  to Seurat object metadata."))
  }
  return(seurat_obj)
}

QCcalulate <- function(obj, ratio=5) {
  meta <- obj@meta.data
  print(nrow(meta))
  
  # Clean the sample field to avoid invisible characters
  meta$sample <- as.character(meta$sample)
  meta$sample <- trimws(meta$sample)
  
  meta$is_outlier <- FALSE
  print(table(meta$is_outlier))

  samples <- unique(meta$sample)

  for (sam in samples) {
    
    idx <- which(meta$sample == sam)
    sub <- meta[idx, , drop = FALSE]

    nFeature_thresh <- median(sub$nFeature_RNA,na.rm = TRUE) + ratio * mad(sub$nFeature_RNA,na.rm = TRUE)
    nCount_thresh   <- median(sub$nCount_RNA,na.rm = TRUE)   + ratio * mad(sub$nCount_RNA,na.rm = TRUE)
    pct_mt_thresh   <- median(sub$pct_counts_mt,na.rm = TRUE)+ 3 * mad(sub$pct_counts_mt,na.rm = TRUE)

    outlier_sub <- (sub$nFeature_RNA > nFeature_thresh) |
                  (sub$nFeature_RNA <= 200) |
                  (sub$nCount_RNA > nCount_thresh) |
                  (sub$pct_counts_mt > pct_mt_thresh) |
                  (sub$log10GenesPerUMI < 0.8) |
                  (sub$log10GenesPerUMI == 1) |
                  (sub$pct_counts_hb > 1)
    meta$is_outlier[idx] <- outlier_sub
    print(length(outlier_sub))
  }
  print(table(meta$is_outlier))
  obj@meta.data <- meta
  return(obj)
}

decontX_sample <- function(sample_name, input_dir, output_dir) {
  batch = ifelse(str_ends(sample_name,'-CX'),'FFPE',
                ifelse(str_ends(sample_name,'_1114'),'20241114_FFPE',
                  ifelse(str_ends(sample_name,'_0619'),'20250619_FFPE','snRNAseq')))
  if (batch == '20250619_FFPE') {
    dirs <- file.path(input_dir, gsub('_0619','',sample_name))
  } else {
    dirs <- file.path(input_dir, sample_name)
  }
  if (dir.exists(paste0(dirs,'/filtered_feature_bc_matrix'))) {
    counts_filtered <- Read10X(paste0(dirs,'/filtered_feature_bc_matrix'))
    #counts_raw <- Read10X(paste0(dirs,'/raw_feature_bc_matrix'))
  } else if (dir.exists(paste0(dirs,'/sample_filtered_feature_bc_matrix'))) {
    counts_filtered <- Read10X(paste0(dirs,'/sample_filtered_feature_bc_matrix'))
    #counts_raw <- Read10X(paste0(dirs,'/sample_raw_feature_bc_matrix'))
  } else {
    counts_filtered <- Read10X(paste0(dirs,'/count/sample_filtered_feature_bc_matrix'))
    #counts_raw <- Read10X(paste0(dirs,'/count/sample_raw_feature_bc_matrix'))
  }

  #counts_raw = counts_raw[rownames(counts_filtered),]
  sample_name = gsub('-CX', '', sample_name)

  # seurat_object = CreateSeuratObject(counts = counts_filtered, 
  #                                   project = sample_name,
  #                                   min.cells = 3,
  #                                   min.features = 0)

  sce = SingleCellExperiment(list(counts = counts_filtered))
  #sce.raw = SingleCellExperiment(list(counts = counts_raw))
  # decontx_result_bg <- celda::decontX(x = sce,
  #                                 background = sce.raw,
  #                                 maxIter = 1000,
  #                                 varGenes = 2000)
  decontx_result <- celda::decontX(x = sce,
                                maxIter = 1000,
                                varGenes = 2000)
  
  pdf(file = paste0(output_dir, '/figs/', sample_name, '_decontX_diagnostics.pdf'), width = 8, height = 6)
  #p = celda::plotDecontXContamination(decontx_result_bg) + ggtitle(paste0(sample_name, ' decontX with background'))
  #print(p)
  p2 = celda::plotDecontXContamination(decontx_result) + ggtitle(paste0(sample_name, ' decontX without background'))
  print(p2)
  dev.off()
  seurat_object = CreateSeuratObject(counts = round(decontXcounts(decontx_result)), 
                                    project = sample_name,
                                    min.cells = 3,
                                    min.features = 0)
  #seurat_object$decontX_contamination_bg <- colData(decontx_result_bg)$decontX_contamination
  # Store decontX contamination values
  seurat_object$decontX_contamination <- colData(decontx_result)$decontX_contamination

  # Create the assay and add it to the Seurat object (assay name: "decontX")
  #seurat_object[["decontX_bg"]] <- Seurat::CreateAssayObject(counts = round(decontXcounts(decontx_result_bg)))
  #seurat_object[["decontX"]] <- Seurat::CreateAssayObject(counts = round(decontXcounts(decontx_result)))
  seurat_object$batch = batch
  seurat_object$sample = sample_name
  
  output_path <- paste0(output_dir, '/', sample_name, '.qs')
  qs_save(seurat_object,file = output_path)
}

#the decontX section was run on the GPU server, so the paths differ
if (FALSE) {
  input_dir1 = "/data/home/chenwenhui/scRNAseq/rawdata_all"
  samples1 = list.dirs(input_dir1, full.names = FALSE, recursive = FALSE)
  output_dir = "/data/home/chenwenhui/scRNAseq/decontX_results_all"
  figs_dir = "/data/home/chenwenhui/scRNAseq/decontX_results_all/figs"
  dir.create(figs_dir, showWarnings = FALSE)
  dir.create(output_dir, showWarnings = FALSE)

  input_dir2 <- '/data/home/chenwenhui/scRNAseq/FFPE_20250619/AM20250604-F1/outs/per_sample_outs'
  samples2 <- list.dirs(input_dir2, full.names = FALSE, recursive = FALSE)
  samples2 <- paste0(samples2, '_0619')

  # Run decontX_sample in parallel for samples1
  library(future)
  library(future.apply)

  nworkers <- 12
  future::plan(future::multisession, workers = nworkers)

  results_samples1 <- future.apply::future_lapply(samples1, function(s) {
    tryCatch({
      decontX_sample(s, input_dir1, output_dir)
      list(sample = s, status = "OK")
    }, error = function(e) {
      list(sample = s, status = "ERROR", message = conditionMessage(e))
    })
  }, future.seed = TRUE)

  # Restore sequential plan and release background workers
  future::plan(sequential)

  # Return the list of parallel results, with each element containing sample and status, and message when failed
  results_samples1
  # Run decontX_sample in parallel for samples2
  future::plan(future::multisession, workers = nworkers)

  results_samples2 <- future.apply::future_lapply(samples2, function(s) {
    tryCatch({
      decontX_sample(s, input_dir2, output_dir)
      list(sample = s, status = "OK")
    }, error = function(e) {
      list(sample = s, status = "ERROR", message = conditionMessage(e))
    })
  }, future.seed = TRUE)

  # Restore sequential plan and release background workers
  future::plan(sequential)

  # Return the list of parallel results
  results_samples2

  gc()
}

