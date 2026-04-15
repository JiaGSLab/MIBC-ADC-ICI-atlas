# ==================================================
# Script: DEG analysis in MP4-high malignant cells
# Author: Wenhui Chen
# Purpose: Performs MAST differential expression analysis between CR and NR within pre-treatment MP4-high malignant cells.
# Input: Main Seurat or summary objects already generated in the project, plus any file paths referenced in the script.
# Output: Intermediate result files and/or figure files written to output_dir.
# Note: Review file paths before running. Keep original analysis logic unchanged.
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
}

malignant.obj = qs_read(file.path(output_dir, 'malignant_obj_intergration.qs'))

mp4_pre <- subset(malignant.obj, subset = (MP4_Status_99 == "MP4_High") &
                                          (Time_point == "Pre"))
Idents(mp4_pre) = "Outcome"
future::plan(sequential)
DEG.MP4 = FindMarkers(mp4_pre,
                      ident.1= "NR",
                      ident.2 = "CR",
                      test.use="MAST",
                      logfc.threshold=0.05,
                      min.pct = 0.02,
                      latent.vars=c("sample"),
                      random.seed=123)

qs_save(DEG.MP4,file.path(output_dir,"DEG.Pre_MP4_MAST.qs"))