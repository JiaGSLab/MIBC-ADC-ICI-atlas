# ==================================================
# Script: CellChat analysis across clinical groups
# Author: Wenhui Chen
# Purpose: Builds CellChat objects for clinical groups and compares intercellular communication patterns.
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
    library(foreach)
    library(doParallel)
    library(furrr)
    library(clustree)
    library(future.apply)
    library(gplots)
    library(ggalluvial)
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
    library(ggpubr)
    library(hdf5r)
    library(qs2)
    library(rsvd)
    library(scStarPlot)
    library(scop)
    library(scproportion)
    library(CellChat)}
  )
  options(future.globals.maxSize = 600 * 1024^3)
  set.seed(123)
  output_dir <- "~/scRNAseq/GPU_IT/QC/cellchat"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  setwd(output_dir)
}

cellchat.analysis=function(obj,output,smooth=TRUE,trim=NULL) {
  cellchat <- createCellChat(object = obj, group.by ="subcluster", assay = "RNA")
  CellChatDB <- CellChatDB.human
  CellChatDB.use <- subsetDB(CellChatDB)
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  #future::plan("multisession", workers = 16)
  future::plan(sequential)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  if (smooth) {
    cellchat <- smoothData(cellchat, adj = PPI.human)
    if (is.null(trim)) {
      cellchat <- computeCommunProb(cellchat, type = "triMean", trim = NULL,raw.use = FALSE) 
    } else {
      cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = trim,raw.use = FALSE)
    }
  } else {
    if (is.null(trim)) {
      cellchat <- computeCommunProb(cellchat, type = "triMean", trim = NULL,raw.use = TRUE)      
    } else {
      cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = trim,raw.use = TRUE)
    }
  }

  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  future::plan(sequential)
  qs_save(cellchat,output)
  gc()
  return(cellchat)
}

seu.combined.IT = qs_read('/home/data/chenwenhui/scRNAseq/GPU_IT/QC/seu_combined_IT_filtered_obj_intergration.qs')
malignant.obj = qs_read("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/malignant_obj_intergration.qs")
seu.combined.IT$subcluster = seu.combined.IT$celltype.new
seu.combined.IT$subcluster = as.character(seu.combined.IT$subcluster)
seu.combined.IT = subset(seu.combined.IT, subset = major == 'Epithelial',invert= T)
seu.combined.IT@meta.data[(malignant.obj@meta.data$Basal_pEMT>0.5) & (malignant.obj@meta.data$Uro_Differentiation<0.5),'subcluster'] = 'Basal_pEMT'
seu.combined.IT@meta.data[(malignant.obj@meta.data$Basal_pEMT<0.5) & (malignant.obj@meta.data$Uro_Differentiation>0.5),'subcluster'] = 'Uro_Differentiation'
seu.combined.IT$subcluster = ifelse(seu.combined.IT$subcluster=='Epithelial','other malignant',seu.combined.IT$subcluster)
#seu.combined.IT = subset(seu.combined.IT, subset = subcluster %in% c('m17_Nutri_ISG','t14_dnT_HAVCR2',"b04_B_MKI67","b05_Pla_SLC1A5","b06_Pla_CD38","b07_Pla_IGHG3","m07_Macro_MKI67","s17_Endo_MKI67","s18_Schwann_SOX10","t05_CD8_MKI67"),invert= T)
seu.combined.IT = subset(seu.combined.IT, subset = celltype.new %in% c('dnT cell','Schwann'),invert= T)
seu.combined.IT$subcluster = ifelse(seu.combined.IT$subcelltype=='t03_CD8_LAG3 (Tex)','t03_CD8_LAG3',seu.combined.IT$subcluster)
seu.combined.IT$subcluster = gsub(" \\(.*\\)", "", seu.combined.IT$subcluster)
seu.combined.IT$samples = as.character(seu.combined.IT$sample)
Pre_CR = subset(seu.combined.IT,subset = Group == 'Pre_CR')
Pre_NR = subset(seu.combined.IT,subset = Group == 'Pre_NR')
Post_NR = subset(seu.combined.IT,subset = Group == 'Post_NR')
Post_CR = subset(seu.combined.IT,subset = Group == 'Post_CR')

# cellchat.Pre_CR = cellchat.analysis(Pre_CR,'cellchat.Pre_CR.qs',smooth=TRUE,trim=0.1)
# gc()
# cellchat.Pre_NR = cellchat.analysis(Pre_NR,'cellchat.Pre_NR.qs',smooth=TRUE,trim=0.1)
# gc()
# cellchat.Post_NR = cellchat.analysis(Post_NR,'cellchat.Post_NR.qs',smooth=TRUE,trim=0.1)
# gc()
# cellchat.Post_CR = cellchat.analysis(Post_CR,'cellchat.Post_CR.qs',smooth=TRUE,trim=0.1)
# gc()
cellchat.Pre_CR = cellchat.analysis(Pre_CR,'cellchat.Pre_CR_NULL.qs',smooth=FALSE,trim=0.1)
gc()
cellchat.Pre_NR = cellchat.analysis(Pre_NR,'cellchat.Pre_NR_NULL.qs',smooth=FALSE,trim=0.1)
gc()
cellchat.Post_NR = cellchat.analysis(Post_NR,'cellchat.Post_NR_NULL.qs',smooth=FALSE,trim=0.1)
gc()
cellchat.Post_CR = cellchat.analysis(Post_CR,'cellchat.Post_CR_NULL.qs',smooth=FALSE,trim=0.1)
gc()

if (FALSE) {
  ##################################
  cellchat.Pre_CR = qs_read('cellchat.Pre_CR_NULL.qs')

  #cellchat.Pre_CR <- netAnalysis_computeCentrality(cellchat.Pre_CR, slot.name = "netP")
  cellchat.Pre_NR = qs_read('cellchat.Pre_NR_NULL.qs')
  cellchat.Post_NR = qs_read('cellchat.Post_NR_NULL.qs')
  cellchat.Post_CR = qs_read('cellchat.Post_CR_NULL.qs')
  #cellchat.Pre_NR <- netAnalysis_computeCentrality(cellchat.Pre_NR, slot.name = "netP")
  #cellchat.Post_NR <- netAnalysis_computeCentrality(cellchat.Post_NR, slot.name = "netP")
  #cellchat.Post_CR <- netAnalysis_computeCentrality(cellchat.Post_CR, slot.name = "netP")


  #####################
  #Pre_NR vs Pre_CR
  object.list = list(Pre_NR = cellchat.Pre_NR, Pre_CR = cellchat.Pre_CR)
  cellchat <- mergeCellChat(object.list, add.names = names(object.list))
  gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
  gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
  gg1 + gg2
  pdf('interaction_and_strength_celltype.pdf',width=24,height=12)
  gg1 <- netVisual_heatmap(cellchat)
  #> Do heatmap based on a merged object
  gg2 <- netVisual_heatmap(cellchat, measure = "weight")
  #> Do heatmap based on a merged object
  gg1 + gg2
  dev.off()


  num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
  weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
  gg <- list()
  for (i in 1:length(object.list)) {
    gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax,label.size=3)
  }
  #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  patchwork::wrap_plots(plots = gg)
  ggsave(filename = "income_outcome_scatter_celltype.pdf",width = 12,height = 6)


  gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = 21, stacked = T, do.stat = TRUE)
  gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = 21, stacked = F, do.stat = TRUE)

  gg1 + gg2
  ggsave('rankNET_celltype.pdf',width = 8,height = 14)

  pdf('signaling_patterns_celltype.pdf',width = 12,height = 18)
  i = 1
  pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
  pathway.union = c('ADGRG','PDL2','PD-L1','CD160','CCL','CXCL','CD86','CD40','CD80')
  ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
  ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
  draw(ht1 + ht2)
  dev.off()

  pdf('signaling_patterns_incoming_celltype.pdf',width = 10,height = 7)
  i = 1
  pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
  pathway.union = c('ADGRG','PDL2','PD-L1','CD160','CCL','CXCL','CD86','CD40','CD80')
  ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
  ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
  draw(ht1 + ht2)
  dev.off()

  pdf('netVisual_bubble_increase_celltype.pdf',width = 20,height = 30)
  gg1 <- netVisual_bubble(cellchat, sources.use =NULL, targets.use = 21,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in CR", angle.x = 45, remove.isolate = T)
  gg1
  dev.off()

  pdf('netVisual_bubble_decrease_celltype.pdf',width = 20,height = 25)
  gg1 <- netVisual_bubble(cellchat, sources.use =NULL, targets.use = 21,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in CR", angle.x = 45, remove.isolate = T)
  gg1
  dev.off()


  curated_inhibitory_pathways <- c("PD-L1", "PDL2", "VISTA", "BTLA", "TGFb", 
                                  "SPP1", "CD39", "VEGF", "CCL", "CXCL",'ADGRG','CD160','CD86','CD40')

  # Use pairLR.use or signaling parameters during plotting
  netVisual_bubble(cellchat, 
                  signaling = curated_inhibitory_pathways,
                  sources.use = NULL, 
                  targets.use = 21,
                  comparison = c(1, 2), 
                  angle.x = 45,max.dataset = 1)
  ggsave('curated_inhibitory_pathways_NR_celltype.pdf',width = 12,height = 8)
  netVisual_bubble(cellchat, 
                  signaling = curated_inhibitory_pathways,
                  sources.use = NULL, 
                  targets.use = 21,
                  comparison = c(1, 2), 
                  angle.x = 45,max.dataset = 2)
  ggsave('curated_inhibitory_pathways_CR_celltype.pdf',width = 12,height = 8)

  pos.dataset = "Pre_NR"
  features.name = paste0(pos.dataset, ".merged")
  cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 
  net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
  net.up <- subsetCommunication(cellchat, net = net, datasets = "Pre_NR",ligand.logFC = 0.05, receptor.logFC = NULL)
  net.down <- subsetCommunication(cellchat, net = net, datasets = "Pre_CR",ligand.logFC = -0.05, receptor.logFC = NULL)
  gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
  gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
  netVisual_chord_gene(object.list[[1]],
                      sources.use = NULL, 
                      targets.use = 4, 
                      slot.name = 'net', 
                      net = net.up, 
                      lab.cex = 0.8, 
                      small.gap = 3.5, 
                      title.name = paste0("Up-regulated signaling in ", names(object.list)[1]),
                      color.use = colors_90)
  netVisual_chord_gene(object.listNULL,
                      sources.use = c(7:24), 
                      targets.use = 46, 
                      slot.name = 'net', 
                      net = net.down, 
                      lab.cex = 0.8, 
                      small.gap = 3.5, 
                      title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
  computeEnrichmentScore(net.up, species = 'human', variable.both = TRUE)

  pdf('netVisual_bubble_pre_nr_basal_pEMT.pdf',width = 40,height = 60)
  netVisual_bubble(cellchat, 
                  signaling = NULL,
                  sources.use = c(7:57), 
                  targets.use = 7,
                  comparison = c(1, 2), 
                  angle.x = 45,max.dataset = 1)
  dev.off()
  pdf('netVisual_bubble_pre_cr_pEMT.pdf',width = 40,height = 60)
  netVisual_bubble(cellchat, 
                  signaling = NULL,
                  sources.use = c(7:57), 
                  targets.use = 7,
                  comparison = c(1, 2), 
                  angle.x = 45,max.dataset = 2)
  dev.off()
}
