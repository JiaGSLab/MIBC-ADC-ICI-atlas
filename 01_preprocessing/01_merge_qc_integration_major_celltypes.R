# ==================================================
# Script: Merge, QC, integration, and major cell type assignment
# Author: Wenhui Chen
# Purpose: Merges per-sample objects, performs QC and Harmony integration, assigns major cell types, and generates broad cell class objects.
# Input: decontX per-sample qs files, sample metadata, and helper functions from setup_scRNAseq_environment.R.
# Output: Integrated combined object, major cell type objects, and QC or overview plots.
# Note: Run scripts in workflow order. Keep original path definitions unless updating the project layout.
# ==================================================
# Set parallel execution plan
future::plan(future::multisession, workers = 51)
decontX.output.dir = "/home/data/chenwenhui/scRNAseq/GPU_IT/decontX_results_all"
samples_all <- list.files(decontX.output.dir, full.names = FALSE, pattern = "*.qs")
seu.list <- future.apply::future_lapply(samples_all, function(x) {
  seu = qs_read(file.path(decontX.output.dir, x))
  seu = subset(seu, subset=(nFeature_RNA>200)&(nCount_RNA>200))
  colnames(seu) <- paste0(seu$sample, '_', colnames(seu))
  seu$pct_counts_mt <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu$pct_counts_hb <- PercentageFeatureSet(seu, pattern = "^HB[^(PES)]")
  seu$pct_counts_rp <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
  seu$log10GenesPerUMI <- log10(seu$nFeature_RNA) / log10(seu$nCount_RNA)
  #seu[['decontX_bg']] <- NULL  # remove the assay containing background decontamination results
  return(seu)
}, future.seed = TRUE)

# Restore sequential plan
future::plan(sequential)
gc()
names(seu.list) <- gsub('.qs','',samples_all)
seu.combined <- merge(seu.list[[1]], y = seu.list[-1])
seu.combined <- QCcalulate(seu.combined)
seu.combined <- subset(seu.combined,subset = is_outlier==FALSE)

seu.combined$orig.ident = seu.combined$sample
# Select columns from sample_info that need to be added, excluding columns already present in meta.data
cols_to_add <- setdiff(colnames(sample_info), colnames(seu.combined@meta.data))

# Construct a metadata table aligned to cell order, with cell names as row names
meta_to_add <- sample_info |> 
  as.data.frame() |> 
  dplyr::select(dplyr::any_of(cols_to_add))

meta_matched <- meta_to_add[match(seu.combined$orig.ident, sample_info$orig.ident), , drop = FALSE]
rownames(meta_matched) <- colnames(seu.combined)

# Add the matched metadata to the Seurat object
seu.combined <- Seurat::AddMetaData(seu.combined, metadata = meta_matched)

# Return the updated object
seu.combined
qs_save(seu.combined, file = file.path(output_dir, 'seu_combined_decontX.qs'))
seu.combined.IT <- subset(seu.combined, subset = Therapy == 'IT')
seu.combined.IC <- subset(seu.combined, subset = Therapy == 'IC')
seu.combined.IT$Group = paste0(seu.combined.IT$Time_point,'_',seu.combined.IT$Outcome)
group_sample_counts <- seu.combined.IT@meta.data |>
  as.data.frame() |>
  tibble::as_tibble(rownames = "cell") |>
  dplyr::select(Group, sample) |>
  dplyr::distinct() |>
  dplyr::count(Group, name = "n_samples") |>
  dplyr::arrange(desc(n_samples))

group_sample_counts
qs_save(seu.combined.IT, file = file.path(output_dir, 'seu_combined_decontX_IT.qs'))
qs_save(seu.combined.IC, file = file.path(output_dir, 'seu_combined_decontX_IC.qs'))
rm(seu.combined.IT,seu.combined.IC,seu.combined)

seu.combined.IT <- qs_read(file.path(output_dir, 'seu_combined_decontX_IT.qs'))
DefaultAssay(seu.combined.IT) <- 'RNA'
seu.combined.IT <-  subset(seu.combined.IT, subset = batch=='snRNAseq',invert=TRUE)
seu.combined.IT = JoinLayers(seu.combined.IT)
seu.combined.IT[['RNA']] = split(seu.combined.IT[['RNA']],f=seu.combined.IT$sample)
# seu.combined.IT <- subset(seu.combined.IT, subset = (pct_counts_mt < 20) &
#                                                     (pct_counts_hb < 1) &
#                                                     (nFeature_RNA > 200) &
#                                                     (log10GenesPerUMI > 0.8) &
#                                                     (log10GenesPerUMI < 1) &
#                                                     (nFeature_RNA < 7000))
seu.combined.IT = NormalizeData(seu.combined.IT)
s.genes = intersect(cc.genes$s.genes, rownames(seu.combined.IT))
g2m.genes = intersect(cc.genes$g2m.genes, rownames(seu.combined.IT))
cc.genes = list(s.genes = s.genes, g2m.genes = g2m.genes)
seu.combined.IT = CellCycleScoring(seu.combined.IT,
                                s.features = cc.genes$s.genes,
                                g2m.features = cc.genes$g2m.genes)

qs_save(seu.combined.IT,file=file.path(output_dir, 'seu_combined_decontX_IT_remove_snRNAseq.qs'))
#seu.combined.IT = qs_read(file.path(output_dir, 'seu_combined_decontX_IT_remove_snRNAseq.qs'))
seu.combined.IT = seu.combined.IT %>%
  FindVariableFeatures(.,selection.method = 'vst',nfeatures = 2500)
hvg = VariableFeatures(seu.combined.IT)
MT.genes = intersect(rownames(seu.combined.IT)[grepl('^MT-',rownames(seu.combined.IT))], rownames(seu.combined.IT))
HB.genes = intersect(rownames(seu.combined.IT)[grepl('^HB[^(PES)]',rownames(seu.combined.IT))], rownames(seu.combined.IT))
RP.genes = intersect(rownames(seu.combined.IT)[grepl('^RP[SL]',rownames(seu.combined.IT))], rownames(seu.combined.IT))
blacklist = c(MT.genes,HB.genes,RP.genes)
hvg = setdiff(hvg,blacklist)
VariableFeatures(seu.combined.IT) = hvg
seu.combined.IT = seu.combined.IT %>%
  ScaleData(.,vars.to.regress = c("S.Score", "G2M.Score",'pct_counts_mt')) %>%
  RunPCA(.)
ElbowPlot(seu.combined.IT,ndims = 50)
npcs = 50
seu.combined.IT = seu.combined.IT %>%
  IntegrateLayers(.,method = HarmonyIntegration,
                    new.reduction = "harmony",orig.reduction = "pca",max.iter.harmony = 20L,verbose=T) %>%
  RunUMAP(.,reduction = "harmony",dims = 1:npcs) %>%
  FindNeighbors(., reduction = "harmony", dims = 1:npcs,verbose = F) %>%
  FindClusters(., resolution = c(1),algorithm=4,verbose = F,random.seed=1234)
seu.combined.IT = JoinLayers(seu.combined.IT)

# Create DotPlot and DimPlot and save them into the same PNG
png(file.path(output_dir, "DotPlot_DimPlot_combined.png"), width = 3600, height = 1800, res = 150)

# Create DotPlot
dotplot <- DotPlot(seu.combined.IT, features = markers,
        cols = 'RdBu',dot.scale = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())

# Create DimPlot
dimplot <- DimPlot(seu.combined.IT, cols = colors_91, raster = F, alpha = 0.2, label = T, label.size = 5) + NoLegend()

# Combine the two plots with patchwork
combined_plot <- dotplot + dimplot +
  plot_layout(ncol = 2, widths = c(2, 1)) +
  plot_annotation(tag_levels = 'A')

print(combined_plot)
dev.off()

# Also display in the console
#print(combined_plot)
# Create a combined QC figure: FeaturePlots for nFeature_RNA, nCount_RNA, and pct_counts_mt, together with a DimPlot of is_outlier
png(file.path(output_dir, "QC_combined.png"), width = 3600, height = 3600, res = 150)

# Create three FeaturePlots
feature_plot1 <- FeaturePlot(seu.combined.IT, features = 'nFeature_RNA', cols = c('#87CEEB','#FFB6C1'), order = T, alpha = 0.2,raster=F) + 
  ggtitle('nFeature_RNA')

feature_plot2 <- FeaturePlot(seu.combined.IT, features = 'nCount_RNA', cols = c('#87CEEB','#FFB6C1'), order = T, alpha = 0.2,raster=F) + 
  ggtitle('nCount_RNA')

feature_plot3 <- FeaturePlot(seu.combined.IT, features = 'pct_counts_mt', cols = c('#87CEEB','#FFB6C1'), order = T, alpha = 0.2,raster=F) + 
  ggtitle('pct_counts_mt')

# Create a DimPlot of is_outlier using light red and light blue
#outlier_colors <- c('FALSE' = '#87CEEB', 'TRUE' = '#FFB6C1')  # light blue and light red

# Ensure is_outlier is stored as a factor
#seu.combined.IT$is_outlier <- factor(seu.combined.IT$is_outlier, levels = c('FALSE', 'TRUE'))

dimplot_outlier <- DimPlot(seu.combined.IT, group.by = 'is_outlier', cols = c('#87CEEB','#FFB6C1'), 
                          raster = F, alpha = 0.2, label = F,split.by = 'is_outlier') + 
  ggtitle('Outlier Cells') 

# Combine the four plots into a 2x2 layout with patchwork
combined_qc_plot <-
  (feature_plot1 | feature_plot2) /
  (feature_plot3 | dimplot_outlier) +
  plot_annotation(tag_levels = 'A')



print(combined_qc_plot)
dev.off()

qs_save(seu.combined.IT,file=file.path(output_dir, 'seu_combined_decontX_IT_remove_snRNAseq_integration.qs'))

markers.cluster = cosg(seu.combined.IT,
                      groups = 'all',mu = 100,expressed_pct = 0.3)
View(markers.cluster$names)

seu.decontX.IT = seu.combined.IT
doublet = c(25,26,31,35,37,40,44,46,47,49,50,51)
Myeloid.cluster=c(10,22,36,39)
Mast.cluster=c(30)
Neutrophils.cluster=c(27)
T.cluster=c(7,9)
B.cluster=c(15)
Plasma.cluster=c(16)
SMC.cluster=c(23)
Pericyte.cluster=c(20)
Fibroblast.cluster=c(3)
Endothelials.cluster=c(12)
Epithelials.cluster =c(1,2,4,5,6,8,11,13,14,17,18,19,21,24,28,29,32,33,34,38,41,42,43,45,48,52)
seu.decontX.IT = subset(seu.decontX.IT,idents=doublet,invert=TRUE)
Idents(seu.decontX.IT, cells = WhichCells(seu.decontX.IT, idents = Epithelials.cluster)) <- "Epithelial"
Idents(seu.decontX.IT, cells = WhichCells(seu.decontX.IT, idents = Myeloid.cluster)) <- "Myeloid"
Idents(seu.decontX.IT, cells = WhichCells(seu.decontX.IT, idents = Mast.cluster)) <- "Mast"
Idents(seu.decontX.IT, cells = WhichCells(seu.decontX.IT, idents = Neutrophils.cluster)) <- "Neutrophil"
Idents(seu.decontX.IT, cells = WhichCells(seu.decontX.IT, idents = T.cluster)) <- "T"
Idents(seu.decontX.IT, cells = WhichCells(seu.decontX.IT, idents = B.cluster)) <- "B"
Idents(seu.decontX.IT, cells = WhichCells(seu.decontX.IT, idents = Plasma.cluster)) <- "Plasma"
Idents(seu.decontX.IT, cells = WhichCells(seu.decontX.IT, idents = SMC.cluster)) <- "SMC"
Idents(seu.decontX.IT, cells = WhichCells(seu.decontX.IT, idents = Pericyte.cluster)) <- "Pericyte"
Idents(seu.decontX.IT, cells = WhichCells(seu.decontX.IT, idents = Fibroblast.cluster)) <- "Fibroblast"
Idents(seu.decontX.IT, cells = WhichCells(seu.decontX.IT, idents = Endothelials.cluster)) <- "Endothelial"
seu.decontX.IT$celltype = Idents(seu.decontX.IT)
qs_save(seu.decontX.IT,file=file.path(output_dir, 'seu_combined_decontX_IT_intergration_celltype.qs'))
png(file.path(output_dir, "celltype_dimplot.png"), width = 1800, height = 1800, res = 300)
DimPlot(seu.decontX.IT,cols = colors_91,raster =F,alpha=0.2,label =T,label.size = 5) + NoLegend()
dev.off()

Epi.obj = subset(seu.decontX.IT,idents = c('Epithelial'))
Myeloid.obj = subset(seu.decontX.IT,idents = c('Myeloid','Mast','Neutrophil'))
Tcell.obj = subset(seu.decontX.IT,idents = c('T'))
Bcell.obj = subset(seu.decontX.IT,idents = c('B','Plasma'))
Stromal.obj = subset(seu.decontX.IT,idents = c('SMC','Pericyte','Fibroblast','Endothelial'))
#rm(seu.decontX.IT)
gc()
qs_save(Epi.obj,file=file.path(output_dir, 'Epi_obj.qs'))
qs_save(Myeloid.obj,file=file.path(output_dir, 'Myeloid_obj.qs'))
qs_save(Tcell.obj,file=file.path(output_dir, 'Tcell_obj.qs'))
qs_save(Bcell.obj,file=file.path(output_dir, 'Bcell_obj.qs'))
qs_save(Stromal.obj,file=file.path(output_dir, 'Stromal_obj.qs'))

rm(seu.decontX.IT,seu.combined.IT)
gc()

subcelltype_analysis <- function(obj,output_dir,prefix,nfeatures = 2500,npcs = 30,res = seq(0.1,1.5,0.1),loess.span=0.3) {
  meta.raw = obj@meta.data
  meta.raw = meta.raw[,!grepl('^RNA_snn_res',colnames(meta.raw))]
  obj = CreateSeuratObject(counts = obj@assays$RNA$counts,project = prefix,min.cells = 10,min.features = 200)
  meta = obj@meta.data
  add.cols = setdiff(colnames(meta.raw),colnames(meta))
  meta[add.cols] = meta.raw[row.names(meta),add.cols]
  obj@meta.data = meta
  obj$orig.ident = obj$sample
  #obj@meta.data = obj@meta.data[,!grepl('^RNA_snn_res',colnames(obj@meta.data))]
  obj[['RNA']] = split(obj[['RNA']],f=obj$sample)
  obj = obj %>%
    NormalizeData(.,verbose = T) %>%
    FindVariableFeatures(.,selection.method = 'vst',nfeatures = nfeatures,loess.span = loess.span)
  hvg = VariableFeatures(obj)
  MT.genes = rownames(obj)[grepl('^MT-',rownames(obj))]
  HB.genes = rownames(obj)[grepl('^HB[^(PES)]',rownames(obj))]
  RP.genes = rownames(obj)[grepl('^RP[SL]',rownames(obj))]
  blacklist = c(MT.genes,HB.genes,RP.genes)
  hvg = setdiff(hvg,blacklist)
  VariableFeatures(obj) = hvg
  obj = obj %>%
    ScaleData(.,vars.to.regress = c("S.Score", "G2M.Score",'pct_counts_mt')) %>%
    RunPCA(.,verbose = T) %>%
    IntegrateLayers(.,method = HarmonyIntegration,
                    new.reduction = "harmony",orig.reduction = "pca",max.iter.harmony = 20L,verbose=T) %>%
    RunUMAP(.,reduction = "harmony",dims = 1:npcs) %>%
    FindNeighbors(., reduction = "harmony", dims = 1:npcs,verbose = F) %>%
    FindClusters(., resolution = res,algorithm=4,verbose = F,random.seed=1234)
  obj = JoinLayers(obj)
  qs_save(obj,file=file.path(output_dir, paste0(prefix,'_obj_intergration.qs')))
  return(obj)
}


# Define the object list and corresponding prefixes
obj.list = list(Epi.obj, Myeloid.obj, Tcell.obj, Bcell.obj, Stromal.obj)
prefix.list = c('Epi', 'Myeloid', 'Tcell', 'Bcell', 'Stromal')

for (i in 1:length(obj.list)) {
  obj.list[[i]] = subcelltype_analysis(obj.list[[i]],output_dir,prefix.list[[i]])
}


