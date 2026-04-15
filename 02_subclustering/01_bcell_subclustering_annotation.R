# ==================================================
# Script: B cell subclustering and annotation
# Author: Wenhui Chen
# Purpose: Performs B cell subclustering, manual annotation, and related visualization.
# Input: Bcell_obj_intergration.qs and shared variables created in the preprocessing workflow.
# Output: Annotated B cell object and B cell subtype plots.
# Note: Run scripts in workflow order. Keep original path definitions unless updating the project layout.
# ==================================================

Bcell.obj = qs_read(file.path(output_dir, 'Bcell_obj_intergration_new.qs'))
obj = Bcell.obj
prefix = 'Bcell'
nfeatures = 1000
npcs = 10
res = seq(0.1,1.5,0.1)
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
    FindVariableFeatures(.,selection.method = 'vst',nfeatures = nfeatures)
hvg = VariableFeatures(obj)
MT.genes = rownames(obj)[grepl('^MT-',rownames(obj))]
HB.genes = rownames(obj)[grepl('^HB[^(PES)]',rownames(obj))]
RP.genes = rownames(obj)[grepl('^RP[SL]',rownames(obj))]
blacklist = c(MT.genes,HB.genes,RP.genes)
hvg = setdiff(hvg,blacklist)
VariableFeatures(obj) = hvg
obj = obj %>%
    ScaleData(.,vars.to.regress = c("S.Score", "G2M.Score",'pct_counts_mt')) %>%
    RunPCA(.,verbose = T)
viz_results <- visualize_pc_selection(obj, max_pcs = 50)

obj = obj %>%
    IntegrateLayers(.,method = HarmonyIntegration,
                    new.reduction = "harmony",orig.reduction = "pca",max.iter.harmony = 20L,verbose=T) %>%
    RunUMAP(.,reduction = "harmony",dims = 1:npcs) %>%
    FindNeighbors(., reduction = "harmony", dims = 1:npcs,verbose = F) %>%
    FindClusters(., resolution = res,algorithm=4,verbose = F,random.seed=1234)
DimPlot(obj,reduction = 'umap',group.by = 'RNA_snn_res.0.8',cols = colors_90,raster = F,label=T,label.size=5) + NoLegend() + ggtitle('Bcell')
obj = JoinLayers(obj)
Idents(obj) = "RNA_snn_res.0.8"
DotPlot(obj,
  features = c('MS4A1','CD19','MZB1','JCHAIN',"CD79A","CD79B","MKI67"),
  cols = 'RdBu') +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
B.obj = subset(obj,idents = c(2,5,7,8,9))
Pla.obj = subset(obj,idents = c(1,3,4,6,11))
markers.B = cosg(B.obj,
            groups='all',
            assay = 'RNA',
            mu=100,
            remove_lowly_expressed = TRUE,
            expressed_pct = 0.3,
            n_genes_user=50)
View(markers.B$names)
markers.Pla = cosg(Pla.obj,
            groups='all',
            assay = 'RNA',
            mu=100,
            remove_lowly_expressed = TRUE,
            expressed_pct = 0.3,
            n_genes_user=50)
View(markers.Pla$names)
DotPlot(Pla.obj,
  features = markers.Pla$names[1:20,'6'],
  cols = 'RdBu') +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
Bcell.obj = obj
Bcell.obj = subset(Bcell.obj,idents = c(10),invert=T)
Idents(Bcell.obj, cells = WhichCells(Bcell.obj, idents = c(2))) <- "b01_B_CD82"
Idents(Bcell.obj, cells = WhichCells(Bcell.obj, idents = c(5))) <- "b02_B_CIRBP"
Idents(Bcell.obj, cells = WhichCells(Bcell.obj, idents = c(7))) <- "b03_B_IGHD"
Idents(Bcell.obj, cells = WhichCells(Bcell.obj, idents = c(8))) <- "b04_B_ISG"
Idents(Bcell.obj, cells = WhichCells(Bcell.obj, idents = c(9))) <- "b05_B_MKI67"
Idents(Bcell.obj, cells = WhichCells(Bcell.obj, idents = c(1))) <- "b06_Pla_CCNL2"
Idents(Bcell.obj, cells = WhichCells(Bcell.obj, idents = c(3,4))) <- "b07_Pla_B2M"
Idents(Bcell.obj, cells = WhichCells(Bcell.obj, idents = c(6))) <- "b08_Pla_FTH1"
Idents(Bcell.obj, cells = WhichCells(Bcell.obj, idents = c(11))) <- "b09_Pla_MKI67"
Bcell.obj$subcelltype = Idents(Bcell.obj)
Idents(Bcell.obj) = 'subcelltype'
qs_save(Bcell.obj,file.path(output_dir, 'Bcell_obj_intergration_new.qs'))
DimPlot(Bcell.obj,reduction = 'umap',cols = colors_90,raster = F,alpha=0.7) + ggtitle('Bcell') + theme(aspect.ratio = 1)
ggsave(file.path(output_dir, "Figures",'Bcell_subcelltype_dimplot_new.pdf'),width = 8,height =6)
plot_density_dimplot_new(seurat_object = Bcell.obj,
                    group.by = 'subcelltype',
                    split.by = 'Group',
                    ncol = 2,
                    title = 'Bcell',
                    color_palette = 'magma',group_colors = colors_90)
ggsave(file.path(output_dir,"Figures", 'Bcell_subcelltype_density_new.png'),width = 14,height =12,dpi = 300)




