# ==================================================
# Script: Myeloid subclustering and annotation
# Author: Wenhui Chen
# Purpose: Performs myeloid subclustering, manual annotation, and related visualization.
# Input: Myeloid_obj_intergration.qs and shared variables created in the preprocessing workflow.
# Output: Annotated myeloid object, updated epithelial object, and subtype plots.
# Note: Run scripts in workflow order. Keep original path definitions unless updating the project layout.
# ==================================================
#################################
#Myeloid
Myeloid.obj = qs_read(file.path(output_dir, 'Myeloid_obj_intergration.qs'))
obj = Myeloid.obj
prefix = 'Myeloid'
nfeatures = 1500 
npcs = 20
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
recommended_pcs <- calculate_optimal_pcs(obj)

obj = obj %>%
    IntegrateLayers(.,method = HarmonyIntegration,
                    new.reduction = "harmony",orig.reduction = "pca",max.iter.harmony = 20L,verbose=T) %>%
    RunUMAP(.,reduction = "harmony",dims = 1:npcs) %>%
    FindNeighbors(., reduction = "harmony", dims = 1:npcs,verbose = F) %>%
    FindClusters(., resolution = res,algorithm=4,verbose = F,random.seed=1234)
DimPlot(obj,reduction = 'umap',group.by = 'RNA_snn_res.0.8',cols = colors_90,raster = F,label=T,label.size=5) + NoLegend() + ggtitle('Myeloid')
obj = JoinLayers(obj)
Idents(obj) = "RNA_snn_res.0.8"
DotPlot(obj,
  features = c('LYZ','CD68','CD163','CSF3R','CPA3','LILRA4','CLEC9A','CD1C','CCR7','FCGR3A','CD14','FCGR3B',"MKI67"),
  cols = 'RdBu') +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
Macro.obj = subset(obj,idents = c(1,2,4,6,8,10,13,14))
markers.M = cosg(obj,
            groups='all',
            assay = 'RNA',
            mu=100,
            remove_lowly_expressed = TRUE,
            expressed_pct = 0.3,
            n_genes_user=50)
View(markers.M$names)
markers.Macro = cosg(Macro.obj,
            groups='all',
            assay = 'RNA',
            mu=100,
            remove_lowly_expressed = TRUE,
            expressed_pct = 0.3,
            n_genes_user=50)
View(markers.Macro$names)
DotPlot(Macro.obj,
  features = markers.Macro$names[1:20,'14'],
  cols = 'RdBu') +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
Myeloid.obj = obj
Myeloid.obj = subset(Myeloid.obj,idents = c(15),invert=T)
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(1)) ) <- "m01_Macro_B2M"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(2)) ) <- "m02_Macro_LYVE1"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(4)) ) <- "m03_Macro_SPP1"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(8)) ) <- "m04_Macro_CXCL10"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(10)) ) <- "m05_Macro_SAT1"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(14)) ) <- "m06_Macro_MMP9"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(13)) ) <- "m07_Macro_MKI67"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(6)) ) <- "m08_Mono_VCAN"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(12)) ) <- "m09_cDC1_CLEC9A"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(7)) ) <- "m10_cDC2_CD1C"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(9)) ) <- "m11_cDC3_CCR7"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(11)) ) <- "m12_pDC_LILRA4"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(3)) ) <- "m13_Nutri_CXCR1"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(5)) ) <- "m14_MAST_CPA3"
Myeloid.obj$subcelltype = Idents(Myeloid.obj)
Idents(Myeloid.obj) = "subcelltype"
Myeloid.obj$Group = paste0(Myeloid.obj$Time_point,'_',Myeloid.obj$Outcome)
qs_save(Myeloid.obj,file.path(output_dir, 'Myeloid_obj_intergration_new.qs'))
DimPlot(Myeloid.obj,reduction = 'umap',cols = colors_90,raster = F,alpha=0.7) + ggtitle('Myeloid') + theme(aspect.ratio = 1)
ggsave(file.path(output_dir, "Figures",'Myeloid_subcelltype_dimplot_new.pdf'),width = 8,height =6)
plot_density_dimplot_new(seurat_object = Myeloid.obj,
                    group.by = 'subcelltype',
                    split.by = 'Group',
                    ncol = 2,
                    title = 'Myeloid',
                    color_palette = 'magma',group_colors = colors_90)
ggsave(file.path(output_dir,"Figures", 'Myeloid_subcelltype_density_new.png'),width = 14,height =12,dpi = 300)




