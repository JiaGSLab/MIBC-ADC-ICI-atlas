# ==================================================
# Script: T cell subclustering and annotation
# Author: Wenhui Chen
# Purpose: Performs T cell subclustering, manual annotation, and related visualization.
# Input: Tcell_obj_intergration.qs and shared variables created in the preprocessing workflow.
# Output: Annotated T cell object and T cell subtype plots.
# Note: Run scripts in workflow order. Keep original path definitions unless updating the project layout.
# ==================================================
#T cells

Tcell.obj = qs_read(file.path(output_dir, 'Tcell_obj_intergration.qs'))
obj = Tcell.obj
prefix = 'Tcell'
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
ElbowPlot(obj,ndims=50)
obj = obj %>%
    IntegrateLayers(.,method = HarmonyIntegration,
                    new.reduction = "harmony",orig.reduction = "pca",max.iter.harmony = 20L,verbose=T) %>%
    RunUMAP(.,reduction = "harmony",dims = 1:npcs) %>%
    FindNeighbors(., reduction = "harmony", dims = 1:npcs,verbose = F) %>%
    FindClusters(., resolution = res,algorithm=4,verbose = F,random.seed=1234)
DimPlot(obj,reduction = 'umap',group.by = 'RNA_snn_res.1',cols = colors_90,raster = F,label=T,label.size=5) + NoLegend() + ggtitle('Tcell')
obj = JoinLayers(obj)
Idents(obj) = "RNA_snn_res.1"
DotPlot(obj,
  features = c('CD2','CD3D','CD3E','CD4','CD8A','CD8B','TRDC','TRAC'),
  cols = 'RdBu') +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
markers.T = cosg(obj,
              groups='all',
              assay = 'RNA',
              mu=100,
              remove_lowly_expressed = TRUE,
              expressed_pct = 0.3,
              n_genes_user=50)
View(markers.T$names)
CD8T.obj = subset(obj,idents = c(1,3,6,7,12,15,16))
CD4T.obj = subset(obj,idents = c(2,4,5,8,9,10,14))
gdT.obj = subset(obj,idents = c(11,13))
markers.CD8T = cosg(CD8T.obj,
            groups='all',
            assay = 'RNA',
            mu=100,
            remove_lowly_expressed = TRUE,
            expressed_pct = 0.3,
            n_genes_user=50)
View(markers.CD8T$names)
markers.CD4T = cosg(CD4T.obj,
            groups='all',
            assay = 'RNA',
            mu=100,
            remove_lowly_expressed = TRUE,
            expressed_pct = 0.3,
            n_genes_user=50)
View(markers.CD4T$names)
markers.gdT = cosg(gdT.obj,
            groups='all',
            assay = 'RNA',
            mu=100,
            remove_lowly_expressed = TRUE,
            expressed_pct = 0.3,
            n_genes_user=50)
View(markers.gdT$names)
DotPlot(obj,
  features = markers.gdT$names[1:20,'13'],
  cols = 'RdBu') +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
Tcell.obj = obj
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(15))) <- "t01_CD8_IL7R"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(1))) <- "t02_CD8_AKNA"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(3))) <- "t03_CD8_EOMES"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(7))) <- "t04_CD8_ZNF683"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(6))) <- "t05_CD8_HAVCR2"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(12))) <- "t06_CD8_ISG"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(16))) <- "t07_CD8_MKI67"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(8))) <- "t08_CD4_IL7R"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(2))) <- "t09_CD4_TCF7"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(9))) <- "t10_CD4_CXCL13"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(10))) <- "t11_CD4_RORA"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(4))) <- "t12_CD4_ETS1"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(5))) <- "t13_CD4_TNFRSF9"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(14))) <- "t14_CD4_ISG"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(11))) <- "t15_gdT_CX3CR1"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(13))) <- "t16_gdT_FCER1G"
Tcell.obj$subcelltype = Idents(Tcell.obj)
Tcell.obj$Group = paste0(Tcell.obj$Time_point,'_',Tcell.obj$Outcome)
qs_save(Tcell.obj,file.path(output_dir, 'Tcell_obj_intergration_new.qs'))
DimPlot(Tcell.obj,reduction = 'umap',cols = colors_90,raster = F,alpha=0.7) + ggtitle('Tcell') + theme(aspect.ratio = 1)
ggsave(file.path(output_dir,"Figures", 'Tcell_subcelltype_dimplot_new.pdf'),width = 8,height =6)
plot_density_dimplot_new(seurat_object = Tcell.obj,
                    group.by = 'subcelltype',
                    split.by = 'Group',
                    ncol = 2,
                    title = 'Tcell',
                    color_palette = 'magma',group_colors = colors_90)
ggsave(file.path(output_dir,"Figures", 'Tcell_subcelltype_density_new.png'),width = 14,height =12,dpi = 300)


obj <- RunCytoTRACE(obj,species = "human",cores = 10)
pdf(file.path(output_dir, "Figures",'Tcell_CytoTRACE.pdf'),width = 16,height = 8)
CytoTRACEPlot(obj,group.by='RNA_snn_res.1')
dev.off()
obj <- RunSlingshot(obj,
                  group.by = "RNA_snn_res.1",
                  reduction = "umap"
                )
obj <- RunMonocle2(obj,
                  group.by = "RNA_snn_res.1"
                )
names(obj@tools$Monocle2)
trajectory <-obj@tools$Monocle2$trajectory
p1 <- CellDimPlot(
                  obj,
                  group.by = "Monocle2_State",
                  reduction = "DDRTree",
                  label = TRUE,
                  theme_use = "theme_blank"
                )
p1 + trajectory
p1 <- CellDimPlot(
                  obj,
                  group.by = "RNA_snn_res.1",
                  reduction = "DDRTree",
                  label = TRUE,
                  theme_use = "theme_blank"
                )
p1 + trajectory
p1 <- CellDimPlot(
                  obj,
                  group.by = "Group",
                  reduction = "DDRTree",
                  label = TRUE,
                  theme_use = "theme_blank"
                )
p1 + trajectory  
FeatureDimPlot(
    obj,
    features = "Monocle2_Pseudotime",
    reduction = "umap",
    theme_use = "theme_blank"
  )
CD8T.obj = subset(obj,idents=c(1,3,4,7))
markers.CD8T = cosg(CD8T.obj,
              groups='all',
              assay = 'RNA',
              mu=100,
              remove_lowly_expressed = TRUE,
              expressed_pct = 0.3,
              n_genes_user=100)
View(markers.CD8T$names)


