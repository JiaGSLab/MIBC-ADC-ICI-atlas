# ==================================================
# Script: Stromal subclustering and annotation
# Author: Wenhui Chen
# Purpose: Performs stromal subclustering, manual annotation, and related visualization.
# Input: Stromal_obj_intergration.qs and shared variables created in the preprocessing workflow.
# Output: Annotated stromal object and subtype plots.
# Note: Run scripts in workflow order. Keep original path definitions unless updating the project layout.
# ==================================================
##################################

Stromal.obj = qs_read(file.path(output_dir, 'Stromal_obj_intergration.qs'))
obj = Stromal.obj
prefix = 'Stromal'
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
DimPlot(obj,reduction = 'umap',group.by = 'RNA_snn_res.0.8',cols = colors_90,raster = F,label=T,label.size=5) + NoLegend() + ggtitle('Stromal')
obj = JoinLayers(obj)
Idents(obj) = "RNA_snn_res.0.8"
DotPlot(obj,
  features = c('COL1A1','COL3A1','DCN','LUM','VWF','PECAM1','CDH5','ACTA2','TAGLN','MYH11','RGS5','PDGFRB',"MKI67","CPA3"),
  cols = 'RdBu') +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
markers.S = cosg(obj,
            groups='all',
            assay = 'RNA',
            mu=100,
            remove_lowly_expressed = TRUE,
            expressed_pct = 0.3,
            n_genes_user=50)
View(markers.S$names)
Fib.obj = subset(obj,idents=c(1,2,4,9,14,15))
Endo.obj = subset(obj,idents=c(5,7,10,11,13))
SMC_Per=subset(obj,idents=c(3,6,8,12))
markers.Fib = cosg(Fib.obj,
            groups='all',
            assay = 'RNA',
            mu=100,
            remove_lowly_expressed = TRUE,
            expressed_pct = 0.3,
            n_genes_user=50)
View(markers.Fib$names)
markers.Endo = cosg(Endo.obj,
          groups='all',
          assay = 'RNA',
          mu=100,
          remove_lowly_expressed = TRUE,
          expressed_pct = 0.3,
          n_genes_user=50)
View(markers.Endo$names)
markers.SMC_Per = cosg(SMC_Per,
          groups='all',
          assay = 'RNA',
          mu=100,
          remove_lowly_expressed = TRUE,
          expressed_pct = 0.3,
          n_genes_user=50)
View(markers.SMC_Per$names)
DotPlot(Fib.obj,
  features = markers.Fib$names[1:20,'15'],
  cols = 'RdBu') +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(Endo.obj,
  features = markers.Endo$names[1:20,'11'],
  cols = 'RdBu') +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(SMC_Per,
  features = markers.SMC_Per$names[1:20,'12'],
  cols = 'RdBu') +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
Stromal.obj = obj
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(1))) <- "s01_Fib_WNT4"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(2))) <- "s02_Fib_COL14A1"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(4))) <- "s03_Fib_S100A6"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(9))) <- "s04_Fib_POSTN"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(14))) <- "s05_Fib_IL11"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(15))) <- "s06_Fib_CCL19"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(3))) <- "s07_SMC_DES"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(12))) <- "s08_vSMC_MUSTN1"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(6))) <- "s09_Per_PDGFRB"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(5))) <- "s10_Endo_ACKR1"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(7))) <- "s11_Endo_TMSB4X"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(10))) <- "s12_Endo_SEMA3G"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(11))) <- "s13_Endo_RGCC"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(13))) <- "s14_Endo_PROX1"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(8))) <- "s15_vSMC_MKI67"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(16))) <- "s16_Schwann_SOX10"
Stromal.obj$subcelltype = Idents(Stromal.obj)
Stromal.obj$Group = paste0(Stromal.obj$Time_point,'_',Stromal.obj$Outcome)
qs_save(Stromal.obj,file.path(output_dir, 'Stromal_obj_intergration_new.qs'))
DimPlot(Stromal.obj,reduction = 'umap',cols = colors_90,raster = F,alpha=0.7) + ggtitle('Stromal') + theme(aspect.ratio = 1)
ggsave(file.path(output_dir,"Figures", 'Stromal_subcelltype_dimplot_new.pdf'),width = 8,height =6)
plot_density_dimplot_new(seurat_object = Stromal.obj,
                    group.by = 'subcelltype',
                    split.by = 'Group',
                    ncol = 2,
                    title = 'Stromal',
                    color_palette = 'magma',group_colors = colors_90)
ggsave(file.path(output_dir,"Figures", 'Stromal_subcelltype_density_new.png'),width = 14,height =12,dpi = 300)



######
