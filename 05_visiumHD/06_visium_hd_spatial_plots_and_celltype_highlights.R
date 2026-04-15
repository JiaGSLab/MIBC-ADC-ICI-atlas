# ==================================================
# Script: Visium HD spatial plots and celltype highlight views
# Author: Wenhui Chen
# Purpose: Generates celltype highlight plots and other spatial visualization outputs from the annotated Visium HD slide objects.
# Input: Annotated object_slide01, object_slide02, and object_slide03 objects created earlier.
# Output: Spatial plot files for celltype overview and selected highlight combinations.
# Note: Review file paths before running. Keep original analysis logic unchanged.
# ==================================================
##########################
Idents(object_slide01) = 'first_type'
celltype = unique(object_slide01@meta.data$first_type)
celltype = setdiff(celltype,c('unknown','reject'))
celltype.list = lapply(celltype,function(x) {
  rownames(object_slide01@meta.data[object_slide01@meta.data$first_type == x | object_slide01@meta.data$second_type == x,])})
names(celltype.list) = celltype
celltype.list = celltype.list[celltype.list != '']
p <- SpatialDimPlot(object_slide01,
  cells.highlight = celltype.list,
  cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
  ,pt.size.factor = 4,ncol = 5) + NoLegend()
p
ggsave(p, filename = "~/visiumHD/slide01_celltype.png", width = 30, height = 36, dpi = 300)
Endo_KRT5 = celltype.list[c('Endothelial','e03_Epi_KRT5')]
p <- SpatialDimPlot(object_slide01,
  cells.highlight = Endo_KRT5,
  cols.highlight = c("#FFFF00",'#8dd3c7', "grey50")) + NoLegend()
p
ggsave(p, filename = "~/visiumHD/slide01_Endo_KRT5.png", width = 10, height = 10, dpi = 300)

Idents(object_slide02) = 'first_type'
celltype.list = lapply(celltype,function(x) {
  rownames(object_slide02@meta.data[object_slide02@meta.data$first_type == x | object_slide02@meta.data$second_type == x,])})
names(celltype.list) = celltype
celltype.list = celltype.list[celltype.list != '']
#cells <- CellsByIdentities(object_slide02, idents = c('Endothelial','e03_Epi_KRT5'))
p <- SpatialDimPlot(object_slide02,
  cells.highlight = celltype.list,
  cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
  ,pt.size.factor = 4,ncol = 5) + NoLegend()
p
ggsave(p, filename = "~/visiumHD/slide02_celltype.png", width = 30, height = 36, dpi = 300)
Endo_KRT5 = celltype.list[c('Endothelial','e03_Epi_KRT5')]
p <- SpatialDimPlot(object_slide02,
  cells.highlight = Endo_KRT5,
  cols.highlight = c("#FFFF00",'#8dd3c7', "grey50")) + NoLegend()
p
ggsave(p, filename = "~/visiumHD/slide02_Endo_KRT5.png", width = 10, height = 10, dpi = 300)


Idents(object_slide03) = 'first_type'
celltype.list = lapply(celltype,function(x) {
  rownames(object_slide03@meta.data[object_slide03@meta.data$first_type == x | object_slide03@meta.data$second_type == x,])})
names(celltype.list) = celltype
celltype.list = celltype.list[celltype.list != '']
#cells <- CellsByIdentities(object_slide03, idents = c('Endothelial','e03_Epi_KRT5'))
p <- SpatialDimPlot(object_slide03,
  cells.highlight = celltype.list,
  cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
  ,pt.size.factor = 4,ncol = 5) + NoLegend()
p
ggsave(p, filename = "~/visiumHD/slide03_celltype.png", width = 30, height = 36, dpi = 300)
Endo_KRT5 = celltype.list[c('Endothelial','e03_Epi_KRT5')]
p <- SpatialDimPlot(object_slide03,
  cells.highlight = Endo_KRT5,
  cols.highlight = c("#FFFF00",'#8dd3c7', "grey50")) + NoLegend()
p
ggsave(p, filename = "~/visiumHD/slide03_Endo_KRT5.png", width = 10, height = 10, dpi = 300)


#################
