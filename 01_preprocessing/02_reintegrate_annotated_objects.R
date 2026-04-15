# ==================================================
# Script: Reintegrate annotated cell type objects
# Author: Wenhui Chen
# Purpose: Merges annotated compartment-level objects back into a filtered integrated object for downstream analyses.
# Input: Annotated Myeloid, Epi, Stromal, Tcell, and Bcell Seurat objects.
# Output: Final integrated Seurat object with updated celltype and subcelltype labels.
# Note: Run scripts in workflow order. Keep original path definitions unless updating the project layout.
# ==================================================
#Merge Epi, Stromal, Myeloid, T, and B Seurat objects
Myeloid.obj = qs_read(file.path(output_dir, 'Myeloid_obj_intergration_new.qs'))
Epi.obj = qs_read(file.path(output_dir, 'Epi_obj_intergration.qs'))
Stromal.obj = qs_read(file.path(output_dir, 'Stromal_obj_intergration_new.qs'))
Tcell.obj = qs_read(file.path(output_dir, 'Tcell_obj_intergration_new.qs'))
Bcell.obj = qs_read(file.path(output_dir, 'Bcell_obj_intergration_new.qs'))

seu.combined.IT = merge(
    x = Myeloid.obj,
    y = list(Epi.obj, Stromal.obj, Tcell.obj, Bcell.obj),
    add.cell.ids = NULL,   # because prefixes have already been added manually
    merge.data = TRUE
)
seu.combined.IT = JoinLayers(seu.combined.IT)
seu.combined.IT = subcelltype_analysis(obj = seu.combined.IT,
                                        output_dir = output_dir,
                                        prefix = 'seu_combined_IT_filtered',
                                        npcs = 50,
                                        res = c(1))
meta = seu.combined.IT@meta.data
seu.combined.IT = qs_read('/home/data/chenwenhui/scRNAseq/GPU_IT/QC/seu_combined_IT_filtered_obj_intergration.qs')
seu.combined.IT = subset(seu.combined.IT,cells=rownames(meta))
seu.combined.IT@meta.data$subcelltype = meta[rownames(seu.combined.IT@meta.data),'subcelltype']
Idents(seu.combined.IT) = 'subcelltype'
seu.combined.IT = subset(seu.combined.IT,idents=c('s16_Schwann_SOX10'),invert=T)
seu.combined.IT$major = ifelse(seu.combined.IT$celltype == 'Epithelial', 'Epithelial',
                              ifelse(seu.combined.IT$celltype %in% c('Endothelial','Fibroblast','Pericyte','SMC'), 'Stromal',
                              ifelse(seu.combined.IT$celltype %in% c('B','Plasma'), 'B cell',
                              ifelse(seu.combined.IT$celltype %in% c('T'), 'T cell',
                              ifelse(seu.combined.IT$celltype %in% c('Myeloid','Mast','Neutrophil'), 'Myeloid','unknown')))))
seu.combined.IT$major = ifelse((seu.combined.IT$celltype == 'Epithelial' & seu.combined.IT$copykat_prediction == 'aneuploid'),'malignant',seu.combined.IT$major)
seu.combined.IT$celltype = ifelse((seu.combined.IT$celltype == 'Epithelial' & seu.combined.IT$copykat_prediction == 'aneuploid'),'malignant',seu.combined.IT$celltype)
seu.combined.IT$subcelltype = ifelse(seu.combined.IT$celltype == 'malignant','malignant',seu.combined.IT$subcelltype)
qs_save(seu.combined.IT,file.path(output_dir, 'seu_combined_IT_filtered_obj_intergration.qs'))

B.subcelltype = unique(seu.combined.IT$subcelltype)[grepl('_B_',unique(seu.combined.IT$subcelltype))]
Plasma.subcelltype = unique(seu.combined.IT$subcelltype)[grepl('_Pla_',unique(seu.combined.IT$subcelltype))]
Macro.subcelltype = unique(seu.combined.IT$subcelltype)[grepl('_Macro_',unique(seu.combined.IT$subcelltype))]
Mono.subcelltype = unique(seu.combined.IT$subcelltype)[grepl('_Mono_',unique(seu.combined.IT$subcelltype))]
cDC1.subcelltype = unique(seu.combined.IT$subcelltype)[grepl('_cDC1_',unique(seu.combined.IT$subcelltype))]
cDC2.subcelltype = unique(seu.combined.IT$subcelltype)[grepl('_cDC2_',unique(seu.combined.IT$subcelltype))]
cDC3.subcelltype = unique(seu.combined.IT$subcelltype)[grepl('_cDC3_',unique(seu.combined.IT$subcelltype))]
pDC.subcelltype = unique(seu.combined.IT$subcelltype)[grepl('_pDC_',unique(seu.combined.IT$subcelltype))]
Nutri.subcelltype = unique(seu.combined.IT$subcelltype)[grepl('_Nutri_',unique(seu.combined.IT$subcelltype))]
MAST.subcelltype = unique(seu.combined.IT$subcelltype)[grepl('_MAST_',unique(seu.combined.IT$subcelltype))]
Fib.subcelltype = unique(seu.combined.IT$subcelltype)[grepl('_Fib_',unique(seu.combined.IT$subcelltype))]
SMC.subcelltype = unique(seu.combined.IT$subcelltype)[grepl('_SMC_',unique(seu.combined.IT$subcelltype))]
vSMC.subcelltype = unique(seu.combined.IT$subcelltype)[grepl('_vSMC_',unique(seu.combined.IT$subcelltype))]
Per.subcelltype = unique(seu.combined.IT$subcelltype)[grepl('_Per_',unique(seu.combined.IT$subcelltype))]
Endo.subcelltype = unique(seu.combined.IT$subcelltype)[grepl('_Endo_',unique(seu.combined.IT$subcelltype))]
CD8T.subcelltype = unique(seu.combined.IT$subcelltype)[grepl('_CD8_',unique(seu.combined.IT$subcelltype))]
Treg.subcelltype = c("t12_CD4_ETS1",'t13_CD4_TNFRSF9')
CD4T.subcelltype = unique(seu.combined.IT$subcelltype)[grepl('_CD4_',unique(seu.combined.IT$subcelltype))]
CD4T.subcelltype = setdiff(CD4T.subcelltype,Treg.subcelltype)
gdT.subcelltype = unique(seu.combined.IT$subcelltype)[grepl('_gdT_',unique(seu.combined.IT$subcelltype))]


seu.combined.IT$celltype.new = seu.combined.IT$subcelltype
Idents(seu.combined.IT) = 'celltype.new'
Idents(seu.combined.IT,cells=WhichCells(seu.combined.IT,idents = B.subcelltype)) = 'B cell'
Idents(seu.combined.IT,cells=WhichCells(seu.combined.IT,idents = Plasma.subcelltype)) = 'Plasma cell'
Idents(seu.combined.IT,cells=WhichCells(seu.combined.IT,idents = Macro.subcelltype)) = 'Macrophage'
Idents(seu.combined.IT,cells=WhichCells(seu.combined.IT,idents = Mono.subcelltype)) = 'Monocyte'
Idents(seu.combined.IT,cells=WhichCells(seu.combined.IT,idents = cDC1.subcelltype)) = 'cDC1'
Idents(seu.combined.IT,cells=WhichCells(seu.combined.IT,idents = cDC2.subcelltype)) = 'cDC2'
Idents(seu.combined.IT,cells=WhichCells(seu.combined.IT,idents = cDC3.subcelltype)) = 'cDC3'
Idents(seu.combined.IT,cells=WhichCells(seu.combined.IT,idents = pDC.subcelltype)) = 'pDC'
Idents(seu.combined.IT,cells=WhichCells(seu.combined.IT,idents = Nutri.subcelltype)) = 'Neutrophil'
Idents(seu.combined.IT,cells=WhichCells(seu.combined.IT,idents = MAST.subcelltype)) = 'MAST cell'
Idents(seu.combined.IT,cells=WhichCells(seu.combined.IT,idents = Fib.subcelltype)) = 'Fibroblast'
Idents(seu.combined.IT,cells=WhichCells(seu.combined.IT,idents = SMC.subcelltype)) = 'SMC'
Idents(seu.combined.IT,cells=WhichCells(seu.combined.IT,idents = vSMC.subcelltype)) = 'vSMC'
Idents(seu.combined.IT,cells=WhichCells(seu.combined.IT,idents = Per.subcelltype)) = 'Pericyte'
Idents(seu.combined.IT,cells=WhichCells(seu.combined.IT,idents = Endo.subcelltype)) = 'Endothelial'
Idents(seu.combined.IT,cells=WhichCells(seu.combined.IT,idents = CD8T.subcelltype)) = 'CD8T cell'
Idents(seu.combined.IT,cells=WhichCells(seu.combined.IT,idents = CD4T.subcelltype)) = 'CD4T cell'
Idents(seu.combined.IT,cells=WhichCells(seu.combined.IT,idents = Treg.subcelltype)) = 'Treg'
Idents(seu.combined.IT,cells=WhichCells(seu.combined.IT,idents = gdT.subcelltype)) = 'gdT cell'
seu.combined.IT$celltype.new = Idents(seu.combined.IT)
Idents(seu.combined.IT) = 'subcelltype'
qs_save(seu.combined.IT,file.path(output_dir, 'seu_combined_IT_filtered_obj_intergration.qs'))

#####################
