# ==================================================
# Script: Attach RCTD labels to Visium HD objects and merge slides
# Author: Wenhui Chen
# Purpose: Loads per-slide RCTD outputs, adds spot annotations and sample metadata to Visium HD Seurat objects, then merges slides into a combined object.
# Input: Per-slide RCTD result files, Visium HD binned outputs, and sample barcode mapping tables.
# Output: Per-slide Seurat objects with RCTD labels and the merged ST.combined.obj object.
# Note: Review file paths before running. Keep original analysis logic unchanged.
# ==================================================
if (TRUE) {
  suppressPackageStartupMessages({library(arrow)
    library(spacexr)
    library(ggplot2)
    library(Seurat)
    library(dplyr)
    library(forcats)
    library(tidyr)
    library(purrr)
    library(tibble)
    library(qs2)
    })
  options(future.globals.maxSize = 600 * 1024^3)
  set.seed(12356)
  colors_15 <- c("#A52A2A", "#053E7A", "#E46B4F", "#8B8B00", "#B8860B",  
                "#CD661D", "#8B5A2B", "#CDAD00", "#2F4F4F", "#483D8B", 
                "#7A378B", "#4A7088", "#556B2F", "#008B8B", "#CD7054")
  load('~/scRNAseq/colors_90.RData') #this is a color list and can be commented out if not needed
  output_dir <- "~/visiumHD"
  Ref.sample = c('ZZJ-1','DYZ-1_1114','LJS-1_1114','LY-1_1114','FRD-1_1114',
                'LK-1_1114','ZY-1_1114','HYS-1_1114','LSY-1_1114','LS-1')
  setwd(output_dir)
}


#############################
#Add RCTD results to the Seurat object
load("~/visiumHD/RCTD_tissue1_celltype.RData")
object_slide01 <- Load10X_Spatial(data.dir = '~/visiumHD/FWSP20252181/FWSP20252181-1/binned_outputs/square_008um')
SpatialFeaturePlot(object_slide01, features = 'nFeature_Spatial',min.cutoff = 50) + theme(legend.position = "right")
RCTD.result = myRCTD@results$results_df[,c('spot_class','first_type','second_type')]
RCTD.result$spot_class = as.character(RCTD.result$spot_class)
RCTD.result$first_type = as.character(RCTD.result$first_type)
RCTD.result$second_type = as.character(RCTD.result$second_type) 
#object_slide01 = subset(object_slide01, subset = nFeature_Spatial >50)
object_slide01$orig.ident = "slide01"
object_slide01$spot_class = 'reject'
object_slide01$first_type = 'unknown'
object_slide01$second_type = 'unknown'
meta = object_slide01@meta.data
RCTD.result  = RCTD.result[intersect(rownames(RCTD.result),rownames(meta)),]
meta[rownames(RCTD.result),'spot_class'] = RCTD.result$spot_class
meta[rownames(RCTD.result),'first_type'] = RCTD.result$first_type
meta[rownames(RCTD.result),'second_type'] = RCTD.result$second_type
meta$first_type = ifelse(meta$spot_class == 'reject', 'unknown', meta$first_type)
meta$second_type = ifelse(meta$spot_class == 'doublet_certain', meta$second_type, 'unknown')
table(meta$spot_class)
table(meta$first_type)
table(meta$second_type)
#object_slide01@meta.data = meta
#Add spatial position information to metadata
tissue_pos <-read_parquet("~/visiumHD/FWSP20252181/FWSP20252181-1/binned_outputs/square_008um/spatial/rename_tissue_positions.parquet")
tissue_pos=as.data.frame(tissue_pos)
rownames(tissue_pos)=tissue_pos$barcode
all_in  = rownames(meta) %in% rownames(tissue_pos)
sum(all_in) == length(rownames(meta))
tissue_pos = tissue_pos[rownames(meta),]
tissue_pos$barcode = NULL
meta = cbind(meta, tissue_pos)
View(meta)


#Add patient information to metadata
meta$sample = 'unknown'
meta$Therapy = 'IT'
meta$Outcome = 'NR'
meta$sample_id = 'pt'
CJW_barcode = read.csv('~/visiumHD/sample_location/CJW_6652.csv',header = T)
colnames(CJW_barcode) = c('Barcode','sample')
rownames(CJW_barcode) = CJW_barcode$Barcode
CJW_barcode = CJW_barcode[intersect(rownames(CJW_barcode),rownames(meta)),]
meta[rownames(CJW_barcode),'sample'] = CJW_barcode$sample
meta[rownames(CJW_barcode),'Outcome'] = 'NR'
meta[rownames(CJW_barcode),'sample_id'] = 'pt22'


LS_barcode = read.csv('~/visiumHD/sample_location/LS_0991.csv',header = T)
colnames(LS_barcode) = c('Barcode','sample')
rownames(LS_barcode) = LS_barcode$Barcode
LS_barcode = LS_barcode[intersect(rownames(LS_barcode),rownames(meta)),]
meta[rownames(LS_barcode),'sample'] = LS_barcode$sample
meta[rownames(LS_barcode),'Outcome'] = 'CR'
meta[rownames(LS_barcode),'sample_id'] = 'pt28'
meta[rownames(LS_barcode),'Therapy'] = 'IC'

FRD_barcode = read.csv('~/visiumHD/sample_location/FRD_5923.csv',header = T)
colnames(FRD_barcode) = c('Barcode','sample')
rownames(FRD_barcode) = FRD_barcode$Barcode
FRD_barcode = FRD_barcode[intersect(rownames(FRD_barcode),rownames(meta)),]
meta[rownames(FRD_barcode),'sample'] = FRD_barcode$sample
meta[rownames(FRD_barcode),'Outcome'] = 'CR'
meta[rownames(FRD_barcode),'sample_id'] = 'pt36'

LK_barcode = read.csv('~/visiumHD/sample_location/LK_2297.csv',header = T)
colnames(LK_barcode) = c('Barcode','sample')
rownames(LK_barcode) = LK_barcode$Barcode
LK_barcode = LK_barcode[intersect(rownames(LK_barcode),rownames(meta)),]
meta[rownames(LK_barcode),'sample'] = LK_barcode$sample
meta[rownames(LK_barcode),'Outcome'] = 'CR'
meta[rownames(LK_barcode),'sample_id'] = 'pt39'

table(meta$sample)

object_slide01@meta.data = meta
object_slide01 = subset(object_slide01, subset = sample != 'unknown')
colnames(object_slide01) = paste0(object_slide01$sample,'_',colnames(object_slide01))
qs_save(object_slide01, file = '~/visiumHD/object_slide01.qs')


#################################
#slide02
load("~/visiumHD/RCTD_tissue2_celltype.RData")
object_slide02 <- Load10X_Spatial(data.dir = "~/visiumHD/FWSP20252181/FWSP20252181-2/binned_outputs/square_008um/")
SpatialFeaturePlot(object_slide02, features = 'nFeature_Spatial',min.cutoff = 50) + theme(legend.position = "right")
RCTD.result = myRCTD@results$results_df[,c('spot_class','first_type','second_type')]
RCTD.result$spot_class = as.character(RCTD.result$spot_class)
RCTD.result$first_type = as.character(RCTD.result$first_type)
RCTD.result$second_type = as.character(RCTD.result$second_type) 
#object_slide02 = subset(object_slide02, subset = nFeature_Spatial >50)
object_slide02$orig.ident = "slide02"
object_slide02$spot_class = 'reject'
object_slide02$first_type = 'unknown'
object_slide02$second_type = 'unknown'
meta = object_slide02@meta.data
RCTD.result = RCTD.result[intersect(rownames(RCTD.result),rownames(meta)),]
meta[rownames(RCTD.result),'spot_class'] = RCTD.result$spot_class
meta[rownames(RCTD.result),'first_type'] = RCTD.result$first_type
meta[rownames(RCTD.result),'second_type'] = RCTD.result$second_type
meta$first_type = ifelse(meta$spot_class == 'reject', 'unknown', meta$first_type)
meta$second_type = ifelse(meta$spot_class == 'doublet_certain', meta$second_type, 'unknown')
table(meta$spot_class)
table(meta$first_type)
table(meta$second_type)
#object_slide02@meta.data = meta
#Add spatial position information to metadata
tissue_pos <-read_parquet("~/visiumHD/FWSP20252181/FWSP20252181-2/binned_outputs/square_008um/spatial/rename_tissue_positions.parquet")
tissue_pos=as.data.frame(tissue_pos)
rownames(tissue_pos)=tissue_pos$barcode
all_in  = rownames(meta) %in% rownames(tissue_pos)
sum(all_in) == length(rownames(meta))
tissue_pos = tissue_pos[rownames(meta),]
tissue_pos$barcode = NULL
meta = cbind(meta, tissue_pos)
View(meta)

#Add patient information to metadata
meta$sample = 'unknown'
meta$Therapy = 'IT'
meta$Outcome = 'NR'
meta$sample_id = 'pt'
LY_barcode = read.csv('~/visiumHD/sample_location/LY_8646.csv',header = T)
colnames(LY_barcode) = c('Barcode','sample')
rownames(LY_barcode) = LY_barcode$Barcode
LY_barcode = LY_barcode[intersect(rownames(LY_barcode),rownames(meta)),]
meta[rownames(LY_barcode),'sample'] = LY_barcode$sample
meta[rownames(LY_barcode),'Outcome'] = 'CR'
meta[rownames(LY_barcode),'sample_id'] = 'pt35'

DYZ_barcode = read.csv('~/visiumHD/sample_location/DYZ_3639.csv',header = T)
colnames(DYZ_barcode) = c('Barcode','sample')
rownames(DYZ_barcode) = DYZ_barcode$Barcode
DYZ_barcode = DYZ_barcode[intersect(rownames(DYZ_barcode),rownames(meta)),]
meta[rownames(DYZ_barcode),'sample'] = DYZ_barcode$sample
meta[rownames(DYZ_barcode),'Outcome'] = 'NR'
meta[rownames(DYZ_barcode),'sample_id'] = 'pt31'

HYS_barcode = read.csv('~/visiumHD/sample_location/HYS_1561.csv',header = T)
colnames(HYS_barcode) = c('Barcode','sample')
rownames(HYS_barcode) = HYS_barcode$Barcode
HYS_barcode = HYS_barcode[intersect(rownames(HYS_barcode),rownames(meta)),]
meta[rownames(HYS_barcode),'sample'] = HYS_barcode$sample
meta[rownames(HYS_barcode),'Outcome'] = 'CR'
meta[rownames(HYS_barcode),'sample_id'] = 'pt42'

LSY_barcode = read.csv('~/visiumHD/sample_location/LSY_0068.csv',header = T)
colnames(LSY_barcode) = c('Barcode','sample')
rownames(LSY_barcode) = LSY_barcode$Barcode
LSY_barcode = LSY_barcode[intersect(rownames(LSY_barcode),rownames(meta)),]
meta[rownames(LSY_barcode),'sample'] = LSY_barcode$sample
meta[rownames(LSY_barcode),'Outcome'] = 'CR'
meta[rownames(LSY_barcode),'sample_id'] = 'pt44'

table(meta$sample)
object_slide02@meta.data = meta
object_slide02 = subset(object_slide02, subset = sample != 'unknown')
colnames(object_slide02) = paste0(object_slide02$sample,'_',colnames(object_slide02))
qs_save(object_slide02, file = '~/visiumHD/object_slide02.qs')  


####################################
#slide03
load("~/visiumHD/RCTD_tissue3_celltype.RData")
object_slide03 <- Load10X_Spatial(data.dir = "~/visiumHD/FWSP20252181/FWSP20252181-3/binned_outputs/square_008um/")
SpatialFeaturePlot(object_slide03, features = 'nFeature_Spatial',min.cutoff = 50) + theme(legend.position = "right")
RCTD.result = myRCTD@results$results_df[,c('spot_class','first_type','second_type')]
RCTD.result$spot_class = as.character(RCTD.result$spot_class)
RCTD.result$first_type = as.character(RCTD.result$first_type)
RCTD.result$second_type = as.character(RCTD.result$second_type) 
#object_slide03 = subset(object_slide03,subset = nFeature_Spatial >50)
object_slide03$orig.ident = "slide03"
object_slide03$spot_class = 'reject'
object_slide03$first_type = 'unknown'
object_slide03$second_type = 'unknown'
meta = object_slide03@meta.data
RCTD.result = RCTD.result[intersect(rownames(RCTD.result),rownames(meta)),]
meta[rownames(RCTD.result),'spot_class'] = RCTD.result$spot_class
meta[rownames(RCTD.result),'first_type'] = RCTD.result$first_type
meta[rownames(RCTD.result),'second_type'] = RCTD.result$second_type
meta$first_type = ifelse(meta$spot_class == 'reject', 'unknown', meta$first_type)
meta$second_type = ifelse(meta$spot_class == 'doublet_certain', meta$second_type, 'unknown')
table(meta$spot_class)
table(meta$first_type)
table(meta$second_type)
#object_slide03@meta.data = meta
#Add spatial position information to metadata
tissue_pos <-read_parquet("~/visiumHD/FWSP20252181/FWSP20252181-3/binned_outputs/square_008um/spatial/rename_tissue_positions.parquet")
tissue_pos=as.data.frame(tissue_pos)
rownames(tissue_pos)=tissue_pos$barcode
all_in  = rownames(meta) %in% rownames(tissue_pos)
sum(all_in) == length(rownames(meta))
tissue_pos = tissue_pos[rownames(meta),]
tissue_pos$barcode = NULL
meta = cbind(meta, tissue_pos)
View(meta)

#Add patient information to metadata
meta$sample = 'unknown'
meta$Therapy = 'IT'
meta$Outcome = 'NR'
meta$sample_id = 'pt'
LJS_barcode = read.csv('~/visiumHD/sample_location/LJS_2573.csv',header = T)
colnames(LJS_barcode) = c('Barcode','sample')
rownames(LJS_barcode) = LJS_barcode$Barcode
LJS_barcode = LJS_barcode[intersect(rownames(LJS_barcode),rownames(meta)),]
meta[rownames(LJS_barcode),'sample'] = LJS_barcode$sample
meta[rownames(LJS_barcode),'Outcome'] = 'CR'
meta[rownames(LJS_barcode),'sample_id'] = 'pt32'

ZZJ_barcode = read.csv('~/visiumHD/sample_location/ZZJ_2362.csv',header = T)
colnames(ZZJ_barcode) = c('Barcode','sample')
rownames(ZZJ_barcode) = ZZJ_barcode$Barcode
ZZJ_barcode = ZZJ_barcode[intersect(rownames(ZZJ_barcode),rownames(meta)),]
meta[rownames(ZZJ_barcode),'sample'] = ZZJ_barcode$sample
meta[rownames(ZZJ_barcode),'Outcome'] = 'NR'
meta[rownames(ZZJ_barcode),'sample_id'] = 'pt23'

ZZJ_barcode = read.csv('~/visiumHD/sample_location/ZZJ_8916.csv',header = T)
colnames(ZZJ_barcode) = c('Barcode','sample')
rownames(ZZJ_barcode) = ZZJ_barcode$Barcode
ZZJ_barcode = ZZJ_barcode[intersect(rownames(ZZJ_barcode),rownames(meta)),]
meta[rownames(ZZJ_barcode),'sample'] = ZZJ_barcode$sample
meta[rownames(ZZJ_barcode),'Outcome'] = 'CR'
meta[rownames(ZZJ_barcode),'sample_id'] = 'pt34'

ZS_barcode = read.csv('~/visiumHD/sample_location/ZS_6423.csv',header = T)
colnames(ZS_barcode) = c('Barcode','sample')
rownames(ZS_barcode) = ZS_barcode$Barcode
ZS_barcode = ZS_barcode[intersect(rownames(ZS_barcode),rownames(meta)),]
meta[rownames(ZS_barcode),'sample'] = ZS_barcode$sample
meta[rownames(ZS_barcode),'Outcome'] = 'NR'
meta[rownames(ZS_barcode),'sample_id'] = 'pt40'

table(meta$sample)
object_slide03@meta.data = meta
object_slide03 = subset(object_slide03, subset = sample != 'unknown')
SpatialFeaturePlot(object_slide03, features = 'nFeature_Spatial',min.cutoff = 50) + theme(legend.position = "right")
colnames(object_slide03) = paste0(object_slide03$sample,'_',colnames(object_slide03))
qs_save(object_slide03, file = '~/visiumHD/object_slide03.qs')

combined.obj = merge(x = object_slide01,y=list(object_slide02,object_slide03),project = "VisiumHD_Project")
qs_save(combined.obj,file = '~/visiumHD/ST.combined.obj.qs')


