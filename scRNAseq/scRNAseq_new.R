# >>> GENERATED FILE NOTE
# scRNAseq_new.R is a minimally organized copy of scRNAseq.R.
# The original scRNAseq.R is unchanged.
# Only figure-producing blocks already consolidated in Fig1.R-Fig5.R are commented below.
# Upstream object construction, integration, scoring and result generation remain executable.
# Each commented block records its original scRNAseq.R line range and destination figure script.
# <<< END GENERATED FILE NOTE

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
  load('~/scRNAseq/colors_90.RData') #是一个颜色list，可以注释掉
  colors_90 = colors_91
  #load('~/scRNAseq/QC/blacklist.RData')  #TCR 1:247, BCR 248:683，TCR和免疫球蛋白，還有幾個是綫粒體基因和紅細胞基因,在高变基因中去掉有利于聚类，不去掉影响也不大，可以注释掉
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
  non_immune.samples = c('LLZ-1','LS-2','LHF-1','ZZS-2') #在分析非肿瘤细胞时，需要先过滤掉这些样本，因为细胞数量太少
  PR.samples = c('pt22','pt33','pt34','pt40','pt55','pt56','pt25','pt1','pt10') #确定有部分响应的患者
  output_dir <- "~/scRNAseq/GPU_IT/QC"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  setwd(output_dir)
  source("/home/data/chenwenhui/scRNAseq/GPU_IT/scripts/function_for_plot.R")
}

#
sample_info = read.csv('/home/data/chenwenhui/sample_info.csv', header = TRUE, stringsAsFactors = FALSE,fill = TRUE)
sample_info = sample_info[sample_info$scFFPE.seq==1,]
sample_info$Time_point = ifelse(sample_info$Time_point == 'pre-treatment','Pre','Post')
sample_info$Outcome = ifelse(sample_info$Outcome == 'response','CR','NR')
sample_info$idx = NULL

calculate_simple_gene_set_score <- function(seurat_obj, gene_set, assay = "RNA", slot = "data", score_name = "sig_score") {
  #判断gene_set是否为list，必须是list输入
  if (!is.list(gene_set)) {
    stop("gene_set 必须是 list 类型。")
  }

  #判断score_name是否是1一个元素或者和gene_set一样长，否则报错
  if (length(score_name) != 1 && length(score_name) != length(gene_set)) {
    stop("score_name 必须是单个元素或与 gene_set 长度相同。")
  }
  # 1. 检查基因集是否存在于数据中
  # 获取当前 Assay 中所有表达的基因
  expressed_genes <- rownames(Seurat::GetAssayData(object = seurat_obj, assay = assay, slot = slot))

  if (length(score_name) ==1) {
      score_name = paste0(score_name,1:length(gene_set),'_AveExp')
  } else {
      score_name = paste0(score_name,'_AveExp')
  }
  
  for (i in 1:length(gene_set)) {
    # 找出共同存在的基因
    valid_genes <- intersect(gene_set[[i]], expressed_genes)

    # 检查是否有基因匹配
    if (length(valid_genes) == 0) {
      print("提供的基因集 (gene_set) 中没有基因在 Seurat 对象的 Assay 数据中被找到。请检查基因名称。")
      next
    }

    if (length(valid_genes) < length(gene_set[[i]])) {
      warning(
        paste0(
          "注意：提供的基因集 (", length(gene_set[[i]]), "个基因) 中只有 ", 
          length(valid_genes), " 个基因存在于数据中。将使用这 ", 
          length(valid_genes), " 个基因进行计算。"
        )
      )
    }

    # 2. 提取表达矩阵
    expression_matrix <- Seurat::GetAssayData(object = seurat_obj, assay = assay, slot = slot)[valid_genes,]
    print(dim(expression_matrix))
    # 3. 计算每个细胞的平均表达量
    # matrixStats::rowMeans2 或 base::rowMeans 可以更快
    # 由于是转置矩阵，我们需要计算每一列的平均值
    # colMeans 计算每一列的平均值 (即每个细胞的平均值)
    gene_set_scores <- colMeans(expression_matrix)

    # 4. 添加到 Seurat 对象的元数据中
    seurat_obj@meta.data[score_name[i]] <- gene_set_scores

    message(paste0("成功计算并添加了评分：", score_name[i], " 到 Seurat 对象的元数据中。"))
  }
  return(seurat_obj)
}

QCcalulate <- function(obj, ratio=5) {
  meta <- obj@meta.data
  print(nrow(meta))
  
  # 清理 sample 字段（防止 invisible characters）
  meta$sample <- as.character(meta$sample)
  meta$sample <- trimws(meta$sample)
  
  meta$is_outlier <- FALSE
  print(table(meta$is_outlier))

  samples <- unique(meta$sample)

  for (sam in samples) {
    
    idx <- which(meta$sample == sam)
    sub <- meta[idx, , drop = FALSE]

    nFeature_thresh <- median(sub$nFeature_RNA,na.rm = TRUE) + ratio * mad(sub$nFeature_RNA,na.rm = TRUE)
    nCount_thresh   <- median(sub$nCount_RNA,na.rm = TRUE)   + ratio * mad(sub$nCount_RNA,na.rm = TRUE)
    pct_mt_thresh   <- median(sub$pct_counts_mt,na.rm = TRUE)+ 3 * mad(sub$pct_counts_mt,na.rm = TRUE)

    outlier_sub <- (sub$nFeature_RNA > nFeature_thresh) |
                  (sub$nFeature_RNA <= 200) |
                  (sub$nCount_RNA > nCount_thresh) |
                  (sub$pct_counts_mt > pct_mt_thresh) |
                  (sub$log10GenesPerUMI < 0.8) |
                  (sub$log10GenesPerUMI == 1) |
                  (sub$pct_counts_hb > 1)
    meta$is_outlier[idx] <- outlier_sub
    print(length(outlier_sub))
  }
  print(table(meta$is_outlier))
  obj@meta.data <- meta
  return(obj)
}

decontX_sample <- function(sample_name, input_dir, output_dir) {
  batch = ifelse(str_ends(sample_name,'-CX'),'FFPE',
                ifelse(str_ends(sample_name,'_1114'),'20241114_FFPE',
                  ifelse(str_ends(sample_name,'_0619'),'20250619_FFPE','snRNAseq')))
  if (batch == '20250619_FFPE') {
    dirs <- file.path(input_dir, gsub('_0619','',sample_name))
  } else {
    dirs <- file.path(input_dir, sample_name)
  }
  if (dir.exists(paste0(dirs,'/filtered_feature_bc_matrix'))) {
    counts_filtered <- Read10X(paste0(dirs,'/filtered_feature_bc_matrix'))
    #counts_raw <- Read10X(paste0(dirs,'/raw_feature_bc_matrix'))
  } else if (dir.exists(paste0(dirs,'/sample_filtered_feature_bc_matrix'))) {
    counts_filtered <- Read10X(paste0(dirs,'/sample_filtered_feature_bc_matrix'))
    #counts_raw <- Read10X(paste0(dirs,'/sample_raw_feature_bc_matrix'))
  } else {
    counts_filtered <- Read10X(paste0(dirs,'/count/sample_filtered_feature_bc_matrix'))
    #counts_raw <- Read10X(paste0(dirs,'/count/sample_raw_feature_bc_matrix'))
  }

  #counts_raw = counts_raw[rownames(counts_filtered),]
  sample_name = gsub('-CX', '', sample_name)

  # seurat_object = CreateSeuratObject(counts = counts_filtered, 
  #                                   project = sample_name,
  #                                   min.cells = 3,
  #                                   min.features = 0)

  sce = SingleCellExperiment(list(counts = counts_filtered))
  #sce.raw = SingleCellExperiment(list(counts = counts_raw))
  # decontx_result_bg <- celda::decontX(x = sce,
  #                                 background = sce.raw,
  #                                 maxIter = 1000,
  #                                 varGenes = 2000)
  decontx_result <- celda::decontX(x = sce,
                                maxIter = 1000,
                                varGenes = 2000)
  
  pdf(file = paste0(output_dir, '/figs/', sample_name, '_decontX_diagnostics.pdf'), width = 8, height = 6)
  #p = celda::plotDecontXContamination(decontx_result_bg) + ggtitle(paste0(sample_name, ' decontX with background'))
  #print(p)
  p2 = celda::plotDecontXContamination(decontx_result) + ggtitle(paste0(sample_name, ' decontX without background'))
  print(p2)
  dev.off()
  seurat_object = CreateSeuratObject(counts = round(decontXcounts(decontx_result)), 
                                    project = sample_name,
                                    min.cells = 3,
                                    min.features = 0)
  #seurat_object$decontX_contamination_bg <- colData(decontx_result_bg)$decontX_contamination
  # 保存 decontX 污染率
  seurat_object$decontX_contamination <- colData(decontx_result)$decontX_contamination

  # 建立 assay 并加入 Seurat object（assay 名称为 "decontX"）
  #seurat_object[["decontX_bg"]] <- Seurat::CreateAssayObject(counts = round(decontXcounts(decontx_result_bg)))
  #seurat_object[["decontX"]] <- Seurat::CreateAssayObject(counts = round(decontXcounts(decontx_result)))
  seurat_object$batch = batch
  seurat_object$sample = sample_name
  
  output_path <- paste0(output_dir, '/', sample_name, '.qs')
  qs_save(seurat_object,file = output_path)
}

#decontX部分是在GPU服务器运行，路径不一致
if (FALSE) {
  input_dir1 = "/data/home/chenwenhui/scRNAseq/rawdata_all"
  samples1 = list.dirs(input_dir1, full.names = FALSE, recursive = FALSE)
  output_dir = "/data/home/chenwenhui/scRNAseq/decontX_results_all"
  figs_dir = "/data/home/chenwenhui/scRNAseq/decontX_results_all/figs"
  dir.create(figs_dir, showWarnings = FALSE)
  dir.create(output_dir, showWarnings = FALSE)

  input_dir2 <- '/data/home/chenwenhui/scRNAseq/FFPE_20250619/AM20250604-F1/outs/per_sample_outs'
  samples2 <- list.dirs(input_dir2, full.names = FALSE, recursive = FALSE)
  samples2 <- paste0(samples2, '_0619')

  # 并行执行 samples1 的 decontX_sample
  library(future)
  library(future.apply)

  nworkers <- 12
  future::plan(future::multisession, workers = nworkers)

  results_samples1 <- future.apply::future_lapply(samples1, function(s) {
    tryCatch({
      decontX_sample(s, input_dir1, output_dir)
      list(sample = s, status = "OK")
    }, error = function(e) {
      list(sample = s, status = "ERROR", message = conditionMessage(e))
    })
  }, future.seed = TRUE)

  # 恢复为顺序计划，释放后台进程
  future::plan(sequential)

  # 返回并行执行的结果列表（每个元素包含 sample 和 status，失败时含 message）
  results_samples1
  # 并行执行 samples2 的 decontX_sample
  future::plan(future::multisession, workers = nworkers)

  results_samples2 <- future.apply::future_lapply(samples2, function(s) {
    tryCatch({
      decontX_sample(s, input_dir2, output_dir)
      list(sample = s, status = "OK")
    }, error = function(e) {
      list(sample = s, status = "ERROR", message = conditionMessage(e))
    })
  }, future.seed = TRUE)

  # 恢复为顺序计划，释放后台进程
  future::plan(sequential)

  # 返回并行执行的结果列表
  results_samples2

  gc()
}

# 设置并行计划
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
  #seu[['decontX_bg']] <- NULL  # 删除含有背景去污的 assay
  return(seu)
}, future.seed = TRUE)

# 恢复为顺序计划
future::plan(sequential)
gc()
names(seu.list) <- gsub('.qs','',samples_all)
seu.combined <- merge(seu.list[[1]], y = seu.list[-1])
seu.combined <- QCcalulate(seu.combined)
seu.combined <- subset(seu.combined,subset = is_outlier==FALSE)

seu.combined$orig.ident = seu.combined$sample
# 从 sample_info 中挑选出需要加入的列（排除已经存在于 meta.data 的列）
cols_to_add <- setdiff(colnames(sample_info), colnames(seu.combined@meta.data))

# 构建按细胞顺序对齐的元数据表（行名为细胞名）
meta_to_add <- sample_info |> 
  as.data.frame() |> 
  dplyr::select(dplyr::any_of(cols_to_add))

meta_matched <- meta_to_add[match(seu.combined$orig.ident, sample_info$orig.ident), , drop = FALSE]
rownames(meta_matched) <- colnames(seu.combined)

# 将匹配好的元数据添加到 Seurat 对象
seu.combined <- Seurat::AddMetaData(seu.combined, metadata = meta_matched)

# 返回更新后的对象
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

# 创建DotPlot和DimPlot并保存为同一个PNG
png(file.path(output_dir, "DotPlot_DimPlot_combined.png"), width = 3600, height = 1800, res = 150)

# 创建DotPlot
dotplot <- DotPlot(seu.combined.IT, features = markers,
        cols = 'RdBu',dot.scale = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())

# 创建DimPlot
dimplot <- DimPlot(seu.combined.IT, cols = colors_91, raster = F, alpha = 0.2, label = T, label.size = 5) + NoLegend()

# 使用patchwork将两个图组合在一起
combined_plot <- dotplot + dimplot +
  plot_layout(ncol = 2, widths = c(2, 1)) +
  plot_annotation(tag_levels = 'A')

print(combined_plot)
dev.off()

# 同时在控制台显示
#print(combined_plot)
# 创建QC组合图：nFeature_RNA、nCount_RNA、pct_counts_mt的FeaturePlot和is_outlier的DimPlot
png(file.path(output_dir, "QC_combined.png"), width = 3600, height = 3600, res = 150)

# 创建三个FeaturePlot
feature_plot1 <- FeaturePlot(seu.combined.IT, features = 'nFeature_RNA', cols = c('#87CEEB','#FFB6C1'), order = T, alpha = 0.2,raster=F) + 
  ggtitle('nFeature_RNA')

feature_plot2 <- FeaturePlot(seu.combined.IT, features = 'nCount_RNA', cols = c('#87CEEB','#FFB6C1'), order = T, alpha = 0.2,raster=F) + 
  ggtitle('nCount_RNA')

feature_plot3 <- FeaturePlot(seu.combined.IT, features = 'pct_counts_mt', cols = c('#87CEEB','#FFB6C1'), order = T, alpha = 0.2,raster=F) + 
  ggtitle('pct_counts_mt')

# 创建is_outlier的DimPlot，使用浅红和浅蓝色
#outlier_colors <- c('FALSE' = '#87CEEB', 'TRUE' = '#FFB6C1')  # 浅蓝和浅红

# 确保is_outlier是因子类型
#seu.combined.IT$is_outlier <- factor(seu.combined.IT$is_outlier, levels = c('FALSE', 'TRUE'))

dimplot_outlier <- DimPlot(seu.combined.IT, group.by = 'is_outlier', cols = c('#87CEEB','#FFB6C1'), 
                          raster = F, alpha = 0.2, label = F,split.by = 'is_outlier') + 
  ggtitle('Outlier Cells') 

# 使用patchwork将四个图组合成2x2布局
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


# 定义对象列表和对应的前缀
obj.list = list(Epi.obj, Myeloid.obj, Tcell.obj, Bcell.obj, Stromal.obj)
prefix.list = c('Epi', 'Myeloid', 'Tcell', 'Bcell', 'Stromal')

for (i in 1:length(obj.list)) {
  obj.list[[i]] = subcelltype_analysis(obj.list[[i]],output_dir,prefix.list[[i]])
}


#Bcell.obj = obj.list[[4]]
Bcell.obj = qs_read(file.path(output_dir, 'Bcell_obj_intergration.qs'))
DimPlot(Bcell.obj,reduction = 'umap',group.by = 'RNA_snn_res.0.8',cols = colors_90,raster = F,label=T,label.size=5) + NoLegend() + ggtitle('Bcell')
Idents(Bcell.obj)='RNA_snn_res.0.8'
markers.B = cosg(Bcell.obj,
                groups='all',
                assay = 'RNA',
                mu=100,
                remove_lowly_expressed = TRUE,
                expressed_pct = 0.3,
                n_genes_user=50)
View(markers.B$names)
DotPlot(Bcell.obj,
  features = c('MS4A1','CD19','MZB1','JCHAIN'),
  cols = 'RdBu') +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(Bcell.obj,
  features = markers.B$names[31:45,'5'],
  cols = 'RdBu') +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
#注释B细胞亚群
Bcell.obj = subset(Bcell.obj,idents = c(6,7,9,10,11,12,13),invert=T)
Idents(Bcell.obj, cells = WhichCells(Bcell.obj, idents = c(1))) <- "b01_B_IER5"
Idents(Bcell.obj, cells = WhichCells(Bcell.obj, idents = c(4))) <- "b02_B_ARHGEF1"
Idents(Bcell.obj, cells = WhichCells(Bcell.obj, idents = c(8))) <- "b03_B_IGHD"
Idents(Bcell.obj, cells = WhichCells(Bcell.obj, idents = c(14))) <- "b04_B_MKI67"
Idents(Bcell.obj, cells = WhichCells(Bcell.obj, idents = c(2))) <- "b05_Pla_SLC1A5"
Idents(Bcell.obj, cells = WhichCells(Bcell.obj, idents = c(3))) <- "b06_Pla_CD38"
Idents(Bcell.obj, cells = WhichCells(Bcell.obj, idents = c(5))) <- "b07_Pla_IGHG3"
Bcell.obj$subcelltype = Idents(Bcell.obj)
Bcell.obj$Group = paste0(Bcell.obj$Time_point,'_',Bcell.obj$Outcome)
#baocun
Bcell.obj = subcelltype_analysis(Bcell.obj,output_dir,'Bcell')
Idents(Bcell.obj) = 'subcelltype'
qs_save(Bcell.obj,file.path(output_dir, 'Bcell_obj_intergration.qs'))
DimPlot(Bcell.obj,reduction = 'umap',cols = colors_90,raster = F,alpha=0.7) + ggtitle('Bcell') + theme(aspect.ratio = 1)
ggsave(file.path(output_dir, 'Bcell_subcelltype_dimplot.pdf'),width = 8,height =6)
plot_density_dimplot_new(seurat_object = Bcell.obj,
                    group.by = 'subcelltype',
                    split.by = 'Group',
                    ncol = 2,
                    title = 'Bcell',
                    color_palette = 'magma',group_colors = colors_90)
ggsave(file.path(output_dir, 'Bcell_subcelltype_density.png'),width = 14,height =12,dpi = 300)

#########################################
#测试一下不同的k对B细胞聚类的影响
if (FALSE) {
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

}


#T细胞
Tcell.obj = qs_read(file.path(output_dir, 'Tcell_obj_intergration.qs'))
#Tcell.obj = obj.list[[3]]
DimPlot(Tcell.obj,reduction = 'umap',group.by = 'RNA_snn_res.0.8',cols = colors_90,raster = F,label=T,label.size=5) + NoLegend() + ggtitle('Tcell')
Idents(Tcell.obj)='RNA_snn_res.0.8'
DotPlot(Tcell.obj,
  features = c('CD2','CD3D','CD3E','CD4','CD8A','CD8B','TRDC','TRAC'),
  cols = 'RdBu') +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
markers.T = cosg(Tcell.obj,
                groups='all',
                assay = 'RNA',
                mu=100,
                remove_lowly_expressed = TRUE,
                expressed_pct = 0.3,
                n_genes_user=50)
View(markers.T$names)
DotPlot(Tcell.obj,
  features = markers.T$names[1:15,'17'],
  cols = 'RdBu') +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(Tcell.obj,
  features = c('FCGR3A','NCAM1'),
  cols = 'RdBu',idents = c(1,2,3,13,14,15,16)) +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
Tcell.obj = subset(Tcell.obj,idents = c(10,12,15),invert=T)

Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(1))) <- "t01_CD8_KLRG1"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(2))) <- "t02_CD8_CCR5 (Teff)"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(3))) <- "t03_CD8_LAG3 (Tex)"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(9))) <- "t04_CD8_IFIT1"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(16))) <- "t05_CD8_MKI67"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(4))) <- "t06_CD4_IL7R (Tn)"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(7))) <- "t07_CD4_S1PR1 (Tcm)"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(8))) <- "t08_CD4_CXCL13 (Tfh)"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(5))) <- "t09_CD4_CCR4 (Treg)"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(6))) <- "t10_CD4_ICOS (Treg)"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(11))) <- "t11_CD4_RORA (Th)"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(13))) <- "t12_NKT_FCGR3A"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(14))) <- "t13_NK_NCAM1"
Idents(Tcell.obj,cells = WhichCells(Tcell.obj, idents = c(17))) <- "t14_dnT_HAVCR2"
Tcell.obj$subcelltype = Idents(Tcell.obj)
Tcell.obj$Group = paste0(Tcell.obj$Time_point,'_',Tcell.obj$Outcome)
qs_save(Tcell.obj,file.path(output_dir, 'Tcell_obj_intergration.qs'))
# >>> COMMENTED FIGURE BLOCK: Fig2.R / Fig2A: T-cell subtype UMAP and density exports [original lines 751-759]
# DimPlot(Tcell.obj,reduction = 'umap',cols = colors_90,raster = F,alpha=0.7) + ggtitle('Tcell') + theme(aspect.ratio = 1)
# ggsave(file.path(output_dir, 'Tcell_subcelltype_dimplot.pdf'),width = 8,height =6)
# plot_density_dimplot_new(seurat_object = Tcell.obj,
#                     group.by = 'subcelltype',
#                     split.by = 'Group',
#                     ncol = 2,
#                     title = 'Tcell',
#                     color_palette = 'magma',group_colors = colors_90)
# ggsave(file.path(output_dir, 'Tcell_subcelltype_density.png'),width = 14,height =12,dpi = 300)
# <<< END COMMENTED FIGURE BLOCK



#########################################
#测试一下不同的k对T细胞聚类的影响
if (FALSE) {
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
}

#################################
#Myeloid
Myeloid.obj = qs_read(file.path(output_dir, 'Myeloid_obj_intergration.qs'))
Myeloid.obj.raw = Myeloid.obj
DimPlot(Myeloid.obj,reduction = 'umap',group.by = 'RNA_snn_res.1.2',cols = colors_90,raster = F,label=T,label.size=5) + NoLegend() + ggtitle('Myeloid')
Idents(Myeloid.obj)='RNA_snn_res.1.2'
DotPlot(Myeloid.obj,
  features = c('LYZ','CD68','CSF3R','CPA3','LILRA4','CLEC9A','CD1C','CCR7','FCGR3A','CD14','FCGR3B'),
  cols = 'RdBu') +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
VlnPlot(Myeloid.obj,
  features = c('LYZ','CD68','CSF3R','CPA3','LILRA4','CLEC9A','CD1C','CCR7','FCGR3A','CD14','FCGR3B'),
  cols = colors_90,
  fill.by = 'ident',
  stack=T,flip=T)
markers.Myeloid = cosg(Myeloid.obj,
                        groups='all',
                        assay = 'RNA',
                        mu=50,
                        remove_lowly_expressed = FALSE,
                        expressed_pct = 0.2,
                        n_genes_user=50)
View(markers.Myeloid$names)
DotPlot(Myeloid.obj,
  features = markers.Myeloid$names[1:15,'11'],
  cols = 'RdBu') +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
Myeloid.obj = subset(Myeloid.obj,idents = c(5,7,8,12,14,18,22,23,24),invert=T)
Myeloid.obj = subcelltype_analysis(Myeloid.obj,output_dir,'Myeloid')
Myeloid.obj = subset(Myeloid.obj,idents = c(22,23),invert=T)
Myeloid.obj = subcelltype_analysis(Myeloid.obj,output_dir,'Myeloid')
Myeloid.obj = subset(Myeloid.obj,idents = c(19),invert=T)
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(1)) ) <- "m01_Macro_APOE"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(2)) ) <- "m02_Macro_TNFRSF25"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(5)) ) <- "m03_Macro_SPP1"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(9)) ) <- "m04_Macro_VASH1"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(10)) ) <- "m05_Macro_HTRA4"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(13)) ) <- "m06_Macro_FOLR2"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(17)) ) <- "m07_Macro_MKI67"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(18)) ) <- "m08_Macro_MMP9"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(16)) ) <- "m09_Macro_ISG"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(14)) ) <- "m10_Macro_NRG1"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(6)) ) <- "m11_Mono_VCAN"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(15)) ) <- "m12_cDC1_CLEC9A"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(7)) ) <- "m13_cDC2_CD1C"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(8)) ) <- "m14_cDC3_CCR7"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(12)) ) <- "m15_pDC_LILRA4"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(4)) ) <- "m16_Nutri_CXCR1"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(11)) ) <- "m17_Nutri_ISG"
Idents(Myeloid.obj,cells = WhichCells(Myeloid.obj, idents = c(3)) ) <- "m18_MAST_CPA3"
Myeloid.obj$subcelltype = Idents(Myeloid.obj)
Myeloid.obj$Group = paste0(Myeloid.obj$Time_point,'_',Myeloid.obj$Outcome)
qs_save(Myeloid.obj,file.path(output_dir, 'Myeloid_obj_intergration.qs'))
DimPlot(Myeloid.obj,reduction = 'umap',cols = colors_90,raster = F,alpha=0.7) + ggtitle('Myeloid') + theme(aspect.ratio = 1)
ggsave(file.path(output_dir, 'Myeloid_subcelltype_dimplot.pdf'),width = 8,height =6)
plot_density_dimplot_new(seurat_object = Myeloid.obj,
                    group.by = 'subcelltype',
                    split.by = 'Group',
                    ncol = 2,
                    title = 'Myeloid',
                    color_palette = 'magma',group_colors = colors_90)
ggsave(file.path(output_dir, 'Myeloid_subcelltype_density.png'),width = 14,height =12,dpi = 300)


#########################################
#测试一下不同的k对B细胞聚类的影响
if (FALSE) {
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

}


######################################
Epi.obj = qs_read(file.path(output_dir, 'Epi_obj_intergration.qs'))
DimPlot(Epi.obj,reduction = 'umap',group.by = 'RNA_snn_res.0.2',cols = colors_90,raster = F,label=T,label.size=5) + NoLegend() + ggtitle('Epi')
DimPlot(Epi.obj,reduction = 'umap',group.by = 'sample',cols = colors_90,raster = F,alpha = 0.2) + NoLegend() + ggtitle('Sample distribution in Epi cluster')
ggsave(file.path(output_dir, 'Epi_sample_distribution.pdf'),width = 6,height =6)
# Idents(Epi.obj)='RNA_snn_res.0.2'
# Epi.markers = cosg(Epi.obj,
#                   groups='all',
#                   assay = 'RNA',
#                   mu=100,
#                   remove_lowly_expressed = TRUE,
#                   expressed_pct = 0.3,
#                   n_genes_user=50)
# View(Epi.markers$names)
# DotPlot(Epi.obj,
#   features = Epi.markers$names[1:15,'6'],
#   cols = 'RdBu') +
# theme(axis.text.x = element_text(angle = 45, hjust = 1))
# #Epi.obj = subset(Epi.obj,idents = c(28,32,37,59),invert=T)
# Idents(Epi.obj,cells = WhichCells(Epi.obj, idents = c(1)) ) <- "e01_Epi_PNCK"
# Idents(Epi.obj,cells = WhichCells(Epi.obj, idents = c(2)) ) <- "e02_Epi_TMSB10"
# Idents(Epi.obj,cells = WhichCells(Epi.obj, idents = c(3)) ) <- "e03_Epi_MKI67"
# Idents(Epi.obj,cells = WhichCells(Epi.obj, idents = c(4)) ) <- "e04_Epi_UPK3A"
# Idents(Epi.obj,cells = WhichCells(Epi.obj, idents = c(5)) ) <- "e05_Epi_KRT5"
# Idents(Epi.obj,cells = WhichCells(Epi.obj, idents = c(6)) ) <- "e06_Epi_PIGR"
# Idents(Epi.obj,cells = WhichCells(Epi.obj, idents = c(7)) ) <- "e07_Epi_IFI"
# Idents(Epi.obj,cells = WhichCells(Epi.obj, idents = c(8)) ) <- "e08_Epi_PTHLH"
Epi.obj$subcelltype = Epi.obj$celltype
Epi.obj$Group = paste0(Epi.obj$Time_point,'_',Epi.obj$Outcome)
qs_save(Epi.obj,file.path(output_dir, 'Epi_obj_intergration.qs'))




##################################
Stromal.obj = qs_read(file.path(output_dir, 'Stromal_obj_intergration.qs'))
DimPlot(Stromal.obj,reduction = 'umap',group.by = 'RNA_snn_res.0.9',cols = colors_90,raster = F,label=T,label.size=5) + NoLegend() + ggtitle('Stromal')
DimPlot(Stromal.obj,reduction = 'umap',group.by = 'sample',cols = colors_90,raster = F,alpha = 0.5) + NoLegend() + ggtitle('Stromal')

Idents(Stromal.obj)='RNA_snn_res.0.9'
DotPlot(Stromal.obj,
  features = c('COL1A1','COL3A1','DCN','LUM','VWF','PECAM1','CDH5','ACTA2','TAGLN','MYH11','RGS5','PDGFRB'),
  cols = 'RdBu') +
theme(axis.text.x = element_text(angle = 45, hjust = 1))

markers.Stromal = cosg(Stromal.obj,
                      groups='all',
                      assay = 'RNA',
                      mu=100,
                      remove_lowly_expressed = TRUE,
                      expressed_pct = 0.3,
                      n_genes_user=50)
View(markers.Stromal$names)
DotPlot(Stromal.obj,
  features = markers.Stromal$names[1:15,'18'],
  cols = 'RdBu') +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
Stromal.obj = subset(Stromal.obj,idents = c(6,9,12,14,18,20),invert=T)
Stromal.obj = subcelltype_analysis(Stromal.obj,output_dir,'Stromal')


Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(1))) <- "s01_Fib_APOE"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(2))) <- "s02_Fib_CXCL14"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(3))) <- "s03_Fib_SFRP1"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(5))) <- "s04_Fib_PLAGL1"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(8))) <- "s05_Fib_POSTN"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(15))) <- "s06_Fib_CCL19"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(17))) <- "s07_Fib_MMP3"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(4))) <- "s08_SMC_SPEG"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(13))) <- "s09_SMC_PCP4"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(7))) <- "s10_Per_GJC1"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(14))) <- "s11_Per_RERGL"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(6))) <- "s12_Endo_HECW2"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(10))) <- "s13_Endo_ACKR1"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(11))) <- "s14_Endo_SELP"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(12))) <- "s15_Endo_SEMA3G"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(16))) <- "s16_Endo_PROX1"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(9))) <- "s17_Endo_MKI67"
Idents(Stromal.obj, cells = WhichCells(Stromal.obj, idents = c(18))) <- "s18_Schwann_SOX10"
Stromal.obj$subcelltype = Idents(Stromal.obj)
Stromal.obj$Group = paste0(Stromal.obj$Time_point,'_',Stromal.obj$Outcome)
qs_save(Stromal.obj,file.path(output_dir, 'Stromal_obj_intergration.qs'))
DimPlot(Stromal.obj,reduction = 'umap',cols = colors_90,raster = F,alpha=0.7) + ggtitle('Stromal') + theme(aspect.ratio = 1)
ggsave(file.path(output_dir, 'Stromal_subcelltype_dimplot.pdf'),width = 8,height =6)
plot_density_dimplot_new(seurat_object = Stromal.obj,
                    group.by = 'subcelltype',
                    split.by = 'Group',
                    ncol = 2,
                    title = 'Stromal',
                    color_palette = 'magma',group_colors = colors_90)
ggsave(file.path(output_dir, 'Stromal_subcelltype_density.png'),width = 14,height =12,dpi = 300)


if (FALSE) {
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

}

######
#合并Epi和Stromal，Myeloid，T以及B seurat object
Myeloid.obj = qs_read(file.path(output_dir, 'Myeloid_obj_intergration_new.qs'))
Epi.obj = qs_read(file.path(output_dir, 'Epi_obj_intergration.qs'))
Stromal.obj = qs_read(file.path(output_dir, 'Stromal_obj_intergration_new.qs'))
Tcell.obj = qs_read(file.path(output_dir, 'Tcell_obj_intergration_new.qs'))
Bcell.obj = qs_read(file.path(output_dir, 'Bcell_obj_intergration_new.qs'))

seu.combined.IT = merge(
    x = Myeloid.obj,
    y = list(Epi.obj, Stromal.obj, Tcell.obj, Bcell.obj),
    add.cell.ids = NULL,   # 因为我们已经手动加 prefix 了
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
#构建pyscenic输入
#Myeloid
dir.create(file.path(output_dir, 'pyscenic'))
dir.create(file.path(output_dir, 'pyscenic','input'))
getmatrix <- function(obj,output) {
  counts = as.matrix(obj@assays$RNA$counts)
  colnames(counts) = colnames(obj)
  rownames(counts) = rownames(obj)
  write.csv(counts,file = output)
}
malignant.obj = qs_read(file.path(output_dir, 'malignant_obj_intergration.qs'))
Idents(malignant.obj) = 'sample'
set.seed(123)
malignant.sub <- subset(x = malignant.obj, downsample = 1500)
getmatrix(Myeloid.obj,file.path(output_dir, 'pyscenic','input','Myeloid_counts.csv'))
getmatrix(Stromal.obj,file.path(output_dir, 'pyscenic','input','Stromal_counts.csv'))
getmatrix(Tcell.obj,file.path(output_dir, 'pyscenic','input','Tcell_counts.csv'))
getmatrix(Bcell.obj,file.path(output_dir, 'pyscenic','input','Bcell_counts.csv'))
getmatrix(malignant.sub,file.path(output_dir, 'pyscenic','input','malignant_counts.csv'))

##################
#加上copykat的结果
copykat.prediction = qs_read(file.path(output_dir, 'copykat_output','copykat_prediction.qs'))
View(copykat.prediction)
rownames(copykat.prediction) = copykat.prediction$cell.names
#判断细胞id是否相同
length(setdiff(colnames(seu.combined.IT),rownames(copykat.prediction)))
seu.combined.IT$copykat_prediction = copykat.prediction[colnames(seu.combined.IT),'copykat.pred']
Epi.obj$copykat_prediction = copykat.prediction[colnames(Epi.obj),'copykat.pred']
malignant.obj = subset(Epi.obj,subset = copykat_prediction == 'aneuploid')
qs_save(malignant.obj,file.path(output_dir, 'malignant_obj_intergration.qs'))

#nmf分析肿瘤细胞保守MP
malignant.obj = qs_read(file.path(output_dir, 'malignant_obj_intergration.qs'))
nmf_output = qs_read(file.path(output_dir, 'NMF_output', 'nmf_results_all_ranks_multi.qs'))
obj = malignant.obj
meta.raw = obj@meta.data
meta.raw = meta.raw[,!grepl('^RNA_snn_res',colnames(meta.raw))]
obj = CreateSeuratObject(counts = obj@assays$RNA$counts,project = "malignant",min.cells = 10,min.features = 200)
meta = obj@meta.data
add.cols = setdiff(colnames(meta.raw),colnames(meta))
meta[add.cols] = meta.raw[row.names(meta),add.cols]
obj@meta.data = meta
obj$orig.ident = obj$sample
#obj@meta.data = obj@meta.data[,!grepl('^RNA_snn_res',colnames(obj@meta.data))]
obj[['RNA']] = split(obj[['RNA']],f=obj$sample)
obj = obj %>%
  NormalizeData(.,verbose = T) %>%
  FindVariableFeatures(.,selection.method = 'vst',nfeatures = 2500)
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
ElbowPlot(obj,ndims = 50)
npcs = 30
obj = obj %>%
    RunUMAP(.,reduction = "pca",dims = 1:npcs,reduction.name = 'umap.raw') %>%
    IntegrateLayers(.,method = HarmonyIntegration,
                    new.reduction = "harmony",orig.reduction = "pca",max.iter.harmony = 20L,verbose=T) %>%
    RunUMAP(.,reduction = "harmony",dims = 1:npcs,reduction.name = 'umap.harmony') %>%
    FindNeighbors(., reduction = "harmony", dims = 1:npcs,verbose = F) %>%
    FindClusters(., resolution = 0.5,algorithm=4,verbose = F,random.seed=1234)
obj = JoinLayers(obj)
malignant.obj = obj
# >>> COMMENTED FIGURE BLOCK: Fig4.R / FigS4A: malignant UMAP exports [original lines 1471-1476]
# DimPlot(malignant.obj,reduction = 'umap.harmony',group.by = 'sample',cols = colors_90,raster = F,alpha=0.4) + ggtitle('Malignant harmony') + theme(aspect.ratio = 1) + NoLegend()
# ggsave(file.path(output_dir, 'Malignant_harmony_dimplot.pdf'),width = 8,height =6)
# DimPlot(malignant.obj,reduction = 'umap.raw',group.by = 'sample',cols = colors_90,raster = F,alpha=0.4) + ggtitle('Malignant without intergration') + theme(aspect.ratio = 1) + NoLegend()
# ggsave(file.path(output_dir, 'Malignant_without_intergration_dimplot.pdf'),width = 8,height =6)
# DimPlot(malignant.obj,reduction = 'umap.raw',group.by = 'Sample_ID',cols = colors_90,raster = F,alpha=0.4) + ggtitle('Malignant Sample_ID') + theme(aspect.ratio = 1) + NoLegend()
# ggsave(file.path(output_dir, 'Malignant_without_intergration_Sample_ID_dimplot.pdf'),width = 8,height =6)
# <<< END COMMENTED FIGURE BLOCK

## --- 提取每个 Rank 的 Top 50 基因程序 ---


# 定义提取函数
extract_top_genes <- function(w_matrix, n_top = 50) {
  # w_matrix 的行是基因，列是 Program
  top_genes <- apply(w_matrix, 2, function(column) {
    # 根据贡献度得分排序，取前 n_top 个基因名
    names(sort(column, decreasing = TRUE))[1:n_top]
  })
  return(top_genes) # 返回一个矩阵：Rows=50, Cols=Rank数量
}

# 转换结构
nmf_programs <- lapply(names(nmf_output), function(sample_id) {
  message("Processing sample: ", sample_id)
  
  sample_data <- nmf_output[[sample_id]]
  
  # 获取该样本下所有运行过的 Rank (如 rank_4, rank_5...)
  all_ranks_top <- lapply(names(sample_data), function(r_name) {
    w_mat <- sample_data[[r_name]]$w
    
    # 提取该 Rank 下的所有 Programs 的 Top 50 基因
    top_mat <- extract_top_genes(w_mat, n_top = 50)
    
    # 为了区分，给程序列名加上 Rank 前缀，例如 "R4_1", "R4_2" ...
    colnames(top_mat) <- paste0(gsub("rank_", "R", r_name), "_", 1:ncol(top_mat))
    return(top_mat)
  })
  
  # 将该样本不同 Rank 的结果横向合并成一个大矩阵
  # 最终矩阵：行是 1:50, 列是该样本所有的 Program (4+5+6+7+8+9 = 39 列)
  sample_combined_mat <- do.call(cbind, all_ranks_top)
  
  return(sample_combined_mat)
})

# 保持列表名称与样本名一致
names(nmf_programs) <- names(nmf_output)

# 检查一下其中一个样本的结果维度
# 预期应该是 50 行，列数等于 (4+5+6+7+8+9) = 39 列
View(nmf_programs[[1]])

intra_stats <- lapply(names(nmf_programs), function(x) {
  z <- nmf_programs[[x]]
  
  # 计算样本内部 Program 之间的基因交集数量
  intra_intersect <- apply(z, 2, function(col_i) {
    apply(z, 2, function(col_j) length(intersect(col_i, col_j)))
  })
  
  # 找到每个 Program 与其他 Program (非自身) 的最大交集数
  # sort(y, decreasing = T)[2] 是为了跳过与自身的重叠（自身重叠永远是 50）
  max_overlaps <- apply(intra_intersect, 2, function(y) sort(y, decreasing = T)[2])
  
  # 统计达到 35 个基因重合的程序数量
  n_pass <- sum(max_overlaps >= 35)
  
  return(data.frame(
    sample = x, 
    total_programs = ncol(z), 
    passed_intra_35 = n_pass,
    avg_max_overlap = mean(max_overlaps) # 顺便看看平均重合水平
  ))
})

intra_summary <- do.call(rbind, intra_stats)
print(intra_summary)

# 检查总存活数
if(sum(intra_summary$passed_intra_35) == 0) {
  message("警告：没有样本能通过 intra_min = 35 的筛选，建议下调该阈值。")
}

nmf_programs_renamed <- lapply(names(nmf_programs), function(s_name) {
  mat <- nmf_programs[[s_name]]
  # 将列名从 "R4_1" 改为 "SampleName_R4_1"
  colnames(mat) <- paste0(s_name, "_", colnames(mat))
  return(mat)
})
names(nmf_programs_renamed) <- names(nmf_programs)
robust_nmf_programs <- function(nmf_programs, intra_min = 35, intra_max = 10, inter_filter = TRUE, inter_min = 10) {
    
    # 1. Intra-sample selection
    intra_intersect <- lapply(nmf_programs, function(z) {
        apply(z, 2, function(x) apply(z, 2, function(y) length(intersect(x, y))))
    })
    
    intra_intersect_max <- lapply(intra_intersect, function(x) {
        apply(x, 2, function(y) sort(y, decreasing = TRUE)[2])
    })
    
    # 筛选满足内部重叠度的程序
    nmf_sel <- lapply(names(nmf_programs), function(x) {
        mat <- nmf_programs[[x]][, intra_intersect_max[[x]] >= intra_min, drop = FALSE]
        return(mat)
    })
    names(nmf_sel) <- names(nmf_programs)
    
    # 2. Inter-sample selection
    nmf_sel_unlist <- do.call(cbind, nmf_sel)
    if (is.null(nmf_sel_unlist) || ncol(nmf_sel_unlist) == 0) return(NULL)
    
    inter_intersect <- apply(nmf_sel_unlist, 2, function(x) {
        apply(nmf_sel_unlist, 2, function(y) length(intersect(x, y)))
    })
    
    final_filter <- NULL 
    
    for(i in names(nmf_sel)) {
        # 匹配该样本的所有列 (现在列名已经是 Sample_R4_1 格式)
        current_sample_cols <- grep(paste0("^", i, "_"), colnames(inter_intersect))
        
        if(length(current_sample_cols) == 0) next
        
        # 提取 其他样本 vs 当前样本 的重合度矩阵
        a <- inter_intersect[-current_sample_cols, current_sample_cols, drop = FALSE]
        
        if (nrow(a) == 0) {
            # 如果只有一个样本，则无法进行 inter-filter，直接保留当前样本内最稳定的
            b <- rep(50, length(current_sample_cols))
            names(b) <- colnames(inter_intersect)[current_sample_cols]
        } else {
            # 获取每个程序在其他样本中的最大重叠
            b <- sort(apply(a, 2, max), decreasing = TRUE)
        }
        
        # 跨样本过滤：必须在至少一个其他样本中找到相似程序 (重叠 >= 10)
        if(inter_filter == TRUE) {
            b <- b[b >= inter_min] 
        }
        
        # 3. 冗余剔除
        if(length(b) >= 1) {
            c <- names(b)[1]
            if(length(b) > 1) {
                for(y in 2:length(b)) {
                    # 如果待选程序与已选程序的重叠 <= 10，则认为是一个独立的 Program
                    if(max(inter_intersect[c, names(b[y])]) <= intra_max) {
                        c <- c(c, names(b[y]))
                    }
                }
            }
            final_filter <- c(final_filter, c)
        }
    }
    return(final_filter)
}

# 运行
final_program_names <- robust_nmf_programs(nmf_programs_renamed)
print(final_program_names)

# 1. 提取所有选定程序的 Top 50 基因矩阵
# nmf_programs_renamed 是之前带样本名前缀的 list
all_programs_mat <- do.call(cbind, nmf_programs_renamed)
robust_programs_mat <- all_programs_mat[, final_program_names]

# 2. 构建注释信息 (Annotation)
# 假设你有一个样本元数据表格叫做 sample_metadata，包含 'sample' 和 'group' 列
# 如果没有，请手动创建一个
anno_df <- data.frame(
  Program = colnames(robust_programs_mat),
  Sample = gsub("_R[0-9].*", "", colnames(robust_programs_mat)) # 提取样本名
)
# 匹配分组信息 (请确保你的 metadata 中有对应的 sample 名)
metadata = malignant.obj@meta.data
anno_df$Group <- metadata$Group[match(anno_df$Sample, metadata$sample)]

# 3. 计算 Jaccard 相似度矩阵
calc_jaccard_mat <- function(mat) {
  n <- ncol(mat)
  j_mat <- matrix(0, n, n)
  colnames(j_mat) <- rownames(j_mat) <- colnames(mat)
  for(i in 1:n) {
    for(j in i:n) {
      it <- length(intersect(mat[,i], mat[,j]))
      un <- length(union(mat[,i], mat[,j]))
      j_mat[i,j] <- j_mat[j,i] <- it/un
    }
  }
  return(j_mat)
}
jaccard_mat <- calc_jaccard_mat(robust_programs_mat)

library(ComplexHeatmap)
library(circlize)

# 设置颜色
# 假设 colors_90 已经存在，我们提取对应样本数的颜色
unique_samples <- unique(anno_df$Sample)
sample_colors <- setNames(colors_90[1:length(unique_samples)], unique_samples)
group_colors <- c("Pre_NR" = "#A52A2A", "Pre_CR" = "#053E7A","Post_NR" = "#E46B4F","Post_CR" = "#BFBF73") # 根据实际修改

# 定义热图顶部注释
column_ha = HeatmapAnnotation(
  Sample = anno_df$Sample,
  Group = anno_df$Group,
  col = list(Sample = sample_colors,Group = group_colors), 
  show_annotation_name = TRUE
)

dist_obj <- as.dist(1 - jaccard_mat)
hc <- hclust(dist_obj, method = "average")
# 绘图
pdf(file.path(output_dir, "Figures","NMF_Programs_Jaccard_Heatmap.pdf"), width = 12, height = 10)
hp <- Heatmap(jaccard_mat, 
        name = "Jaccard",
        column_title = "NMF Programs Clustering (Aligned Diagonal)",
        top_annotation = column_ha,
        
        cluster_rows = hc,         # 强制行聚类顺序
        cluster_columns = hc,      # 强制列聚类顺序
        
        show_row_names = FALSE, 
        show_column_names = FALSE,
        col = colorRamp2(c(0.05, 0.1, 0.2), c("white", "#fddbc7", "#b2182b")), # 调整颜色映射
        border = TRUE)
draw(hp)
dev.off()

ordered_indices <- hc$order

# 2. 生成有序的 Jaccard 矩阵
jaccard_mat_ordered <- jaccard_mat[ordered_indices, ordered_indices]

# 3. 同时也给你的基因程序矩阵排序（方便后面按 ID 提取基因）
robust_programs_mat_ordered <- robust_programs_mat[, ordered_indices]

# 4. 保存一个带 ID 的基因列表到本地，方便你后续对照
# 第一行是序号，方便你手动记下 MP 的范围
robust_programs_ordered_with_id <- rbind(ID = 1:ncol(robust_programs_mat_ordered), 
                                        robust_programs_mat_ordered)
write.csv(robust_programs_ordered_with_id, file.path(output_dir,"Robust_Programs_Ordered_with_ID.csv"))

diag_ids <- 1:ncol(jaccard_mat_ordered)
pdf(file.path(output_dir,"Figures","Jaccard_Heatmap_With_Diagonal_IDs.pdf"), width = 14, height = 13)
ht <- Heatmap(jaccard_mat_ordered, 
        name = "Jaccard",
        column_title = "NMF Programs: IDs on Diagonal (Zoom in to see)",
        
        # 关闭内部自动聚类，使用我们已经排好的顺序
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        
        # 如果你之前算过分群，依然可以加上切分线
        #column_split = cluster_assignments[ordered_indices], 
        #row_split = cluster_assignments[ordered_indices],
        #gap = unit(2, "mm"),
        
        # 核心：在格子内画 ID
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(i == j) {
            grid.text(diag_ids[i], x, y, 
                      gp = gpar(fontsize = 5, col = "black", fontface = "bold"))
          }
        },
        
        # 视觉美化
        show_row_names = FALSE, 
        show_column_names = FALSE,
        border = TRUE,
        col = colorRamp2(c(0.05, 0.1, 0.2), c("white", "#fddbc7", "#b2182b")),
        
        # 顶部注释（注意这里也需要按顺序排）
        top_annotation = HeatmapAnnotation(
          Sample = anno_df$Sample[match(colnames(jaccard_mat_ordered), anno_df$Program)],
          col = list(Sample = sample_colors),
          show_annotation_name = TRUE
        ))
draw(ht)
dev.off()

View(robust_programs_mat_with_id)
MPs = list(MP1 = c(1:23),
          MP2 = c(24:34),
          MP3 = c(37:64),
          MP4 = c(69:92),
          MP5 = c(96:101),
          MP6 = c(111:120),
          MP7 = c(121:134),
          MP8 = c(135:155),
          MP9 = c(169:187),
          MP10 = c(188:194))

# 1. 提取每个 MP 的候选基因 (频率 >= 25%)
# 使用之前排序后的矩阵 robust_programs_mat_ordered
mps_candidate_genes <- lapply(MPs, function(range_idx) {
  # 提取该范围内的程序
  selected_progs <- colnames(robust_programs_mat_ordered)[range_idx]
  member_genes_mat <- robust_programs_mat_ordered[, selected_progs, drop=FALSE]
  
  # 统计基因出现频率
  gene_counts <- table(member_genes_mat)
  # 筛选频率 >= 25% 的基因
  pass_genes <- names(gene_counts[gene_counts >= (length(selected_progs) * 0.25)])
  return(pass_genes)
})

# 2. 在 Seurat 对象中计算 MP 初始评分
# 假设 malignant.obj 是你的单细胞对象
message("Calculating initial Module Scores...")
malignant.obj <- AddModuleScore(
  object = malignant.obj,
  features = mps_candidate_genes,
  name = "MP_Score_", # 结果会存为 MP_Score_1, MP_Score_2 ...
  nbin = 24, 
  seed = 123
)

# 3. 计算相关性并提取最终 Top 30
# 提取表达矩阵 (建议使用 scale.data 或 data 层)
# 为了计算速度和准确性，通常针对 VariableFeatures 或全基因
final_mp_signatures_list <- lapply(seq_along(MPs), function(i) {
  mp_name <- names(MPs)[i]
  score_col <- paste0("MP_Score_", i)
  candidates <- mps_candidate_genes[[i]] # 取得該 MP 的候選基因集
  
  message(paste0("正在精煉 ", mp_name, "，候選基因數：", length(candidates)))
  
  # 取得評分向量
  mp_score_vector <- malignant.obj@meta.data[[score_col]]
  
  # 僅提取候選基因的表達矩陣
  # 注意：確保候選基因都在表達矩陣中
  valid_candidates <- intersect(candidates, rownames(malignant.obj))
  subset_expr <- GetAssayData(malignant.obj, assay = "RNA", layer = "data")[valid_candidates, ]
  
  # 計算候選基因與總分的 Pearson 相關性
  # cor 需要輸入 矩陣(cells x genes) 和 向量(cells)
  cor_values <- cor(t(as.matrix(subset_expr)), mp_score_vector)
  
  # 從候選基因中取出相關性最高的前 30 個
  top_30 <- names(sort(cor_values[,1], decreasing = TRUE)[1:30])
  
  return(top_30)
})

names(final_mp_signatures_list) <- names(MPs)

# 4. 保存结果
qs_save(final_mp_signatures_list, file.path(output_dir, "Final_MP_Top30_Signatures.qs"))
#final_mp_signatures_list = qs_read(file.path(output_dir, "Final_MP_Top30_Signatures.qs"))
# 1. 创建 MP 标注向量
# 初始化为 "None"
mp_annotation <- rep("None", ncol(jaccard_mat_ordered))

# 根据你提供的 MPs 列表填充名称
for(mp_name in names(MPs)) {
  indices <- MPs[[mp_name]]
  mp_annotation[indices] <- mp_name
}

# 2. 构建注释 Data Frame
# 确保 anno_df 的顺序已经通过之前的 match 操作与 jaccard_mat_ordered 对齐
anno_final <- data.frame(
  MP = mp_annotation,
  Patient = anno_df$Sample[match(colnames(jaccard_mat_ordered), anno_df$Program)],
  Group = anno_df$Group[match(colnames(jaccard_mat_ordered), anno_df$Program)] # 确保 metadata 中有 Group
)

# 1. MP 颜色：10 个 MP 用鲜艳颜色，None 用白色
u_mps <- names(MPs)
mp_colors <- setNames(c(brewer.pal(10, "Paired"), "white"), c(u_mps, "None"))

# 2. Patient 颜色：使用你之前的 colors_90
u_patients <- unique(anno_final$Patient)
patient_colors <- setNames(colors_90[1:length(u_patients)], u_patients)

# 3. Group 颜色 (示例，请根据实际修改)
u_groups <- unique(anno_final$Group)
group_colors <- c("Pre_NR" = "#A52A2A", "Pre_CR" = "#053E7A","Post_NR" = "#E46B4F","Post_CR" = "#BFBF73")

col_ha <- HeatmapAnnotation(
  MP = anno_final$MP,
  Group = anno_final$Group,
  Patient = anno_final$Patient,
  col = list(
    MP = mp_colors,
    Patient = patient_colors,
    Group = group_colors
  ),
  # 去掉颜色块的边框
  border = FALSE, 
  # gp = gpar(lty = 0) 会去掉色块之间的缝隙/边线
  #gp = gpar(lty = 0), 
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 10)
)

# >>> COMMENTED FIGURE BLOCK: Fig4.R / Fig4A: annotated MP similarity heatmap export [original lines 1870-1910]
# pdf(file.path(output_dir,'Figures',"Final_MPs_Similarity_Heatmap_Annotated.pdf"), width = 14, height = 12)
# 
# ht <- Heatmap(jaccard_mat_ordered, 
#         name = "Jaccard Index",
#         column_title = "Consensus Meta-Programs (MP1-MP10)",
#         
#         # 基础设置
#         cluster_rows = FALSE, 
#         cluster_columns = FALSE,
#         show_row_names = FALSE, 
#         show_column_names = FALSE,
#         
#         # 顶部注释
#         top_annotation = col_ha,
#         
#         # 物理切分设置
#         column_split = factor(anno_final$MP, levels = c(names(MPs), "None")),
#         row_split = factor(anno_final$MP, levels = c(names(MPs), "None")),
#         cluster_column_slices = FALSE,
#         cluster_row_slices = FALSE,
#         gap = unit(1.5, "mm"), # 保持 Cluster 间的断开感
#         
#         # 颜色映射 (10个基因阈值 0.111)
#         col = colorRamp2(c(0.05, 0.1, 0.3), c("white", "#fddbc7", "#b2182b")),
#         
#         # 热图主体的格子边框（如果你想要热图格点也无边框，可以设为 NA）
#         rect_gp = gpar(col = NA), 
#         
#         # 对角线 ID 标注
#         # layer_fun = function(j, i, x, y, width, height, fill) {
#         #   v_index <- which(i == j)
#         #   if(length(v_index) > 0) {
#         #     grid.text(1:ncol(jaccard_mat_ordered), x[v_index], y[v_index], 
#         #               gp = gpar(fontsize = 3, col = "black"))
#         #   }
#         # }
# )
# 
# # 绘图时关闭图例的边框
# draw(ht, merge_legend = TRUE, annotation_legend_side = "right")
# dev.off()
# <<< END COMMENTED FIGURE BLOCK

websiteLive <- getOption("enrichR.live")
setEnrichrSite("Enrichr")
dbs <- c("GO_Biological_Process_2025","KEGG_2021_Human","MSigDB_Hallmark_2020")
bg_genes <- rownames(malignant.obj)
all_mp_enrichr <- lapply(names(final_mp_signatures_list), function(mp_name) {
  message(paste0("正在通过 enrichR 分析 ", mp_name, "..."))
  
  # 获取该 MP 的 Top 30 基因
  genes <- final_mp_signatures_list[[mp_name]]
  
  # 执行富集，加入 background 参数
  res <- enrichr(genes = genes, 
                databases = dbs, 
                background = bg_genes)
  
  return(res)
})
names(all_mp_enrichr) <- names(final_mp_signatures_list)

# 4. 保存结果
qs_save(all_mp_enrichr, "All_MP_enrichR_Results.qs")

pdf(file.path(output_dir,"Figures","MPs_EnrichR_Annotation_Plots.pdf"), width = 10, height = 7)

for (mp_name in names(all_mp_enrichr)) {
  
  # 提取 Hallmark 数据库结果
  plot_data <- all_mp_enrichr[[mp_name]]$MSigDB_Hallmark_2020
  
  if (!is.null(plot_data) && nrow(plot_data) > 0) {
    
    # 1. 数据预处理
    plot_df <- plot_data %>%
      # 过滤显著条目
      filter(Adjusted.P.value < 0.05) %>%
      # 稳健地计算基因数量：统计分号个数并加 1
      # 如果 Genes 为空或 NA，则设为 0
      mutate(
        Gene_Count = ifelse(is.na(Genes) | Genes == "", 
                            0, 
                            str_count(Genes, ";") + 1)
      ) %>%
      # 选取 Combined.Score 最高的前 10 个
      top_n(10, wt = Combined.Score) %>%
      # 按照得分重新排序，以便在图中从高到低排列
      arrange(Combined.Score) %>%
      mutate(Term = factor(Term, levels = unique(Term)))
    
    # 2. 检查是否有有效行后绘图
    if (nrow(plot_df) > 0) {
      p <- ggplot(plot_df, aes(x = Combined.Score, 
                              y = Term, 
                              color = Adjusted.P.value, 
                              size = Gene_Count)) +
        geom_point() +
        # 颜色映射：显著性高（P值小）为红色，显著性低为蓝色
        scale_color_gradient(low = "#E41A1C", high = "#377EB8") +
        scale_size_continuous(range = c(4, 10)) + 
        theme_bw() +
        theme(
          axis.text.y = element_text(size = 9, face = "bold"),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5)
        ) +
        labs(title = paste0(mp_name, " - Hallmark Enrichment Analysis"),
            subtitle = "Background: Malignant Cell Transcriptome",
            x = "Combined Score", 
            y = "",
            size = "Gene Count",
            color = "Adj. P-value")
      
      print(p)
    }
  }
  #KEGG_2021_Human
  plot_data <- all_mp_enrichr[[mp_name]]$KEGG_2021_Human
  
  if (!is.null(plot_data) && nrow(plot_data) > 0) {
    
    # 1. 数据预处理
    plot_df <- plot_data %>%
      # 过滤显著条目
      filter(Adjusted.P.value < 0.05) %>%
      # 稳健地计算基因数量：统计分号个数并加 1
      # 如果 Genes 为空或 NA，则设为 0
      mutate(
        Gene_Count = ifelse(is.na(Genes) | Genes == "", 
                            0, 
                            str_count(Genes, ";") + 1)
      ) %>%
      # 选取 Combined.Score 最高的前 10 个
      top_n(10, wt = Combined.Score) %>%
      # 按照得分重新排序，以便在图中从高到低排列
      arrange(Combined.Score) %>%
      mutate(Term = factor(Term, levels = unique(Term)))
    
    # 2. 检查是否有有效行后绘图
    if (nrow(plot_df) > 0) {
      p <- ggplot(plot_df, aes(x = Combined.Score, 
                              y = Term, 
                              color = Adjusted.P.value, 
                              size = Gene_Count)) +
        geom_point() +
        # 颜色映射：显著性高（P值小）为红色，显著性低为蓝色
        scale_color_gradient(low = "#E41A1C", high = "#377EB8") +
        scale_size_continuous(range = c(4, 10)) + 
        theme_bw() +
        theme(
          axis.text.y = element_text(size = 9, face = "bold"),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5)
        ) +
        labs(title = paste0(mp_name, " - KEGG Enrichment Analysis"),
        subtitle = "Background: Malignant Cell Transcriptome",
        x = "Combined Score", 
        y = "",
        size = "Gene Count",
        color = "Adj. P-value")
      print(p)
    }
  }
  #GO_Biological_Process_2025
  plot_data <- all_mp_enrichr[[mp_name]]$GO_Biological_Process_2025
  
  if (!is.null(plot_data) && nrow(plot_data) > 0) {
    
    # 1. 数据预处理
    plot_df <- plot_data %>%
      # 过滤显著条目
      filter(Adjusted.P.value < 0.05) %>%
      # 稳健地计算基因数量：统计分号个数并加 1
      # 如果 Genes 为空或 NA，则设为 0
      mutate(
        Gene_Count = ifelse(is.na(Genes) | Genes == "", 
                            0, 
                            str_count(Genes, ";") + 1)
      ) %>%
      # 选取 Combined.Score 最高的前 10 个
      top_n(10, wt = Combined.Score) %>%
      # 按照
      arrange(Combined.Score) %>%
      mutate(Term = factor(Term, levels = unique(Term)))
    
    if (nrow(plot_df) > 0) {
      p <- ggplot(plot_df, aes(x = Combined.Score, 
                              y = Term, 
                              color = Adjusted.P.value, 
                              size = Gene_Count)) +
        geom_point() +
        # 颜色映射：显著性高（P值小）为红色，显著性低为蓝色
        scale_color_gradient(low = "#E41A1C", high = "#377EB8") +
        scale_size_continuous(range = c(4, 10)) + 
        theme_bw() +
        theme(
          axis.text.y = element_text(size = 9, face = "bold"),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5)
        ) +
        labs(title = paste0(mp_name, " - GO Biological Process Enrichment Analysis"),
        subtitle = "Background: Malignant Cell Transcriptome",
        x = "Combined Score", 
        y = "",
        size = "Gene Count",
        color = "Adj. P-value")
      print(p)
    }
  }
}
dev.off()

names(final_mp_signatures_list) <- c(
  "Cycle_G2M", "Surface_Trafficking", "RNA_Splicing", "Basal_pEMT", 
  "Stress_IEG", "IFN_Response", "Hypoxia", "Metabolic", 
  "Uro_Differentiation", "Ep_Inflammation"
)
malignant.obj <- AddModuleScore(
  object = malignant.obj,
  features = final_mp_signatures_list,
  name = "MP_Final_", # 会生成 MP_Final_1 到 10
  nbin = 24, seed = 123
)
new_names <- names(final_mp_signatures_list)
new_names = c(
  "Cycle_G2M", "Surface_Trafficking", "RNA_Splicing", "Basal_pEMT", 
  "Stress_IEG", "IFN_Response", "Hypoxia", "Metabolic", 
  "Uro_Differentiation", "Ep_Inflammation"
)
for(i in 1:length(new_names)){
  colnames(malignant.obj@meta.data)[which(colnames(malignant.obj@meta.data) == paste0("MP_Final_", i))] <- new_names[i]
}

avg_mp_mat <- malignant.obj@meta.data %>%
  group_by(sample, Group) %>%
  summarise(across(all_of(new_names), mean)) %>%
  ungroup()

# 2. 提取 Group 信息并排序（确保热图按组排列）
avg_mp_mat <- avg_mp_mat %>% 
  arrange(factor(Group, levels = c("Pre_NR", "Pre_CR", "Post_NR", "Post_CR")))

# 3. 转为矩阵格式，供 Heatmap 使用
plot_mat <- as.matrix(avg_mp_mat[, new_names])
rownames(plot_mat) <- avg_mp_mat$sample
plot_mat <- t(plot_mat) # 转置为 MPs x Patients

# 4. 准备 Group 注释
group_info <- avg_mp_mat$Group
group_colors <- c("Pre_NR" = "#A52A2A", "Pre_CR" = "#053E7A", 
                  "Post_NR" = "#E46B4F", "Post_CR" = "#BFBF73")

                  # 定义顶部注释
top_anno <- HeatmapAnnotation(
  Group = group_info,
  col = list(Group = group_colors),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  # 去掉边框和缝隙
  gp = gpar(lty = 0),
  border = FALSE
)

# 绘图
pdf(file.path(output_dir,"Figures","MP_Activity_Cross_Groups_Final.pdf"), width = 10, height = 6)

ht <- Heatmap(plot_mat, 
        name = "Mean Score",
        column_title = "Meta-Program Activity across Clinical Groups",
        
        # 聚类设置
        cluster_rows = TRUE,          # MP 之间可以聚类看功能相似性
        cluster_columns = FALSE,      # 样本按 Group 顺序排列，不进行列聚类
        
        # 注释与分割
        top_annotation = top_anno,
        column_split = factor(group_info, levels = c("Pre_NR", "Pre_CR", "Post_NR", "Post_CR")),
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        
        # 隐藏横坐标样本名
        show_column_names = FALSE,
        show_row_names = TRUE,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 10),
        
        # 颜色映射 (Z-score 化处理让对比更明显)
        # 如果原始得分差异较小，可以使用 scale() 处理 plot_mat
        col = colorRamp2(c(min(plot_mat), 0, max(plot_mat)), c("#4575B4", "white", "#D73027")),
        
        # 视觉细节
        rect_gp = gpar(col = "white", lwd = 0.5), # 样本间极细的白线增加精致感
        border = TRUE,
        gap = unit(2, "mm") # 组与组之间的间隙
)

draw(ht, merge_legend = TRUE)
dev.off()


# 1. 准备统计数据
stat_df <- malignant.obj@meta.data %>%
  dplyr::select(sample, Group, all_of(new_names)) %>%
  tidyr::pivot_longer(cols = all_of(new_names), names_to = "MP", values_to = "Score")

# 2. 批量绘制小提琴图并计算 P 值
# 重点比较：Pre_CR vs Pre_NR (预后因素) 和 Pre_NR vs Post_NR (治疗演化)
comparisons_list <- list(c("Pre_CR", "Pre_NR"),c('Post_CR',"Post_NR"), c("Pre_NR", "Post_NR"))

# >>> COMMENTED FIGURE BLOCK: Fig4.R / Fig4B: MP global comparison violin export [original lines 2178-2185]
# pdf(file.path(output_dir,"Figures","MP_Global_Comparison_Violin.pdf"), width = 16, height = 10)
# ggviolin(stat_df, x = "Group", y = "Score", fill = "Group",
#         palette = group_colors, add = "boxplot", add.params = list(fill = "white")) +
#   facet_wrap(~MP, scales = "free_y", ncol = 5) +
#   stat_compare_means(comparisons = comparisons_list, label = "p.signif") +
#   theme_classic() + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# dev.off()
# <<< END COMMENTED FIGURE BLOCK

library(ggridges)

# >>> COMMENTED FIGURE BLOCK: Fig4.R / Fig4B: Basal-pEMT and Uro-differentiation ridge export [original lines 2189-2196]
# ggplot(stat_df %>% filter(MP %in% c("Basal_pEMT", "Uro_Differentiation")), 
#       aes(x = Score, y = Group, fill = Group)) +
#   geom_density_ridges(alpha = 0.7,quantile_lines = TRUE) +
#   scale_fill_manual(values = group_colors) +
#   facet_wrap(~MP, scales = "free") +
#   theme_ridges() +
#   labs(title = "Population-level Shift of MP Activity")
# ggsave(file.path(output_dir,"Figures",'Fig',"Fig4B_MP_Global_Comparison_Ridge_MP4_MP9.pdf"), width = 6, height = 2.5)
# <<< END COMMENTED FIGURE BLOCK

ggplot(stat_df, 
      aes(x = Score, y = Group, fill = Group)) +
  geom_density_ridges(alpha = 0.7) +
  scale_fill_manual(values = group_colors) +
  facet_wrap(~MP, scales = "free") +
  theme_ridges() +
  labs(title = "Population-level Shift of MP Activity")
ggsave(file.path(output_dir,"Figures","MP_Global_Comparison_Ridge.pdf"), width = 20, height = 20)

# 1. 提取 10 个 MP 的评分矩阵
mp_scores <- malignant.obj@meta.data[, new_names]

# 2. 定义分配函数
assign_mp <- function(row) {
  max_val <- max(row)
  if (max_val > 0) {
    return(names(row)[which.max(row)])
  } else {
    return("Unassigned")
  }
}

# 3. 应用到每一行并存回 Metadata
malignant.obj$MP_ID <- apply(mp_scores, 1, assign_mp)

# 看看分配情况
table(malignant.obj$MP_ID)


# 1. 计算比例矩阵
prop_df <- malignant.obj@meta.data %>%
  group_by(sample, Group, MP_ID) %>%
  summarise(n = n()) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup() %>%
  filter(MP_ID != "Unassigned") # 过滤掉未分配的，专注于功能亚群

# 2. 补齐缺失值（如果某个样本完全没有某个 MP，需补 0 以便统计）
prop_df <- prop_df %>%
  complete(nesting(sample,Group), MP_ID, fill = list(n = 0, proportion = 0))

# 3. 绘制样本层面的比例差异分析 (以 Patient 为统计单位)
# 这里我们比较 Pre_NR vs Post_NR
pdf(file.path(output_dir,"Figures","MP_Proportion_sample_Level_Stat.pdf"), width = 14, height = 10)

ggboxplot(prop_df, x = "Group", y = "proportion", fill = "Group",
          palette = group_colors, add = "jitter",
          facet.by = "MP_ID", scales = "free_y") +
  stat_compare_means(comparisons = list(c("Pre_CR", "Pre_NR")), 
                    label = "p.signif", method = "wilcox.test",hide.ns = TRUE,tip.length = 0.01) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Cell Proportion (per Patient)", title = "MP Composition Shift at Patient Level")

dev.off()



# 1. 确定排序基准：提取每个样本中 Basal_pEMT 的比例
prop_df <- malignant.obj@meta.data %>%
  group_by(sample, Group, MP_ID) %>%
  summarise(n = n(), .groups = "drop") %>% # 使用 .groups="drop" 自动取消分组
  group_by(sample) %>% # 重新按病人分组以计算每个病人的百分比
  mutate(proportion = n / sum(n)) %>% 
  ungroup()

# 2. 提取 Basal_pEMT 排序基准
basal_sort_order <- prop_df %>%
  dplyr::filter(MP_ID == "Basal_pEMT") %>%
  dplyr::select(sample, proportion) %>%
  dplyr::rename(basal_prop = proportion)

# 3. 合并排序基准并补齐缺失值
# 有些样本可能完全没有 Basal 细胞，left_join 后会产生 NA，需要补 0
prop_df_sorted <- prop_df %>%
  left_join(basal_sort_order, by = "sample") %>%
  mutate(basal_prop = replace_na(basal_prop, 0))

# 4. 执行组内排序逻辑
group_lvls <- c("Pre_NR", "Pre_CR", "Post_NR", "Post_CR")

prop_df_sorted <- prop_df_sorted %>%
  mutate(Group = factor(Group, levels = group_lvls)) %>%
  # 关键步骤：先按 Group 排序，组内按 basal_prop 降序排
  arrange(Group, desc(basal_prop)) %>%
  # 锁定 Patient 的因子水平，使其在绘图中保持此顺序
  mutate(sample = factor(sample, levels = unique(sample)))

# 定义 MP 的颜色（10个 MP 使用对比鲜明的配色，Unassigned 用浅灰色）
mp_palette <- c(
  "Basal_pEMT"          = "#E41A1C", # 红色突出显示
  "Uro_Differentiation" = "#377EB8", # 蓝色对比
  "Cycle_G2M"           = "#4DAF4A",
  "RNA_Splicing"        = "#984EA3",
  "Hypoxia"             = "#FF7F00",
  "IFN_Response"        = "#FFFF33",
  "Stress_IEG"          = "#A65628",
  "Metabolic" = "#F781BF",
  "Surface_Trafficking" = "#999999",
  "Ep_Inflammation"     = "#66C2A5"
)

pdf(file.path(output_dir,"Figures","MP_Composition_Sorted_by_Basal.pdf"), width = 12, height = 7)

ggplot(prop_df_sorted, aes(x = sample, y = proportion, fill = MP_ID)) +
  geom_bar(stat = "identity", width = 0.8) +
  # 分组分面，确保 Group 之间有物理隔离
  facet_grid(. ~ Group, scales = "free_x", space = "free", switch = "x") +
  scale_fill_manual(values = mp_palette) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),     # 隐藏密集的样本名
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.5, "lines"), # 组间间距
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Tumor MP Composition (Sorted by Basal_pEMT Content)",
    subtitle = "Patients within each group are ordered by descending proportion of Basal_pEMT cells",
    x = "Patients grouped by Clinical Response",
    y = "Relative Proportion of MP States",
    fill = "Meta-Program ID"
  )

dev.off()


# 1. 计算每个患者在 10 个 MP 上的平均得分
patient_avg_df <- malignant.obj@meta.data %>%
  group_by(sample, Group) %>%
  summarise(across(all_of(new_names), mean), .groups = "drop")

# 2. 将数据转换为长格式以便于批量绘图
patient_avg_long <- patient_avg_df %>%
  tidyr::pivot_longer(cols = all_of(new_names), names_to = "MP", values_to = "Average_Score")

# 3. 设定比较组（重点关注预后和治疗演化）
my_comparisons <- list(c("Pre_CR", "Pre_NR"), c("Pre_NR", "Post_NR"))
# 定义主图需要的 MP 名称
main_mps <- c("Basal_pEMT", "Uro_Differentiation")

# 自动获取剩下的 MP 名称作为附图
supp_mps <- setdiff(new_names, main_mps)
# 1. 过滤主图数据
df_main <- patient_avg_long %>% dplyr::filter(MP %in% main_mps)

# >>> COMMENTED FIGURE BLOCK: Fig4.R / Fig4C and FigS4C: patient-level MP boxplot exports [original lines 2347-2388]
# pdf(file.path(output_dir, "Figures", "Fig","Fig4C_Main_MP4_MP9_Boxplots.pdf"), width = 3, height = 2.5)
# 
# p_main <- ggpubr::ggboxplot(df_main, x = "Group", y = "Average_Score",
#                     fill = "Group", palette = group_colors,
#                     add = "jitter", 
#                     facet.by = "MP", scales = "free_y", ncol = 2) +
#   ggpubr::stat_compare_means(comparisons = my_comparisons, 
#                              method = "wilcox.test", 
#                              label = "p.format") +
#   ggplot2::theme_classic() +
#   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
#                  strip.background = ggplot2::element_rect(fill = "grey90", color = NA),
#                  strip.text = ggplot2::element_text(face = "bold", size = 12),
#                  legend.position = "none") +
#   ggplot2::labs(title = "Core Meta-Programs: Basal vs Urothelial",
#                 y = "Mean Module Score per Patient", x = "")
# 
# print(p_main)
# dev.off()
# 
# # 2. 过滤附图数据
# df_supp <- patient_avg_long %>% dplyr::filter(MP %in% supp_mps)
# 
# pdf(file.path(output_dir, "Figures","FigS", "FigS04C_Supplementary_MP_Boxplots.pdf"), width = 6, height = 5)
# 
# p_supp <- ggpubr::ggboxplot(df_supp, x = "Group", y = "Average_Score",
#                     fill = "Group", palette = group_colors,
#                     add = "jitter", 
#                     facet.by = "MP", scales = "free_y", ncol = 4) +
#   ggpubr::stat_compare_means(comparisons = my_comparisons, 
#                              method = "wilcox.test", 
#                              label = "p.format") +
#   ggplot2::theme_classic() +
#   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
#                  strip.background = ggplot2::element_rect(fill = "grey90", color = NA),
#                  strip.text = ggplot2::element_text(size = 10),
#                  legend.position = "none") +
#   ggplot2::labs(title = "Supplementary Meta-Programs Analysis",
#                 y = "Mean Module Score per Patient", x = "")
# 
# print(p_supp)
# dev.off()
# <<< END COMMENTED FIGURE BLOCK

pdf(file.path(output_dir,"Figures","MP_Patient_Level_Significance_Boxplots.pdf"), width = 16, height = 10)

p <- ggboxplot(patient_avg_long, x = "Group", y = "Average_Score",
          fill = "Group", palette = group_colors,
          add = "jitter", # 显示代表每个患者的点
          facet.by = "MP", scales = "free_y", ncol = 5) +
  stat_compare_means(comparisons = my_comparisons, 
                    method = "wilcox.test", 
                    label = "p.format") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "grey90"),
        legend.position = "none") +
  labs(title = "Patient-level Statistical Comparison of MP Activities",
      y = "Mean Module Score per Patient", x = "")

print(p)
dev.off()

do.tissueDist.full(malignant.obj@meta.data,
                  meta.cluster = "MP_ID",
                  fig.output = file.path(output_dir,"Figures","MP_Group_Distribution.pdf"),
                  bar_colors = colors_15,pdf_width = 7
                  )

# 1. 提取评分并合并元数据
score_meta <- cbind(malignant.obj@meta.data[, c("sample", "Group")], 
                    malignant.obj@meta.data[, new_names])

# 2. 计算活跃细胞比例 (Score > 0)
active_prop_df <- score_meta %>%
  pivot_longer(cols = all_of(new_names), names_to = "MP", values_to = "Score") %>%
  group_by(sample, Group, MP) %>%
  summarise(
    Active_Count = sum(Score > 0),
    Total_Count = n(),
    Active_Proportion = Active_Count / Total_Count,
    .groups = "drop"
  )

# 3. 设定比较组
my_comparisons <- list(c("Pre_CR", "Pre_NR"), c("Pre_NR", "Post_NR"))

pdf(file.path(output_dir,"Figures","MP_Active_Cell_Proportion_Significance.pdf"), width = 16, height = 10)

ggboxplot(active_prop_df, x = "Group", y = "Active_Proportion",
          fill = "Group", palette = group_colors,
          add = "jitter", 
          facet.by = "MP", scales = "free_y", ncol = 5) +
  stat_compare_means(comparisons = my_comparisons, 
                    method = "wilcox.test", 
                    label = "p.format") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold")) +
  labs(title = "Statistical Comparison of Active Cell Proportion (Score > 0)",
      y = "Proportion of Cells with Score > 0 per Patient", x = "")

dev.off()


# 1. 提取评分矩阵
score_mat <- malignant.obj@meta.data[, new_names]

# 2. 计算每个 MP 的全局均值 (Threshold)
# 也可以考虑使用中位数 median，取决于分布是否对称
#thresholds <- colMeans(score_mat) 
thresholds <- apply(score_mat, 2, median)

# 3. 统计每个样本中大于均值的细胞比例
active_prop_df <- malignant.obj@meta.data %>%
  dplyr::select(sample, Group, all_of(new_names)) %>%
  pivot_longer(cols = all_of(new_names), names_to = "MP", values_to = "Score") %>%
  group_by(MP) %>%
  # 动态匹配对应 MP 的全局均值
  mutate(Is_High = Score > thresholds[unique(MP)]) %>%
  group_by(sample, Group, MP) %>%
  summarise(
    High_Active_Count = sum(Is_High),
    Total_Count = n(),
    High_Active_Proportion = High_Active_Count / Total_Count,
    .groups = "drop"
  )

# 4. 设定比较组
my_comparisons <- list(c("Pre_CR", "Pre_NR"), c("Pre_NR", "Post_NR"))
pdf(file.path(output_dir,"Figures","MP_High_Active_Proportion_GlobalMean_Threshold.pdf"), width = 16, height = 10)

ggboxplot(active_prop_df, x = "Group", y = "High_Active_Proportion",
          fill = "Group", palette = group_colors,
          add = "jitter", 
          facet.by = "MP", scales = "free_y", ncol = 5) +
  stat_compare_means(comparisons = my_comparisons, 
                    method = "wilcox.test", 
                    label = "p.format") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold")) +
  labs(title = "Proportion of High-Active Cells (Score > Global Mean)",
      y = "Proportion of Cells per Patient", x = "")

dev.off()


# 1. 提取所有惡性細胞的評分和分組信息
cor_data <- malignant.obj@meta.data %>%
  dplyr::select(Group, all_of(new_names))

# 2. 定義一個函數來計算並美化相關性矩陣
get_cor_mat <- function(df, group_name) {
  # 過濾特定組
  sub_df <- df %>% dplyr::filter(Group == group_name) %>% dplyr::select(-Group)
  # 計算相關性
  cor_mat <- cor(sub_df, method = "pearson")
  # 轉為長格式方便 ggplot
  cor_df <- melt(cor_mat)
  cor_df$Group <- group_name
  return(cor_df)
}

# 3. 批量計算四個組的相關性
all_groups <- c("Pre_CR", "Pre_NR", "Post_NR", "Post_CR")
cor_combined <- do.call(rbind, lapply(all_groups, function(g) get_cor_mat(cor_data, g)))

# 4. 繪製對比熱圖
pdf(file.path(output_dir,"Figures","MP_Correlation_Comparison_by_Group.pdf"), width = 16, height = 12)

ggplot(cor_combined, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C", limit = c(-1, 1)) +
  facet_wrap(~Group) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Co-activation of MPs across Groups (Cell-level)",
      x = "", y = "", fill = "Pearson Cor")

dev.off()

malignant.obj@meta.data %>%
  dplyr::filter(Group %in% c("Pre_CR", "Post_NR")) %>%
  ggplot(aes(x = Basal_pEMT, y = Hypoxia)) +
  geom_bin2d(bins = 50) +
  scale_fill_viridis_c() +
  facet_wrap(~Group) +
  geom_smooth(method = "lm", color = "red") +
  theme_minimal() +
  labs(title = "Correlation Shift: MP4 vs MP7")

ggscatter(malignant.obj@meta.data, x = "Basal_pEMT", y = "Metabolic",
          color = "Group", palette = group_colors,
          facet.by = "Group", add = "reg.line", 
          conf.int = TRUE, alpha = 0.05) +
  stat_cor(method = "pearson") +
  labs(title = "Evolution of Program Coupling: MP4 vs MP8")

qs_save(malignant.obj,file.path(output_dir, 'malignant_obj_intergration.qs'))
#测试不同k以及如何分MP
if (FALSE) {
  pdf(file.path(output_dir, "Figures","NMF_MPs_Selection_K7-K14.pdf"), width = 14, height = 12)

  k_values <- 7:14

  for (k in k_values) {
    message(paste0("Processing k = ", k))
    
    # 计算当前的聚类分配
    cluster_ids <- cutree(hc, k = k)
    
    # 统计每个 Cluster 的样本跨度 (Sample Span)
    sample_span <- sapply(1:k, function(i) {
      length(unique(anno_df$Sample[cluster_ids == i]))
    })
    avg_span <- round(mean(sample_span), 2)
    
    # 顶部注释
    col_ha <- HeatmapAnnotation(
      Cluster = as.character(cluster_ids),
      Sample = anno_df$Sample,
      Group = anno_df$Group,
      col = list(
        Sample = sample_colors,
        Cluster = setNames(ggsci::pal_d3("category20")(k), as.character(1:k))
      ),
      show_annotation_name = TRUE
    )
    
    # 绘制热图
    ht <- Heatmap(jaccard_mat, 
            name = "Jaccard",
            column_title = paste0("K = ", k, " | Avg Sample Span: ", avg_span, " samples/cluster"),
            column_title_gp = gpar(fontsize = 16, fontface = "bold"),
            top_annotation = col_ha,
            
            # 物理切分
            cluster_rows = hc,
            cluster_columns = hc,
            row_split = k, 
            column_split = k,
            gap = unit(2, "mm"),
            
            show_row_names = FALSE, 
            show_column_names = FALSE,
            border = TRUE,
            col = colorRamp2(c(0, 0.12, 0.4), c("white", "#FFD700", "#FF4500")),
            
            # 在右侧添加简单的统计文本注释
            right_annotation = rowAnnotation(
              span = anno_text(as.character(sample_span[cluster_ids]), 
                              which = "row", gp = gpar(fontsize = 6))
            ))
    
    # 画图并开启新的一页
    draw(ht, merge_legend = TRUE)
    grid.text(paste0("Current K: ", k, " - Evaluate blocks coherence"), 
              x = unit(0.5, "npc"), y = unit(0.05, "npc"), gp = gpar(fontsize = 10, col = "grey"))
  }

  dev.off()
  message("PDF generated: NMF_MPs_Selection_K8-K13.pdf")

  # 定义 Gavish 逻辑下的 MP 提取
  build_gavish_mps <- function(robust_programs_mat, min_overlap = 10) {
    prog_names <- colnames(robust_programs_mat)
    # 计算全局交集矩阵
    inter_intersect <- apply(robust_programs_mat, 2, function(x) {
      apply(robust_programs_mat, 2, function(y) length(intersect(x, y)))
    })
    
    remaining_progs <- prog_names
    final_mps <- list()
    
    while(length(remaining_progs) > 0) {
      # 找到当前剩余池中，与其他样本重叠最广的 Program 作为 Founder
      sub_intersect <- inter_intersect[remaining_progs, remaining_progs, drop=FALSE]
      if(ncol(sub_intersect) == 1) break
      
      diag(sub_intersect) <- 0
      # 选出一个与其他样本有最高重叠的种子
      founder <- names(which.max(apply(sub_intersect, 2, max)))
      
      # 初始 Cluster
      current_cluster <- founder
      
      # 迭代吸纳：只要与当前 Cluster 中任意成员重叠 >= 10，就加入
      while(TRUE) {
        candidates <- setdiff(remaining_progs, current_cluster)
        if(length(candidates) == 0) break
        
        # 计算候选者与当前 Cluster 成员的最大重叠度
        overlaps <- apply(inter_intersect[candidates, current_cluster, drop=FALSE], 1, max)
        best_cand <- names(which.max(overlaps))
        
        if(max(overlaps) >= min_overlap) {
          current_cluster <- c(current_cluster, best_cand)
        } else {
          break
        }
      }
      
      # 记录 MP (要求必须包含至少 2 个不同样本的程序，确保其为 Meta-program)
      samples_in_cluster <- gsub("_R[0-9].*", "", current_cluster)
      if(length(unique(samples_in_cluster)) >= 2) {
        mp_id <- paste0("MP_", length(final_mps) + 1)
        final_mps[[mp_id]] <- current_cluster
      }
      
      # 从池子中移除已分类的（无论是否组成 MP，为了避免死循环）
      remaining_progs <- setdiff(remaining_progs, current_cluster)
    }
    return(final_mps)
  }

  # 运行得到 MP 列表
  gavish_mps <- build_gavish_mps(robust_programs_mat,min_overlap=15)

  all_included_progs <- unlist(gavish_mps)
  n_included <- length(all_included_progs)
  n_excluded <- ncol(robust_programs_mat) - n_included

  message(sprintf("分析完成：\n总程序数: %d\n入选 MP 的程序数: %d\n被剔除的孤儿程序数: %d", 
                  ncol(robust_programs_mat), n_included, n_excluded))

  # 看看每个 MP 覆盖了多少个样本
  mp_sample_coverage <- sapply(gavish_mps, function(x) length(unique(gsub("_R[0-9].*", "", x))))
  print(table(mp_sample_coverage))

  # 按 MP 顺序重新排列矩阵
  ordered_progs <- unlist(gavish_mps)
  sub_jaccard <- jaccard_mat[ordered_progs, ordered_progs]

  # 定义物理分割线
  split_factor <- rep(names(gavish_mps), sapply(gavish_mps, length))
  split_factor <- factor(split_factor, levels = names(gavish_mps))

  pdf(file.path(output_dir,"Figures", "Gavish_Iterative_MPs_Heatmap.pdf"), width = 12, height = 11)
  Heatmap(sub_jaccard, 
          name = "Jaccard",
          column_title = "Meta-Programs Defined by Gavish Iteration (Overlap >= 10)",
          cluster_rows = FALSE, 
          cluster_columns = FALSE,
          row_split = split_factor,
          column_split = split_factor,
          gap = unit(2, "mm"), # 这里的间隙就是被算法“断开”的差异较大程序
          
          top_annotation = column_ha[ordered_progs, ], # 沿用之前的 colors_90 注释
          show_row_names = FALSE, show_column_names = FALSE,
          col = colorRamp2(c(0, 0.111, 0.4), c("white", "#FFD700", "#FF4500")),
          border = TRUE)
  dev.off()
}

malignant.obj@meta.data %>%
  filter(Group %in% c("Pre_CR", "Pre_NR")) %>%
  ggplot(aes(x = Basal_pEMT, y = Uro_Differentiation)) +
  geom_bin2d(bins = 50) +
  scale_fill_viridis_c() +
  facet_wrap(~Group) +
  geom_smooth(method = "lm", color = "red") +
  theme_minimal() +
  labs(title = "Correlation Shift: MP4 vs MP9")




####################################
#看NMF过滤基因的表达量阈值
calc_gene_q98_thresholds <- function(
  seu,
  sample_col = "sample",
  assay = "RNA",
  layer = "data",              # 建议用 NormalizeData 后的 data（log1p(CP10K)）
  min_cells = 200,
  prob = 0.98,
  top_ns = c(2000, 2500, 3000),
  verbose = TRUE
) {
  # 1) 确认 sample 列
  md <- seu@meta.data
  if (!sample_col %in% colnames(md)) {
    if ("orig.ident" %in% colnames(md)) {
      if (isTRUE(verbose)) message("sample_col 不存在，改用 'orig.ident' 作为分组列。")
      sample_col <- "orig.ident"
    } else {
      stop("在 meta.data 中找不到 sample_col，也找不到 orig.ident。请指定正确的分组列名。")
    }
  }

  # 2) 按 sample 拆分
  seu_list <- Seurat::SplitObject(seu, split.by = sample_col)

  # 3) 逐 sample 计算 q98
  res_list <- lapply(names(seu_list), function(sid) {
    obj <- seu_list[[sid]]
    n_cells <- ncol(obj)
    if (n_cells < min_cells) {
      if (isTRUE(verbose)) message(sprintf("Skip %s: nCells=%d < %d", sid, n_cells, min_cells))
      return(NULL)
    }

    mat <- SeuratObject::GetAssayData(obj, assay = assay, layer = layer)  # features x cells
    n_genes <- nrow(mat)

    # 推荐：稀疏矩阵用 sparseMatrixStats::rowQuantiles（非常关键，速度和内存都好很多）
    q <- NULL
    if (inherits(mat, "dgCMatrix") && requireNamespace("sparseMatrixStats", quietly = TRUE)) {
      q <- sparseMatrixStats::rowQuantiles(mat, probs = prob, na.rm = TRUE)
    } else if (requireNamespace("matrixStats", quietly = TRUE)) {
      # mat 若不是稀疏矩阵也可用 matrixStats；若是稀疏矩阵但没 sparseMatrixStats，这里会较慢/可能强制转 dense
      q <- matrixStats::rowQuantiles(as.matrix(mat), probs = prob, na.rm = TRUE)
    } else {
      stop("建议安装 sparseMatrixStats（稀疏矩阵更快）：install.packages('sparseMatrixStats')")
    }

    names(q) <- rownames(mat)

    # 从大到小排序：第 N 名对应的值就是“topN 基因阈值”
    q_sorted <- sort(q, decreasing = TRUE, na.last = TRUE)

    out <- data.frame(
      sample = sid,
      nCells = n_cells,
      nGenes = n_genes,
      stringsAsFactors = FALSE
    )

    for (N in top_ns) {
      coln <- paste0("q", as.integer(prob * 100), "_top", N)
      out[[coln]] <- if (length(q_sorted) >= N) unname(q_sorted[N]) else NA_real_
    }

    out
  })

  res <- do.call(rbind, res_list)
  if (is.null(res) || nrow(res) == 0) {
    warning("没有任何 sample 通过 min_cells 过滤。")
    return(res)
  }

  rownames(res) <- res$sample
  res
}

df_thr <- calc_gene_q98_thresholds(
  seu = malignant.obj,
  sample_col = "sample",   # 如果你的列不是 sample，请改成实际列名
  assay = "RNA",
  layer = "data",
  min_cells = 200,
  prob = 0.98,
  top_ns = c(2000, 2500, 3000)
)

View(df_thr)
write.csv(df_thr,file.path(output_dir,'malignant_expression_thresholds.csv'))

############################################################
#分析progeny和scMetabolism的结果
malignant.obj = qs_read(file.path(output_dir, 'malignant_obj_intergration.qs'))
scMetabolism.obj = qs_read(file.path(output_dir, 'seu_combined_IT_filtered_obj_intergration_metabolism.qs'))
progeny.obj = qs_read(file.path(output_dir, 'seu_combined_IT_filtered_obj_intergration_progeny.qs'))
Metabolism.counts = scMetabolism.obj@assays$METABOLISM$score
colnames(Metabolism.counts) = gsub("\\.","-",colnames(Metabolism.counts))
Metabolism.counts=t(Metabolism.counts)
View(Metabolism.counts)
progeny.counts = t(progeny.obj@assays$progeny$data)
View(progeny.counts)
colnames(progeny.counts)[4] <- "Hypoxia_P"
meta = malignant.obj@meta.data
meta[,colnames(Metabolism.counts)] = Metabolism.counts[rownames(meta),]
meta[,colnames(progeny.counts)] = progeny.counts[rownames(meta),]
metabolism.cols = colnames(Metabolism.counts)
progeny.cols = colnames(progeny.counts)
malignant.obj@meta.data = meta
MP.cols = unique(malignant.obj$MP_ID)
qs_save(malignant.obj,file.path(output_dir, 'malignant_obj_intergration.qs'))

###########################
#差异分析
DEG.Pre_NR_vs_Pre_CR  = qs_read(file.path(output_dir,"muscat","DEG_output_Pre_NR-Pre_CR.qs"))
DEG.Pre_NR_vs_Post_NR = qs_read(file.path(output_dir,"muscat","DEG_output_Pre_NR-Post_NR.qs"))
View(DEG.Pre_NR_vs_Pre_CR[[9]]$celltype_de$celltype_de$de_output_tidy)
DE_genes = DEG.Pre_NR_vs_Pre_CR[[10]]$celltype_de$celltype_de$de_output_tidy  %>%
          filter(contrast == 'Pre_NR-Pre_CR') %>% 
          filter(p_adj <= 0.05 & abs(logFC) >= 1) %>% 
          arrange(p_adj) %>% pull(gene) %>% unique()
violin_plot = make_DEgene_violin_plot(sce = sce, gene_oi = gene_oi, celltype_oi = celltype_oi, group_id = group_id, sample_id = sample_id, celltype_id = celltype_id)

# 定义相关性计算函数
get_cor_matrix <- function(obj, group_name, cols_x, cols_y) {
  # 提取特定组的数据
  data <- obj@meta.data %>% filter(Group == group_name)
  
  # 计算 X (MPs) 与 Y (Metabolism/Progeny) 的相关性
  cor_mat <- cor(data[, cols_x], data[, cols_y], method = "pearson")
  return(cor_mat)
}

# 1. 计算 Metabolism 的组间相关性
metab_cor_list <- lapply(unique(malignant.obj$Group), function(g) {
  get_cor_matrix(malignant.obj, g, MP.cols, metabolism.cols)
})
names(metab_cor_list) <- unique(malignant.obj$Group)

# 2. 计算 PROGENY 的组间相关性
progeny_cor_list <- lapply(unique(malignant.obj$Group), function(g) {
  get_cor_matrix(malignant.obj, g, MP.cols, progeny.cols)
})
names(progeny_cor_list) <- unique(malignant.obj$Group)

pdf(file.path(output_dir, 'Figures','Metabolism_Heatmap.pdf'),width = 12,height = 7)
for (i in 1:length(metab_cor_list)){
  ht_metab <- Heatmap(metab_cor_list[[i]], 
        name = "Pearson Cor",
        column_title = paste0("Metabolism Correlation with ",names(metab_cor_list)[i]),
        col = colorRamp2(c(-0.5, 0, 0.5), c("#377EB8", "white", "#E41A1C")),
        cluster_rows = F, cluster_columns = F,
        show_column_names = TRUE, column_names_gp = gpar(fontsize = 8))
  draw(ht_metab)
}
dev.off()

pdf(file.path(output_dir, 'Figures','Progeny_Heatmap.pdf'),width = 7,height = 5)
for (i in 1:length(progeny_cor_list)){
  ht_progeny <- Heatmap(progeny_cor_list[[i]], 
                        name = "Pearson Cor",
                        column_title = paste0("Progeny Correlation with ",names(progeny_cor_list)[i]),
                        col = colorRamp2(c(-0.5, 0, 0.5), c("#377EB8", "white", "#E41A1C")),
                        cluster_rows = F, cluster_columns = F,
                        show_column_names = TRUE, column_names_gp = gpar(fontsize = 8))
  draw(ht_progeny)
}
dev.off()
pdf(file.path(output_dir, 'Figures','Progeny_DotPlot.pdf'),width = 12,height = 12)
VlnPlot(malignant.obj,features = progeny.cols,group.by = 'MP_ID',split.by = 'Group',cols = colors_90,stack=T,flip=T)
dev.off()

# 1. 筛选 Pre_CR 和 Pre_NR 组的数据
plot_data <- malignant.obj@meta.data %>%
  dplyr::filter(Group %in% c("Pre_CR", "Pre_NR")) %>%
  # 选择需要的列：MP4 和 四个通路
  # 请确保列名与你的 progeny.cols 中一致（例如 "EGFR", "Estrogen" 等）
  dplyr::select(Group, Basal_pEMT, EGFR, Estrogen, Hypoxia_P, p53, WNT) 

# 2. 转换为长格式以便于 facet 分面展示
plot_long <- plot_data %>%
  pivot_longer(cols = c("EGFR", "Estrogen", "Hypoxia_P", "p53","WNT"), 
              names_to = "Pathway", values_to = "Pathway_Score")
group_colors_subset <- c("Pre_CR" = "#053E7A", "Pre_NR" = "#A52A2A")



pdf(file.path(output_dir, 'Figures', 'MP4_vs_Pathways_Bin2d.pdf'), width = 6, height = 9)

# 建议使用 ggplot2 核心语法以获得更好的控制力
p <- ggplot(plot_long, aes(x = Basal_pEMT, y = Pathway_Score)) +
  # 使用 geom_bin2d 代替散点，bins 调整格子的细腻程度（建议 50-80）
  geom_bin2d(bins = 100) + 
  
  # 使用 viridis 颜色梯度：黄色代表高密度，紫色代表低密度
  scale_fill_viridis(option = "magma", direction = -1, name = "Cell Count") +
  
  # 添加回归线（注意：在分箱图之上添加线需要明确颜色，防止被底色遮盖）
  geom_smooth(method = "lm", color = "cyan", size = 1, fill = "lightgrey", alpha = 0.2) +
  
  # 按照 Pathway 和 Group 进行分面
  facet_grid(Pathway ~ Group, scales = "free_y") +
  
  # 自动计算相关系数
  stat_cor(method = "pearson", 
           label.x.npc = "left", 
           label.y.npc = "top", 
           size = 3.5, 
           color = "white") + # 白色字体在深色背景下更清晰
  
  # 样式美化
  theme_minimal() +
  labs(title = "Density Distribution: MP4 and Signaling Pathways Coupling",
       subtitle = "Density-based visualization to avoid overplotting",
       x = "MP4 (Basal_pEMT) Score",
       y = "PROGENY Pathway Activity") +
  theme(strip.background = element_rect(fill = "grey20"),
        strip.text = element_text(color = "white", face = "bold"),
        panel.grid = element_blank(), # 去掉网格线让密度感更强
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(p)

dev.off()


plot_data <- malignant.obj@meta.data %>%
  dplyr::filter(Group %in% c("Pre_CR", "Pre_NR")) %>%
  # 选择需要的列：MP4 和 四个通路
  # 请确保列名与你的 progeny.cols 中一致（例如 "EGFR", "Estrogen" 等）
  dplyr::select(Group, Basal_pEMT,'Metabolism of xenobiotics by cytochrome P450',
                'Drug metabolism - cytochrome P450', 'Drug metabolism - other enzymes') 

# 2. 转换为长格式以便于 facet 分面展示
plot_long <- plot_data %>%
  pivot_longer(cols = c("Metabolism of xenobiotics by cytochrome P450",
              "Drug metabolism - cytochrome P450", "Drug metabolism - other enzymes"), 
              names_to = "Pathway", values_to = "Pathway_Score")
group_colors_subset <- c("Pre_CR" = "#053E7A", "Pre_NR" = "#A52A2A")

pdf(file.path(output_dir, 'Figures','MP4_vs_metabolism_Pathways_PreCorrelation.pdf'), width = 6, height = 7)

p <- ggplot(plot_long, aes(x = Basal_pEMT, y = Pathway_Score)) +
  # 使用 geom_bin2d 代替散点，bins 调整格子的细腻程度（建议 50-80）
  geom_bin2d(bins = 100) + 
  
  # 使用 viridis 颜色梯度：黄色代表高密度，紫色代表低密度
  scale_fill_viridis(option = "magma", direction = -1, name = "Cell Count") +
  
  # 添加回归线（注意：在分箱图之上添加线需要明确颜色，防止被底色遮盖）
  geom_smooth(method = "lm", color = "cyan", size = 1, fill = "lightgrey", alpha = 0.2) +
  
  # 按照 Pathway 和 Group 进行分面
  facet_grid(Pathway ~ Group, scales = "free_y") +
  
  # 自动计算相关系数
  stat_cor(method = "pearson", 
           label.x.npc = "left", 
           label.y.npc = "top", 
           size = 3.5, 
           color = "white") + # 白色字体在深色背景下更清晰
  
  # 样式美化
  theme_minimal() +
  labs(title = "Density Distribution: MP4 and Signaling Pathways Coupling",
       subtitle = "Density-based visualization to avoid overplotting",
       x = "MP4 (Basal_pEMT) Score",
       y = "PROGENY Pathway Activity") +
  theme(strip.background = element_rect(fill = "grey20"),
        strip.text = element_text(color = "white", face = "bold"),
        panel.grid = element_blank(), # 去掉网格线让密度感更强
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(p)

dev.off()

plot_df <- data.frame(
  MP4 = malignant.obj$Basal_pEMT,
  MP9 = malignant.obj$Uro_Differentiation,
  ERBB2 = FetchData(malignant.obj, vars = "ERBB2")[,1],
  Group = factor(malignant.obj$Group, levels = c("Pre_CR", "Pre_NR", "Post_CR", "Post_NR"))
)


# 2. 绘制 MP4 vs ERBB2 密度图
p_mp4 <- ggplot(plot_df, aes(x = MP4, y = ERBB2)) +
  geom_bin2d(bins = 80) + 
  scale_fill_viridis_c(option = "magma", trans = "log10") + 
  facet_wrap(~Group, ncol = 4) +
  geom_smooth(method = "lm", color = "cyan", size = 0.8, linetype = "dashed") +
  stat_cor(method = "pearson", color = "white", size = 3.5, 
          label.x.npc = "left", label.y.npc = "top") + # 显式显示 R 和 P
  theme_bw() +
  labs(title = "Negative Coupling: Basal Program (MP4) vs ERBB2",
      x = "MP4 (Basal_pEMT) Score", y = "ERBB2 Expression", fill = "Log10(Count)") +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "#A52A2A"),
        strip.text = element_text(color = "white", face = "bold"))

# 3. 绘制 MP9 vs ERBB2 密度图
p_mp9 <- ggplot(plot_df, aes(x = MP9, y = ERBB2)) +
  geom_bin2d(bins = 80) + 
  scale_fill_viridis_c(option = "magma", trans = "log10") + 
  facet_wrap(~Group, ncol = 4) +
  geom_smooth(method = "lm", color = "cyan", size = 0.8, linetype = "dashed") +
  stat_cor(method = "pearson", color = "white", size = 3.5, 
          label.x.npc = "left", label.y.npc = "top") + # 显式显示 R 和 P
  theme_bw() +
  labs(title = "Positive Coupling: Luminal Program (MP9) vs ERBB2",
      x = "MP9 (Uro_Differentiation) Score", y = "ERBB2 Expression", fill = "Log10(Count)") +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "#053E7A"),
        strip.text = element_text(color = "white", face = "bold"))

# 保存图片
pdf(file.path(output_dir,"Figures","MP_vs_ERBB2_Density.pdf"), width = 8, height = 6)
print(p_mp4 / p_mp9)
dev.off()

plot_df <- data.frame(
  MP4 = malignant.obj$Basal_pEMT,
  MP9 = malignant.obj$Uro_Differentiation,
  PDL1 = FetchData(malignant.obj, vars = "CD274")[,1],
  Group = factor(malignant.obj$Group, levels = c("Pre_CR", "Pre_NR", "Post_CR", "Post_NR"))
)
plot_df <- plot_df %>% filter(PDL1 > 0)

# 2. 绘制 MP4 vs ERBB2 密度图
p_mp4 <- ggplot(plot_df, aes(x = MP4, y = PDL1)) +
  geom_bin2d(bins = 80) + 
  scale_fill_viridis_c(option = "magma", trans = "log10") + 
  facet_wrap(~Group, ncol = 4) +
  geom_smooth(method = "lm", color = "cyan", size = 0.8, linetype = "dashed") +
  stat_cor(method = "pearson", color = "white", size = 3.5, 
          label.x.npc = "left", label.y.npc = "top") + # 显式显示 R 和 P
  theme_bw() +
  labs(title = "Negative Coupling: Basal Program (MP4) vs PD-L1",
      x = "MP4 (Basal_pEMT) Score", y = "PD-L1 Expression", fill = "Log10(Count)") +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "#A52A2A"),
        strip.text = element_text(color = "white", face = "bold"))

# 3. 绘制 MP9 vs ERBB2 密度图
p_mp9 <- ggplot(plot_df, aes(x = MP9, y = PDL1)) +
  geom_bin2d(bins = 80) + 
  scale_fill_viridis_c(option = "magma", trans = "log10") + 
  facet_wrap(~Group, ncol = 4) +
  geom_smooth(method = "lm", color = "cyan", size = 0.8, linetype = "dashed") +
  stat_cor(method = "pearson", color = "white", size = 3.5, 
          label.x.npc = "left", label.y.npc = "top") + # 显式显示 R 和 P
  theme_bw() +
  labs(title = "Positive Coupling: Luminal Program (MP9) vs PD-L1",
      x = "MP9 (Uro_Differentiation) Score", y = "PD-L1 Expression", fill = "Log10(Count)") +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "#053E7A"),
        strip.text = element_text(color = "white", face = "bold"))

# 保存图片
pdf(file.path(output_dir,"Figures","MP_vs_PDL1_Density.pdf"), width = 8, height = 6)
print(p_mp4 / p_mp9)
dev.off()

# >>> COMMENTED FIGURE BLOCK: Fig1.R and Fig2.R: subtype composition, RoE and fraction plots [original lines 3076-3163]
# ######################################################################
# #画subcelltype比例图
# Myeloid.obj = qs_read(file.path(output_dir, 'Myeloid_obj_intergration.qs'))
# Myeloid.obj$Group = factor(Myeloid.obj$Group,levels = c('Pre_NR','Pre_CR','Post_NR','Post_CR'))
# #CellStatPlot(Myeloid.obj,stat.by = 'subcelltype',group.by = 'Group',palcolor = colors_90)
# CellStatPlot(Myeloid.obj,stat.by = 'subcelltype',group.by = 'Group',palcolor = colors_90,plot_type = 'dot')
# ggsave(file.path(output_dir, 'Figures','Myeloid_subcelltype_ratio.pdf'),width = 6,height =7)
# do.tissueDist(Myeloid.obj@meta.data,
#               fig.output = file.path(output_dir, 'Figures','Myeloid_subcelltype_Roe.pdf'))
# 
# do.tissueDist.full(
#     Myeloid.obj@meta.data,
#     meta.cluster = 'subcelltype',
#     Group.info = 'Group',
#     bar_colors = colors_15, # 传入自定义颜色
#     pdf_width = 10,         # 增加宽度以容纳 18 个亚群
#     fig.output = file.path(output_dir, 'Figures', 'Myeloid_subcelltype_Roe_with_bar.pdf')
# )
# 
# Stromal.obj = qs_read(file.path(output_dir, 'Stromal_obj_intergration.qs'))
# Stromal.obj$Group = factor(Stromal.obj$Group,levels = c('Pre_NR','Pre_CR','Post_NR','Post_CR'))
# CellStatPlot(Stromal.obj,stat.by = 'subcelltype',group.by = 'Group',palcolor = colors_90,plot_type = 'dot')
# ggsave(file.path(output_dir, 'Figures','Stromal_subcelltype_ratio.pdf'),width = 6,height =7)
# do.tissueDist(Stromal.obj@meta.data,
#               fig.output = file.path(output_dir, 'Figures','Stromal_subcelltype_Roe.pdf'))
# do.tissueDist.full(
#     Stromal.obj@meta.data,
#     meta.cluster = 'subcelltype',
#     Group.info = 'Group',
#     bar_colors = colors_15, # 传入自定义颜色
#     pdf_width = 10,         # 增加宽度以容纳 18 个亚群
#     fig.output = file.path(output_dir, 'Figures', 'Stromal_subcelltype_Roe_with_bar.pdf')
# )
# 
# 
# Tcell.obj = qs_read(file.path(output_dir, 'Tcell_obj_intergration.qs'))
# Tcell.obj$Group = factor(Tcell.obj$Group,levels = c('Pre_NR','Pre_CR','Post_NR','Post_CR'))
# CellStatPlot(Tcell.obj,stat.by = 'subcelltype',group.by = 'Group',palcolor = colors_90,plot_type = 'dot')
# ggsave(file.path(output_dir, 'Figures','Tcell_subcelltype_ratio.pdf'),width = 6,height =7)
# do.tissueDist(Tcell.obj@meta.data,
#               fig.output = file.path(output_dir, 'Figures','Tcell_subcelltype_Roe.pdf'))
# do.tissueDist.full(
#     Tcell.obj@meta.data,
#     meta.cluster = 'subcelltype',
#     Group.info = 'Group',
#     bar_colors = colors_15, # 传入自定义颜色
#     pdf_width = 8,         # 增加宽度以容纳 18 个亚群
#     fig.output = file.path(output_dir, 'Figures', 'Tcell_subcelltype_Roe_with_bar.pdf')
# )
# 
# Bcell.obj = qs_read(file.path(output_dir, 'Bcell_obj_intergration.qs'))
# Bcell.obj$Group = factor(Bcell.obj$Group,levels = c('Pre_NR','Pre_CR','Post_NR','Post_CR'))
# CellStatPlot(Bcell.obj,stat.by = 'subcelltype',group.by = 'Group',palcolor = colors_90,plot_type = 'dot')
# ggsave(file.path(output_dir, 'Figures','Bcell_subcelltype_ratio.pdf'),width = 6,height =5)
# do.tissueDist(Bcell.obj@meta.data,
#               fig.output = file.path(output_dir, 'Figures','Bcell_subcelltype_Roe.pdf'),pdf.height = 4)
# do.tissueDist.full(
#     Bcell.obj@meta.data,
#     meta.cluster = 'subcelltype',
#     Group.info = 'Group',
#     bar_colors = colors_15, # 传入自定义颜色
#     pdf_width = 3,         # 增加宽度以容纳 18 个亚群
#     pdf_height = 3.5,
#     top_bar_height_cm = 2,
#     fig.output = file.path(output_dir, 'Figures', 'Bcell_subcelltype_Roe_with_bar.pdf')
# )
# 
# 
# plot_cell_fraction_expert(Myeloid.obj,
#                           celltype_col = 'subcelltype',
#                           comparisons = list(c('Pre_NR','Pre_CR'),c('Pre_NR','Post_NR')),
#                           fig_output = file.path(output_dir, 'Figures','Myeloid_subcelltype_fraction_boxplot.pdf'),
#                           min_fraction = 0.02)
# plot_cell_fraction_expert(Stromal.obj,
#                           celltype_col = 'subcelltype',
#                           comparisons = list(c('Pre_NR','Pre_CR'),c('Pre_NR','Post_NR')),
#                           fig_output = file.path(output_dir, 'Figures','Stromal_subcelltype_fraction_boxplot.pdf'),
#                           min_fraction = 0.02)
# plot_cell_fraction_expert(Tcell.obj,
#                           celltype_col = 'subcelltype',
#                           comparisons = list(c('Pre_NR','Pre_CR'),c('Pre_NR','Post_NR')),
#                           fig_output = file.path(output_dir, 'Figures','Tcell_subcelltype_fraction_boxplot.pdf'),
#                           min_fraction = 0.02)
# plot_cell_fraction_expert(Bcell.obj,
#                           celltype_col = 'subcelltype',
#                           comparisons = list(c('Pre_NR','Pre_CR'),c('Pre_NR','Post_NR')),
#                           fig_output = file.path(output_dir, 'Figures','Bcell_subcelltype_fraction_boxplot.pdf'),
#                           min_fraction = 0.02)
# <<< END COMMENTED FIGURE BLOCK


################################################
# >>> COMMENTED FIGURE BLOCK: Fig2.R / Fig2B-C: T-cell markers and signature heatmaps [original lines 3167-3211]
# #T细胞分析
# markers.T = cosg(Tcell.obj,
#                 groups='all',
#                 assay = 'RNA',
#                 mu=100,
#                 remove_lowly_expressed = TRUE,
#                 expressed_pct = 0.3,
#                 n_genes_user=50)
# View(markers.T$names)
# markers.Tcell = c("KLRG1","GZMK",'CCL5',"EOMES","CCR5","SAMD3","ZNF683","HAVCR2","LAG3","IFIT1","MX1",
#                   "MKI67","TOP2A","IL7R","TCF7","CCR7","S1PR1","CD40LG","CXCL13","FKBP5","PDCD1",
#                 "RTKN2","FOXP3","CCR4","CCR8","IL2RA","ICOS","RORA","KLRB1","GPR183","CX3CR1","FCGR3A",
#                 "NCAM1","FCER1G","GNLY","MYB","PRDM8")
# 
# p=jjDotPlot(object = Tcell.obj,
#           gene = markers.Tcell,
#           id = 'subcelltype',
#           #gene.order = rev(markers.Tcell),
#           #cluster.order = colnames(markers.Tcell$names),
#           ytree = F) + theme(axis.text.x = element_text(angle = 45,vjust = 0.2))
# p
# ggsave(file.path(output_dir, 'Figures','Tcell_markers.pdf'),plot=p,width = 12,height =4)
# 
# Tcell.signiture = read.csv("~/scRNAseq/T_cell_gene_list.csv",header = T)
# #dataframe to list
# Tcell.signiture = as.list(Tcell.signiture)
# Tcell.signiture = lapply(Tcell.signiture,function(x) {
#   x = gsub("\\.","-",x)
#   x = intersect(x,rownames(Tcell.obj))
# })
# 
# Tcell.obj  = AddModuleScore(Tcell.obj,features = Tcell.signiture,name = 'Tcell_signiture')
# colnames(Tcell.obj@meta.data)[grepl("^Tcell_signiture",colnames(Tcell.obj@meta.data))]  = paste0(names(Tcell.signiture),'_score')
# features = paste0(names(Tcell.signiture),'_score')
# ht1 <- GroupHeatmap(Tcell.obj,features = features,group.by = 'subcelltype')
# ht1$plot
# ggsave(file.path(output_dir, 'Figures','Tcell_signiture_heatmap_Groupheatmap.pdf'),plot=ht1$plot,width = 8,height =8)
# do.SignatureHeatmap(Tcell.obj,features = features,fig_output = file.path(output_dir, 'Figures','Tcell_signiture_heatmap.pdf'),plot_margin = c(1,1,1,1))
# 
# Tcell.obj = calculate_simple_gene_set_score(Tcell.obj,gene_set = Tcell.signiture, score_name = names(Tcell.signiture))
# features = paste0(names(Tcell.signiture),'_AveExp')
# ht1 <- GroupHeatmap(Tcell.obj,features = features,group.by = 'subcelltype')
# ht1$plot
# ggsave(file.path(output_dir, 'Figures','Tcell_signiture_heatmap_Groupheatmap_AveExp.pdf'),plot=ht1$plot,width = 8,height =8)
# do.SignatureHeatmap(Tcell.obj,features = features,fig_output = file.path(output_dir, 'Figures','Tcell_signiture_heatmap_AveExp.pdf'),plot_margin = c(1,1,1,1))
# <<< END COMMENTED FIGURE BLOCK

#monocle2轨迹分析
run_monocle2 <- function(seurat,rm=T){
  #outfile='cds_by_RNA_snn_res.0.1.Rdata'
  library(monocle) 
  expr_matrix=seurat@assays$RNA$counts#使用counts表达值
  sample_sheet<-seurat@meta.data#将实验信息赋值新变量
  gene_annotation=data.frame(gene_short_name=rownames(seurat))#构建一个含有基因名字的数据框
  rownames(gene_annotation)=rownames(seurat)#将上述数据框的行名赋值基因名字
  
  pd <- new("AnnotatedDataFrame", data = sample_sheet)#将实验信息变量转化为monocel可以接收的对象
  fd <- new("AnnotatedDataFrame", data = gene_annotation)#将基因注释变量转化为monocle可以接收的对象
  cds <- newCellDataSet(expr_matrix, phenoData = pd, featureData = fd,
                        expressionFamily=negbinomial.size())#创建一个monocle的对象
  
  cds #cellData对象；monocle独有
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  
  ## 拟时序只需要3个步骤： 
  # 第一步: 挑选合适的基因. 有多个方法，
  # 例如提供已知的基因集，或者自己做差异分析后拿到的统计学显著的差异基因列表 
  if(rm){
    diff_test_res <- differentialGeneTest(cds,
                                          fullModelFormulaStr = " ~ subcelltype + sample ",
                                          reducedModelFormulaStr = " ~ sample",
                                          relative_expr=TRUE, cores=20) 
  }else{
    diff_test_res <- differentialGeneTest(cds,
                                          fullModelFormulaStr = " ~ subcelltype ",
                                          #reducedModelFormulaStr = " ~ orig.ident",
                                          relative_expr=TRUE, cores=1) 
  }
  ordering_genes = VariableFeatures(seurat)
  ordering_genes <- row.names (base::subset(diff_test_res, qval < 0.01))
  ordering_genes
  cds <- setOrderingFilter(cds, ordering_genes)
  plot_ordering_genes(cds) 
  
  # 第二步: 降维。降维的目的是为了更好的展示数据。函数里提供了很多种方法,
  # 不同方法的最后展示的图都不太一样, 其中“DDRTree”是Monocle2使用的默认方法
  if(rm){
    cds <- reduceDimension(cds, 
                          max_components = 2, 
                          num_dim = 6, 
                          reduction_method = 'DDRTree', 
                          residualModelFormulaStr = "~sample", #去除样本影响
                          verbose = F)
  }
  else{
    cds <- reduceDimension(cds, 
                          max_components = 2, 
                          num_dim = 6, 
                          reduction_method = 'DDRTree', 
                          #residualModelFormulaStr = "~orig.ident", #去除样本影响
                          verbose = F)
  }

  # 第三步: 对细胞进行排序
  cds <- orderCells(cds)
  #save(cds,file = outfile)
  return(cds)
  
}
Tcell.obj = qs_read(file.path(output_dir, 'Tcell_obj_intergration.qs'))
table(Tcell.obj$subcelltype)
CD8T = subset(Tcell.obj,subset = subcelltype %in% c("t01_CD8_KLRG1",'t02_CD8_CCR5 (Teff)','t03_CD8_LAG3 (Tex)','t04_CD8_IFIT1'))
CD8T = subcelltype_analysis(CD8T,output_dir = output_dir,prefix = 'CD8T')
CD8T = qs_read(file.path(output_dir, 'CD8T_obj_intergration.qs'))
DimPlot(CD8T,group.by = 'subcelltype')
CD8.cds = run_monocle2(CD8T)
# >>> COMMENTED FIGURE BLOCK: Fig2.R / Fig2D: CD8 T-cell subtype fraction plot [original lines 3283-3287]
# plot_cell_fraction_expert(CD8T,
#                           celltype_col = 'subcelltype',
#                           comparisons = list(c('Pre_NR','Pre_CR'),c('Pre_NR','Post_NR')),
#                           fig_output = file.path(output_dir, 'Figures','CD8T_subcelltype_fraction_boxplot.pdf'),
#                           min_fraction = 0.02,pdf_width=6,pdf_height=4)
# <<< END COMMENTED FIGURE BLOCK
subcelltype_analysis_SCT <- function(obj,output_dir,prefix,npcs = 30,res = seq(0.1,1.5,0.1)) {
  meta.raw = obj@meta.data
  meta.raw = meta.raw[,!grepl('^RNA_snn_res',colnames(meta.raw))]
  meta.raw = meta.raw[,!grepl('^SCT_snn_res',colnames(meta.raw))]
  obj = CreateSeuratObject(counts = obj@assays$RNA$counts,project = prefix,min.cells = 10,min.features = 200)
  meta = obj@meta.data
  add.cols = setdiff(colnames(meta.raw),colnames(meta))
  meta[add.cols] = meta.raw[row.names(meta),add.cols]
  obj@meta.data = meta
  obj$orig.ident = obj$sample

  # --- 新增过滤功能：过滤掉细胞数少于 20 的样本 ---
  sample_table <- table(obj$sample)
  keep_samples <- names(sample_table)[sample_table >= 20]
  remove_samples <- names(sample_table)[sample_table < 20]
  
  if (length(remove_samples) > 0) {
    message("！！！警告：以下样本细胞数少于 20，将被剔除：", paste(remove_samples, collapse = ", "))
    # 执行过滤
    obj <- subset(obj, subset = sample %in% keep_samples)
  } else {
    message("所有样本细胞数均大于等于 20，无需额外过滤。")
  }
  #obj@meta.data = obj@meta.data[,!grepl('^RNA_snn_res',colnames(obj@meta.data))]
  obj[['RNA']] = split(obj[['RNA']],f=obj$sample)
  future::plan(sequential)
  obj = NormalizeData(obj,verbose = T)
  obj = SCTransform(obj, 
                    vars.to.regress = c("S.Score", "G2M.Score",'nCount_RNA','pct_counts_mt'),
                    verbose = FALSE,
                    vst.flavor = "v2",
                    min_cells=0)

  hvg = VariableFeatures(obj)
  MT.genes = rownames(obj)[grepl('^MT-',rownames(obj))]
  HB.genes = rownames(obj)[grepl('^HB[^(PES)]',rownames(obj))]
  RP.genes = rownames(obj)[grepl('^RP[SL]',rownames(obj))]
  blacklist = c(MT.genes,HB.genes,RP.genes)
  hvg = setdiff(hvg,blacklist)
  VariableFeatures(obj) = hvg
  scaled_genes <- rownames(GetAssayData(obj, assay = "SCT", layer = "scale.data"))
  missing_genes <- setdiff(hvg, scaled_genes)
  all_present <- length(missing_genes) == 0
  if (all_present) {
    message("所有 HVG 基因均已标准化。")
  } else {
    message("！！！警告：以下 HVG 基因未被标准化：", paste(missing_genes, collapse = ", "))
  }
  obj = obj %>%
    RunPCA(.,verbose = T) %>%
    IntegrateLayers(.,method = HarmonyIntegration,
                    new.reduction = "harmony",orig.reduction = "pca",assay='SCT',max.iter.harmony = 20L,verbose=T) %>%
    RunUMAP(.,reduction = "harmony",dims = 1:npcs) %>%
    FindNeighbors(., reduction = "harmony", dims = 1:npcs,verbose = F) %>%
    FindClusters(., resolution = res,algorithm=4,verbose = F,random.seed=1234)
  DefaultAssay(obj) = "RNA"
  obj = JoinLayers(obj)
  qs_save(obj,file=file.path(output_dir, paste0(prefix,'_obj_intergration.qs')))
  return(obj)
}

CD8T = qs_read(file.path(output_dir, 'CD8T_obj_intergration.qs'))
CD8T_SCT = subcelltype_analysis_SCT(CD8T,output_dir = output_dir,prefix = 'CD8T_SCT')
plot_density_dimplot_new(seurat_object = CD8T_SCT,
                    group.by = 'SCT_snn_res.0.3',
                    split.by = 'Group',
                    ncol = 2,
                    title = 'CD8 Tcell',
                    color_palette = 'magma',group_colors = colors_90)
ggsave(file.path(output_dir, 'CD8_Tcell_subcelltype_density_SCT.png'),width = 14,height =12,dpi = 300)
DimPlot(CD8T_SCT,group.by = 'subcelltype')
DimPlot(CD8T,group.by = 'RNA_snn_res.0.3',cols = colors_15,label = T)
Idents(CD8T) = 'RNA_snn_res.0.3'
markers.CD8T = cosg(CD8T,
                      groups='all',
                      assay = 'RNA',
                      mu=100,
                      remove_lowly_expressed = TRUE,
                      expressed_pct = 0.3,
                      n_genes_user=50)
View(markers.CD8T$names)
seurat_to_monocle3 = function(seurat_obj) {
  # 假設你的 Seurat 物件名為 seurat_obj

  # 1. 提取表達矩陣 (使用 counts 矩陣)
  expression_matrix <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")

  # 2. 提取細胞元數據 (包含 celltype)
  cell_metadata <- seurat_obj@meta.data

  # 3. 提取基因元數據
  gene_annotation <- data.frame(
    gene_short_name = rownames(expression_matrix), 
    row.names = rownames(expression_matrix)
  )

  # 4. 構建 Monocle 3 的 cds 物件
  cds <- new_cell_data_set(
    expression_matrix,
    cell_metadata = cell_metadata,
    gene_metadata = gene_annotation
  )

  # 5. 【關鍵】直接讀取 Seurat 的 UMAP 座標
  # 注意：Monocle 3 的單元格順序必須與 Seurat 一致
  reducedDims(cds)$UMAP <- seurat_obj@reductions$umap@cell.embeddings

  # 6. 【關鍵】讀取 Seurat 的聚類分群或 CellType
  # 雖然 metadata 已經帶入，但 Monocle 3 有特定的聚類槽位，可以手動同步
  cds@clusters$UMAP$clusters <- seurat_obj$subcelltype

  # 7. (選做) 如果之後要算擬時序，建議執行基本的預處理流程，但跳過降維
  # 因為我們已經手動指定了 UMAP 座標
  #cds <- monocle3::estimate_size_factors(cds)
  # 如果需要保留 Monocle3 的分群邏輯，可以執行：
  # cds <- monocle3::cluster_cells(cds, reduction_method = "UMAP")
  return(cds)
}
cds <- seurat_to_monocle3(CD8T)
cds <- cluster_cells(cds, reduction_method = "UMAP")

# 3. 现在再运行 learn_graph 应该就不会报错了
# cds <- learn_graph(cds, use_partition = TRUE, close_loop = FALSE,
#                             learn_graph_control = list(rann.k = 25,minimal_branch_len = 10))
cds <- learn_graph(cds, use_partition = TRUE, close_loop = T,
                            learn_graph_control = list(rann.k = 25,minimal_branch_len = 10))
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
          label_leaves = FALSE,  label_branch_points = FALSE)
ggsave(file.path(output_dir, 'CD8_Tcell_subcelltype_pseudotime_UMAP2.pdf'),width = 5,height =5)
colData(cds)$pseudotime <- pseudotime(cds)
ggsave(file.path(output_dir, 'CD8_Tcell_subcelltype_pseudotime_UMAP.pdf'),width = 14,height =12)
DimPlot(CD8T,group.by = 'subcelltype',cols=colors_15)
ggsave(file.path(output_dir, 'CD8_Tcell_subcelltype_UMAP.pdf'))
plot_df <- as.data.frame(SummarizedExperiment::colData(cds))
plot_df <- plot_df[base::is.finite(plot_df$pseudotime), ]
plot_df = plot_df[plot_df$Group %in% c("Pre_NR","Pre_CR"),]

# 2. 计算每个 Group 的总平均比例
# 这将用于绘制水平虚线
group_proportions <- plot_df %>%
  dplyr::count(Group) %>%
  dplyr::mutate(proportion = n / sum(n))

# 3. 顶部：百分比堆叠柱状图 + 平均比例虚线
p_top <- ggplot2::ggplot(plot_df, ggplot2::aes(x = pseudotime, fill = Group)) +
  # 百分比堆叠柱状图
  ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(count)), 
                          position = "fill", 
                          bins = 50,          
                          color = "white",    
                          linewidth = 0.1) +
  
  # 添加水平虚线：使用 geom_hline
  # 遍历 group_proportions 计算出的每个组别的比例
  ggplot2::geom_hline(data = group_proportions, 
                      ggplot2::aes(yintercept = proportion, color = Group), # 用 color 映射组别
                      linetype = "dashed",    # 虚线
                      linewidth = 0.8) +      # 虚线粗细
  
  viridis::scale_fill_viridis(discrete = TRUE, option = "D") +
  viridis::scale_color_viridis(discrete = TRUE, option = "D") + # 为虚线也设置相同的色盲友好色
  
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none",
                 plot.margin = ggplot2::margin(0, 5, 0, 5))

# 4. 底部：密度分布图 (与之前相同)
p_bottom <- ggplot2::ggplot(plot_df, ggplot2::aes(x = pseudotime, fill = Group)) +
  ggplot2::geom_density(alpha = 0.6) +
  viridis::scale_fill_viridis(discrete = TRUE, option = "D") +
  ggplot2::theme_classic() +
  ggplot2::labs(x = "Pseudotime (Naive -> Exhaustion)", 
                y = "Density") +
  ggplot2::theme(legend.position = "bottom",
                 plot.margin = ggplot2::margin(0, 5, 5, 5))

# 5. 使用 patchwork 组合
final_plot <- p_top / p_bottom + 
  patchwork::plot_layout(heights = c(1, 3)) & 
  patchwork::plot_annotation(title = "CD8T Compositional Dynamics along Pseudotime")

# 显示图片
base::print(final_plot)
ggsave(file.path(output_dir,"Figures","CD8T_cell_Compositional_Dynamics_along_Pseudotime2.pdf"),plot=final_plot,width=8,height=6)

plot_df <- as.data.frame(colData(cds))
# 过滤掉 pseudotime 为无穷大的细胞（未联通的轨迹）
plot_df <- plot_df[is.finite(plot_df$pseudotime), ]

# 绘制山脊图 (Ridgeline Plot) 或 密度图
ggplot2::ggplot(plot_df, ggplot2::aes(x = pseudotime, y = Group, fill = Group)) +
  # 使用 ggridges 包可以做出更像 Monocle 2 的效果
  ggridges::geom_density_ridges(alpha = 0.7, scale = 1.5, color = "white") +
  # 显式调用色盲友好色：viridis
  viridis::scale_fill_viridis(discrete = TRUE, option = "D") + 
  ggplot2::theme_classic() +
  ggplot2::labs(title = "CD8+ T Exhaustion Progression",
                x = "Pseudotime (Progression to Terminal Exhaustion)",
                y = "Experimental Group") +
  ggplot2::theme(legend.position = "none")

Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=10)
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
                  pull(gene_short_name) %>% as.character()
Track_genes_sig=c('GZMK','HAVCR2','GZMB','KLRG1','LAG3','ENTPD1','NR4A2',"TNFRSF9")
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="subcelltype", 
                              min_expr=0.5, ncol = 2)


plot_data <- data.frame(
  pseudotime = monocle3::pseudotime(cds),
  group = SummarizedExperiment::colData(cds)$Group, # 假設為 R vs NR
  cell_type = SummarizedExperiment::colData(cds)$subcelltype
)

# 過濾出 Tex 細胞進行對比
tex_data <- plot_data[(plot_data$cell_type == "t03_CD8_LAG3 (Tex)") & (plot_data$group %in% c("Pre_NR","Pre_CR")), ]

# 2. 使用色盲友好顏色繪製密度圖
# viridis 顏色：Response(紫色/深藍)，Non-Responder(黃色/綠色)
ggplot2::ggplot(tex_data, ggplot2::aes(x = pseudotime, fill = group)) +
  ggplot2::geom_density(alpha = 0.6) +
  viridis::scale_fill_viridis(discrete = TRUE, option = "D") + 
  ggplot2::theme_classic() +
  ggplot2::labs(title = "Distribution of Tex Cells along Pseudotime",
                x = "Pseudotime (State of Exhaustion)",
                y = "Density")
keep_cells <- which(colData(cds)$Group %in% c("Pre_NR","Pre_CR"))

cds_pre <- cds[, keep_cells]

plot_cells(cds,
          color_cells_by = "pseudotime",
          label_cell_groups = FALSE,
          label_leaves = FALSE, # 顯示葉節點
          label_branch_points = FALSE,
          graph_label_size = 4) +
  viridis::scale_color_viridis(discrete = FALSE) +
  ggplot2::facet_wrap(~Group)


plot_df <- as.data.frame(SummarizedExperiment::colData(cds)) %>%
  dplyr::filter(base::is.finite(pseudotime)) %>%
  dplyr::mutate(
    subcelltype = base::sub(" \\(.*\\)", "", as.character(subcelltype)),
    Group = base::as.factor(Group)
  )

# 2. 定义需要比较的组别（针对你的 Pre_NR 和 Pre_CR）
# 确保这些名字与你 Group 列中的实际取值完全一致
my_comparisons <- list(c("Pre_NR", "Pre_CR"))

# 3. 绘制分面小提琴图
ggplot2::ggplot(plot_df, ggplot2::aes(x = Group, y = pseudotime, fill = Group)) +
  # 绘制小提琴图，trim = FALSE 保留完整分布形状
  ggplot2::geom_violin(trim = FALSE, alpha = 0.7, color = "black", linewidth = 0.3) +
  # 叠加箱线图，方便查看中位数和四分位数
  ggplot2::geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = "black") +
  
  # 分面：按 subcelltype 分开，scales = "free_y" 可以让不同亚群的 y 轴独立（如果需要）
  ggplot2::facet_wrap(~subcelltype, scales = "free_y") +
  
  # 统计分析：添加 p-value
  # label 可以选 "p.format" (显示具体数值) 或 "p.signif" (显示星号)
  ggpubr::stat_compare_means(comparisons = my_comparisons, 
                             method = "wilcox.test", 
                             label = "p.signif") +
  scale_fill_manual(values = colors_15) +
  

  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    panel.grid = ggplot2::element_blank(), # 去掉网格线更符合高水平杂志要求
    strip.background = ggplot2::element_rect(fill = "white"), # 分面标签背景转白
    strip.text = ggplot2::element_text(face = "bold")
  ) +
  ggplot2::labs(x = "Patient Group", y = "Pseudotime (Exhaustion Score)")
ggsave(file.path(output_dir,"Figures","CD8T_cell_Compositional_Dynamics_along_Pseudotime_violin.pdf"))






