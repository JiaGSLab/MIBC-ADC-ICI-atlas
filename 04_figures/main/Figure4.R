# ==================================================
# Script: Figure 4 plotting
# Author: Wenhui Chen
# Purpose: Generates spatial and transcriptomic plots used in Figure 4.
# Input: Main Seurat or summary objects already generated in the project, plus any file paths referenced in the script.
# Output: Intermediate result files and/or figure files written to output_dir.
# Note: Review file paths before running. Keep original analysis logic unchanged.
# ==================================================
########################
#Fig. 4F quantification of KRT5+ and UPK2+ cell proportions
df = read.csv('/home/data/chenwenhui/scRNAseq/GPU_IT/QC/Figures/Fig/Fig4F.csv')
# 1. Calculate proportions and select the required columns
plot_df <- df %>%
  dplyr::mutate(
    KRT5_Proportion = (KRT5p / total.cell) * 100,
    UPK2_Proportion = (UPK2p / total.cell) * 100
  ) %>%
  dplyr::select(Group, KRT5_Proportion, UPK2_Proportion)

# 2. Convert to long format for faceted plotting
plot_df_long <- plot_df %>%
  tidyr::pivot_longer(
    cols = c("KRT5_Proportion", "UPK2_Proportion"),
    names_to = "CellType",
    values_to = "Percentage"
  )

# 3. Ensure Group is a factor and place CR first
plot_df_long$Group <- factor(plot_df_long$Group, levels = c("CR", "NR"))
# Define a color-blind-friendly palette (CR: vermilion, NR: sky blue)
cb_palette <- c("CR" = "#E41A1C", "NR" = "#377EB8")

p <- ggpubr::ggbarplot(
  plot_df_long, 
  x = "Group", 
  y = "Percentage",
  fill = "Group", 
  palette = cb_palette,
  add = "mean_se",          # Add mean and standard error
  add.params = list(width = 0.5),
  position = position_dodge(0.8),
  facet.by = "CellType",     # Facet by cell type
  scales = "free_y",
  error.plot = "upper_errorbar"
) +
  # Overlay raw data points to show heterogeneity among FOVs
  ggplot2::geom_jitter(
    width = 0.2, 
    size = 1.5,
    shape = 16, 
    color = "black"
  ) +
  # Add statistical testing (Wilcoxon rank-sum test)
  ggpubr::stat_compare_means(
    method = "wilcox.test", 
    label = "p.format", 
    label.x = 1.5
  ) +
  ggplot2::theme_classic() +
  ggplot2::labs(
    title = "Comparison of Cell Proportions within FOVs",
    y = "Positive Cell Proportion (%)",
    x = ""
  ) +
  ggplot2::theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "none" # because the group labels are already shown on the x-axis
  )

print(p)
ggsave(filename = file.path(output_dir,"Figures","Fig","Fig4F_mIHC_cell_proportion.pdf"),width=3,height = 2.5,plot=p)


####################################################
#Visium HD scoring
HD.combined.obj = qs_read('/home/data/chenwenhui/visiumHD/ST.combined.obj.qs')
HD.combined.obj$Outcome = ifelse(HD.combined.obj$sample=='ZZJ_8916','CR',HD.combined.obj$Outcome)
HD.combined.obj = subset(HD.combined.obj,subset=Therapy=="IT")
HD.combined.obj = NormalizeData(HD.combined.obj)
HD.combined.obj$pct_counts_mt <- PercentageFeatureSet(HD.combined.obj, pattern = "^MT-")
VlnPlot(HD.combined.obj,features = "nCount_Spatial",group.by="sample_id",pt.size=0,raster = T)
VlnPlot(HD.combined.obj,features = "nFeature_Spatial",group.by="sample_id",pt.size=0,raster = T)
VlnPlot(HD.combined.obj,features = "pct_counts_mt",group.by="sample_id",pt.size=0,raster = T)
HD.combined.obj = subset(HD.combined.obj,subset=((nCount_Spatial>50) &
                                                (nCount_Spatial<4000) &
                                                (nFeature_Spatial>50) &
                                                (nFeature_Spatial<2000) &
                                                (pct_counts_mt<25)))
HD.combined.obj = subset(HD.combined.obj,subset=orig.ident %in% c('slide02','slide03'))
MP.genelist = qs_read("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/Final_MP_Top30_Signatures.qs")
HD.combined.obj <- AddModuleScore(HD.combined.obj, features = MP.genelist,name="MPscore")

new_names = c(
  "Cycle_G2M", "Surface_Trafficking", "RNA_Splicing", "Basal_pEMT", 
  "Stress_IEG", "IFN_Response", "Hypoxia", "Metabolic", 
  "Uro_Differentiation", "Ep_Inflammation"
)
HD.combined.obj@meta.data = HD.combined.obj@meta.data[,setdiff(colnames(HD.combined.obj@meta.data),new_names)]
for(i in 1:length(new_names)){
  colnames(HD.combined.obj@meta.data)[which(colnames(HD.combined.obj@meta.data) == paste0("MPscore", i))] <- new_names[i]
}
#Select spots mapped as epithelial by RCTD
epi.clusters = c('e01_Epi_PSCA','e02_Epi_MKI67','e03_Epi_KRT5','e04_Epi_NDRG1','e05_Epi_FOS',
                'e06_Epi_S100A9','e07_Epi_KRT13','e08_Epi_CXCL10','e09_Epi_PIGR','e10_Epi_MMP7',"e11_Epi_KRT16")
HD.epi = subset(HD.combined.obj,subset=((first_type %in% epi.clusters) |
                                        (second_type %in% epi.clusters)))
HD.epi@meta.data = HD.epi@meta.data[,setdiff(colnames(HD.epi@meta.data),new_names)]
HD.epi <- AddModuleScore(HD.epi, features = MP.genelist,name="MPscore")
new_names = c(
  "Cycle_G2M", "Surface_Trafficking", "RNA_Splicing", "Basal_pEMT", 
  "Stress_IEG", "IFN_Response", "Hypoxia", "Metabolic", 
  "Uro_Differentiation", "Ep_Inflammation"
)
for(i in 1:length(new_names)){
  colnames(HD.epi@meta.data)[which(colnames(HD.epi@meta.data) == paste0("MPscore", i))] <- new_names[i]
}
hd_sample_avg <- HD.epi@meta.data %>%
  dplyr::group_by(sample_id, Outcome) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(new_names), mean, na.rm = TRUE), .groups = "drop")

# 2. Convert to long format for plotting
hd_long <- hd_sample_avg %>%
  tidyr::pivot_longer(cols = dplyr::all_of(new_names), names_to = "MP", values_to = "Score")

# 3. Set group order and separate main vs supplementary MPs
hd_long$Outcome <- factor(hd_long$Outcome, levels = c("CR", "NR"))
main_mps <- c("Basal_pEMT", "Uro_Differentiation")
supp_mps <- setdiff(new_names, main_mps)

# Define color-blind-friendly colors
group_cols <- c("CR" = "#E41A1C", "NR" = "#377EB8")
# Filter data for the main figure
df_main <- hd_long %>% dplyr::filter(MP %in% main_mps)
library(coin)

calc_exact_p <- function(df, mp_name, alternative_type) {
  sub_df <- df %>% dplyr::filter(MP == mp_name)
  sub_df$Outcome <- factor(sub_df$Outcome, levels = c("CR", "NR"))
  
  test_res <- coin::independence_test(
    Score ~ Outcome, 
    data = sub_df, 
    distribution = "exact", 
    alternative = alternative_type
  )
  return(as.numeric(coin::pvalue(test_res)))
}

# 2. Calculate P values for the core MPs
# For Uro_Diff, CR is expected to be higher than NR -> alternative = "greater"
p_uro <- calc_exact_p(hd_long, "Uro_Differentiation", "greater")
# For Basal_pEMT, NR is expected to be higher than CR -> alternative = "less" (because Outcome is ordered as CR, NR)
p_basal <- calc_exact_p(hd_long, "Basal_pEMT", "less")

# Construct a data frame for annotations
p_label_df <- data.frame(
  MP = c("Uro_Differentiation", "Basal_pEMT"),
  p_val = c(p_uro, p_basal),
  y_pos = c(max(hd_long$Score[hd_long$MP == "Uro_Differentiation"]) * 1.1, 
            max(hd_long$Score[hd_long$MP == "Basal_pEMT"]) * 1.1)
) %>%
  mutate(label = sprintf("Permutation p = %.4f", p_val))


pdf(file.path(output_dir,"Figures","Fig","Fig4G_MP4_MP9_visiumHD_static1.pdf"), width=3,height = 2.5)

# Modify the main plotting function to include the alternative parameter
p_main <- ggpubr::ggboxplot(
  df_main, x = "Outcome", y = "Score",
  fill = "Outcome", palette = group_cols,
  add = "jitter", 
  add.params = list(shape = 16, size = 1.5, alpha = 1),
  facet.by = "MP", scales = "free_y", ncol = 2
) +
  # Specify different one-sided directions for different MPs
  # Note: this should be adjusted for each specific MP, or set directly in stat_compare_means
  ggpubr::stat_compare_means(
    method = "wilcox.test", 
    method.args = list(alternative = "greater"), # assuming CR > NR
    label = "p.format", 
    label.x = 1.5
  ) +
  ggplot2::theme_classic() +
  ggplot2::labs(subtitle = "One-sided Wilcoxon Test (Alternative: Greater)") +
  NoLegend()


# Note: for Basal_pEMT, if NR > CR, you should set alternative = "less" when the x-axis order is CR, NR

print(p_main)
dev.off()

# Filter data for supplementary plots
df_supp <- hd_long %>% dplyr::filter(MP %in% supp_mps)

pdf(file.path(output_dir,"Figures","FigS","FigS4_Supplementary_Figure_VisiumHD_MPs.pdf"), width = 6, height = 5)

p_supp <- ggpubr::ggboxplot(
  df_supp, x = "Outcome", y = "Score",
  fill = "Outcome", palette = group_cols,
  add = "jitter", 
  facet.by = "MP", scales = "free_y", ncol = 4
) +
  ggpubr::stat_compare_means(method = "wilcox.test", label = "p.format") +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(size = 11),
    legend.position = "bottom"
  ) +
  ggplot2::labs(
    title = "Supplementary Meta-Programs across 12 Samples",
    y = "Mean Score per Sample", x = ""
  )

print(p_supp)
dev.off()

############################
##
#Spatial mapping of cell types from RCTD
RCTD.obj = subset(HD.combined.obj,subset=first_type=="unknown",invert =T)
RCTD.obj$celltype = ifelse(RCTD.obj$first_type %in% epi.clusters,"Epithelial",RCTD.obj$first_type)
RCTD.obj$celltype = factor(RCTD.obj$celltype,levels = c("Epithelial","Endothelial","Pericyte","SMC","Fibroblast","Macrophage",
                                                        "Monocyte","Mast cell","DC","Neutrophil","B","Plasma","CD8 T","CD4 T",
                                                        "Treg","γδT"))
slides = c('slide01','slide02','slide03')
for(slide in slides){
  RCTD.obj.slide = subset(RCTD.obj,subset=orig.ident==slide)
  p = SpatialPlot(RCTD.obj.slide,group.by = 'celltype',pt.size.factor = 4,image.alpha = 0.2) +
    scale_fill_manual(values = colors_90) +
    theme(legend.position = "none",aspect.ratio = 1)
  ggsave(file.path(output_dir,"Figures","Fig",paste0("Fig4G_",slide,"_celltype_spatialplot.png")),plot=p,width=4,height = 4,dpi=600,bg = "transparent")
}
p = SpatialPlot(RCTD.obj.slide,group.by = 'celltype',pt.size.factor = 4,image.alpha = 0.2) +
    scale_fill_manual(values = colors_90) +
    theme(legend.position = "right",aspect.ratio = 1)
ggsave(file.path(output_dir,"Figures","Fig",paste0("Fig4G_",slide,"_celltype_spatialplot.pdf")),plot=p,width=4,height = 4)

#Visium HD spatial FeaturePlot
# 1. Prepare parameters
# Get all unique sample IDs
all_samples <- unique(HD.combined.obj$sample_id)

for (samp_id in all_samples) {
  sub_obj <- subset(HD.combined.obj, subset = sample_id == samp_id)
  group_info <- unique(sub_obj$Outcome)
  p = SpatialFeaturePlot(sub_obj,features = 'Uro_Differentiation',pt.size.factor = 5,image.alpha = 0.2) +
    scale_fill_viridis_c(option = "magma", limits = c(-0.1, 0.2)) + 
    # scale_fill_gradient2(
    #   low = "#56B4E9", mid = "white", high = "#d73027", 
    #   midpoint = 0, 
    #   limits = c(-0.2, 0.2) 
    #   ) +
    theme(legend.position = "none")
  ggsave(file.path(output_dir,"Figures","Fig",paste0("Fig4G_",samp_id,"_",group_info,"_Uro_Differentiation_spatialplot.png")),plot=p,width=2,height = 2,dpi=600,bg = "transparent")
  p = SpatialFeaturePlot(sub_obj,features = 'Basal_pEMT',pt.size.factor = 5,image.alpha = 0.2) +
    scale_fill_viridis_c(option = "magma", limits = c(-0.1, 1.5)) + 
    # scale_fill_gradient2(
    #   low = "#56B4E9", mid = "white", high = "#d73027", 
    #   midpoint = 0, 
    #   limits = c(-0.2, 2) 
    #   ) +
    theme(legend.position = "none")
  ggsave(file.path(output_dir,"Figures","Fig",paste0("Fig4G_",samp_id,"_",group_info,"_Basal_pEMT_spatialplot.png")),plot=p,width=2,height = 2,dpi=600,bg = "transparent")
}
p = SpatialFeaturePlot(sub_obj,features = 'Basal_pEMT',pt.size.factor = 2,image.alpha = 0.2) +
    scale_fill_viridis_c(option = "magma", limits = c(-0.1, 1.5)) + 
    # scale_fill_gradient2(
    #   low = "#56B4E9", mid = "white", high = "#d73027", 
    #   midpoint = 0, 
    #   limits = c(-0.2, 2) 
    #   ) +
    theme()
  ggsave(file.path(output_dir,"Figures","Fig",paste0("Fig4G_",samp_id,"_",group_info,"_Basal_pEMT_spatialplot.pdf")),plot=p,width=4,height = 4)


#######################################################
# 1. Extract plotting data (kept unchanged)
plot_df <- data.frame(
  MP4 = malignant.obj$Basal_pEMT,
  MP9 = malignant.obj$Uro_Differentiation,
  ERBB2 = Seurat::FetchData(malignant.obj, vars = "ERBB2")[,1]
)

# 2. Calculate the global slope
fit4 <- stats::lm(ERBB2 ~ MP4, data = plot_df)
slope4 <- stats::coef(fit4)[2]

fit9 <- stats::lm(ERBB2 ~ MP9, data = plot_df)
slope9 <- stats::coef(fit9)[2]
p_mp4_global <- ggplot2::ggplot(plot_df, ggplot2::aes(x = MP4, y = ERBB2)) +
  # Rasterized density layer
  ggplot2::geom_bin2d(bins = 100) + 
  ggplot2::scale_fill_viridis_c(option = "magma", trans = "log10") + 
  # Fitted line: cyan dashed line
  ggplot2::geom_smooth(method = "lm", color = "cyan", linewidth = 1, linetype = "dashed") +
  # Annotate R and P values
  ggpubr::stat_cor(method = "pearson", color = "black", size = 4, 
                   label.x.npc = "left", label.y.npc = "top") + 
  # Annotate slope and manually lower the y-position to avoid overlap
  ggplot2::annotate("text", x = min(plot_df$MP4), y = max(plot_df$ERBB2) * 0.85, 
                    label = paste("Slope =", round(slope4, 3)), 
                    color = "black", hjust = 0, fontface = "italic", size = 4) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "Global Correlation: MP4 vs ERBB2",
                x = "MP4 (Basal_pEMT) Score", y = "ERBB2 Expression") +
  ggplot2::theme(panel.grid = element_blank(),
                 plot.title = ggplot2::element_text(face = "bold"),aspect.ratio = 1)

# 4. Plot MP9 vs ERBB2 (positive correlation)
p_mp9_global <- ggplot2::ggplot(plot_df, ggplot2::aes(x = MP9, y = ERBB2)) +
  ggplot2::geom_bin2d(bins = 100) + 
  ggplot2::scale_fill_viridis_c(option = "magma", trans = "log10") + 
  ggplot2::geom_smooth(method = "lm", color = "cyan", linewidth = 1, linetype = "dashed") +
  ggpubr::stat_cor(method = "pearson", color = "black", size = 4, 
                   label.x.npc = "left", label.y.npc = "top") + 
  ggplot2::annotate("text", x = min(plot_df$MP9), y = max(plot_df$ERBB2) * 0.85, 
                    label = paste("Slope =", round(slope9, 3)), 
                    color = "black", hjust = 0, fontface = "italic", size = 4) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "Global Correlation: MP9 vs ERBB2",
                x = "MP9 (Uro_Differentiation) Score", y = "ERBB2 Expression") +
  ggplot2::theme(panel.grid = element_blank(),
                 plot.title = ggplot2::element_text(face = "bold"),aspect.ratio = 1)

# Combined display
combined_plot <- p_mp4_global / p_mp9_global
ggplot2::ggsave(
  filename = file.path(output_dir, "Figures","Fig", "Fig4I_MP_vs_ERBB2_Global_Density.pdf"),
  plot = combined_plot, 
  width = 4, 
  height = 8,
  device = grDevices::cairo_pdf
)


#########################################################
#Fig4H Use TCGA bulk RNA-seq to test whether pEMT_Basal score is negatively correlated with HER2 expression
library(GSVA)
library(SeuratExtend)
bulk.TCGA = readRDS("/home/data/chenwenhui/TCGA/data/BLCA/01_BLCA_tpm_filt.rds")
MP.genelist = qs_read("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/Final_MP_Top30_Signatures.qs")
library(robustbase)
library(org.Hs.eg.db)

MP4.gene = MP.genelist$MP4
MP9.gene = MP.genelist$MP9
ensembleid = rownames(bulk.TCGA)
ensembleid = gsub("\\..*","",ensembleid)
symbol = EnsemblToGenesymbol(ensembleid,spe="human",keep.seq=T)
bulk.log <- log2(bulk.TCGA + 1)

# 2. Extract log-transformed ERBB2 expression
# Handle duplicated symbols by averaging all corresponding rows
erbb2_idx <- which(symbol == "ERBB2")
if(length(erbb2_idx) == 0) stop("ERBB2 not found in symbol vector!")
erbb2_expr <- colMeans(bulk.log[erbb2_idx, , drop = FALSE])

# 3. Calculate MP4 and MP9 signature scores based on mean log-transformed expression
# Only calculate scores using genes present in the matrix
mp4_score <- colMeans(bulk.log[symbol %in% MP4.gene, , drop = FALSE])
mp9_score <- colMeans(bulk.log[symbol %in% MP9.gene, , drop = FALSE])

# 4. Construct the plotting data frame
plot_df <- data.frame(
  Sample = colnames(bulk.log),
  ERBB2 = as.numeric(erbb2_expr),
  MP4_Score = as.numeric(mp4_score),
  MP9_Score = as.numeric(mp9_score)
)
threshold_25 <- stats::quantile(plot_df$ERBB2, 0.25)

# 2. Filter out the lowest 25% of samples and retain the HER2-active group
plot_df <- plot_df[plot_df$ERBB2 > threshold_25, ]

# 5. Plot correlations after log transformation
# MP4 correlation plot
p1 <- ggplot(plot_df, aes(x = MP4_Score, y = ERBB2)) +
  geom_point(color = "#A52A2A", alpha = 0.5, size = 1.5,shape=16) +
  geom_smooth(method = "glm", color = "black", fill = "gray80", alpha = 0.2) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_classic() +
  labs(title = "",
       x = "MP4 Signature Score",
       y = "ERBB2 Expression") +
  theme(aspect.ratio=1)

# MP9 correlation plot
p2 <- ggplot(plot_df, aes(x = MP9_Score, y = ERBB2)) +
  geom_point(color = "#053E7A", alpha = 0.5, size = 1.5,shape=16) +
  geom_smooth(method = "glm", color = "black", fill = "gray80", alpha = 0.2) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_classic() +
  labs(title = "",
       x = "MP9 Signature Score",
       y = " ERBB2 Expression") + 
  theme(aspect.ratio=1)


# Combined display
cowplot::plot_grid(p1, p2, ncol = 2)


gene_sets <- list(
  MP4_Basal_pEMT = MP4.gene,
  MP9_Differentiation = MP9.gene
)
clean_ids <- gsub("\\..*", "", rownames(bulk.TCGA))

# 2. Build the Ensembl-to-symbol mapping table
gene_map <- AnnotationDbi::select(org.Hs.eg.db, 
                                  keys = clean_ids, 
                                  columns = "SYMBOL", 
                                  keytype = "ENSEMBL") %>%
  dplyr::filter(!is.na(SYMBOL)) # Remove unmatched entries

# 3. Convert the matrix and handle duplicated symbols by averaging
# This step converts the Ensembl-based matrix to a symbol-based matrix and merges duplicates
bulk_symbol <- bulk.TCGA %>%
  as.data.frame() %>%
  mutate(ENSEMBL = clean_ids) %>%
  inner_join(gene_map, by = "ENSEMBL") %>%
  dplyr::select(-ENSEMBL) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(across(everything(), mean)) %>% # Average duplicated symbols
  column_to_rownames("SYMBOL") %>%
  as.matrix()
ssgsea_params <- ssgseaParam(bulk_symbol, gene_sets) 

# 2. Calculate scores
ssgsea_scores <- gsva(ssgsea_params)

# 3. Extract and filter the top 75%
bulk_log <- log2(bulk_symbol + 1)

# 2. Extract log-transformed HER2 (ERBB2) expression
# The extracted result is a vector at this stage
erbb2_log <- bulk_log["ERBB2", ]

# 3. Reconstruct the analysis data frame
plot_df_log <- data.frame(
  Sample = colnames(bulk_log),
  ERBB2 = as.numeric(erbb2_log),
  MP4_ssGSEA = as.numeric(ssgsea_scores["MP4_Basal_pEMT", ]),
  MP9_ssGSEA = as.numeric(ssgsea_scores["MP9_Differentiation", ])
)

# 4. Apply the 25% filtering step again, based on log-transformed values
threshold_25 <- stats::quantile(plot_df_log$ERBB2, 0.25)
plot_df_final <- plot_df_log[plot_df_log$ERBB2 > threshold_25, ]
#plot_df_final = plot_df_log

p1 <- ggplot(plot_df_final, aes(x = MP4_ssGSEA, y = ERBB2)) +
  geom_point(color = "#A52A2A", alpha = 0.5, size = 1.5,shape=16) +
  geom_smooth(method = "glm", color = "black", fill = "gray80", alpha = 0.2) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_classic() +
  labs(title = "",
       x = "MP4 Signature Score",
       y = "ERBB2 Expression") +
  theme(aspect.ratio=1)

# MP9 correlation plot
p2 <- ggplot(plot_df_final, aes(x = MP9_ssGSEA, y = ERBB2)) +
  geom_point(color = "#053E7A", alpha = 0.5, size = 1.5,shape=16) +
  geom_smooth(method = "glm", color = "black", fill = "gray80", alpha = 0.2) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_classic() +
  labs(title = "",
       x = "MP9 Signature Score",
       y = " ERBB2 Expression") + 
  theme(aspect.ratio=1)

cowplot::plot_grid(p1, p2, ncol = 2)

ggsave(file.path(output_dir,"Figures","Fig","Fig4J_TCGA_BLCA_bulk.pdf"),width=6,height=3)

########################################################
#Imvigor210shuj
patient.data = fread('/home/data/chenwenhui/scRNAseq/public/blca_iatlas_imvigor210_2017/data_clinical_patient.txt',
                          header=T,skip=4)
sample.data = fread('/home/data/chenwenhui/scRNAseq/public/blca_iatlas_imvigor210_2017/data_clinical_sample.txt',
                          header=T,skip=4)
View(sample.data)
exp_tpm <- fread("/home/data/chenwenhui/scRNAseq/public/blca_iatlas_imvigor210_2017/data_mrna_seq_tpm.txt", header=T,data.table = F)
rownames(exp_tpm) = exp_tpm$Hugo_Symbol
exp_tpm = exp_tpm[,2:ncol(exp_tpm)]
View(exp_tpm)
exp_tpm =as.matrix(exp_tpm)
all(colnames(exp_tpm)==sample.data$PATIENT_ID)

MP.genelist = qs_read("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/Final_MP_Top30_Signatures.qs")
MP4.gene = MP.genelist$MP4
MP9.gene = MP.genelist$MP9
mp_list <- list(
  MP4_Basal_pEMT = MP4.gene,
  MP9_Differentiation = MP9.gene
)
ssgsea_param <- ssgseaParam(
  exprData = exp_tpm, 
  geneSets = mp_list, 
  minSize = 5,           # Filter out gene sets with too few members
  maxSize = 500, 
  normalize = TRUE       # Internal normalization is usually recommended for ssGSEA
)
mp_scores_res <- gsva(ssgsea_param)
mp_scores_df <- as.data.frame(t(mp_scores_res))
colnames(mp_scores_df) <- c("MP4_Basal_pEMT", "MP9_Uro_Diff")
mp_scores_df$Sample_ID <- rownames(mp_scores_df)
View(mp_scores_df)
her2_exp <- data.frame(Sample_ID = colnames(exp_tpm), 
                       HER2 = as.numeric(exp_tpm["ERBB2", ]))
plot_df <- left_join(mp_scores_df, her2_exp, by = "Sample_ID")
p1 <- ggplot(plot_df, aes(x = MP4_Basal_pEMT, y = HER2)) +
  geom_point(color = "#A52A2A", alpha = 0.5, size = 1.5,shape=16) +
  geom_smooth(method = "glm", color = "black", fill = "gray80", alpha = 0.2) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_classic() +
  labs(title = "",
       x = "MP4 Signature Score",
       y = "ERBB2 Expression") +
  theme(aspect.ratio=1)

# MP9 correlation plot
p2 <- ggplot(plot_df, aes(x = MP9_Uro_Diff, y = HER2)) +
  geom_point(color = "#053E7A", alpha = 0.5, size = 1.5,shape=16) +
  geom_smooth(method = "glm", color = "black", fill = "gray80", alpha = 0.2) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_classic() +
  labs(title = "",
       x = "MP9 Signature Score",
       y = " ERBB2 Expression") + 
  theme(aspect.ratio=1)

cowplot::plot_grid(p1, p2, ncol = 2)

ggsave(file.path(output_dir,"Figures","FigS","Fig4J_IMvigor210_MPvsHER2.pdf"),width=6,height=3)

#################################################
#Relationship of MP4 and MP9 scores with response category
sample.data$SAMPLE_ID = as.character(sample.data$SAMPLE_ID)
plot_data <- mp_scores_df %>%
  dplyr::inner_join(sample.data, by = c("Sample_ID" = "SAMPLE_ID")) %>%
  # Remove samples without response information
  dplyr::filter(!is.na(RESPONSE) & RESPONSE != "") %>%
  # Convert RESPONSE to a factor and order it as CR -> PR -> SD -> PD
  dplyr::mutate(RESPONSE = factor(RESPONSE, levels = c("Complete Response", 
                                                      "Partial Response", 
                                                      "Stable Disease", 
                                                      "Progressive Disease")))

# Check the number of observations after cleaning
table(plot_data$RESPONSE)
# Define colors using a color-blind-friendly palette according to response logic
# You could also use built-in ggpubr palettes such as "npg" or "jco" here
res_cols <- c("#207122", "#62A765", "#E1C239", "#D55E00") # green-to-red gradient

# 2. Plot MP4 (Basal_pEMT)
p_mp4 <- ggpubr::ggboxplot(
  plot_data, x = "RESPONSE", y = "MP4_Basal_pEMT",
  fill = "RESPONSE", palette = res_cols,
  add = "jitter", # Overlay scatter points to show the distribution
  add.params = list(size = 1.5, alpha = 0.4, shape = 16) # solid circles that are friendly for AI editing
) +
  # Key point: use complete response as the reference group for comparison
  ggpubr::stat_compare_means(
    method = "wilcox.test", 
    label = "p.format", # display significance labels; "p.format" can also be used
    ref.group = "Complete Response"
  ) +
  ggplot2::labs(title = "MP4 (Basal_pEMT) vs Atezolizumab Response",
                y = "ssGSEA Score", x = "") +
  ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 legend.position = "none")

# 3. Plot MP9 (Uro_Diff)
p_mp9 <- ggpubr::ggboxplot(
  plot_data, x = "RESPONSE", y = "MP9_Uro_Diff",
  fill = "RESPONSE", palette = res_cols,
  add = "jitter", 
  add.params = list(size = 1.5, alpha = 0.4, shape = 16)
) +
  ggpubr::stat_compare_means(
    method = "wilcox.test", 
    label = "p.format",
    ref.group = "Complete Response"
  ) +
  ggplot2::labs(title = "MP9 (Uro_Diff) vs Atezolizumab Response",
                y = "ssGSEA Score", x = "") +
  ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 legend.position = "none")

# Combined display
library(patchwork)
p_final <- p_mp4 | p_mp9
print(p_final)
ggsave(file.path(output_dir,"Figures","FigS","Fig4J_IMvigor210_MPwithRESPONSE.pdf"),width=6,height=3,plot=p_final)



#####################################################
bulk.TCGA = readRDS("/home/data/chenwenhui/TCGA/data/ESCA/01_ESCA_tpm_filt.rds")
MP.genelist = qs_read("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/Final_MP_Top30_Signatures.qs")

MP4.gene = MP.genelist$MP4
MP9.gene = MP.genelist$MP9
ensembleid = rownames(bulk.TCGA)
ensembleid = gsub("\\..*","",ensembleid)
symbol = EnsemblToGenesymbol(ensembleid,spe="human",keep.seq=T)
bulk.log <- log2(bulk.TCGA + 1)

# 2. Extract log-transformed ERBB2 expression
# Handle duplicated symbols by averaging all corresponding rows
erbb2_idx <- which(symbol == "ERBB2")
if(length(erbb2_idx) == 0) stop("ERBB2 not found in symbol vector!")
erbb2_expr <- colMeans(bulk.log[erbb2_idx, , drop = FALSE])

# 3. Calculate MP4 and MP9 signature scores based on mean log-transformed expression
# Only calculate scores using genes present in the matrix
mp4_score <- colMeans(bulk.log[symbol %in% MP4.gene, , drop = FALSE])
mp9_score <- colMeans(bulk.log[symbol %in% MP9.gene, , drop = FALSE])

# 4. Construct the plotting data frame
plot_df <- data.frame(
  Sample = colnames(bulk.log),
  ERBB2 = as.numeric(erbb2_expr),
  MP4_Score = as.numeric(mp4_score),
  MP9_Score = as.numeric(mp9_score)
)
threshold_25 <- stats::quantile(plot_df$ERBB2, 0.25)

# 2. Filter out the lowest 25% of samples and retain the HER2-active group
plot_df <- plot_df[plot_df$ERBB2 > threshold_25, ]

# 5. Plot correlations after log transformation
# MP4 correlation plot
p1 <- ggplot(plot_df, aes(x = MP4_Score, y = ERBB2)) +
  geom_point(color = "#A52A2A", alpha = 0.5, size = 1.5,shape=16) +
  geom_smooth(method = "glm", color = "black", fill = "gray80", alpha = 0.2) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_classic() +
  labs(title = "",
       x = "MP4 Signature Score",
       y = "ERBB2 Expression") +
  theme(aspect.ratio=1)

# MP9 correlation plot
p2 <- ggplot(plot_df, aes(x = MP9_Score, y = ERBB2)) +
  geom_point(color = "#053E7A", alpha = 0.5, size = 1.5,shape=16) +
  geom_smooth(method = "glm", color = "black", fill = "gray80", alpha = 0.2) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_classic() +
  labs(title = "",
       x = "MP9 Signature Score",
       y = " ERBB2 Expression") + 
  theme(aspect.ratio=1)


# Combined display
cowplot::plot_grid(p1, p2, ncol = 2)


gene_sets <- list(
  MP4_Basal_pEMT = MP4.gene,
  MP9_Differentiation = MP9.gene
)
clean_ids <- gsub("\\..*", "", rownames(bulk.TCGA))

# 2. Build the Ensembl-to-symbol mapping table
gene_map <- AnnotationDbi::select(org.Hs.eg.db, 
                                  keys = clean_ids, 
                                  columns = "SYMBOL", 
                                  keytype = "ENSEMBL") %>%
  dplyr::filter(!is.na(SYMBOL)) # Remove unmatched entries

# 3. Convert the matrix and handle duplicated symbols by averaging
# This step converts the Ensembl-based matrix to a symbol-based matrix and merges duplicates
bulk_symbol <- bulk.TCGA %>%
  as.data.frame() %>%
  mutate(ENSEMBL = clean_ids) %>%
  inner_join(gene_map, by = "ENSEMBL") %>%
  dplyr::select(-ENSEMBL) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(across(everything(), mean)) %>% # Average duplicated symbols
  column_to_rownames("SYMBOL") %>%
  as.matrix()
ssgsea_params <- ssgseaParam(bulk_symbol, gene_sets) 

# 2. Calculate scores
ssgsea_scores <- gsva(ssgsea_params)

# 3. Extract and filter the top 75%
bulk_log <- log2(bulk_symbol + 1)

# 2. Extract log-transformed HER2 (ERBB2) expression
# The extracted result is a vector at this stage
erbb2_log <- bulk_log["ERBB2", ]

# 3. Reconstruct the analysis data frame
plot_df_log <- data.frame(
  Sample = colnames(bulk_log),
  ERBB2 = as.numeric(erbb2_log),
  MP4_ssGSEA = as.numeric(ssgsea_scores["MP4_Basal_pEMT", ]),
  MP9_ssGSEA = as.numeric(ssgsea_scores["MP9_Differentiation", ])
)

# 4. Apply the 25% filtering step again, based on log-transformed values
# threshold_25 <- stats::quantile(plot_df_log$ERBB2, 0.25)
# plot_df_final <- plot_df_log[plot_df_log$ERBB2 > threshold_25, ]
plot_df_final = plot_df_log

p1 <- ggplot(plot_df_final, aes(x = MP4_ssGSEA, y = ERBB2)) +
  geom_point(color = "#A52A2A", alpha = 0.5, size = 1.5,shape=16) +
  geom_smooth(method = "glm", color = "black", fill = "gray80", alpha = 0.2) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_classic() +
  labs(title = "",
       x = "MP4 Signature Score",
       y = "ERBB2 Expression") +
  theme(aspect.ratio=1)

# MP9 correlation plot
p2 <- ggplot(plot_df_final, aes(x = MP9_ssGSEA, y = ERBB2)) +
  geom_point(color = "#053E7A", alpha = 0.5, size = 1.5,shape=16) +
  geom_smooth(method = "glm", color = "black", fill = "gray80", alpha = 0.2) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_classic() +
  labs(title = "",
       x = "MP9 Signature Score",
       y = " ERBB2 Expression") + 
  theme(aspect.ratio=1)

cowplot::plot_grid(p1, p2, ncol = 2)

ggsave(file.path(output_dir,"Figures","FigS","Fig4J_TCGA_ESCA_bulk.pdf"),width=6,height=3)


##############################################################
bulk.TCGA = readRDS("/home/data/chenwenhui/TCGA/data/STAD/01_STAD_tpm_filt.rds")
MP.genelist = qs_read("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/Final_MP_Top30_Signatures.qs")

MP4.gene = MP.genelist$MP4
MP9.gene = MP.genelist$MP9
ensembleid = rownames(bulk.TCGA)
ensembleid = gsub("\\..*","",ensembleid)
symbol = EnsemblToGenesymbol(ensembleid,spe="human",keep.seq=T)
bulk.log <- log2(bulk.TCGA + 1)

# 2. Extract log-transformed ERBB2 expression
# Handle duplicated symbols by averaging all corresponding rows
erbb2_idx <- which(symbol == "ERBB2")
if(length(erbb2_idx) == 0) stop("ERBB2 not found in symbol vector!")
erbb2_expr <- colMeans(bulk.log[erbb2_idx, , drop = FALSE])

# 3. Calculate MP4 and MP9 signature scores based on mean log-transformed expression
# Only calculate scores using genes present in the matrix
mp4_score <- colMeans(bulk.log[symbol %in% MP4.gene, , drop = FALSE])
mp9_score <- colMeans(bulk.log[symbol %in% MP9.gene, , drop = FALSE])

# 4. Construct the plotting data frame
plot_df <- data.frame(
  Sample = colnames(bulk.log),
  ERBB2 = as.numeric(erbb2_expr),
  MP4_Score = as.numeric(mp4_score),
  MP9_Score = as.numeric(mp9_score)
)
threshold_25 <- stats::quantile(plot_df$ERBB2, 0.25)

# 2. Filter out the lowest 25% of samples and retain the HER2-active group
plot_df <- plot_df[plot_df$ERBB2 > threshold_25, ]

# 5. Plot correlations after log transformation
# MP4 correlation plot
p1 <- ggplot(plot_df, aes(x = MP4_Score, y = ERBB2)) +
  geom_point(color = "#A52A2A", alpha = 0.5, size = 1.5,shape=16) +
  geom_smooth(method = "glm", color = "black", fill = "gray80", alpha = 0.2) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_classic() +
  labs(title = "",
       x = "MP4 Signature Score",
       y = "ERBB2 Expression") +
  theme(aspect.ratio=1)

# MP9 correlation plot
p2 <- ggplot(plot_df, aes(x = MP9_Score, y = ERBB2)) +
  geom_point(color = "#053E7A", alpha = 0.5, size = 1.5,shape=16) +
  geom_smooth(method = "glm", color = "black", fill = "gray80", alpha = 0.2) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_classic() +
  labs(title = "",
       x = "MP9 Signature Score",
       y = " ERBB2 Expression") + 
  theme(aspect.ratio=1)


# Combined display
cowplot::plot_grid(p1, p2, ncol = 2)


gene_sets <- list(
  MP4_Basal_pEMT = MP4.gene,
  MP9_Differentiation = MP9.gene
)
clean_ids <- gsub("\\..*", "", rownames(bulk.TCGA))

# 2. Build the Ensembl-to-symbol mapping table
gene_map <- AnnotationDbi::select(org.Hs.eg.db, 
                                  keys = clean_ids, 
                                  columns = "SYMBOL", 
                                  keytype = "ENSEMBL") %>%
  dplyr::filter(!is.na(SYMBOL)) # Remove unmatched entries

# 3. Convert the matrix and handle duplicated symbols by averaging
# This step converts the Ensembl-based matrix to a symbol-based matrix and merges duplicates
bulk_symbol <- bulk.TCGA %>%
  as.data.frame() %>%
  mutate(ENSEMBL = clean_ids) %>%
  inner_join(gene_map, by = "ENSEMBL") %>%
  dplyr::select(-ENSEMBL) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(across(everything(), mean)) %>% # Average duplicated symbols
  column_to_rownames("SYMBOL") %>%
  as.matrix()
ssgsea_params <- ssgseaParam(bulk_symbol, gene_sets) 

# 2. Calculate scores
ssgsea_scores <- gsva(ssgsea_params)

# 3. Extract and filter the top 75%
bulk_log <- log2(bulk_symbol + 1)

# 2. Extract log-transformed HER2 (ERBB2) expression
# The extracted result is a vector at this stage
erbb2_log <- bulk_log["ERBB2", ]

# 3. Reconstruct the analysis data frame
plot_df_log <- data.frame(
  Sample = colnames(bulk_log),
  ERBB2 = as.numeric(erbb2_log),
  MP4_ssGSEA = as.numeric(ssgsea_scores["MP4_Basal_pEMT", ]),
  MP9_ssGSEA = as.numeric(ssgsea_scores["MP9_Differentiation", ])
)

# 4. Apply the 25% filtering step again, based on log-transformed values
# threshold_25 <- stats::quantile(plot_df_log$ERBB2, 0.25)
# plot_df_final <- plot_df_log[plot_df_log$ERBB2 > threshold_25, ]
plot_df_final = plot_df_log

p1 <- ggplot(plot_df_final, aes(x = MP4_ssGSEA, y = ERBB2)) +
  geom_point(color = "#A52A2A", alpha = 0.5, size = 1.5,shape=16) +
  geom_smooth(method = "glm", color = "black", fill = "gray80", alpha = 0.2) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_classic() +
  labs(title = "",
       x = "MP4 Signature Score",
       y = "ERBB2 Expression") +
  theme(aspect.ratio=1)

# MP9 correlation plot
p2 <- ggplot(plot_df_final, aes(x = MP9_ssGSEA, y = ERBB2)) +
  geom_point(color = "#053E7A", alpha = 0.5, size = 1.5,shape=16) +
  geom_smooth(method = "glm", color = "black", fill = "gray80", alpha = 0.2) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_classic() +
  labs(title = "",
       x = "MP9 Signature Score",
       y = " ERBB2 Expression") + 
  theme(aspect.ratio=1)

cowplot::plot_grid(p1, p2, ncol = 2)

ggsave(file.path(output_dir,"Figures","FigS","Fig4J_TCGA_STAD_bulk.pdf"),width=6,height=3)



#######################################
#Cell-based experiments
#Fig4M EMT qPCR RT4
qpcr.data = read.csv('/home/data/chenwenhui/scRNAseq/GPU_IT/QC/cell_experiment/rt4-qPCR.csv')
View(qpcr.data)
df_norm <- qpcr.data %>%
  dplyr::group_by(Marker) %>%
  dplyr::mutate(Ctrl_Mean = mean(Expression[Group == "Control"])) %>%
  dplyr::mutate(Fold_Change = Expression / Ctrl_Mean) %>%
  dplyr::ungroup()
df_norm$Marker <- factor(df_norm$Marker, levels=c("E-cadherin", "N-cadherin", "Vimentin", "HER2"))
df_norm$Group <- factor(df_norm$Group, levels=c("Control", "EMT"))
dodge_w <- 0.5
p_final <- ggplot2::ggplot(df_norm, aes(x = Marker, y = Fold_Change, fill = Group)) +
  # Layer 1: bar plot showing means
  # Explicitly set position_dodge(dodge_w)
  ggplot2::geom_bar(stat = "summary", fun = mean, position = ggplot2::position_dodge(width=dodge_w), 
                  color = "black", width = 0.5, linewidth = 0.5) +
  
  # Layer 2: error bars (standard error)
  # Use the same dodge_w offset
  ggplot2::geom_errorbar(stat = "summary", fun.data = mean_se, 
                        position = ggplot2::position_dodge(width=dodge_w), 
                        width = 0.25, color = "black", linewidth = 0.5) +
  ggplot2::geom_point(aes(group = Group), 
                      position = ggplot2::position_dodge(width = dodge_w), 
                      size = 2, shape = 16, fill = "black", color = "black", 
                      stroke = 0.8, alpha = 1) +

  # Add between-group P values using t-test annotations
  ggpubr::stat_compare_means(aes(group = Group), method = "t.test", 
                            label = "p.format", hide.ns = FALSE,
                            label.y = max(df_norm$Fold_Change) * 1.2) +
  
  # Color-blind-friendly palette (Control: orange, EMT: blue)
  ggplot2::scale_fill_manual(values = c("Control" = "#A52A2A", "EMT" = "#053E7A")) +
  
  ggplot2::theme_classic() +
  ggplot2::labs(title = "RT4", 
                y = "Relative mRNA Fold Change", x = "") +
  ggplot2::theme(
    axis.text.x = element_text(face = "bold", size = 12, color = "black"),
    axis.title.y = element_text(face = "bold"),
    legend.title = ggplot2::element_blank(),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
p_final
# 5. Save as PDF
ggplot2::ggsave(file.path(output_dir,"Figures","Fig","Fig4M_qPCR_RT4_EMT_Markers.pdf"), plot = p_final, width = 4, height = 3)

#T24
qpcr.data = read.csv('/home/data/chenwenhui/scRNAseq/GPU_IT/QC/cell_experiment/t24-qPCR.csv')
df_norm <- qpcr.data %>%
  dplyr::group_by(Marker) %>%
  dplyr::mutate(Ctrl_Mean = mean(Expression[Group == "Control"])) %>%
  dplyr::mutate(Fold_Change = Expression / Ctrl_Mean) %>%
  dplyr::ungroup()
df_norm$Marker <- factor(df_norm$Marker, levels=c("E-cadherin", "N-cadherin", "Vimentin", "HER2"))
df_norm$Group <- factor(df_norm$Group, levels=c("Control", "EMT"))
dodge_w <- 0.5
p_final <- ggplot2::ggplot(df_norm, aes(x = Marker, y = Fold_Change, fill = Group)) +
  # Layer 1: bar plot showing means
  # Explicitly set position_dodge(dodge_w)
  ggplot2::geom_bar(stat = "summary", fun = mean, position = ggplot2::position_dodge(width=dodge_w), 
                  color = "black", width = 0.5, linewidth = 0.5) +
  
  # Layer 2: error bars (standard error)
  # Use the same dodge_w offset
  ggplot2::geom_errorbar(stat = "summary", fun.data = mean_se, 
                        position = ggplot2::position_dodge(width=dodge_w), 
                        width = 0.25, color = "black", linewidth = 0.5) +
  ggplot2::geom_point(aes(group = Group), 
                      position = ggplot2::position_dodge(width = dodge_w), 
                      size = 2, shape = 16, fill = "black", color = "black", 
                      stroke = 0.8, alpha = 1) +

  # Add between-group P values using t-test annotations
  ggpubr::stat_compare_means(aes(group = Group), method = "t.test", 
                            label = "p.format", hide.ns = FALSE,
                            label.y = max(df_norm$Fold_Change) * 1.2) +
  
  # Color-blind-friendly palette (Control: orange, EMT: blue)
  ggplot2::scale_fill_manual(values = c("Control" = "#A52A2A", "EMT" = "#053E7A")) +
  
  ggplot2::theme_classic() +
  ggplot2::labs(title = "T24", 
                y = "Relative mRNA Fold Change", x = "") +
  ggplot2::theme(
    axis.text.x = element_text(face = "bold", size = 12, color = "black"),
    axis.title.y = element_text(face = "bold"),
    legend.title = ggplot2::element_blank(),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
p_final
# 5. Save as PDF
ggplot2::ggsave(file.path(output_dir,"Figures","Fig","Fig4M_qPCR_T24_EMT_Markers.pdf"), plot = p_final, width = 4, height = 3)


################################
# Fig. 4N CCK-8 cell-viability curves
#RT4
library(drc)
OD.data =  read.csv('/home/data/chenwenhui/scRNAseq/GPU_IT/QC/cell_experiment/RT4_CCK8_data2.csv')
baseline.data = read.csv('/home/data/chenwenhui/scRNAseq/GPU_IT/QC/cell_experiment/RT4_CCK8_data2_baseline.csv')
View(OD.data)
View(baseline.data)
mean.baseline = mean(baseline.data$baseline_OD)
control_0 <- c(OD.data$Control1[1], OD.data$Control2[1], OD.data$Control3[1])
treated_0 <- c(OD.data$Treated1[1], OD.data$Treated2[1], OD.data$Treated3[1])
# Calculate cell viability (%)
df_viability <- OD.data %>%
  mutate(
    # Normalize each replicate separately: (OD - blank) / (OD at 0 concentration - blank) × 100
    Control1_viability = (Control1 - mean.baseline) / (control_0[1] - mean.baseline) * 100,
    Control2_viability = (Control2 - mean.baseline) / (control_0[2] - mean.baseline) * 100,
    Control3_viability = (Control3 - mean.baseline) / (control_0[3] - mean.baseline) * 100,
    Treated1_viability = (Treated1 - mean.baseline) / (treated_0[1] - mean.baseline) * 100,
    Treated2_viability = (Treated2 - mean.baseline) / (treated_0[2] - mean.baseline) * 100,
    Treated3_viability = (Treated3 - mean.baseline) / (treated_0[3] - mean.baseline) * 100
  ) %>%
  dplyr::select(inducer_con, tidyselect::ends_with("_viability"))
#print(df_viability)

df_long <- df_viability %>%
  tidyr::pivot_longer(
    cols = -inducer_con,
    names_to = "Sample",
    values_to = "Viability"
  ) %>%
  dplyr::mutate(
    Group = ifelse(stringr::str_detect(Sample, "Control"), "Control", "Treated"),
    Group = factor(Group, levels = c("Control", "Treated"))
  )

df_summary <- df_long %>%
  dplyr::group_by(inducer_con, Group) %>%
  dplyr::summarise(
    Mean = mean(Viability),
    SD = sd(Viability),
    SE = sd(Viability) / sqrt(dplyr::n()),
    n = dplyr::n(),
    .groups = "drop"
  )

min_conc <- min(df_long$inducer_con[df_long$inducer_con > 0])
df_long_fit <- df_long %>%
  dplyr::mutate(conc_fit = ifelse(inducer_con == 0, min_conc / 100, inducer_con))

# Fit 4PL models by group
# LL.4() is the 4-parameter log-logistic model: Y = Bottom + (Top - Bottom) / (1 + (X/IC50)^HillSlope)
fit_4pl <- drc::drm(Viability ~ conc_fit, 
                    curveid = Group,
                    data = df_long_fit, 
                    fct = drc::LL.4(names = c("HillSlope", "Bottom", "Top", "IC50")))

ic50_values <- drc::ED(fit_4pl, 50, interval = "delta")
fit_combined <- drc::drm(Viability ~ conc_fit, 
                         data = df_long_fit, 
                         fct = drc::LL.4())
comparison <- anova(fit_combined, fit_4pl)
print(comparison)

ic50_comparison <- drc::compParm(fit_4pl, "IC50", "-")
print(ic50_comparison)

conc_seq <- exp(seq(log(min_conc/100), log(max(df_long$inducer_con)), length.out = 100))

pred_control <- predict(fit_4pl, newdata = data.frame(conc_fit = conc_seq, Group = "Control"), 
                        interval = "confidence")
pred_treated <- predict(fit_4pl, newdata = data.frame(conc_fit = conc_seq, Group = "Treated"), 
                        interval = "confidence")
 
df_pred <- tibble::tibble(
  conc = rep(conc_seq, 2),
  Group = factor(rep(c("Control", "Treated"), each = length(conc_seq)), 
                 levels = c("Control", "Treated")),
  Predicted = c(pred_control[, "Prediction"], pred_treated[, "Prediction"])
)
 
# if presentconfidence interval
if (ncol(pred_control) > 1) {
  df_pred$Lower <- c(pred_control[, "Lower"], pred_treated[, "Lower"])
  df_pred$Upper <- c(pred_control[, "Upper"], pred_treated[, "Upper"])
}
 
# ============================================
# 7. Draw the plot with 4PL fitted curves
# ============================================
 
# Handle zero concentration values for plotting
df_summary_plot <- df_summary %>%
  dplyr::mutate(conc_plot = ifelse(inducer_con == 0, min_conc / 100, inducer_con))
 
df_long_plot <- df_long %>%
  dplyr::mutate(conc_plot = ifelse(inducer_con == 0, min_conc / 100, inducer_con))
 
# Extract IC50 values for annotation
ic50_control <- ic50_values["e:Control:50", "Estimate"]
ic50_treated <- ic50_values["e:Treated:50", "Estimate"]
 
p_4pl <- ggplot2::ggplot() +
  # Layer 1: 4PL fitted curves
  ggplot2::geom_line(data = df_pred, ggplot2::aes(x = conc, y = Predicted, color = Group), 
                     linewidth = 1.2) +
  
  # Layer 2: confidence interval (if present)
  # ggplot2::geom_ribbon(data = df_pred, ggplot2::aes(x = conc, ymin = Lower, ymax = Upper, fill = Group), 
  #             alpha = 0.2) +
  
  # Layer 3: mean data points with error bars
  ggplot2::geom_errorbar(data = df_summary_plot, 
                         ggplot2::aes(x = conc_plot, ymin = Mean - SE, ymax = Mean + SE, color = Group),
                         width = 0.1, linewidth = 0.8) +
  ggplot2::geom_point(data = df_summary_plot, 
                      ggplot2::aes(x = conc_plot, y = Mean, color = Group, shape = Group),
                      size = 4, fill = "white", stroke = 1.5) +
  
  # Layer 4: scatter points for individual replicates
  ggplot2::geom_point(data = df_long_plot, 
                      ggplot2::aes(x = conc_plot, y = Viability, color = Group),
                      alpha = 0.3, size = 2) +
  
  # Vertical lines marking IC50
  ggplot2::geom_vline(xintercept = ic50_control, color = "#A52A2A", linetype = "dashed", alpha = 0.7) +
  ggplot2::geom_vline(xintercept = ic50_treated, color = "#053E7A", linetype = "dashed", alpha = 0.7) +
  
  # Log-scale axis
  ggplot2::scale_x_log10(
    breaks = c(min_conc/100, 0.01, 0.1, 1),
    labels = c("0", "0.01", "0.1", "1")
  ) +
  
  # Color palette
  ggplot2::scale_color_manual(values = c("Control" = "#A52A2A", "Treated" = "#053E7A")) +
  ggplot2::scale_fill_manual(values = c("Control" = "#A52A2A", "Treated" = "#053E7A")) +
  ggplot2::scale_shape_manual(values = c("Control" = 16, "Treated" = 22)) +
  
  # Theme
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::labs(
    title = "RT4",
    x = "Drug Concentration (mg/mL)",
    y = "Cell Viability (%)",
    color = "Group",
    fill = "Group",
    shape = "Group"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = ggplot2::element_text(face = "bold"),
    axis.text = ggplot2::element_text(color = "black"),
    legend.position = "top"
  ) +
  
  # Y-axis range
  ggplot2::ylim(0, 120) +
  
  # Add IC50 annotation
  ggplot2::annotate("text", x = 0.001, y = 75, 
                    label = paste0("IC50 = ", round(ic50_control, 3), " mg/mL"),
                    color = "#A52A2A", size = 3.5, hjust = 0) +
  ggplot2::annotate("text", 0.001, y = 65, 
                    label = paste0("IC50 = ", round(ic50_treated, 3), " mg/mL"),
                    color = "#053E7A", size = 3.5, hjust = 0) +
  ggplot2::annotate("text", x = 0.001, y = 55, 
                    label = paste0("p = ", comparison$`p value`[2]),
                    color = "black", size = 3.5, hjust = 0)

 
print(p_4pl)
 
# Save the figure
ggplot2::ggsave(file.path(output_dir,"Figures","Fig","Fig4N_RT4_cell_viability_curve.pdf"), p_4pl, width = 4, height = 3)


#T24
OD.data =  read.csv('/home/data/chenwenhui/scRNAseq/GPU_IT/QC/cell_experiment/T24_CCK8_data2.csv')
baseline.data = read.csv('/home/data/chenwenhui/scRNAseq/GPU_IT/QC/cell_experiment/T24_CCK8_data2_baseline.csv')
View(OD.data)
View(baseline.data)
mean.baseline = mean(baseline.data$baseline_OD)
control_0 <- c(OD.data$Control1[1], OD.data$Control2[1], OD.data$Control3[1])
treated_0 <- c(OD.data$Treated1[1], OD.data$Treated2[1], OD.data$Treated3[1])
# Calculate cell viability (%)
df_viability <- OD.data %>%
  mutate(
    # Normalize each replicate separately: (OD - blank) / (OD at 0 concentration - blank) × 100
    Control1_viability = (Control1 - mean.baseline) / (control_0[1] - mean.baseline) * 100,
    Control2_viability = (Control2 - mean.baseline) / (control_0[2] - mean.baseline) * 100,
    Control3_viability = (Control3 - mean.baseline) / (control_0[3] - mean.baseline) * 100,
    Treated1_viability = (Treated1 - mean.baseline) / (treated_0[1] - mean.baseline) * 100,
    Treated2_viability = (Treated2 - mean.baseline) / (treated_0[2] - mean.baseline) * 100,
    Treated3_viability = (Treated3 - mean.baseline) / (treated_0[3] - mean.baseline) * 100
  ) %>%
  dplyr::select(inducer_con, tidyselect::ends_with("_viability"))
#print(df_viability)

df_long <- df_viability %>%
  tidyr::pivot_longer(
    cols = -inducer_con,
    names_to = "Sample",
    values_to = "Viability"
  ) %>%
  dplyr::mutate(
    Group = ifelse(stringr::str_detect(Sample, "Control"), "Control", "Treated"),
    Group = factor(Group, levels = c("Control", "Treated"))
  )

df_summary <- df_long %>%
  dplyr::group_by(inducer_con, Group) %>%
  dplyr::summarise(
    Mean = mean(Viability),
    SD = sd(Viability),
    SE = sd(Viability) / sqrt(dplyr::n()),
    n = dplyr::n(),
    .groups = "drop"
  )

min_conc <- min(df_long$inducer_con[df_long$inducer_con > 0])
df_long_fit <- df_long %>%
  dplyr::mutate(conc_fit = ifelse(inducer_con == 0, min_conc / 100, inducer_con))

# Fit 4PL models by group
# LL.4() is the 4-parameter log-logistic model: Y = Bottom + (Top - Bottom) / (1 + (X/IC50)^HillSlope)
fit_4pl <- drc::drm(Viability ~ conc_fit, 
                    curveid = Group,
                    data = df_long_fit, 
                    fct = drc::LL.4(names = c("HillSlope", "Bottom", "Top", "IC50")))

ic50_values <- drc::ED(fit_4pl, 50, interval = "delta")
fit_combined <- drc::drm(Viability ~ conc_fit, 
                         data = df_long_fit, 
                         fct = drc::LL.4())
comparison <- anova(fit_combined, fit_4pl)
print(comparison)

ic50_comparison <- drc::compParm(fit_4pl, "IC50", "-")
print(ic50_comparison)

conc_seq <- exp(seq(log(min_conc/100), log(max(df_long$inducer_con)), length.out = 100))

pred_control <- predict(fit_4pl, newdata = data.frame(conc_fit = conc_seq, Group = "Control"), 
                        interval = "confidence")
pred_treated <- predict(fit_4pl, newdata = data.frame(conc_fit = conc_seq, Group = "Treated"), 
                        interval = "confidence")
 
df_pred <- tibble::tibble(
  conc = rep(conc_seq, 2),
  Group = factor(rep(c("Control", "Treated"), each = length(conc_seq)), 
                 levels = c("Control", "Treated")),
  Predicted = c(pred_control[, "Prediction"], pred_treated[, "Prediction"])
)
 
# if presentconfidence interval
if (ncol(pred_control) > 1) {
  df_pred$Lower <- c(pred_control[, "Lower"], pred_treated[, "Lower"])
  df_pred$Upper <- c(pred_control[, "Upper"], pred_treated[, "Upper"])
}
 
# ============================================
# 7. Draw the plot with 4PL fitted curves
# ============================================
 
# Handle zero concentration values for plotting
df_summary_plot <- df_summary %>%
  dplyr::mutate(conc_plot = ifelse(inducer_con == 0, min_conc / 100, inducer_con))
 
df_long_plot <- df_long %>%
  dplyr::mutate(conc_plot = ifelse(inducer_con == 0, min_conc / 100, inducer_con))
 
# Extract IC50 values for annotation
ic50_control <- ic50_values["e:Control:50", "Estimate"]
ic50_treated <- ic50_values["e:Treated:50", "Estimate"]
 
p_4pl <- ggplot2::ggplot() +
  # Layer 1: 4PL fitted curves
  ggplot2::geom_line(data = df_pred, ggplot2::aes(x = conc, y = Predicted, color = Group), 
                     linewidth = 1.2) +
  
  # Layer 2: confidence interval (if present)
  # ggplot2::geom_ribbon(data = df_pred, ggplot2::aes(x = conc, ymin = Lower, ymax = Upper, fill = Group), 
  #             alpha = 0.2) +
  
  # Layer 3: mean data points with error bars
  ggplot2::geom_errorbar(data = df_summary_plot, 
                         ggplot2::aes(x = conc_plot, ymin = Mean - SE, ymax = Mean + SE, color = Group),
                         width = 0.1, linewidth = 0.8) +
  ggplot2::geom_point(data = df_summary_plot, 
                      ggplot2::aes(x = conc_plot, y = Mean, color = Group, shape = Group),
                      size = 4, fill = "white", stroke = 1.5) +
  
  # Layer 4: scatter points for individual replicates
  ggplot2::geom_point(data = df_long_plot, 
                      ggplot2::aes(x = conc_plot, y = Viability, color = Group),
                      alpha = 0.3, size = 2) +
  
  # Vertical lines marking IC50
  ggplot2::geom_vline(xintercept = ic50_control, color = "#A52A2A", linetype = "dashed", alpha = 0.7) +
  ggplot2::geom_vline(xintercept = ic50_treated, color = "#053E7A", linetype = "dashed", alpha = 0.7) +
  
  # Log-scale axis
  ggplot2::scale_x_log10(
    breaks = c(min_conc/100, 0.01, 0.1, 1),
    labels = c("0", "0.01", "0.1", "1")
  ) +
  
  # Color palette
  ggplot2::scale_color_manual(values = c("Control" = "#A52A2A", "Treated" = "#053E7A")) +
  ggplot2::scale_fill_manual(values = c("Control" = "#A52A2A", "Treated" = "#053E7A")) +
  ggplot2::scale_shape_manual(values = c("Control" = 16, "Treated" = 22)) +
  
  # Theme
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::labs(
    title = "T24",
    x = "Drug Concentration (mg/mL)",
    y = "Cell Viability (%)",
    color = "Group",
    fill = "Group",
    shape = "Group"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = ggplot2::element_text(face = "bold"),
    axis.text = ggplot2::element_text(color = "black"),
    legend.position = "top"
  ) +
  
  # Y-axis range
  ggplot2::ylim(0, 120) +
  
  # Add IC50 annotation
  ggplot2::annotate("text", x = 0.001, y = 75, 
                    label = paste0("IC50 = ", round(ic50_control, 3), " mg/mL"),
                    color = "#A52A2A", size = 3.5, hjust = 0) +
  ggplot2::annotate("text", 0.001, y = 65, 
                    label = paste0("IC50 = ", round(ic50_treated, 3), " mg/mL"),
                    color = "#053E7A", size = 3.5, hjust = 0) +
  ggplot2::annotate("text", x = 0.001, y = 55, 
                    label = paste0("p = ", comparison$`p value`[2]),
                    color = "black", size = 3.5, hjust = 0)

 
print(p_4pl)
 
# Save the figure
ggplot2::ggsave(file.path(output_dir,"Figures","Fig","Fig4N_T24_cell_viability_curve.pdf"), p_4pl, width = 4, height = 3)


