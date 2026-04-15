# ==================================================
# Script: Supplementary Figure S1 plotting
# Author: Wenhui Chen
# Purpose: Generates plots used in Supplementary Figure S1.
# Input: Main Seurat or summary objects already generated in the project, plus any file paths referenced in the script.
# Output: Intermediate result files and/or figure files written to output_dir.
# Note: Review file paths before running. Keep original analysis logic unchanged.
# ==================================================
# Load required packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(patchwork)
library(stringr)
library(scales)
library(stats)

output_dir <- "~/scRNAseq/GPU_IT/QC"
source("/home/data/chenwenhui/scRNAseq/GPU_IT/scripts/function_for_plot.R", echo = TRUE)
if (!dir.exists(file.path(output_dir, "Figures", "FigS"))) {
  dir.create(file.path(output_dir, "Figures", "FigS"), recursive = TRUE)
}

##############
#Fig. S1A relationship between pathology information and response
# 1. Data preprocessing
# Assume df is the data frame you already loaded and pre-cleaned
df <- df_raw %>%
  dplyr::distinct(Sample_ID, .keep_all = TRUE) %>%
  dplyr::mutate(
    Outcome = ifelse(Outcome == 'response', 'CR', 'NR'),
    Outcome = factor(Outcome, levels = c("CR", "NR")),
    Age = as.numeric(Age) 
  )

# 2. Define a color-blind-friendly palette and a Cancer Cell-like theme
# Orange-blue palette: D55E00 (CR), 0072B2 (NR)
color_palette <- c("CR" = "#D55E00", "NR" = "#0072B2") 

theme_cancer_cell <- function() {
  ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(color = "black", angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(color = "black"),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
      strip.background = ggplot2::element_blank()
    )
}

# 3. Function A: continuous-variable analysis (Wilcoxon test + boxplot)
plot_continuous <- function(data, var_name, label_name) {
  sub_df <- data %>% dplyr::filter(!is.na(.data[[var_name]]), !is.na(Outcome))
  
  # Perform Wilcoxon rank-sum test
  test_res <- stats::wilcox.test(as.formula(paste(var_name, "~ Outcome")), data = sub_df)
  p_val <- test_res$p.value
  p_label <- if (p_val < 0.001) "P < 0.001" else paste0("P = ", round(p_val, 3))
  
  ggplot2::ggplot(sub_df, ggplot2::aes(x = Outcome, y = .data[[var_name]], fill = Outcome)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.5, size = 0.6) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.4, size = 1, color = "black") +
    ggplot2::scale_fill_manual(values = color_palette) +
    ggplot2::labs(title = label_name, x = NULL, y = label_name) +
    ggplot2::annotate("text", x = 1.5, y = max(sub_df[[var_name]], na.rm = TRUE) * 1.05, 
                      label = p_label, fontface = "italic", size = 3) +
    theme_cancer_cell() +
    ggplot2::theme(legend.position = "none") # Do not show legend inside subplots
}

# 4. Function B: categorical-variable analysis (Fisher/chi-squared test + stacked plot)
plot_categorical <- function(data, factor_col) {
  sub_df <- data %>%
    dplyr::filter(!is.na(.data[[factor_col]]), .data[[factor_col]] != "", !is.na(Outcome))
  
  tab <- table(sub_df[[factor_col]], sub_df$Outcome)
  p_val <- if (any(stats::chisq.test(tab)$expected < 5)) stats::fisher.test(tab)$p.value else stats::chisq.test(tab)$p.value
  p_label <- if (p_val < 0.001) "P < 0.001" else paste0("P = ", round(p_val, 3))
  
  plot_df <- sub_df %>%
    dplyr::count(.data[[factor_col]], Outcome) %>%
    dplyr::group_by(.data[[factor_col]]) %>%
    dplyr::mutate(prop = n / sum(n))
  
  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[factor_col]], y = prop, fill = Outcome)) +
    ggplot2::geom_col(position = "stack", color = "black", width = 0.7) +
    ggplot2::scale_y_continuous(labels = scales::percent, expand = ggplot2::expansion(mult = c(0, 0.1))) +
    ggplot2::scale_fill_manual(values = color_palette) +
    ggplot2::labs(title = factor_col, x = NULL, y = "Percentage") +
    ggplot2::annotate("text", x = Inf, y = Inf, label = p_label, 
                      hjust = 1.1, vjust = 1.5, fontface = "italic", size = 3) +
    theme_cancer_cell()
}

# 5. Generate and combine subplots (5 columns × 1 row)
p1_age <- plot_continuous(df, "Age", "Age (Years)")
p2_sex <- plot_categorical(df, "Sex")
p3_grd <- plot_categorical(df, "Grade")
p4_inv <- plot_categorical(df, "Invasive")
p5_dif <- plot_categorical(df, "Differentiation")

# Combine plots and place the legend at the bottom
combined_fig <- (p1_age | p2_sex | p3_grd | p4_inv | p5_dif) +
  patchwork::plot_layout(guides = "collect") & 
  ggplot2::theme(legend.position = "bottom")

# 6. Save
# Because there are 5 columns, the width is recommended to be at least 3.5 times the height
ggplot2::ggsave(
  file.path(output_dir, "Figures", "FigS", "FigS01A_pathology_outcome.pdf"), 
  combined_fig, 
  width = 7, 
  height = 3.5
)

##############################################
#Stacked bar plot of major cell types
seu.combined.IT = qs_read("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/seu_combined_IT_filtered_obj_intergration.qs")
umap_df = data.frame(seu.combined.IT$major,
                    seu.combined.IT$Group)
colnames(umap_df)=c('major',"Group")
#umap_df$major = ifelse(umap_df$major %in% c("Endothelial","Fibroblast","Mural cells"),"Stromal",umap_df$major)
umap_df$major=factor(umap_df$major,levels = c('T cell','B cell', 'Myeloid','Fibroblast','Endothelial','Mural cells','Epithelial','malignant'))
p = plot_stacked_bar(umap_df,celltype_col = "major",group_col="Group",
                    output_name = file.path(output_dir,"Figures","FigS","FigS01B_major_celltype_stacked_bar.pdf"),colors=colors_15,pdf_height = 3.5)
p

####################################


##############################################
#sample,Time_point,Outcome UMAP
seu.combined.IT$sampleid = paste0(seu.combined.IT$Sample_ID,"-",seu.combined.IT$Time_point)
umap_df = data.frame(Embeddings(seu.combined.IT, "umap"),
                    seu.combined.IT$sampleid,
                    seu.combined.IT$Time_point,
                    seu.combined.IT$Outcome)
colnames(umap_df)=c('UMAP_1','UMAP_2','sample',"Time_point","Outcome")
umap_df$sample = factor(umap_df$sample,levels = sort(unique(umap_df$sample)))
umap_df$Time_point = factor(umap_df$Time_point,levels = c("Pre","Post"))
umap_df$Outcome = factor(umap_df$Outcome,levels = c("NR","CR"))

p1<- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = sample)) +
       geom_point(size = 0.15,shape = 16, stroke = 0) +
       scale_color_manual(values = c(colors_91[1:54],colors_15[5:6])) +
       theme_classic() + # Use a clean theme
       theme( plot.background = element_blank(), # Remove background
              panel.grid.major = element_blank(), # Remove major grid lines
              panel.grid.minor = element_blank(), # Remove minor grid lines
              axis.title.x = element_blank(), # Remove x-axis title
              axis.title.y = element_blank(), # Remove y-axis title
              axis.text = element_blank(), # Remove axis tick labels
              axis.ticks = element_blank(), # Remove axis tick marks
              axis.line = element_blank(),
              legend.position ="none",
              aspect.ratio = 1)
p1
ggsave(filename ="./Figures/FigS/FigS01C_sample_Umap.png",height = 4,width = 4,plot = p1,dpi = 600,bg="white")


p2 <- ggplot(data = umap_df, aes(x = UMAP_1, y = UMAP_2, color = sample)) +
       geom_point(size = 0, shape = 16,stroke = 0) +
       theme_void()+
       scale_color_manual(values = c(colors_91[1:54],colors_15[5:6]), name ='') +
       theme( aspect.ratio = 1,
              legend.position ='bottom',
              plot.margin = margin(0,0,0,0) ) +
       guides(color = guide_legend( ncol = 4, override.aes = list(size = 2, alpha = 1) )) +
       theme(legend.text = element_text(size = 5),
              legend.spacing.y = unit(0,'cm'),
              legend.key.height = unit(0,"cm"),
              legend.box.spacing = unit(0,'cm'))

ggsave("./Figures/FigS/FigS01C_sample_legend.pdf", plot = p2,width = 10, height = 6, dpi = 500 )


p1<- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = Time_point)) +
       geom_point(size = 0.15,shape = 16, stroke = 0) +
       scale_color_manual(values = colors_15) +
       theme_classic() + # Use a clean theme
       theme( plot.background = element_blank(), # Remove background
              panel.grid.major = element_blank(), # Remove major grid lines
              panel.grid.minor = element_blank(), # Remove minor grid lines
              axis.title.x = element_blank(), # Remove x-axis title
              axis.title.y = element_blank(), # Remove y-axis title
              axis.text = element_blank(), # Remove axis tick labels
              axis.ticks = element_blank(), # Remove axis tick marks
              axis.line = element_blank(),
              legend.position ="none",
              aspect.ratio = 1)
p1
ggsave(filename ="./Figures/FigS/FigS01C_timepoint_Umap.png",height = 4,width = 4,plot = p1,dpi = 600,bg="white")

p1<- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = Outcome)) +
       geom_point(size = 0.15,shape = 16, stroke = 0) +
       scale_color_manual(values = colors_15[3:4]) +
       theme_classic() + # Use a clean theme
       theme( plot.background = element_blank(), # Remove background
              panel.grid.major = element_blank(), # Remove major grid lines
              panel.grid.minor = element_blank(), # Remove minor grid lines
              axis.title.x = element_blank(), # Remove x-axis title
              axis.title.y = element_blank(), # Remove y-axis title
              axis.text = element_blank(), # Remove axis tick labels
              axis.ticks = element_blank(), # Remove axis tick marks
              axis.line = element_blank(),
              legend.position ="none",
              aspect.ratio = 1)
p1
ggsave(filename ="./Figures/FigS/FigS01C_Outcome_Umap.png",height = 4,width = 4,plot = p1,dpi = 600,bg="white")

##########################
#marker featureplot
pdf(file = file.path(output_dir,"Figures","FigS","FigS01F_marker_FeaturePlot.pdf"),width = 6,height = 6)
FeaturePlot(seu.combined.IT,features = c("KRT19","CD3D","CD79A","CD14","CSF3R","COL1A1","VWF","ACTA2",'RGS5'),
            cols=c('lightgrey',"red"),ncol = 3,raster = T) & NoLegend() & NoAxes()
dev.off()

###############################################
#major Ro/e
do.tissueDist.full(seu.combined.IT@meta.data,
                  meta.cluster = "major",
                  fig.output = file.path(output_dir,"Figures","FigS","FigS01D_major_lineage_Distribution.pdf"),
                  bar_colors = colors_15,pdf_width = 7
                  )


#######################################
#subcelltype Ro/e
do.tissueDist(seu.combined.IT@meta.data,
              fig.output = file.path(output_dir,"Figures","FigS","FigS01E_subcelltype_Distribution.pdf"),pdf.height = 10,pdf.width = 3.5)
