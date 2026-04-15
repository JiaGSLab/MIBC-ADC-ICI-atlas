# ==================================================
# Script: Supplementary Figure S2 plotting
# Author: Wenhui Chen
# Purpose: Generates plots used in Supplementary Figure S2.
# Input: Main Seurat or summary objects already generated in the project, plus any file paths referenced in the script.
# Output: Intermediate result files and/or figure files written to output_dir.
# Note: Review file paths before running. Keep original analysis logic unchanged.
# ==================================================
Tcell.obj = qs_read("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/Tcell_obj_intergration_new.qs")
table(Idents(Tcell.obj))

#sample,Time_point,Outcome UMAP

umap_df = data.frame(Embeddings(Tcell.obj, "umap"),
                    Tcell.obj$sample,
                    Tcell.obj$Time_point,
                    Tcell.obj$Outcome)
colnames(umap_df)=c('UMAP_1','UMAP_2','sample',"Time_point","Outcome")
umap_df$sample = factor(umap_df$sample,levels = sort(unique(umap_df$sample)))
umap_df$Time_point = factor(umap_df$Time_point,levels = c("Pre","Post"))
umap_df$Outcome = factor(umap_df$Outcome,levels = c("NR","CR"))

p1<- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = sample)) +
       geom_point(size = 1,shape = 16, stroke = 0) +
       scale_color_manual(values = c(colors_91)) +
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
ggsave(filename ="./Figures/FigS/FigS02A_sample_Umap.png",height = 4,width = 4,plot = p1,dpi = 600,bg="white")


p1<- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = Time_point)) +
       geom_point(size = 1,shape = 16, stroke = 0) +
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
ggsave(filename ="./Figures/FigS/FigS02A_timepoint_Umap.png",height = 4,width = 4,plot = p1,dpi = 600,bg="white")

p1<- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = Outcome)) +
       geom_point(size = 1,shape = 16, stroke = 0) +
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
ggsave(filename ="./Figures/FigS/FigS02A_Outcome_Umap.png",height = 4,width = 4,plot = p1,dpi = 600,bg="white")

##############
#CD4 CD8A TRAC TRDC FOXP3 FeaturePlot
pdf(file = file.path(output_dir,"Figures","FigS","FigS02A_Tcell_marker_FeaturePlot.pdf"),width = 4,height = 4)
FeaturePlot(Tcell.obj,features = c("CD4","FOXP3","CD8A","TRDC"),
            cols=c('lightgrey',"red"),ncol = 2,raster = T,pt.size=3) & NoLegend() & NoAxes()
dev.off()

pdf(file = file.path(output_dir,"Figures","FigS","colorbar.pdf"),width = 4,height = 4)
FeaturePlot(Tcell.obj,features = c("CD4"),
            cols=c('lightgrey',"red"),raster = T,pt.size=3)
dev.off()

####################
#Compare with the CD8 T-cell reference from https://spica.unil.ch/refs
CD8T.ref = readRDS("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/CD8T_human_ref_v1.rds")
CD8T.subcelltype = c("t01_CD8_IL7R","t02_CD8_AKNA","t03_CD8_EOMES","t04_CD8_ZNF683","t05_CD8_HAVCR2","t06_CD8_ISG","t07_CD8_MKI67")
CD8T.obj = subset(Tcell.obj,idents=CD8T.subcelltype)
# 1. Get the shared highly variable genes between the two objects
common_genes <- intersect(
  Seurat::VariableFeatures(CD8T.obj), 
  Seurat::VariableFeatures(CD8T.ref)
)

# 2. Calculate average expression for user-defined subclusters
# Use Idents() as the grouping variable
user_avg <- Seurat::AverageExpression(
  CD8T.obj, 
  features = common_genes, 
  group.by = "ident", 
  slot = "data"
)$RNA

# 3. Calculate average expression for the reference subclusters
ref_avg <- Seurat::AverageExpression(
  CD8T.ref, 
  features = common_genes, 
  group.by = "ident", 
  slot = "data"
)$RNA
# 4. Calculate the Spearman correlation matrix
# Rows are user-defined subclusters and columns are reference subclusters
user_matrix <- as.matrix(user_avg[common_genes, ])
ref_matrix <- as.matrix(ref_avg[common_genes, ])
cor_mat <- stats::cor(user_matrix, ref_matrix, method = "spearman")

z_mat <- t(scale(t(cor_mat)))

# Define a color-blind-friendly palette: dark blue (low) -> white (mid) -> dark red (high)
# Set the color midpoint to the median rather than 0
col_fun <- circlize::colorRamp2(
  c(min(z_mat), 0, max(z_mat)), 
  c("#4575B4", "#FFFFBF", "#D73027")
)

# 7. Draw the heatmap
hp <- ComplexHeatmap::Heatmap(
  z_mat,
  name = "z-score scaled\nsimilarity",
  col = col_fun,
  cluster_rows = TRUE,          # Cluster your subclusters
  cluster_columns = TRUE,       # Cluster the reference subclusters
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  rect_gp = grid::gpar(col = "white", lwd = 0.5), # White borders around cells
  column_names_rot = 45,        # Rotate bottom labels by 45 degrees
  row_names_side = "left",      # Place your subcluster labels on the left
  column_title = "Reference Subtypes (CD8T.ref)",
  row_title = "Your Clusters (CD8T.obj)",
  cell_fun = function(j, i, x, y, width, height, fill) {
    # If you want to display raw P values or correlation coefficients inside the tiles, add that here
  }
)

# Render the heatmap
ComplexHeatmap::draw(hp)
pdf(file=file.path(output_dir, 'Figures',"FigS", 'FigS02C_Tcell_similarity_with_ref.pdf'),width=5,height=3)
ComplexHeatmap::draw(hp)
dev.off()

###########################
#Ro/e
do.tissueDist.full(
    Tcell.obj@meta.data,
    meta.cluster = 'subcelltype',
    Group.info = 'Group',
    bar_colors = colors_15, # Pass custom colors
    pdf_width = 8,         # Increase width to accommodate 18 subclusters
    fig.output = file.path(output_dir, 'Figures',"FigS", 'FigS02B_Tcell_subcelltype_Roe_with_bar.pdf')
)

###############
#Cell proportions
get_prop_df <- function(obj, cluster_col = NULL) {
  if (is.null(cluster_col)) {
    obj$current_clusters <- Seurat::Idents(obj)
    cluster_col <- "current_clusters"
  }
  
  df <- obj@meta.data %>%
    dplyr::group_by(sample, Time_point, Outcome, !!rlang::sym(cluster_col)) %>%
    dplyr::summarise(count = n(), .groups = "drop") %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(percentage = (count / sum(count)) * 100) %>%
    dplyr::ungroup()
  
  # Ensure factor order
  df$Time_point <- factor(df$Time_point, levels = c("Pre", "Post"))
  df$Outcome <- factor(df$Outcome, levels = c("CR", "NR"))
  return(df)
}

# --- Function 2: permutation test with automatic direction ---
auto_perm_p <- function(data, formula, grp_col) {
  grps <- levels(factor(data[[grp_col]]))
  # Calculate medians for the two groups
  val1 <- data$percentage[data[[grp_col]] == grps[1]]
  val2 <- data$percentage[data[[grp_col]] == grps[2]]
  
  # Check sample size; do not test when n < 2
  if(length(val1) < 2 | length(val2) < 2) return(NA)
  
  direction <- ifelse(median(val2) >= median(val1), "less", "greater")
  
  test_res <- coin::independence_test(formula, data = data, 
                                      distribution = "exact", 
                                      alternative = direction)
  return(coin::pvalue(test_res))
}

run_full_stats <- function(prop_df, cluster_col) {
  clusters <- unique(prop_df[[cluster_col]])
  
  results_list <- lapply(clusters, function(cl) {
    # Extract data for the current subcluster
    sub_df <- prop_df[prop_df[[cluster_col]] == cl, ]
    
    # Define the calculation logic and ensure numeric output
    calc_p <- function(data_subset, formula, grp_col) {
      p_raw <- auto_perm_p(data_subset, formula, grp_col)
      return(as.numeric(p_raw)) # Key fix: force conversion to numeric
    }
    
    # Construct a one-row data frame
    data.frame(
      cluster = as.character(cl),
      p_CR_Time = calc_p(dplyr::filter(sub_df, Outcome == "CR"), percentage ~ Time_point, "Time_point"),
      p_NR_Time = calc_p(dplyr::filter(sub_df, Outcome == "NR"), percentage ~ Time_point, "Time_point"),
      p_Pre_Out = calc_p(dplyr::filter(sub_df, Time_point == "Pre"), percentage ~ Outcome, "Outcome"),
      p_Post_Out = calc_p(dplyr::filter(sub_df, Time_point == "Post"), percentage ~ Outcome, "Outcome"),
      stringsAsFactors = FALSE
    )
  })
  
  # Combine results using bind_rows
  final_stats <- dplyr::bind_rows(results_list)
  return(final_stats)
}

plot_clusters_with_p <- function(prop_df, stats_df, cluster_col, title_prefix = "") {
  cb_colors <- c("Pre" = "#E69F00", "Post" = "#56B4E9") # Color-blind-friendly colors
  
  # Calculate the maximum y value of each panel to determine annotation height
  y_max_df <- prop_df %>% 
    dplyr::group_by(!!rlang::sym(cluster_col)) %>% 
    dplyr::summarise(y_limit = max(percentage) * 1.15)
  
  p <- ggplot2::ggplot(prop_df, ggplot2::aes(x = Outcome, y = percentage, fill = Time_point)) +
    ggplot2::stat_summary(fun = "mean", geom = "bar", position = position_dodge(0.8), 
                          width = 0.7, color = "black", alpha = 0.7, size = 0.3) +
    ggplot2::stat_summary(fun.data = "mean_se", geom = "errorbar", 
                          position = position_dodge(0.8), width = 0.2, size = 0.4) +
    ggplot2::geom_point(ggplot2::aes(color = Time_point), position = position_dodge(0.8), 
                        size = 1.2, show.legend = FALSE) +
    ggplot2::facet_wrap(as.formula(paste0("~", cluster_col)), scales = "free_y") +
    ggplot2::scale_fill_manual(values = cb_colors) +
    ggplot2::scale_color_manual(values = cb_colors) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = title_prefix, y = "Cell Proportion (%)") +
    ggplot2::theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))

  # Note: batch annotation of P values in faceted plots is more complicated; it is often better to draw significant types separately or use this loop as a template
  # This shows an example workflow for annotating Pre vs Post within the CR group
  p <- p + ggsignif::geom_signif(
    data = prop_df,
    ggplot2::aes(xmin = 0.8, xmax = 1.2, 
                annotations = sprintf("p=%.3f", stats_df$p_CR_Time[match(get(cluster_col), stats_df$cluster)])),
    manual = TRUE, tip_length = 0.01, vjust = -0.2, size = 0.3, textsize = 2.5
  )
  
  return(p)
}

# 1. Extract data (assuming you want to use the current Idents)
prop_t <- get_prop_df(Tcell.obj)
# 2. Calculate statistics
stats_t <- run_full_stats(prop_t, "current_clusters")
# 3. Plot
plot_t <- plot_clusters_with_p(prop_t, stats_t, "current_clusters", "All T Cell Subsets")
print(plot_t)
ggsave(filename ="./Figures/FigS/FigS02C_tcell_proportion_histogram.pdf",height = 8,width = 8,plot = plot_t,dpi = 600,bg="white")

# 1. Extract data
prop_cd8 <- get_prop_df(CD8T.obj)
# 2. Calculate statistics
stats_cd8 <- run_full_stats(prop_cd8, "current_clusters")
# 3. Plot
plot_cd8 <- plot_clusters_with_p(prop_cd8, stats_cd8, "current_clusters", "CD8+ T Cell Sub-clusters")
print(plot_cd8)
ggsave(filename ="./Figures/FigS/FigS02C_CD8T_proportion_histogram.pdf",height = 8,width = 8,plot = plot_cd8,dpi = 600,bg="white")
