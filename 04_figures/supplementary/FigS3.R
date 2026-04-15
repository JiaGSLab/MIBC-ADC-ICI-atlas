# ==================================================
# Script: Supplementary Figure S3 plotting
# Author: Wenhui Chen
# Purpose: Generates plots used in Supplementary Figure S3.
# Input: Main Seurat or summary objects already generated in the project, plus any file paths referenced in the script.
# Output: Intermediate result files and/or figure files written to output_dir.
# Note: Review file paths before running. Keep original analysis logic unchanged.
# ==================================================
#############################
#Proportions were calculated in Fig3.R; here plot G2M score versus Tex proportion
meta_info <- seu.combined.IT@meta.data[, c("sample", "Group")]

# 2. Remove duplicates to obtain the unique group label for each sample
# This step is important to ensure that each sample maps to only one Group
sample_group_map <- base::unique(meta_info)
# 1. Extract the target data and convert it to a data frame
df_plot <- data.frame(
  Sample = rownames(proportion_result$proportion_matrix),
  Proportion = proportion_result$proportion_matrix[, "t05_CD8_HAVCR2"],
  MP_Score = mp_stats_result$mean_matrix[rownames(proportion_result$proportion_matrix), "Cycle_G2M"]
)

# 2. Merge group information (assuming sample_group_map uses Sample as row names or contains a Sample column)
# If sample_group_map uses sample IDs as row names:
df_plot$Group <- sample_group_map[df_plot$Sample, "Group"]

# 3. Remove samples containing NA values, if present
df_plot <- na.omit(df_plot)

# Calculate the slope of the global linear regression
fit <- stats::lm(MP_Score ~ Proportion, data = df_plot)
slope_val <- stats::coef(fit)[2]
p_val <- summary(fit)$coefficients[2, 4]

# Construct annotation text
slope_label <- paste0("Slope = ", round(slope_val, 3))

# Define color-blind-friendly colors (Okabe-Ito Palette)

cb_palette <- c("Pre_CR" = "#053E7A", 
                "Pre_NR" = "#A52A2A", 
                "Post_NR" = "#E46B4F",
                "Post_CR" = "#8B8B00")

p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Proportion, y = MP_Score)) +
  # 1. Draw scatter points, with color indicating Group
  ggrastr::geom_point_rast(
      ggplot2::aes(color = Group), 
      size = 4, 
      alpha = 0.8,
      raster.dpi = 300 # 300 DPI is sufficient for AI-based figure editing
    ) +
  # 2. Draw the global fitted line (linear regression)
  ggplot2::geom_smooth(method = "lm", color = "black", linetype = "dashed", fill = "grey80", alpha = 0.2) +
  # 3. Set color mapping
  ggplot2::scale_color_manual(values = cb_palette) +
  # 4. Add slope labels and correlation statistics (R and P added automatically with ggpubr)
  ggpubr::stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  ggplot2::annotate("text", x = min(df_plot$Proportion), y = max(df_plot$MP_Score) * 0.95, 
                    label = slope_label, hjust = 0, fontface = "italic", size = 4) +
  # 5. Theme and labels
  ggplot2::theme_classic() +
  ggplot2::labs(
    title = "Correlation: CD8_HAVCR2 vs Cycle_G2M",
    x = "Proportion of t05_CD8_HAVCR2 (%)",
    y = "Cycle_G2M MP Median Score",
    color = "Patient Group"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    axis.text = ggplot2::element_text(color = "black"),
    legend.position = "right"
  )

print(p)
ggsave(filename = file.path(output_dir, 'Figures',"FigS",'FigS04B_Tex_G2M_proportion_scatter.pdf'),plot=p,width = 4,height =3)


############
df_plot <- data.frame(
  Sample = rownames(proportion_result$proportion_matrix),
  Proportion = proportion_result$proportion_matrix[, "t05_CD8_HAVCR2"],
  MP_Score = mp_stats_result$mean_matrix[rownames(proportion_result$proportion_matrix), "IFN_Response"]
)

# 2. Merge group information (assuming sample_group_map uses Sample as row names or contains a Sample column)
# If sample_group_map uses sample IDs as row names:
df_plot$Group <- sample_group_map[df_plot$Sample, "Group"]

# 3. Remove samples containing NA values, if present
df_plot <- na.omit(df_plot)

# Calculate the slope of the global linear regression
fit <- stats::lm(MP_Score ~ Proportion, data = df_plot)
slope_val <- stats::coef(fit)[2]
p_val <- summary(fit)$coefficients[2, 4]

# Construct annotation text
slope_label <- paste0("Slope = ", round(slope_val, 3))

# Define color-blind-friendly colors (Okabe-Ito Palette)

cb_palette <- c("Pre_CR" = "#053E7A", 
                "Pre_NR" = "#A52A2A", 
                "Post_NR" = "#E46B4F",
                "Post_CR" = "#8B8B00")

p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Proportion, y = MP_Score)) +
  # 1. Draw scatter points, with color indicating Group
  ggrastr::geom_point_rast(
      ggplot2::aes(color = Group), 
      size = 4, 
      alpha = 0.8,
      raster.dpi = 300 # 300 DPI is sufficient for AI-based figure editing
    ) +
  # 2. Draw the global fitted line (linear regression)
  ggplot2::geom_smooth(method = "lm", color = "black", linetype = "dashed", fill = "grey80", alpha = 0.2) +
  # 3. Set color mapping
  ggplot2::scale_color_manual(values = cb_palette) +
  # 4. Add slope labels and correlation statistics (R and P added automatically with ggpubr)
  ggpubr::stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  ggplot2::annotate("text", x = min(df_plot$Proportion), y = max(df_plot$MP_Score) * 0.95, 
                    label = slope_label, hjust = 0, fontface = "italic", size = 4) +
  # 5. Theme and labels
  ggplot2::theme_classic() +
  ggplot2::labs(
    title = "Correlation: CD8_HAVCR2 vs IFN_response",
    x = "Proportion of t05_CD8_HAVCR2 (%)",
    y = "IFN_response MP Median Score",
    color = "Patient Group"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    axis.text = ggplot2::element_text(color = "black"),
    legend.position = "right"
  )

print(p)
ggsave(filename = file.path(output_dir, 'Figures',"FigS",'FigS04B_Tex_IFN_proportion_scatter.pdf'),plot=p,width = 4,height =3)


##################################
#Correlation of ISG score with PDCD1 expression and CXCL10
malignant.obj = qs_read(file.path(output_dir, 'malignant_obj_intergration.qs'))
raw_data <- Seurat::FetchData(
  object = malignant.obj,
  vars = c("sample", "Group", "IFN_Response", "CD274", "CXCL10"),
  layer = "data"
)

# 2. Aggregate by sample: calculate median IFN score and average gene expression
sample_df <- raw_data %>%
  dplyr::group_by(sample, Group) %>%
  dplyr::summarise(
    Median_IFN = stats::median(IFN_Response, na.rm = TRUE),
    Avg_CD274 = mean(CD274, na.rm = TRUE),
    Avg_CXCL10 = mean(CXCL10, na.rm = TRUE),
    .groups = "drop"
  )

# Ensure Group order (for example, CR first and NR second)
sample_df$Group <- factor(sample_df$Group)

plot_sample_level <- function(df, y_var, y_label) {
  # 1. Calculate the slope of the global linear regression
  fit_formula <- as.formula(paste(y_var, "~ Median_IFN"))
  fit <- stats::lm(fit_formula, data = df)
  slope_val <- stats::coef(fit)[2]
  
  # 2. Define color-blind-friendly colors (CR: vermilion, NR: sky blue)
  # You can also fine-tune the exact colors to match previous plots
  cb_colors <- c("Pre_CR" = "#053E7A", 
                "Pre_NR" = "#A52A2A", 
                "Post_NR" = "#E46B4F",
                "Post_CR" = "#8B8B00") 
  
  ggplot2::ggplot(df, ggplot2::aes(x = Median_IFN, y = .data[[y_var]])) +
    # Rasterized scatter-point layer: colors indicate Group
    ggrastr::geom_point_rast(
      ggplot2::aes(color = Group), 
      size = 4, 
      alpha = 0.8,
      raster.dpi = 300 # 300 DPI is sufficient for AI-based figure editing
    ) +
    # Global fitted curve (black dashed line; fitted across all samples to show the overall trend)
    ggplot2::geom_smooth(
      ggplot2::aes(group = 1), 
      method = "lm", 
      color = "black", 
      linetype = "dashed", 
      fill = "grey80", 
      alpha = 0.2
    ) +
    # Add correlation coefficient and P value
    ggpubr::stat_cor(method = "pearson", label.x.npc = "left") +
    # Annotate the slope
    ggplot2::annotate(
      "text", 
      x = min(df$Median_IFN), 
      y = max(df[[y_var]]) * 0.95, 
      label = paste("Slope =", round(slope_val, 3)), 
      hjust = 0, fontface = "italic", size = 4
    ) +
    # Adjust styling
    ggplot2::scale_color_manual(values = cb_colors) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste("Sample-level Correlation:", y_label),
      x = "Median IFN Response Score",
      y = paste("Average", y_label, "Expression"),
      color = "Group"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text = ggplot2::element_text(color = "black"),
      legend.position = "right"
    )
}

# Generate the two plots
p1 <- plot_sample_level(sample_df, "Avg_CD274", "CD274")
p2 <- plot_sample_level(sample_df, "Avg_CXCL10", "CXCL10")

# Combine the plots
combined_p <- p1 / p2
ggsave(filename = file.path(output_dir, 'Figures',"FigS",'FigS04C_IFN_CD274_CXCL10_proportion_scatter.pdf'),plot=combined_p,width = 4,height =6)
