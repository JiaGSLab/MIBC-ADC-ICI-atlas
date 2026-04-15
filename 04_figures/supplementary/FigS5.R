# ==================================================
# Script: Supplementary Figure S5 plotting
# Author: Wenhui Chen
# Purpose: Generates plots used in Supplementary Figure S5.
# Input: Main Seurat or summary objects already generated in the project, plus any file paths referenced in the script.
# Output: Intermediate result files and/or figure files written to output_dir.
# Note: Review file paths before running. Keep original analysis logic unchanged.
# ==================================================
##################################################
#Fig. S5A classify cells into MP4-high and MP4-low using GMM
library(mclust)
malignant.obj = qs_read(file.path(output_dir, 'malignant_obj_intergration.qs'))
scores <- malignant.obj$Basal_pEMT

# 2. Fit a Gaussian mixture model (GMM) with exactly two components
# Call explicitly
gmm_fit <- mclust::Mclust(scores, G = 2, verbose = FALSE)

# 3. Extract fitted parameters: means, variances, and mixing proportions
high_comp <- which.max(gmm_fit$parameters$mean)
low_comp  <- which.min(gmm_fit$parameters$mean)

mu <- gmm_fit$parameters$mean
sigma <- sqrt(gmm_fit$parameters$variance$sigmasq)
if(length(sigma) == 1) sigma <- rep(sigma, 2)
lambda <- gmm_fit$parameters$pro

# 2. Define the posterior probability function for belonging to the High component
calc_post_prob <- function(x) {
  p_high <- lambda[high_comp] * dnorm(x, mu[high_comp], sigma[high_comp])
  p_low  <- lambda[low_comp] * dnorm(x, mu[low_comp], sigma[low_comp])
  return(p_high / (p_high + p_low))
}

# 3. Use uniroot to find the x value where the probability equals 0.99
# Search for this point between the two component means
target_prob <- 0.99
threshold_99 <- uniroot(function(x) calc_post_prob(x) - target_prob, 
                        interval = c(min(mu), max(scores)))$root

message("The MP4-high threshold at 99% confidence is: ", round(threshold_99, 3))

# 4. Update group assignment
malignant.obj$MP4_Status_99 <- ifelse(scores >= threshold_99, "MP4_High", "MP4_Low")
malignant.obj$MP4_Status_99 <- factor(malignant.obj$MP4_Status_99, levels = c("MP4_Low", "MP4_High"))

# 6. Construct Gaussian curve data for plotting
# 5. Prepare curve data
x_range <- seq(min(scores) - 0.2, max(scores) + 0.2, length.out = 500)
curve_df <- data.frame(
  x = x_range, 
  c_low = lambda[low_comp] * dnorm(x_range, mu[low_comp], sigma[low_comp]),
  c_high = lambda[high_comp] * dnorm(x_range, mu[high_comp], sigma[high_comp])
)

# 6. Plot using a color-blind-friendly palette
p_99 <- ggplot2::ggplot() +
  ggplot2::geom_density(data = data.frame(s = scores), aes(x = s), 
                        fill = "grey95", color = "grey80", alpha = 0.5) +
  # Draw the two Gaussian components
  ggplot2::geom_line(data = curve_df, aes(x = x, y = c_low), 
                     color = "#56B4E9", size = 1, linetype = "dashed") +
  ggplot2::geom_line(data = curve_df, aes(x = x, y = c_high), 
                     color = "#D55E00", size = 1, linetype = "dashed") +
  # Mark the 0.95 threshold line (bold solid line)
  ggplot2::geom_vline(xintercept = threshold_99, color = "#D55E00", size = 1) +
  # Annotation text
  ggplot2::annotate("text", x = threshold_99 + 0.05, y = 0.5, 
                    label = paste0("Threshold (P>0.99): ", round(threshold_99, 3)),
                    color = "#D55E00", hjust = 0, fontface = "bold") +
  ggplot2::theme_classic() +
  ggplot2::labs(
    title = "Refined MP4 Partitioning (Confidence > 99%)",
    x = "MP4 Score", y = "Density"
  )

print(p_99)
ggsave(file.path(output_dir,"Figures","FigS","FigS5A_MP4_GMM_distribution.pdf"),plot=p_99,width=3,height = 2.5)

##################################
#Bar plot for MP4-high proportion
library(coin)
plot_df <- malignant.obj@meta.data %>%
  dplyr::group_by(sample, Outcome, Time_point) %>%
  dplyr::summarise(
    Total_Cells = n(),
    High_Count = sum(MP4_Status_99 == "MP4_High"),
    Percentage = (High_Count / Total_Cells) * 100,
    .groups = "drop"
  )

# 2. Set factor order
plot_df$Outcome <- factor(plot_df$Outcome, levels = c("CR", "NR"))
plot_df$Time_point <- factor(plot_df$Time_point, levels = c("Pre", "Post"))
# Define the exact permutation test function
get_perm_p <- function(sub_df, group_var) {
  # Use the coin package for exact independence testing
  formula_str <- as.formula(paste("Percentage ~", group_var))
  test_res <- coin::independence_test(formula_str, data = sub_df, distribution = "exact")
  return(as.numeric(coin::pvalue(test_res)))
}

# Calculate P values
p_cr_pre_post <- get_perm_p(plot_df %>% filter(Outcome == "CR"), "Time_point")
p_nr_pre_post <- get_perm_p(plot_df %>% filter(Outcome == "NR"), "Time_point")
p_pre_cr_nr   <- get_perm_p(plot_df %>% filter(Time_point == "Pre"), "Outcome")
p_post_cr_nr  <- get_perm_p(plot_df %>% filter(Time_point == "Post"), "Outcome")

# Define colors: Pre in yellow and Post in blue
tp_cols <- c("Pre" = "#D5A04A", "Post" = "#85C1E9")

p <- ggpubr::ggbarplot(
  plot_df, x = "Outcome", y = "Percentage",
  fill = "Time_point", palette = tp_cols,
  add = "mean_se", # Add mean and standard error
  position = position_dodge(0.8),
  add.params = list(width = 0.5)
) +
  # Overlay sample points to show inter-sample variability
  ggplot2::geom_jitter(
    aes(group = Time_point), 
    position = position_dodge(0.8), alpha = 0.7, size = 2, shape = 16
  ) +
  # Manually add permutation-test P-value annotations
  # Annotate Pre vs Post within the CR group
  ggplot2::annotate("text", x = 1, y = max(plot_df$Percentage) * 1.15, 
                    label = paste0("p = ", round(p_cr_pre_post, 3)), size = 4) +
  # Annotate Pre vs Post within the NR group
  ggplot2::annotate("text", x = 2, y = max(plot_df$Percentage) * 1.15, 
                    label = paste0("p = ", round(p_nr_pre_post, 3)), size = 4) +
  # Set axes and theme
  ggplot2::theme_classic() +
  ggplot2::labs(
    y = "% of MP4-High in Tumor cells",
    x = "",
    fill = "Time Point"
  ) +
  ggplot2::theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(size = 13),
    legend.position = "bottom"
  )

# Save the figure
pdf(file.path(output_dir,"Figures","FigS","FigS5B_MP4_high_proportion.pdf"), width = 3, height = 2.4)
print(p)
dev.off()


###################################################
#Pathway heatmap
mp4_cells <- subset(malignant.obj, subset = MP4_Status_99 == "MP4_High")

# 2. Define the feature columns according to the provided indices
kegg_cols <- colnames(mp4_cells@meta.data)[46:130]
progeny_cols <- colnames(mp4_cells@meta.data)[c(131:143, 154)]
mp_cols <- colnames(mp4_cells@meta.data)[144:153]

# 3. Function to calculate group-wise mean values
get_group_means <- function(obj, features) {
  dat <- obj@meta.data %>%
    dplyr::group_by(Group) %>%
    # Use the new anonymous-function syntax to remove warnings
    dplyr::summarise(dplyr::across(dplyr::all_of(features), \(x) mean(x, na.rm = TRUE))) %>%
    tibble::column_to_rownames("Group") %>%
    t() # Transpose so that rows are pathways and columns are groups
  return(dat)
}
mp_means <- get_group_means(mp4_cells, mp_cols)

# Row Z-score normalization
mp_z <- t(scale(t(mp_means)))

pdf(file.path(output_dir,"Figures","FigS","FigS5C_MP_heatmap.pdf"),width=4,height=5)
ComplexHeatmap::Heatmap(mp_z, 
        col = colorRamp2(c(-2, 0, 2), c("#56B4E9", "white", "#D55E00")), # Match the col argument exactly
        name = "Z-score",
        column_title = "MP Program Activity (MP4-High Cells)",
        cluster_columns = FALSE, 
        row_names_gp = gpar(fontsize = 10),
        column_names_side = "bottom",
        column_names_rot = 45,
        column_order = c("Pre_NR", "Pre_CR", "Post_NR", "Post_CR"))
dev.off()

progeny_means <- get_group_means(mp4_cells, progeny_cols)

# Row Z-score normalization
progeny_z <- t(scale(t(progeny_means)))

pdf(file.path(output_dir,"Figures","FigS","FigS5C_progeny_heatmap.pdf"),width=3.6,height=6)
ComplexHeatmap::Heatmap(progeny_z, 
        col = colorRamp2(c(-2, 0, 2), c("#56B4E9", "white", "#D55E00")), # Match the col argument exactly
        name = "Z-score",
        column_title = "Progeny pathway Activity",
        cluster_columns = FALSE, 
        row_names_gp = gpar(fontsize = 10),
        column_names_side = "bottom",
        column_names_rot = 45,
        column_order = c("Pre_NR", "Pre_CR", "Post_NR", "Post_CR"))
dev.off()


kegg_means <- get_group_means(mp4_cells, kegg_cols)

# Row Z-score normalization
kegg_z <- t(scale(t(kegg_means)))

pdf(file.path(output_dir,"Figures","FigS","FigS5C_kegg_heatmap.pdf"),width=7,height=16)
ComplexHeatmap::Heatmap(kegg_z, 
        col = colorRamp2(c(-2, 0, 2), c("#56B4E9", "white", "#D55E00")), 
        name = "Z-score",
        column_title = "KEGG pathway Activity",
        cluster_columns = FALSE, 
        row_names_gp = gpar(fontsize = 10),
        column_names_side = "bottom",
        column_names_rot = 45,
        column_order = c("Pre_NR", "Pre_CR", "Post_NR", "Post_CR"))
dev.off()