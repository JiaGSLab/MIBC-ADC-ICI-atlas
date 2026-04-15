# ==================================================
# Script: Cell type proportion analysis
# Author: Wenhui Chen
# Purpose: Calculates cell type proportions by sample and group and generates comparison plots.
# Input: Main Seurat or summary objects already generated in the project, plus any file paths referenced in the script.
# Output: Intermediate result files and/or figure files written to output_dir.
# Note: Review file paths before running. Keep original analysis logic unchanged.
# ==================================================
output_dir <- "~/scRNAseq/GPU_IT/QC"
# Read seu.decontX.IT and calculate the count and proportion of each cell type in each sample
#seu.decontX.IT <- qs_read(file.path(output_dir, 'seu_combined_decontX_IT_intergration_celltype.qs'))
seu.decontX.IT$Group = paste0(seu.decontX.IT$Time_point, '_', seu.decontX.IT$Outcome)
# Calculate the cell count and proportion of each cell type in each sample
celltype_prop_data <- seu.decontX.IT@meta.data |>
  as.data.frame() |>
  tibble::as_tibble(rownames = "cell") |>
  dplyr::select(sample, Group, celltype) |>
  dplyr::group_by(sample, Group, celltype) |>
  dplyr::summarise(n_cells = n(), .groups = 'drop') |>
  dplyr::group_by(sample) |>
  dplyr::mutate(total_cells = sum(n_cells),
                proportion = n_cells / total_cells) |>
  dplyr::ungroup()

View(celltype_prop_data)

# Compare whether there are significant differences between Pre_CR and Pre_NR, and between Post_NR and Pre_NR within Group
# Extract the groups to be compared
compare_groups <- celltype_prop_data |>
  dplyr::filter(Group %in% c('Pre_CR', 'Pre_NR', 'Post_NR')) |>
  dplyr::mutate(comparison_group = case_when(
    Group == 'Pre_CR' ~ 'Pre_CR_vs_Pre_NR',
    Group == 'Pre_NR' ~ 'Pre_CR_vs_Pre_NR',
    Group == 'Post_NR' ~ 'Post_NR_vs_Pre_NR'
  ))

# Perform Wilcoxon rank-sum test
wilcox_results <- data.frame()

# Perform the test separately for each cell type and comparison group
for(ct in unique(compare_groups$celltype)) {
  for(comp in c('Pre_CR_vs_Pre_NR', 'Post_NR_vs_Pre_NR')) {
    
    # Select the corresponding Group based on the comparison
    if(comp == 'Pre_CR_vs_Pre_NR') {
      groups_to_compare <- c('Pre_CR', 'Pre_NR')
    } else {
      groups_to_compare <- c('Post_NR', 'Pre_NR')
    }
    
    # Filter the data
    test_data <- compare_groups |>
      dplyr::filter(celltype == ct, Group %in% groups_to_compare)
    
    # Check whether there are enough data points for testing
    if(nrow(test_data) >= 2 && length(unique(test_data$Group)) == 2) {
      
      # Perform Wilcoxon test
      test_result <- wilcox.test(proportion ~ Group, data = test_data)
      
      # Add results to the result data frame
      wilcox_results <- rbind(wilcox_results, data.frame(
        celltype = ct,
        comparison_group = comp,
        p_value = test_result$p.value
      ))
    }
  }
}

# Adjust p values and add significance labels
if(nrow(wilcox_results) > 0) {
  wilcox_results <- wilcox_results |>
    dplyr::mutate(p_value_adj = p.adjust(p_value, method = 'BH'),
                  significance = case_when(
                    p_value_adj < 0.001 ~ '***',
                    p_value_adj < 0.01 ~ '**',
                    p_value_adj < 0.05 ~ '*',
                    TRUE ~ 'ns'
                  ))
}

print(wilcox_results)

# Draw violin plots, one facet per cell type
library(ggplot2)
library(ggpubr)

# Prepare data for plotting
plot_data <- compare_groups |>
  dplyr::filter(comparison_group %in% c('Pre_CR_vs_Pre_NR', 'Post_NR_vs_Pre_NR'))

# Create violin plot
p <- ggplot(plot_data, aes(x = Group, y = proportion, fill = Group)) +
  geom_violin(scale = 'width', trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = 'white', alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 0.6) +
  facet_wrap(~ celltype, scales = 'free_y', ncol = 4) +
  scale_fill_manual(values = c('Pre_CR' = '#2E8B57', 'Pre_NR' = '#CD5C5C', 'Post_NR' = '#FFA500')) +
  labs(title = 'Cell Type Proportions Across Groups',
       x = 'Group',
       y = 'Proportion') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'bottom')

# Add significance labels
p <- p + stat_compare_means(
  method = 'wilcox.test',
  comparisons = list(c('Pre_CR', 'Pre_NR'), c('Post_NR', 'Pre_NR')),
  label = 'p.signif',
  hide.ns = TRUE
)

# Save the figure
pdf(file = file.path(output_dir, 'celltype_proportions_comparison.pdf'), 
    width = 16, height = 12)
print(p)

# Show an individual plot for each cell type
for(ct in unique(plot_data$celltype)) {
  p_single <- ggplot(dplyr::filter(plot_data, celltype == ct), 
                     aes(x = Group, y = proportion, fill = Group)) +
    geom_violin(scale = 'width', trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = 'white', alpha = 0.5, outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.2), size = 1, alpha = 0.8) +
    scale_fill_manual(values = c('Pre_CR' = '#2E8B57', 'Pre_NR' = '#CD5C5C', 'Post_NR' = '#FFA500')) +
    labs(title = paste('Cell Type:', ct),
         x = 'Group',
         y = 'Proportion') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Add significance labels
  p_single <- p_single + stat_compare_means(
    method = 'wilcox.test',
    comparisons = list(c('Pre_CR', 'Pre_NR'), c('Post_NR', 'Pre_NR')),
    label = 'p.signif',
    hide.ns = TRUE
  )
  
  print(p_single)
}

dev.off()

# Create stacked bar plots of cell type composition, with sample on the x-axis, cell type proportion on the y-axis, and one facet per Group
library(ggplot2)
library(dplyr)

# Prepare data for the stacked bar plot, ordered by Epithelial proportion
stacked_bar_data <- celltype_prop_data |>
  dplyr::filter(Group %in% c('Pre_CR', 'Pre_NR', 'Post_CR', 'Post_NR')) |>
  dplyr::group_by(sample) |>
  dplyr::mutate(epithelial_prop = proportion[celltype == 'Epithelial']) |>
  dplyr::ungroup() |>
  dplyr::arrange(Group, desc(epithelial_prop), sample, desc(proportion))

# Set the sample order as a factor
sample_order <- stacked_bar_data |>
  dplyr::distinct(sample, Group, epithelial_prop) |>
  dplyr::arrange(Group, desc(epithelial_prop)) |>
  dplyr::pull(sample)

stacked_bar_data$sample <- factor(stacked_bar_data$sample, levels = sample_order)

# Create stacked bar plot
stacked_plot <- ggplot(stacked_bar_data, aes(x = sample, y = proportion, fill = celltype)) +
  geom_col(position = 'stack') +
  facet_wrap(~ Group, scales = 'free_x', ncol = 2) +
  scale_fill_manual(values = colors_91[1:length(unique(stacked_bar_data$celltype))]) +
  labs(title = 'Cell Type Composition by Sample and Group (Sorted by Epithelial Proportion)',
       x = 'Sample',
       y = 'Proportion',
       fill = 'Cell Type') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'right',
        legend.key.size = unit(0.5, 'cm'))

# Save stacked bar plot
pdf(file = file.path(output_dir, 'celltype_stacked_barplot.pdf'), 
    width = 16, height = 10)
print(stacked_plot)
dev.off()

# Draw a stacked bar plot after removing Epithelial cells, ordered by T-cell proportion
non_epithelial_data <- celltype_prop_data |>
  dplyr::filter(Group %in% c('Pre_CR', 'Pre_NR', 'Post_CR', 'Post_NR'),
                celltype != 'Epithelial') |>
  dplyr::group_by(sample) |>
  dplyr::mutate(total_non_epithelial = sum(proportion),
                proportion_non_epithelial = proportion / total_non_epithelial) |>
  dplyr::ungroup() |>
  dplyr::group_by(sample) |>
  dplyr::mutate(t_cell_prop = proportion_non_epithelial[celltype == 'T']) |>
  dplyr::ungroup() |>
  dplyr::arrange(Group, desc(t_cell_prop), sample, desc(proportion_non_epithelial))

# Set the sample order factor after removing Epithelial cells
sample_order_non_epi <- non_epithelial_data |>
  dplyr::distinct(sample, Group, t_cell_prop) |>
  dplyr::arrange(Group, desc(t_cell_prop)) |>
  dplyr::pull(sample)

non_epithelial_data$sample <- factor(non_epithelial_data$sample, levels = sample_order_non_epi)

# Create the stacked bar plot after removing Epithelial cells
non_epithelial_plot <- ggplot(non_epithelial_data, aes(x = sample, y = proportion_non_epithelial, fill = celltype)) +
  geom_col(position = 'stack') +
  facet_wrap(~ Group, scales = 'free_x', ncol = 2) +
  scale_fill_manual(values = colors_91[1:length(unique(non_epithelial_data$celltype))]) +
  labs(title = 'Non-Epithelial Cell Type Composition (Sorted by T Cell Proportion)',
       x = 'Sample',
       y = 'Proportion (Non-Epithelial Cells)',
       fill = 'Cell Type') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'right',
        legend.key.size = unit(0.5, 'cm'))

# Save the stacked bar plot after removing Epithelial cells
pdf(file = file.path(output_dir, 'non_epithelial_stacked_barplot.pdf'), 
    width = 16, height = 10)
print(non_epithelial_plot)
dev.off()

# Save statistical results
write.csv(celltype_prop_data, file.path(output_dir, 'celltype_proportions_data.csv'), row.names = FALSE)
write.csv(wilcox_results, file.path(output_dir, 'wilcox_test_results.csv'), row.names = FALSE)

cat('Analysis completed. Results saved to:\\n')
cat('- celltype_proportions_comparison.pdf: cell type proportion comparison plot\\n')
cat('- celltype_stacked_barplot.pdf: stacked bar plot of cell types\\n')
cat('- celltype_proportions_data.csv: cell type proportion data\\n')
cat('- wilcox_test_results.csv: Wilcoxon test results\\n')