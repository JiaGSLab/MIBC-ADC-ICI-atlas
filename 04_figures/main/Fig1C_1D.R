# ==================================================
# Script: Figure 1C-1D plotting
# Author: Wenhui Chen
# Purpose: Generates UMAP and paired proportion plots used in Figure 1.
# Input: Main Seurat or summary objects already generated in the project, plus any file paths referenced in the script.
# Output: Intermediate result files and/or figure files written to output_dir.
# Note: Review file paths before running. Keep original analysis logic unchanged.
# ==================================================
#Fig. 1C UMAP refinement
seu.combined.IT = qs_read("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/seu_combined_IT_filtered_obj_intergration.qs")
umap_df = data.frame(Embeddings(seu.combined.IT, "umap"),
                     seu.combined.IT$major,
                     seu.combined.IT$subcelltype,
                     seu.combined.IT$Group)
colnames(umap_df)=c('UMAP_1','UMAP_2','major','subcelltype',"Group")
umap_df$major = ifelse(umap_df$major %in% c("Endothelial","Fibroblast","Mural cells"),"Stromal",umap_df$major)
umap_df$major=factor(umap_df$major,levels = c('T cell','B cell', 'Myeloid','Stromal','Epithelial','malignant'))
subcl=sort(unique(umap_df$subcelltype))
subcl=subcl[c(41:56,1:9,11:24,26:40,10,25)]
umap_df$subcelltype=factor(umap_df$subcelltype,levels = subcl)
p <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = major)) +
       geom_point(size = 0.05,shape = 16, stroke = 0) +
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
p
ggsave(filename ="./Figures/Fig1C_major_Umap.png",height = 2,width = 2,plot = p,dpi = 600,bg="white")

p1<- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = subcelltype)) +
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
ggsave(filename ="./Figures/Fig1C_subcelltype_Umap.png",height = 4,width = 4,plot = p1,dpi = 600,bg="white")


p2 <- ggplot(data = umap_df, aes(x = UMAP_1, y = UMAP_2, color = subcelltype)) +
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

ggsave("./Figures/Fig1C_subcelltype_legend.pdf", plot = p2,width = 10, height = 6, dpi = 500 )


#Fig1D
meta_df <- seu.combined.IT@meta.data

prop_data <- meta_df %>%
       dplyr::group_by(sample, Sample_ID, Outcome, celltype,Time_point) %>% # modify this if your column names differ
       dplyr::summarise(count = dplyr::n(), .groups = 'drop') %>%
       dplyr::group_by(sample) %>%
       dplyr::mutate(proportion = count / sum(count)) %>%
       dplyr::ungroup()

# 2. Ensure the time-point order (Pre -> Post)
# If you do not have a Time column, you can infer it from the sample string, for example:
# prop_data <- prop_data %>% 
#   dplyr::mutate(Time = ifelse(grepl("Pre", sample), "Pre", "Post"))
prop_data$Time <- factor(prop_data$Time_point, levels = c("Pre", "Post"))

# 3. Define color-blind-friendly colors (CR: orange, NR: blue)
color_palette <- c("CR" = "#D55E00", "NR" = "#0072B2")

# 4. Plot paired line charts
# Use facet_wrap to display all cell types
p_paired_all <- ggplot2::ggplot(prop_data, ggplot2::aes(x = Time, y = proportion, group = Sample_ID)) +
       # Draw paired lines, with color indicating response
       ggplot2::geom_line(ggplot2::aes(color = Outcome), linewidth = 1, alpha = 0.7) +
       # Plot sample points
       ggplot2::geom_point(ggplot2::aes(color = Outcome), size = 2.5) +
       # Facet the plot, one panel per cell type
       ggplot2::facet_wrap(~celltype, scales = "free_y", ncol = 4) +
       # Format the y-axis as percentages
       ggplot2::scale_y_continuous(labels = scales::percent_format()) +
       # Use the specified color palette
       ggplot2::scale_color_manual(values = color_palette) +
       # Style the plot
       ggplot2::labs(
       title = "Cell Type Proportion Dynamics (Paired Samples)",
       subtitle = paste0("n = 4 pairs (3 NR, 1 CR)"),
       x = NULL,
       y = "Proportion of Total Cells"
       ) +
       ggplot2::theme_bw(base_size = 12) +
       ggplot2::theme(
       panel.grid = ggplot2::element_blank(),
       strip.background = ggplot2::element_rect(fill = "white", color = "black"),
       strip.text = ggplot2::element_text(face = "bold"),
       axis.text.x = ggplot2::element_text(face = "bold", color = "black"),
       legend.position = "bottom",
       plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
       )

# 5. Save the figure
ggplot2::ggsave(
       filename = file.path(output_dir, "Figures", "Fig", "Fig1D_CellType_Paired_Proportion.pdf"),
       plot = p_paired_all,
       width = 10, 
       height = 8
)

##########################################
#Bar plot showing major-cell-type proportions
library(coin)
plot_df <- seu.combined.IT@meta.data %>%
       dplyr::group_by(sampleid, Time_point, Outcome, major) %>%
       dplyr::summarise(count = n(), .groups = "drop") %>%
       dplyr::group_by(sampleid) %>%
       dplyr::mutate(percentage = (count / sum(count)) * 100) %>%
       dplyr::ungroup()
# Harmonize factor levels
plot_df$Time_point <- factor(plot_df$Time_point, levels = c("Pre", "Post"))
plot_df$Outcome <- factor(plot_df$Outcome, levels = c("CR", "NR"))
# Function for one-sided permutation testing with automatic direction
auto_perm_p <- function(data, formula, grp_col) {
       # Extract the two groups and calculate medians
       grps <- levels(factor(data[[grp_col]]))
       med1 <- median(data$percentage[data[[grp_col]] == grps[1]], na.rm = TRUE)
       med2 <- median(data$percentage[data[[grp_col]] == grps[2]], na.rm = TRUE)
       
       # Determine the test direction based on median values
       # If med2 > med1, test for a significant increase; otherwise test for a significant decrease
       #direction <- ifelse(med2 >= med1, "greater", "less")
       direction <- ifelse(med2 >= med1, "less", "greater")
       res <- coin::independence_test(formula, data = data, 
                                   distribution = "exact", 
                                   alternative = direction)
       return(coin::pvalue(res))
       }

all_majors <- unique(plot_df$major)
stats_list <- list()

for (m in all_majors) {
       sub_m <- dplyr::filter(plot_df, major == m)
       
       # 1. CR group: Pre vs Post
       p1 <- auto_perm_p(dplyr::filter(sub_m, Outcome == "CR"), percentage ~ Time_point, "Time_point")
       # 2. NR group: Pre vs Post
       p2 <- auto_perm_p(dplyr::filter(sub_m, Outcome == "NR"), percentage ~ Time_point, "Time_point")
       # 3. Pre-treatment stage: CR vs NR
       p3 <- auto_perm_p(dplyr::filter(sub_m, Time_point == "Pre"), percentage ~ Outcome, "Outcome")
       # 4. Post-treatment stage: CR vs NR
       p4 <- auto_perm_p(dplyr::filter(sub_m, Time_point == "Post"), percentage ~ Outcome, "Outcome")
       
       stats_list[[m]] <- data.frame(major = m, p_CR_Time = p1, p_NR_Time = p2, p_Pre_Out = p3, p_Post_Out = p4)
}

stats_df <- dplyr::bind_rows(stats_list)

# Color-blind-friendly palette
# Pre: #D55E00 (Vermillion), Post: #0072B2 (Blue)
cb_colors <- c("Pre" = "#D55E00", "Post" = "#0072B2")

# Base plot
p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Outcome, y = percentage, fill = Time_point)) +
  # Bar plot
  ggplot2::stat_summary(fun = "mean", geom = "bar", position = position_dodge(0.8), 
                        width = 0.7, color = "black", alpha = 0.7) +
  # Error bars (SEM)
  ggplot2::stat_summary(fun.data = "mean_se", geom = "errorbar", 
                        position = position_dodge(0.8), width = 0.2) +
  # Scatter points
  ggplot2::geom_point(ggplot2::aes(color = Time_point), position = position_dodge(0.8), 
                      size = 1.2, show.legend = FALSE) +
  ggplot2::scale_fill_manual(values = cb_colors) +
  ggplot2::scale_color_manual(values = cb_colors) +
  ggplot2::facet_wrap(~major, scales = "free_y") +
  ggplot2::theme_classic() +
  ggplot2::theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))

# Iteratively add four groups of P-value annotations using manual control with ggsignif
# Annotation 1: Pre vs Post within the CR group (vertical comparison)
p <- p + ggsignif::geom_signif(
  data = plot_df,
  ggplot2::aes(xmin = 0.8, xmax = 1.2, annotations = sprintf("p=%.3f", stats_df$p_CR_Time[match(major, stats_df$major)])),
  manual = TRUE, y_position = 55, tip_length = 0.01, vjust = -0.5, size = 0.3, textsize = 2.5
)

# Annotation 2: Pre vs Post within the NR group (vertical comparison)
p <- p + ggsignif::geom_signif(
  data = plot_df,
  ggplot2::aes(xmin = 1.8, xmax = 2.2, annotations = sprintf("p=%.3f", stats_df$p_NR_Time[match(major, stats_df$major)])),
  manual = TRUE, y_position = 45, tip_length = 0.01, vjust = -0.5, size = 0.3, textsize = 2.5
)

# Annotation 3: CR vs NR at the pre-treatment stage (horizontal comparison)
p <- p + ggsignif::geom_signif(
  data = plot_df,
  ggplot2::aes(xmin = 0.8, xmax = 1.8, annotations = sprintf("p=%.3f", stats_df$p_Pre_Out[match(major, stats_df$major)])),
  manual = TRUE, y_position = 65, tip_length = 0.01, vjust = -0.5, size = 0.3, textsize = 2.5
)

print(p)
ggsave(filename ="./Figures/Fig/Fig1E_major_proportion_histogram.pdf",height = 6,width = 6,plot = p,dpi = 600,bg="white")
