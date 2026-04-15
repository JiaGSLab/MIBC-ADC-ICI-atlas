# ==================================================
# Script: Figure 5 plotting
# Author: Wenhui Chen
# Purpose: Generates malignant cell program and spatial analysis plots used in Figure 5.
# Input: Main Seurat or summary objects already generated in the project, plus any file paths referenced in the script.
# Output: Intermediate result files and/or figure files written to output_dir.
# Note: Review file paths before running. Keep original analysis logic unchanged.
# ==================================================
#################################################
#Fig. 5A MP4 and MP9 are negatively correlated, suggesting that tumour cells switch between two states
malignant.obj = qs_read(file.path(output_dir, 'malignant_obj_intergration.qs'))

plot_df <- malignant.obj@meta.data[, c("Basal_pEMT", "Uro_Differentiation", "Outcome", "Time_point")]

# Then use dplyr for transformation
plot_df <- plot_df %>%
  dplyr::mutate(
    Outcome = factor(Outcome, levels = c("CR", "NR")),
    Time_point = factor(Time_point, levels = c("Pre", "Post"))
  )

# 2. Draw 2D density plot
p5a_density <- ggplot(plot_df, aes(x = Basal_pEMT, y = Uro_Differentiation)) +
  # Use filled 2D density polygons
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  # Use a perceptually uniform color palette
  scale_fill_viridis_c(option = "magma", name = "Cell Density") +
  # Facet by Outcome and Time_point in a 2×2 layout
  facet_grid(Outcome ~ Time_point) + 
  # Set axis ranges to ensure comparable scales across all four panels
  coord_cartesian(xlim = c(min(plot_df$Basal_pEMT), max(plot_df$Basal_pEMT)),
                  ylim = c(min(plot_df$Uro_Differentiation), max(plot_df$Uro_Differentiation))) +
  theme_bw() +
  labs(
    x = "Basal_pEMTScore",
    y = "Uro_Differentiation Score"
  ) +
  theme(
    strip.background = element_rect(fill = "#fff5f0",),
    strip.text = element_text(color = "black",size=6),
    panel.grid = element_blank(),
    aspect.ratio = 1
  )

print(p5a_density)
ggsave(file.path(output_dir,"Figures","Fig","Fig5A_tumor_2D_density_plot.pdf"),plot=p5a_density,width = 3,height = 2.6)

################################################
#Quantify MP4 and MP9 scores in perivascular epithelial cells
library(Rfast)
slides = c('slide01','slide02','slide03')
RCTD.obj = subset(RCTD.obj,subset=orig.ident=="slide01",invert =T) #from Figure4.R
epi.clusters = c('e01_Epi_PSCA','e02_Epi_MKI67','e03_Epi_KRT5','e04_Epi_NDRG1','e05_Epi_FOS',
                'e06_Epi_S100A9','e07_Epi_KRT13','e08_Epi_CXCL10','e09_Epi_PIGR','e10_Epi_MMP7',"e11_Epi_KRT16")
RCTD.obj$first_type  = ifelse(RCTD.obj$first_type %in% epi.clusters, "Epithelial", RCTD.obj$first_type)
RCTD.obj$second_type  = ifelse(RCTD.obj$second_type %in% epi.clusters, "Epithelial", RCTD.obj$second_type)
sample_list <- unique(RCTD.obj$sample_id)
calculate_niche_stats_500 <- function(seurat_obj, s_id) {
  # Extract the single-sample subset
  sub_obj <- subset(seurat_obj, subset = sample_id == s_id)
  coords <- GetTissueCoordinates(sub_obj)
  
  # Identify cell types based on the RCTD result column names
  is_endo <- sub_obj$first_type == "Endothelial" | sub_obj$second_type == "Endothelial"
  is_epi <- sub_obj$first_type == "Epithelial" | sub_obj$second_type == "Epithelial"
  
  endo_coords <- coords[is_endo, c("x", "y")]
  epi_coords <- coords[is_epi, c("x", "y")]
  
  if(nrow(endo_coords) == 0 | nrow(epi_coords) == 0) return(NULL)
  
  # Calculate the distance from each epithelial spot to the nearest vascular spot
  min_dists <- apply(epi_coords, 1, function(p) {
    sqrt(min(rowSums(sweep(endo_coords, 2, p)^2)))
  })
  
  # Construct the base data frame
  res <- data.frame(
    sample = s_id,
    dist = min_dists,
    Basal_pEMT = sub_obj$Basal_pEMT[is_epi],
    Uro_Diff = sub_obj$Uro_Differentiation[is_epi]
  )
  
  # Key change: retain only data within 500 µm
  res <- res %>% dplyr::filter(dist <= 500)
  
  # Bin the data in 50-µm intervals
  res_bin <- res %>%
    dplyr::mutate(bin = cut(dist, breaks = seq(0, 500, by = 50), 
                            labels = FALSE, include.lowest = TRUE)) %>%
    dplyr::filter(!is.na(bin)) %>%
    dplyr::mutate(dist_bin = (bin - 0.5) * 50) %>%
    dplyr::group_by(dist_bin) %>%
    dplyr::summarise(
      mean_Basal = mean(Basal_pEMT, na.rm = TRUE),
      mean_Uro = mean(Uro_Diff, na.rm = TRUE),
      sample = s_id,
      .groups = 'drop'
    )
  
  return(res_bin)
}

# 2. Run iteratively
all_sample_stats_500 <- lapply(sample_list, function(s) calculate_niche_stats_500(RCTD.obj, s)) %>% 
  dplyr::bind_rows()

summary_main <- all_sample_stats_500 %>%
  dplyr::group_by(dist_bin) %>%
  dplyr::summarise(
    # Calculate means
    m_Basal = mean(mean_Basal, na.rm = TRUE),
    m_Uro = mean(mean_Uro, na.rm = TRUE),
    # Calculate standard error (SE) for the shaded region
    se_Basal = sd(mean_Basal, na.rm = TRUE) / sqrt(n()),
    se_Uro = sd(mean_Uro, na.rm = TRUE) / sqrt(n())
  )

# 2. Draw the main plot (mean lines with shaded error bands)
p5f_main <- ggplot2::ggplot(summary_main, aes(x = dist_bin)) +
  # MP4 shaded band and mean line
  ggplot2::geom_ribbon(aes(ymin = m_Basal - se_Basal, ymax = m_Basal + se_Basal), 
                       fill = "#D55E00", alpha = 0.1) +
  ggplot2::geom_line(aes(y = m_Basal, color = "Basal_pEMT (MP4)"), size = 1.5) +
  
  # MP9 shaded band and mean line (scaled by ratio to align the dual y-axes)
  ggplot2::geom_ribbon(aes(ymin = (m_Uro - se_Uro) * ratio, ymax = (m_Uro + se_Uro) * ratio), 
                       fill = "#56B4E9", alpha = 0.1) +
  ggplot2::geom_line(aes(y = m_Uro * ratio, color = "Uro_Diff (MP9)"), size = 1.5) +
  ggplot2::geom_point(data = summary_500, aes(x = dist_bin, y = grand_mean_Basal), 
                      color = "#D55E00", size = 3) +
  ggplot2::geom_point(data = summary_500, aes(x = dist_bin, y = grand_mean_Uro * ratio), 
                      color = "#56B4E9", size = 3) +
  
  # Keep the dual y-axis setting
  ggplot2::scale_y_continuous(
    name = "Basal_pEMT score",
    sec.axis = sec_axis(~ . / ratio, name = "Uro_Differentiation score")
  ) +
  ggplot2::scale_color_manual(values = c("Basal_pEMT (MP4)" = "#D55E00", "Uro_Diff (MP9)" = "#56B4E9")) +
  ggplot2::theme_classic() +
  ggplot2::labs(x = "Distance to Vasculature (µm)") +
  theme(aspect.ratio = 0.7)
print(p5f_main)
ggsave(file.path(output_dir,"Figures","Fig","Fig5B_VisiumHD_MPscore_gradient.pdf"),plot=p5f_main,width = 6,height = 4.8)

p5f_final <- ggplot2::ggplot() +
  # --- Layer 1: thin lines for individual samples (background) ---
  # Basal_pEMT thin lines
  ggplot2::geom_line(data = all_sample_stats_500, 
                     aes(x = dist_bin, y = mean_Basal, group = sample), 
                     color = "#D55E00", alpha = 0.2, size = 0.5) +
  # Uro_Diff thin lines (scaled by ratio)
  ggplot2::geom_line(data = all_sample_stats_500, 
                     aes(x = dist_bin, y = mean_Uro * ratio, group = sample), 
                     color = "#56B4E9", alpha = 0.2, size = 0.5) +
  
  # --- Layer 2: overall trend thick lines (foreground) ---
  # Basal_pEMT mean line
  ggplot2::geom_line(data = summary_500, 
                     aes(x = dist_bin, y = grand_mean_Basal, color = "Basal_pEMT (MP4)"), 
                     size = 1.5) +
  # Uro_Diff mean line (scaled by ratio)
  ggplot2::geom_line(data = summary_500, 
                     aes(x = dist_bin, y = grand_mean_Uro * ratio, color = "Uro_Diff (MP9)"), 
                     size = 1.5) +
  
  # --- Layer 3: point markers ---
  ggplot2::geom_point(data = summary_500, aes(x = dist_bin, y = grand_mean_Basal), 
                      color = "#D55E00", size = 3) +
  ggplot2::geom_point(data = summary_500, aes(x = dist_bin, y = grand_mean_Uro * ratio), 
                      color = "#56B4E9", size = 3) +
  
  # --- Axis settings ---
  ggplot2::scale_y_continuous(
    name = "Basal_pEMT score",
    sec.axis = sec_axis(~ . / ratio, name = "Uro_Differentiation score")
  ) +
  
  # Colors and legend
  ggplot2::scale_color_manual(values = c("Basal_pEMT (MP4)" = "#D55E00", 
                                       "Uro_Diff (MP9)" = "#56B4E9")) +
  
  # Styling
  ggplot2::theme_classic() +
  ggplot2::labs(
    x = "Distance to Vasculature (µm)",
    color = "Programs"
  ) +
  ggplot2::theme(
    axis.title.y.left = element_text(color = "#D55E00", face = "bold", size = 12),
    axis.title.y.right = element_text(color = "#56B4E9", face = "bold", size = 12),
    axis.title.x = element_text(face = "bold"),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )

print(p5f_final)
ggsave(file.path(output_dir,"Figures","FigS","FigS5C_VisiumHD_MPscore_gradient_final.pdf"),plot=p5f_final,width = 6,height = 4.8)

##############################################
sub_obj <- subset(RCTD.obj, subset = sample_id == "pt32")
coords <- GetTissueCoordinates(sub_obj)
library(grDevices)

# 2. Prepare a two-dimensional color matrix (3×3)
# Horizontal axis: Basal_pEMT (MP4); vertical axis: Uro_Diff (MP9)
coords_rot <- data.frame(
  new_x = coords$y,
  new_y = -coords$x,
  Spot = rownames(coords)
)
bottom_left  <- "#e8e8e8" # low / low (gray)
bottom_right <- "#E69F00" # high MP4 / low MP9 (orange)
top_left     <- "#56B4E9" # low MP4 / high MP9 (blue)
top_right    <- "#7030A0" # high / high (purple)

# 2. Manually construct a 4×4 color grid
# First compute the endpoint gradients along the horizontal axis (MP4 dimension)
row_low_mp9  <- colorRampPalette(c(bottom_left, bottom_right))(4)
row_high_mp9 <- colorRampPalette(c(top_left, top_right))(4)

# Then interpolate vertically (MP9 dimension) at each MP4 level
color_grid <- matrix(NA, nrow = 4, ncol = 4)
for (i in 1:4) {
  color_grid[, i] <- colorRampPalette(c(row_low_mp9[i], row_high_mp9[i]))(4)
}

# 3. Map the grid to the data frame generated by expand.grid
# Note: in expand.grid(MP4_lvl = 1:4, MP9_lvl = 1:4), MP4 varies first
bivariate_matrix_4x4 <- expand.grid(MP4_lvl = 1:4, MP9_lvl = 1:4) %>%
  dplyr::mutate(
    # Use mapply to extract colors precisely from the color matrix
    color = mapply(function(m4, m9) color_grid[m9, m4], MP4_lvl, MP9_lvl)
  )

# 4. Merge and order the data
get_4lvl <- function(x) cut(x, breaks = quantile(x, probs = seq(0, 1, 0.25), na.rm = TRUE), labels = FALSE, include.lowest = TRUE)

p5e_main <- ggplot2::ggplot() +
  # Tumour/epithelial layer (ordered by Basal_pEMT so high-score points are drawn on top)
  ggplot2::geom_point(data = dplyr::filter(plot_data, !is_endo), 
                      aes(x = new_x, y = new_y, color = color), 
                      size = 0.2, alpha = 0.7) +
  # Vascular layer (drawn on top, color-blind-friendly)
  ggplot2::geom_point(data = dplyr::filter(plot_data, is_endo), 
                      aes(x = new_x, y = new_y), 
                      color = "#000000", size = 0.4) +
  ggplot2::scale_color_identity() +
  ggplot2::theme_void() + 
  ggplot2::coord_fixed()

# Save the main plot
ggplot2::ggsave(file.path(output_dir,"Figures","Fig","Fig5B_pt32_Spatial_Rotated.png"), plot = p5e_main, 
                width = 7, height = 7, dpi = 600, bg = "white")
p_legend <- ggplot2::ggplot(bivariate_matrix_4x4, aes(x = MP4_lvl, y = MP9_lvl)) +
  ggplot2::geom_tile(aes(fill = color), color = "white", size = 0.3) +
  ggplot2::scale_fill_identity() +
  ggplot2::labs(x = "Basal_pEMT (MP4) →", 
                y = "Uro_Diff (MP9) →") +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(size = 10, face = "bold"),
    panel.grid = ggplot2::element_blank(),
    aspect.ratio = 1
  )

# Save as vector PDF
ggplot2::ggsave(file.path(output_dir,"Figures","Fig","Fig5B_Bivariate_Legend.pdf"), plot = p_legend, 
                width = 3, height = 3)
#############################################
#Differences in MP4-high cells between pre-treatment CR and NR samples
meta_sub <- malignant.obj@meta.data %>%
  dplyr::filter(Time_point == "Pre" & MP4_Status_99 == "MP4_High")

# 2. Define column indices as requested
kegg_cols <- colnames(meta_sub)[46:130]
progeny_cols <- colnames(meta_sub)[c(131:143, 154)]
mp_cols <- colnames(meta_sub)[144:153]
all_target_cols <- c(kegg_cols, progeny_cols, mp_cols)

# 3. Iterate through features: standardization followed by differential testing
diff_results <- lapply(all_target_cols, function(x) {
  val <- meta_sub[[x]]
  
  # Exclude constant columns or columns with all NA values
  if(all(is.na(val)) || sd(val, na.rm = TRUE) == 0) return(NULL)
  
  # Key step: perform Z-score normalization
  # Bring all scores back to the same scale (mean = 0, SD = 1)
  val_std <- scale(val) %>% as.numeric()
  
  # Perform group-wise testing
  cr_indices <- which(meta_sub$Outcome == "CR")
  nr_indices <- which(meta_sub$Outcome == "NR")
  
  test <- stats::wilcox.test(val_std[nr_indices], val_std[cr_indices])
  
  # Calculate the standardized mean difference
  std_mean_diff <- mean(val_std[nr_indices], na.rm = TRUE) - 
                   mean(val_std[cr_indices], na.rm = TRUE)
  
  type <- dplyr::case_when(
    x %in% kegg_cols ~ "KEGG",
    x %in% progeny_cols ~ "PROGENy",
    x %in% mp_cols ~ "MP"
  )
  
  data.frame(
    Feature = x,
    Category = type,
    EffectSize = std_mean_diff, # The x-axis is standardized at this stage
    p_val = test$p.value
  )
}) %>% dplyr::bind_rows()

# 4. Required FDR correction
diff_results <- diff_results %>%
  dplyr::mutate(
    p_adj = p.adjust(p_val, method = "BH"),
    minus_log10p_adj = -log10(p_adj),
    Clean_Name = Feature %>% gsub("KEGG_|PROGENy_", "", .) %>% gsub("_", " ")
  )
label_list <- c("Hypoxia", "Hypoxia_P", "TGFb", "IFN_Response", "Stress_IEG","EGFR",
                "Metabolism of xenobiotics by cytochrome P450","Drug metabolism - cytochrome P450","Drug metabolism - other enzymes", "Glutathione Metabolism")

plot_data <- diff_results %>%
  dplyr::arrange(desc(EffectSize)) %>%
  dplyr::mutate(Rank = row_number())

# 2. Prepare polished labels
# Only label top features with very large effect sizes (NR-high) or very small effect sizes (CR-high)
top_n <- 10
label_data <- plot_data %>%
  dplyr::filter(Rank <= top_n | Rank > (max(Rank) - top_n)) %>%
  # You can also filter by keywords of interest
  dplyr::filter(p_adj < 0.05) 

# 3. Draw the ranked effect size plot
p5c_rank <- ggplot(plot_data, aes(x = Rank, y = EffectSize)) +
  # Background line: EffectSize = 0
  geom_hline(yintercept = 0, linetype = "solid", color = "grey80") +
  # Significance threshold lines, for example |effect size| > 0.2
  geom_hline(yintercept = c(-0.2, 0.2), linetype = "dashed", color = "grey60") +   
  # Draw points with shape mapped to Category and color mapped to Outcome
  geom_point(aes(shape = Category, color = EffectSize > 0), size = 2, alpha = 0.8) +  
  # Colors: orange (#D55E00) indicates NR-high and blue (#56B4E9) indicates CR-high
  scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "#56B4E9"),
                    labels = c("TRUE" = "Higher in NR", "FALSE" = "Higher in CR"),
                    name = "Outcome Bias") +  
  # Shape mapping as requested
  scale_shape_manual(values = c("KEGG" = 15, "MP" = 16, "PROGENy" = 17)) +
  
  # Styling
  theme_classic() +
  labs(
    x = "Rank of Features",
    y = "Effect Size"
  ) +
  theme(
    legend.position = "right",
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14)
  )

print(p5c_rank)
ggsave(file.path(output_dir,"Figures","Fig","Fig5C_MP4_high_pathway_rankplot.pdf"),width = 4,height=4,plot = p5c_rank)


################################
#Differentially expressed genes between MP4-high cells in pre-treatment NR and CR samples
mp4_pre <- subset(malignant.obj, subset = (MP4_Status_99 == "MP4_High") &
                                          (Time_point == "Pre"))
Idents(mp4_pre) = "Outcome"
future::plan(sequential)
DEG.MP4 = FindMarkers(mp4_pre,
                      ident.1= "NR",
                      ident.2 = "CR",
                      test.use="MAST",
                      logfc.threshold=0.05,
                      min.pct = 0.02,
                      latent.vars=c("sample"),
                      random.seed=123)

qs_save(DEG.MP4,"DEG.Pre_MP4_MAST.qs")
DEG.MP4 = qs_read(file.path(output_dir,"DEG.Pre_MP4_MAST.qs"))
View(DEG.MP4)

DEG.MP4_w = FindMarkers(mp4_pre,
                      ident.1= "NR",
                      ident.2 = "CR",
                      logfc.threshold=0.1,
                      min.pct = 0.05,
                      random.seed=123)
DEG.MP4_w$diff.pct = DEG.MP4_w$pct.1-DEG.MP4_w$pct.2
View(DEG.MP4_w)
plot_df <- DEG.MP4_w %>%
  dplyr::mutate(
    # Calculate the difference in expression proportion (assuming pct.1 is NR and pct.2 is CR)
    diff_pct = pct.1 - pct.2,
    # Label categories
    group = dplyr::case_when(
      avg_log2FC > 1 & diff_pct > 0.1 ~ "NR_high",
      avg_log2FC < -1 & diff_pct < -0.1 ~ "CR_high",
      TRUE ~ "NS"
    )
  )
plot_df$gene = rownames(plot_df)
library(GSEABase)
hsa00980.genes <- geneIds(getGmt("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/hsa00980.gmt")[[1]])
hsa00982.genes = geneIds(getGmt("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/hsa00982.gmt")[[1]])
hsa00983.genes = geneIds(getGmt("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/hsa00983.gmt")[[1]])
glutathione.genes = geneIds(getGmt("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/GLUTATHIONE_METABOLISM.gmt")[[1]])
drug.meta.genes = unique(c(hsa00980.genes,hsa00982.genes,hsa00983.genes,glutathione.genes))
drug.meta.genes = intersect(drug.meta.genes,rownames(plot_df[plot_df$group == "NR_high",]))
hypoxia.hallmark = geneIds(getGmt("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/HALLMARK_HYPOXIA.v2026.1.Hs.gmt")[[1]])
hypoxia.harris = geneIds(getGmt("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/HARRIS_HYPOXIA.v2026.1.Hs.gmt")[[1]])
hypoxia.genes = unique(c(hypoxia.hallmark,hypoxia.harris))
hypoxia.genes = intersect(hypoxia.genes,rownames(plot_df[plot_df$group == "CR_high",]))
target_genes <- unique(c("TNFSF10","IGF1R","IGF1",drug.meta.genes,hypoxia.genes))
label_df <- plot_df %>% 
  dplyr::filter(gene %in% target_genes) %>%
  dplyr::mutate(group = factor(group, levels = c("NR_high", "CR_high", "NS")))

p5d_scatter <- ggplot2::ggplot(plot_df, aes(x = diff_pct, y = avg_log2FC)) +
  # Draw background points
  ggplot2::geom_point(aes(color = group), alpha = 0.4, size = 1.5,shape = 16) +
  
  # Set threshold lines
  ggplot2::geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey50") +
  ggplot2::geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  ggplot2::geom_vline(xintercept = c(0), linetype = "solid", color = "black") +
  ggplot2::geom_hline(yintercept = c(0), linetype = "solid", color = "black") +
  ggplot2::geom_point(data = label_df, 
                      aes(fill = group),    # This determines the fill color of the circles
                      shape = 21, 
                      color = "black",      # This determines the color of the black outline
                      size = 3, 
                      stroke = 0.6) +
  
  # Highlight the colors of NR_high and CR_high
  # Orange (#D55E00) denotes NR-driven genes, while blue (#56B4E9) denotes CR-associated genes
  ggplot2::scale_color_manual(values = c("NR_high" = "#D55E00", "CR_high" = "#56B4E9", "NS" = "grey85")) +
  ggplot2::scale_fill_manual(values = c("NR_high" = "#D55E00", "CR_high" = "#56B4E9", "NS" = "grey85")) +
  
  # Automatically label core genes
  ggrepel::geom_text_repel(
    data = dplyr::filter(plot_df, gene %in% target_genes),
    aes(label = gene),
    size = 2,
    box.padding = 0.2, point.padding = 0.3,
    max.overlaps = 15, segment.color = "grey50"
  ) +
  
  # Styling
  ggplot2::theme_classic() +
  ggplot2::labs(
    x = "Difference in Detection Rate (NR - CR)",
    y = "Average Log2 Fold Change (NR / CR)"
  ) +
  ggplot2::theme(
    legend.position = "none",
    axis.title = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    aspect.ratio = 0.7
  )

print(p5d_scatter)
ggsave(file.path(output_dir,"Figures","Fig","Fig5D_MP4_high_DEG.pdf"),plot=p5d_scatter,width = 4,height=3)


tumor.pre = subset(malignant.obj,subset=Time_point == "Pre")
Idents(tumor.pre) = "Outcome"
DEG.pre = FindMarkers(tumor.pre,
                      ident.1= "NR",
                      ident.2 = "CR",
                      logfc.threshold=0.1,
                      min.pct = 0.05,
                      random.seed=123)
DEG.pre$diff.pct = DEG.pre$pct.1-DEG.pre$pct.2
View(DEG.pre)
plot_df <- DEG.pre %>%
  dplyr::mutate(
    # Calculate the difference in expression proportion (assuming pct.1 is NR and pct.2 is CR)
    diff_pct = pct.1 - pct.2,
    # Label categories
    group = dplyr::case_when(
      avg_log2FC > 1 & diff_pct > 0.1 ~ "NR_high",
      avg_log2FC < -1 & diff_pct < -0.1 ~ "CR_high",
      TRUE ~ "NS"
    )
  )
plot_df$gene = rownames(plot_df)
hsa00980.genes <- geneIds(getGmt("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/hsa00980.gmt")[[1]])
hsa00982.genes = geneIds(getGmt("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/hsa00982.gmt")[[1]])
hsa00983.genes = geneIds(getGmt("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/hsa00983.gmt")[[1]])
glutathione.genes = geneIds(getGmt("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/GLUTATHIONE_METABOLISM.gmt")[[1]])
drug.meta.genes = unique(c(hsa00980.genes,hsa00982.genes,hsa00983.genes,glutathione.genes))
drug.meta.genes = intersect(drug.meta.genes,rownames(plot_df[plot_df$group == "NR_high",]))
hypoxia.hallmark = geneIds(getGmt("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/HALLMARK_HYPOXIA.v2026.1.Hs.gmt")[[1]])
hypoxia.harris = geneIds(getGmt("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/HARRIS_HYPOXIA.v2026.1.Hs.gmt")[[1]])
hypoxia.genes = unique(c(hypoxia.hallmark,hypoxia.harris))
hypoxia.genes = intersect(hypoxia.genes,rownames(plot_df[plot_df$group == "CR_high",]))
target_genes <- unique(c("IGF1R","IGF1",drug.meta.genes,hypoxia.genes))
label_df <- plot_df %>% 
  dplyr::filter(gene %in% target_genes) %>%
  dplyr::mutate(group = factor(group, levels = c("NR_high", "CR_high", "NS")))

p5d_scatter <- ggplot2::ggplot(plot_df, aes(x = diff_pct, y = avg_log2FC)) +
  # Draw background points
  ggplot2::geom_point(aes(color = group), alpha = 0.4, size = 1.5,shape = 16) +
  
  # Set threshold lines
  ggplot2::geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey50") +
  ggplot2::geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  ggplot2::geom_vline(xintercept = c(0), linetype = "solid", color = "black") +
  ggplot2::geom_hline(yintercept = c(0), linetype = "solid", color = "black") +
  ggplot2::geom_point(data = label_df, 
                      aes(fill = group),    # This determines the fill color of the circles
                      shape = 21, 
                      color = "black",      # This determines the color of the black outline
                      size = 3, 
                      stroke = 0.6) +
  
  # Highlight the colors of NR_high and CR_high
  # Orange (#D55E00) denotes NR-driven genes, while blue (#56B4E9) denotes CR-associated genes
  ggplot2::scale_color_manual(values = c("NR_high" = "#D55E00", "CR_high" = "#56B4E9", "NS" = "grey85")) +
  ggplot2::scale_fill_manual(values = c("NR_high" = "#D55E00", "CR_high" = "#56B4E9", "NS" = "grey85")) +
  
  # Automatically label core genes
  ggrepel::geom_text_repel(
    data = dplyr::filter(plot_df, gene %in% target_genes),
    aes(label = gene),
    size = 2,
    box.padding = 0.2, point.padding = 0.3,
    max.overlaps = 15, segment.color = "grey50"
  ) +
  
  # Styling
  ggplot2::theme_classic() +
  ggplot2::labs(
    x = "Difference in Detection Rate (NR - CR)",
    y = "Average Log2 Fold Change (NR / CR)"
  ) +
  ggplot2::theme(
    legend.position = "none",
    axis.title = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    aspect.ratio = 0.7
  )

print(p5d_scatter)
ggsave(file.path(output_dir,"Figures","FigS","Fig5D_tumor_DEG.pdf"),plot=p5d_scatter,width = 4,height=3)


##########################################


############################################################
#Fig. 5D show expression differences of genes from hsa00980, hsa00982, and hsa00983 between MP4-high and MP4-low cells in pre-treatment NR samples
library(GSEABase)
hsa00980.genes <- geneIds(getGmt("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/hsa00980.gmt")[[1]])
hsa00982.genes = geneIds(getGmt("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/hsa00982.gmt")[[1]])
hsa00983.genes = geneIds(getGmt("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/hsa00983.gmt")[[1]])
drug.meta.genes = unique(c(hsa00980.genes,hsa00982.genes,hsa00983.genes))
malignant.obj = qs_read(file.path(output_dir, 'malignant_obj_intergration.qs'))
Pre_NR.obj = subset(malignant.obj,subset=Group=="Pre_NR")
drug.meta.genes = intersect(drug.meta.genes,rownames(Pre_NR.obj))

median_mp4 <- median(Pre_NR.obj$Basal_pEMT)
Pre_NR.obj$MP4_status <- ifelse(Pre_NR.obj$Basal_pEMT > median_mp4, "MP4_High", "MP4_Low")

# 3. Run differential analysis for drug.meta.genes
# Set logfc.threshold = 0 explicitly to retain all genes for volcano plotting
Idents(Pre_NR.obj) <- "MP4_status"
de_results <- FindMarkers(Pre_NR.obj, 
                          ident.1 = "MP4_High", 
                          ident.2 = "MP4_Low", 
                          features = drug.meta.genes, # Only analyze genes from drug-metabolism pathways
                          logfc.threshold = 0, 
                          min.pct = 0.01, # slightly filter low-abundance genes
                          only.pos = FALSE)

# Convert to a data frame and add gene names
de_results <- de_results %>% 
  mutate(abs_logFC = abs(avg_log2FC))
log2FC_cut <- 1 # can be adjusted according to the data distribution
p_cut <- 0.01
de_results$gene <- rownames(de_results)

# 5. Label upregulated and downregulated genes
de_results$change <- ifelse(de_results$p_val_adj < p_cut & de_results$avg_log2FC > log2FC_cut, "Up (in MP4 High)",
                    ifelse(de_results$p_val_adj < p_cut & de_results$avg_log2FC < -log2FC_cut, "Down", "Stable"))

de_results$log10p <- -log10(de_results$p_val_adj)
# Handle Inf values and values greater than 100
de_results$log10p[de_results$log10p > 200] <- 200

# 2. Plot
p_volcano <- ggplot(de_results, aes(x = avg_log2FC, y = log10p, color = change)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Up (in MP4 High)" = "#A52A2A", 
                                "Down" = "#053E7A", 
                                "Stable" = "gray80")) +
  # Adjust reference lines according to the upper limit of 100
  geom_vline(xintercept = c(-log2FC_cut, log2FC_cut), lty = 4, col = "black", lwd = 0.6) +
  geom_hline(yintercept = -log10(p_cut), lty = 4, col = "black", lwd = 0.6) +
  
  # Label gene names (note that the y-axis is now mapped to log10p)
  geom_text_repel(data = subset(de_results, log10p > 20 & abs_logFC > 1),
                  aes(label = gene), 
                  size = 3, 
                  box.padding = 0.5, 
                  max.overlaps = 20) +
  
  # Key fix: remove the top and right borders
  theme_classic() + # Use the classic theme, which keeps only the left and bottom axes by default
  theme(
    axis.line = element_line(color = "black"),
    legend.position = "right",
    aspact.ratio =1
    # If you want to keep a theme_bw style while removing the top border:
    # panel.border = element_blank(), 
    # axis.line = element_line(colour = "black")
  ) +
  labs(title = "Metabolic Gene Expression (Capped at 10^100)",
       subtitle = "Pre-treatment Non-Responders (Tumor Cells)",
       x = "log2(Fold Change)", 
       y = "-log10(Adjusted P-value)")

print(p_volcano)
ggsave(file.path(output_dir,"Figures","Fig5D_volcano_metabolic_genes.pdf"), plot = p_volcano, width = 4, height = 3)


################
#Mouse experiment: tumour changes before treatment and one week after treatment
df_raw = read.csv('/home/data/chenwenhui/scRNAseq/GPU_IT/QC/mouse_experiment/Tumor_volume.csv')
print(head(df_raw))
print(table(df_raw$Group, df_raw$Genotype))
df_wide <- df_raw %>%
  tidyr::pivot_wider(
    id_cols = c(ID, Genotype, Group),
    names_from = Time,
    values_from = Tumor_volume
  ) %>%
  dplyr::mutate(
    Growth_Rate = OneWeek / Pre,
    Delta_Volume = OneWeek - Pre,
    Percent_Change = (OneWeek - Pre) / Pre * 100
  )
 
# Set factor order
df_wide$Group <- factor(df_wide$Group, levels = c("Control1", "Control2", "Treated1", "Treated2"))
df_wide$Genotype <- factor(df_wide$Genotype, levels = c("HER2_OE", "Igf1r_KD"))
 
df_raw$Group <- factor(df_raw$Group, levels = c("Control1", "Control2", "Treated1", "Treated2"))
df_raw$Genotype <- factor(df_raw$Genotype, levels = c("HER2_OE", "Igf1r_KD"))
df_raw$Time <- factor(df_raw$Time, levels = c("Pre", "OneWeek"))
 
# Genotype color palette
genotype_colors <- c("HER2_OE" = "#E64B35", "Igf1r_KD" = "#4DBBD5")
 
# ============================================
# 3. Construct paired line plot with color indicating genotype
# ============================================
 
# Calculate group means for the mean trend lines
df_summary <- df_raw %>%
  dplyr::group_by(Group, Genotype, Time) %>%
  dplyr::summarise(
    Mean = mean(Tumor_volume),
    SE = sd(Tumor_volume) / sqrt(dplyr::n()),
    .groups = "drop"
  )
 
# Create x-coordinates for paired connecting lines
# Each Group occupies one position, with Pre and OneWeek placed on the left and right, respectively
group_positions <- c("Control1" = 1, "Control2" = 2, "Treated1" = 3, "Treated2" = 4)
x_offset <- 0.15
 
df_raw <- df_raw %>%
  dplyr::mutate(
    x_base = group_positions[as.character(Group)],
    x_pos = ifelse(Time == "Pre", x_base - x_offset, x_base + x_offset)
  )
 
df_summary <- df_summary %>%
  dplyr::mutate(
    x_base = group_positions[as.character(Group)],
    x_pos = ifelse(Time == "Pre", x_base - x_offset, x_base + x_offset)
  )
 
# Create wide-format data for connecting lines
df_plot_wide <- df_raw %>%
  dplyr::select(ID, Genotype, Group, Time, Tumor_volume, x_base) %>%
  tidyr::pivot_wider(
    id_cols = c(ID, Genotype, Group, x_base),
    names_from = Time,
    values_from = Tumor_volume
  ) %>%
  dplyr::mutate(
    x_pre = x_base - x_offset,
    x_post = x_base + x_offset
  )
 
# Wide-format mean table
df_mean_wide <- df_summary %>%
  tidyr::pivot_wider(
    id_cols = c(Group, Genotype, x_base),
    names_from = Time,
    values_from = c(Mean, SE)
  ) %>%
  dplyr::mutate(
    x_pre = x_base - x_offset,
    x_post = x_base + x_offset,
    Percent_Change = (Mean_OneWeek - Mean_Pre) / Mean_Pre * 100
  )
 
# ============================================
# Construct the plot
# ============================================
 
p <- ggplot2::ggplot() +
  # Layer 1: individual connecting lines
  ggplot2::geom_segment(
    data = df_plot_wide,
    ggplot2::aes(x = x_pre, xend = x_post, y = Pre, yend = OneWeek, color = Genotype),
    linewidth = 0.8, alpha = 0.5
  ) +
  
  # Layer 2: individual scatter points
  ggplot2::geom_point(
    data = df_plot_wide,
    ggplot2::aes(x = x_pre, y = Pre, color = Genotype),
    size = 2.5, alpha = 0.6
  ) +
  ggplot2::geom_point(
    data = df_plot_wide,
    ggplot2::aes(x = x_post, y = OneWeek, color = Genotype),
    size = 2.5, alpha = 0.6
  ) +
  
  # Layer 3: mean trend lines (thick lines)
  ggplot2::geom_segment(
    data = df_mean_wide,
    ggplot2::aes(x = x_pre, xend = x_post, y = Mean_Pre, yend = Mean_OneWeek, color = Genotype),
    linewidth = 2.5, alpha = 1
  ) +
  
  # Layer 4: mean points with error bars
  ggplot2::geom_errorbar(
    data = df_mean_wide,
    ggplot2::aes(x = x_pre, ymin = Mean_Pre - SE_Pre, ymax = Mean_Pre + SE_Pre, color = Genotype),
    width = 0.08, linewidth = 0.8
  ) +
  ggplot2::geom_errorbar(
    data = df_mean_wide,
    ggplot2::aes(x = x_post, ymin = Mean_OneWeek - SE_OneWeek, ymax = Mean_OneWeek + SE_OneWeek, color = Genotype),
    width = 0.08, linewidth = 0.8
  ) +
  ggplot2::geom_point(
    data = df_mean_wide,
    ggplot2::aes(x = x_pre, y = Mean_Pre, color = Genotype),
    size = 4, shape = 16
  ) +
  ggplot2::geom_point(
    data = df_mean_wide,
    ggplot2::aes(x = x_post, y = Mean_OneWeek, color = Genotype),
    size = 4, shape = 16
  ) +
  
  # Layer 5: percentage-change annotations
  ggplot2::geom_text(
    data = df_mean_wide,
    ggplot2::aes(
      x = x_base + 0.35, 
      y = (Mean_Pre + Mean_OneWeek) / 2,
      label = sprintf("%+.0f%%", Percent_Change),
      color = Genotype
    ),
    size = 3.5, fontface = "bold", show.legend = FALSE
  ) +
  
  # Color palette
  ggplot2::scale_color_manual(
    values = genotype_colors,
    labels = c("HER2_OE" = "HER2 OE", "Igf1r_KD" = "Igf1r KD")
  ) +
  
  # X axis
  ggplot2::scale_x_continuous(
    breaks = c(1, 2, 3, 4),
    labels = c("Control1", "Control2", "Treated1", "Treated2"),
    limits = c(0.5, 4.8)
  ) +
  
  # Y axis
  ggplot2::scale_y_continuous(limits = c(0, 650)) +
  
  # Theme
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::labs(
    title = "Tumor Volume Changes: Pre vs One Week Post-Treatment",
    x = "",
    y = expression("Tumor Volume (mm"^3*")"),
    color = "Genotype"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title.y = ggplot2::element_text(face = "bold", size = 12),
    axis.text.x = ggplot2::element_text(face = "bold", size = 11),
    axis.text.y = ggplot2::element_text(size = 10),
    legend.position = "top",
    legend.title = ggplot2::element_text(face = "bold")
  )
 
# Add Pre/1W labels
for (i in 1:4) {
  p <- p + 
    ggplot2::annotate("text", x = i - x_offset, y = -30, label = "Pre", size = 3, color = "gray40") +
    ggplot2::annotate("text", x = i + x_offset, y = -30, label = "1W", size = 3, color = "gray40")
}
 
# Add between-group comparison (growth rate comparison between Treated1 and Treated2)
gr_t1 <- df_wide$Growth_Rate[df_wide$Group == "Treated1"]
gr_t2 <- df_wide$Growth_Rate[df_wide$Group == "Treated2"]
t_result <- t.test(gr_t1, gr_t2)
p_val <- t_result$p.value
 
if (p_val < 0.001) {
  p_text <- "p < 0.001"
} else if (p_val < 0.01) {
  p_text <- sprintf("p = %.3f", p_val)
} else {
  p_text <- sprintf("p = %.2f", p_val)
}
 
p <- p +
  # comparison line
  ggplot2::annotate("segment", x = 3, xend = 4, y = 600, yend = 600, linewidth = 0.8) +
  ggplot2::annotate("segment", x = 3, xend = 3, y = 590, yend = 600, linewidth = 0.8) +
  ggplot2::annotate("segment", x = 4, xend = 4, y = 590, yend = 600, linewidth = 0.8) +
  ggplot2::annotate("text", x = 3.5, y = 620, label = p_text, size = 4, fontface = "bold")
 
print(p)
 
# Save
ggplot2::ggsave(file.path(output_dir,"Figures","Fig","Fig5G_Tumor_Volume_Changes.pdf"), p, width = 8, height = 5)
