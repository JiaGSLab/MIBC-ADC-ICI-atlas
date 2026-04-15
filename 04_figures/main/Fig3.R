# ==================================================
# Script: Figure 3 plotting
# Author: Wenhui Chen
# Purpose: Generates subtype proportion matrices, correlation analyses and related plots for Figure 3.
# Input: Main Seurat or summary objects already generated in the project, plus any file paths referenced in the script.
# Output: Intermediate result files and/or figure files written to output_dir.
# Note: Review file paths before running. Keep original analysis logic unchanged.
# ==================================================
seu.combined.IT = qs_read("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/seu_combined_IT_filtered_obj_intergration.qs")
malignant.obj = qs_read(file.path(output_dir, 'malignant_obj_intergration.qs'))
#Fig. 3A subcluster correlation heatmap
#Calculate, for each sample, the proportion of each subcluster within its major lineage (myeloid, T cell, B cell, fibroblast, mural cell, endothelial cell); for tumour cells, calculate the mean or median score of each MP
seu.combined.IT$major = ifelse(seu.combined.IT$celltype.new %in% c('SMC','vSMC','Pericyte'),'Mural cells',
                              ifelse(seu.combined.IT$celltype.new=='Endothelial','Endothelial',
                                    ifelse(seu.combined.IT$celltype.new=='Fibroblast','Fibroblast',seu.combined.IT$major)))

calculate_subcelltype_proportion <- function(seu_obj, exclude_major = c("Epithelial","malignant")) {
  # Retrieve metadata
  meta <- seu_obj@meta.data
  
  # Exclude cells from specified major lineages
  meta_filtered <- meta[!(meta$major %in% exclude_major), ]
  
  # Group by sample and major lineage to calculate the total number of cells per lineage
  major_counts <- meta_filtered %>%
    dplyr::group_by(sample, major) %>%
    dplyr::summarise(total_cells = n(), .groups = "drop")
  
  # Group by sample, major lineage, and subcelltype to calculate cell counts per subcelltype
  subcelltype_counts <- meta_filtered %>%
    dplyr::group_by(sample, major, subcelltype) %>%
    dplyr::summarise(cell_count = n(), .groups = "drop")
  
  # Merge data and calculate proportions
  proportion_df <- subcelltype_counts %>%
    dplyr::left_join(major_counts, by = c("sample", "major")) %>%
    dplyr::mutate(proportion = cell_count / total_cells) %>%
    dplyr::select(sample, subcelltype, major, proportion)
  
  # Convert to wide format: rows are samples, columns are subcelltypes, and values are proportions
  proportion_matrix <- proportion_df %>%
    tidyr::pivot_wider(
      id_cols = sample,
      names_from = subcelltype,
      values_from = proportion,
      values_fill = 0
    ) %>%
    tibble::column_to_rownames("sample") %>%
    as.matrix()
  
  return(list(
    proportion_matrix = proportion_matrix,
    proportion_df = proportion_df
  ))
}

# Calculate the proportion matrix
proportion_result <- calculate_subcelltype_proportion(seu.combined.IT)

# Print results
print("Proportion matrix (rows are samples, columns are subcelltypes, values are proportions within the major lineage):")
print(proportion_result$proportion_matrix)


# Save results
write.csv(proportion_result$proportion_matrix, 
          file.path(output_dir, "subcelltype_proportion_matrix.csv"), 
          row.names = FALSE)

write.csv(proportion_result$proportion_df, 
          file.path(output_dir, "subcelltype_proportion_details.csv"), 
          row.names = FALSE)

# View the summary of results
cat("\nSummary of results:\n")
cat("Number of samples:", nrow(proportion_result$proportion_matrix), "\n")
cat("Number of subcelltypes:", ncol(proportion_result$proportion_matrix) - 1, "\n")
cat("Excluded major lineage: Epithelial\n")
cat("Results saved to:", file.path(output_dir, "subcelltype_proportion_matrix.csv"), "\n")

MP.cols = as.character(unique(malignant.obj$MP_ID))
################################################
# Calculate the mean, median, and 75th percentile for each MP column in malignant.obj for each sample
calculate_mp_statistics <- function(obj, mp_cols) {
  # Retrieve metadata
  meta <- obj@meta.data
  
  # Ensure mp_cols are present
  available_cols <- intersect(mp_cols, colnames(meta))
  if (length(available_cols) == 0) {
    stop("No MP.cols columns were found in the metadata")
  }
  
  cat("Available MP columns:", paste(available_cols, collapse = ", "), "\n")
  
  # Calculate summary statistics grouped by sample
  stats_df <- meta %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(
      across(all_of(available_cols), 
            list(mean = ~mean(., na.rm = TRUE),
                  median = ~median(., na.rm = TRUE),
                  q75 = ~quantile(., probs = 0.75, na.rm = TRUE)),
            .names = "{col}_{fn}"))
  
  # Split into three matrices
  mean_matrix <- stats_df %>%
    dplyr::select(sample, ends_with("_mean")) %>%
    tidyr::pivot_longer(cols = -sample, names_to = "mp", values_to = "mean") %>%
    tidyr::pivot_wider(id_cols = sample, names_from = mp, values_from = mean) %>%
    tibble::column_to_rownames("sample") %>%
    as.matrix()
  
  median_matrix <- stats_df %>%
    dplyr::select(sample, ends_with("_median")) %>%
    tidyr::pivot_longer(cols = -sample, names_to = "mp", values_to = "median") %>%
    tidyr::pivot_wider(id_cols = sample, names_from = mp, values_from = median) %>%
    tibble::column_to_rownames("sample") %>%
    as.matrix()
  
  q75_matrix <- stats_df %>%
    dplyr::select(sample, ends_with("_q75")) %>%
    tidyr::pivot_longer(cols = -sample, names_to = "mp", values_to = "q75") %>%
    tidyr::pivot_wider(id_cols = sample, names_from = mp, values_from = q75) %>%
    tibble::column_to_rownames("sample") %>%
    as.matrix()
  
  # Clean column names by removing _mean, _median, and _q75 suffixes
  colnames(mean_matrix) <- gsub("_mean", "", colnames(mean_matrix))
  colnames(median_matrix) <- gsub("_median", "", colnames(median_matrix))
  colnames(q75_matrix) <- gsub("_q75", "", colnames(q75_matrix))
  
  return(list(
    mean_matrix = mean_matrix,
    median_matrix = median_matrix,
    q75_matrix = q75_matrix
  ))
}

cat("Calculating MP statistics...\n")
cat("MP columns:", paste(MP.cols, collapse = ", "), "\n")

mp_stats_result <- calculate_mp_statistics(malignant.obj, MP.cols)

# Print results
cat("\n=== Mean matrix ===\n")
print(mp_stats_result$mean_matrix)
cat("\n=== Median matrix ===\n")
print(mp_stats_result$median_matrix)
cat("\n=== 75th percentile matrix ===\n")
print(mp_stats_result$q75_matrix)

# Save results
write.csv(mp_stats_result$mean_matrix, 
          file.path(output_dir, "malignant_mp_mean_matrix.csv"))

write.csv(mp_stats_result$median_matrix, 
          file.path(output_dir, "malignant_mp_median_matrix.csv"))

write.csv(mp_stats_result$q75_matrix, 
          file.path(output_dir, "malignant_mp_q75_matrix.csv"))

# View the summary of results
cat("\nSummary of results:\n")
cat("Number of samples:", nrow(mp_stats_result$mean_matrix), "\n")
cat("Number of MP columns: ", ncol(mp_stats_result$mean_matrix), "\n")
cat("Results saved to:\n")
cat("  - ", file.path(output_dir, "malignant_mp_mean_matrix.csv"), "\n")
cat("  - ", file.path(output_dir, "malignant_mp_median_matrix.csv"), "\n")
cat("  - ", file.path(output_dir, "malignant_mp_q75_matrix.csv"), "\n")

proportion.matrix = proportion_result$proportion_matrix
non_zero_counts <- apply(proportion.matrix, 2, function(x) sum(x > 0))
keep_cols <- names(non_zero_counts[non_zero_counts >= 10])
proportion.matrix <- proportion.matrix[, keep_cols]
MP.median.matrix = mp_stats_result$median_matrix
all(rownames(proportion.matrix) == rownames(MP.median.matrix))
#all.matrix = cbind(proportion.matrix,MP.median.matrix)
all.matrix = proportion.matrix
View(all.matrix)
# Assume all.matrix has been prepared: rows are samples and columns are cell types
# Ensure the matrix is numeric
all.matrix <- as.matrix(all.matrix)

# 1. Initialize the correlation matrix and P-value matrix
cell_types <- colnames(all.matrix)
n_types <- length(cell_types)

cor_matrix <- matrix(NA, nrow = n_types, ncol = n_types, 
                           dimnames = list(cell_types, cell_types))
p_matrix <- matrix(NA, nrow = n_types, ncol = n_types, 
                         dimnames = list(cell_types, cell_types))

# 2. Use nested loops to calculate pairwise correlations
for (i in 1:n_types) {
  for (j in i:n_types) {
    type_a <- cell_types[i]
    type_b <- cell_types[j]
    
    # Extract the two columns of data
    vec_a <- all.matrix[, type_a]
    vec_b <- all.matrix[, type_b]
    
    # Key step: exclude samples where either value is 0
    keep_indices <- base::which(vec_a != 0 & vec_b != 0)
    
    # Do not calculate correlation if too few samples meet the condition
    if (base::length(keep_indices) >= 10) {
      test_res <- stats::cor.test(vec_a[keep_indices], 
                                  vec_b[keep_indices], 
                                  method = "pearson")
      
      # Fill the symmetric matrix
      cor_matrix[i, j] <- cor_matrix[j, i] <- test_res$estimate
      p_matrix[i, j] <- p_matrix[j, i] <- test_res$p.value
    }
  }
}

cor_matrix[is.na(cor_matrix)] <- 0
View(cor_matrix)
View(p_matrix)
dist_mat <- as.dist(1 - abs(cor_matrix))
fit <- hclust(dist_mat, method = "ward.D2")
cluster_order <- fit$order
cor_sorted <- cor_matrix[cluster_order, cluster_order]
p_sorted <- p_matrix[cluster_order, cluster_order]

pdf(file.path(output_dir, "Figures", "FigS","FigS3B_subcelltype_correlation_heatmap.pdf"),width=7,height=7)
corrplot::corrplot(
  cor_sorted, 
  method = "circle", 
  type = "lower",           # Display only the lower triangle
  order = "original",       # Use the manually defined order
  
  # Significance settings
  p.mat = p_sorted, 
  sig.level = 0.05, 
  insig = "blank",          # Leave blank if not significant
  
  # Style customization
  col = grDevices::colorRampPalette(c("#2166ac", "white", "#b2182b"))(200),
  tl.cex = 0.7,             # Label size
  tl.col = "black", 
  cl.pos = "b",             # Place legend at the bottom
  mar = c(0.2, 0.2, 0.2, 0.2),
  diag = FALSE              # Do not display the diagonal to reduce visual clutter
)
dev.off()

####################
#Lollipop plot
target_cluster <- "t05_CD8_HAVCR2"
if (!target_cluster %in% colnames(cor_matrix)) {
  stop("The target subcluster is not present in the matrix column names. Please check the name.")
}

# 2. Construct a data frame
plot_df <- data.frame(
  subcelltype = rownames(cor_matrix),
  cor_value = cor_matrix[, target_cluster],
  p_value = p_matrix[, target_cluster]
)

# 3. Filter for positive correlations (cor > 0), significant associations (p < 0.05), and exclude the target itself
plot_df_filtered <- plot_df %>%
  dplyr::filter(cor_value > 0 & p_value < 0.05 & subcelltype != target_cluster) %>%
  dplyr::arrange(desc(cor_value)) # Sort by correlation coefficient from high to low

# Ensure subcluster names are ordered by numeric value rather than alphabetically in the plot
plot_df_filtered$subcelltype <- factor(plot_df_filtered$subcelltype, 
                                       levels = plot_df_filtered$subcelltype)
main_color <- "#fc8d59" 

p=ggplot2::ggplot(plot_df_filtered, ggplot2::aes(x = subcelltype, y = cor_value)) +
  # Draw vertical segments
  ggplot2::geom_segment(
    ggplot2::aes(x = subcelltype, xend = subcelltype, y = 0, yend = cor_value),
    color = "grey80", size = 1.2
  ) +
  # Draw circle markers
  ggplot2::geom_point(size = 8, color = main_color) +
  # Add correlation coefficients inside the circles
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.2f", cor_value)), 
    color = "white", size = 3, fontface = "bold"
  ) +
  # Set axes and theme
  ggplot2::theme_classic() +
  ggplot2::labs(
    title = paste("Correlated with", target_cluster),
    y = "Correlation Coefficient (R)",
    x = NULL
  ) +
  # Rotate x-axis labels to avoid overlap
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.title.y = ggplot2::element_text(size = 12),
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
  ) +
  # Limit y-axis range (correlation coefficient cannot exceed 1)
  ggplot2::scale_y_continuous(limits = c(0, 1.1), expand = c(0, 0))
ggsave(file.path(output_dir, 'Figures',"Fig",'Fig3C_Tex__Lollipop_chart.pdf'),plot=p,width = 4.5,height =3)
#########################################################

meta_info <- seu.combined.IT@meta.data[, c("sample", "Group")]

# 2. Remove duplicates to obtain the unique group label for each sample
# This step is important to ensure that each sample maps to only one Group
sample_group_map <- base::unique(meta_info)

# 3. Merge group information into all.matrix
# Assume all.matrix is a matrix or data frame with Sample ID as row names
all_df <- base::as.data.frame(all.matrix)
all_df$sample <- base::rownames(all_df)

# Use an explicit merge call to perform a left join
all_df <- base::merge(all_df, sample_group_map, by = "sample", all.x = TRUE)

# Check the merged result
base::table(all_df$Group)
p <- ggplot2::ggplot(all_df, aes(x = t04_CD8_ZNF683, y = t05_CD8_HAVCR2, color = Group)) +
  ggplot2::geom_point(size = 3, alpha = 0.8) +
  # 1. Explicitly set the regression-line color to match the points
  ggplot2::geom_smooth(method = "lm", se = TRUE, alpha = 0.1) +
  
  # 2. Key fix: change fill to color
  ggplot2::scale_color_manual(values = c("Pre_CR" = "#053E7A", 
                                        "Pre_NR" = "#A52A2A", 
                                        "Post_NR" = "#E46B4F",
                                        "Post_CR" = "#8B8B00")) +
  
  # 3. Optional: add scale_fill_manual if you also want the confidence-interval shading to use the same color mapping
  ggplot2::scale_fill_manual(values = c("Pre_CR" = "#053E7A", 
                                        "Pre_NR" = "#A52A2A", 
                                        "Post_NR" = "#E46B4F",
                                        "Post_CR" = "#8B8B00")) +
  
  ggpubr::stat_cor(method = "pearson") +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "Group-specific Correlation: TRM vs Terminal Tex",
    x = "ZNF683+ CD8T Proportion",
    y = "HAVCR2+ CD8T Proportion"
  )

print(p)
ggsave(file.path(output_dir, 'Figures','Fig3B_Trm_Tex_proportion_scatter.pdf'),plot=p,width = 5,height =4)


#Fig. 3C pseudotime
Tcell.obj = qs_read("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/Tcell_obj_intergration_new.qs")
subcl=sort(as.character(unique(Tcell.obj$subcelltype)))
CD8T = subset(Tcell.obj,idents=subcl[1:7])
subcelltype_analysis <- function(obj,output_dir,prefix,nfeatures = 1500,npcs = 20,res = seq(0.1,1.5,0.1),loess.span=0.3,k.parm=20,n.neighbors=30,min.dist=0.3) {
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
    RunUMAP(.,reduction = "harmony",dims = 1:npcs,n.neighbors=n.neighbors,min.dist=min.dist) %>%
    FindNeighbors(., reduction = "harmony",k.param=k.parm, dims = 1:npcs,verbose = F) %>%
    FindClusters(., resolution = res,algorithm=4,verbose = F,random.seed=1234)
  obj = JoinLayers(obj)
  qs_save(obj,file=file.path(output_dir, paste0(prefix,'_obj_intergration.qs')))
  return(obj)
}
CD8T = subcelltype_analysis(CD8T,output_dir,"CD8T",n.neighbors = 20,min.dist=0.1)
Idents(CD8T) = "subcelltype"
p=DimPlot(CD8T,reduction = 'umap',,cols = colors_90,raster = F,label=F) + theme(aspect.ratio = 1)
ggsave(file.path(output_dir, 'Figures','Fig3C_CD8T_umap.pdf'),plot=p,width = 5,height =3) 
seurat_to_monocle3 = function(seurat_obj) {
  expression_matrix <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
  cell_metadata <- seurat_obj@meta.data
  gene_annotation <- data.frame(
    gene_short_name = rownames(expression_matrix), 
    row.names = rownames(expression_matrix)
  )

  cds <- new_cell_data_set(
    expression_matrix,
    cell_metadata = cell_metadata,
    gene_metadata = gene_annotation
  )
  reducedDims(cds)$UMAP <- seurat_obj@reductions$umap@cell.embeddings
  cds@clusters$UMAP$clusters <- seurat_obj$subcelltype
  return(cds)
}
cds <- seurat_to_monocle3(CD8T)
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE, close_loop = FALSE,
                            learn_graph_control = list(rann.k = 25,minimal_branch_len = 10))
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
          label_leaves = FALSE,  label_branch_points = FALSE)
pdf(file.path(output_dir,'Figures',"Fig3C_CD8T_slingshot_linage.pdf"),width = 6,height=4)
CD8T <- RunSlingshot(CD8T,
                  group.by = "subcelltype",
                  reduction = "umap",
                  start="t04_CD8_ZNF683"
                )
dev.off()
CD8T <- RunDynamicFeatures(
                          CD8T,
                          lineages = c("Lineage1", "Lineage3"),
                          n_candidates = 1000,
                          minfreq = 5,
                          cores = 20
                        )
pdf(file.path(output_dir,'Figures',"Fig3D_CD8T_slingshot_celldimplot.pdf"),width = 20,height=20)
CellDimPlot(CD8T,group.by="subcelltype",reduction='UMAP',add_density = T,pt.size=0.01,alpha=0.3,lineages=c("Lineage1", "Lineage3"),split.by="Group") & NoLegend()
dev.off()
pdf(file.path(output_dir,'Figures',"Fig3D_CD8T_slingshot_DynamicPlot.pdf"),width = 12,height=5)
DynamicPlot(
  CD8T,
  lineages = c("Lineage1", "Lineage3"),
  group.by = 'subcelltype',
  features = c(
    "GZMK","EOMES","HAVCR2", "PDCD1", "LAG3", "TIGIT", "CTLA4", "ENTPD1","TOX"
  ),
  add_point = F,
  pt.size=0,
  compare_lineages = TRUE,
  compare_features = FALSE,
  nrow = 2
)
dev.off()

pdf(file.path(output_dir,'Figures',"Fig3D_CD8T_slingshot_DynamicPlot_Group.pdf"),width = 20,height=20)
DynamicPlot(
  CD8T,
  lineages = c("Lineage1", "Lineage3"),
  group.by = "Group",
  features = c(
    "TCF7","IL7R","GZMK","EOMES","HAVCR2", "PDCD1", "LAG3", "TIGIT", "CTLA4", "ENTPD1","TOX","ITGAE"
  ),
  point_palcolor = c("Pre_CR" = "#053E7A", 
                    "Pre_NR" = "#A52A2A", 
                    "Post_NR" = "#E46B4F",
                    "Post_CR" = "#8B8B00"),
  compare_lineages = TRUE,
  compare_features = FALSE,
  aspect.ratio = 1
)
dev.off()

pdf(file.path(output_dir,'Figures',"Fig3D_CD8T_slingshot_DynamicPlot_Group_Lineage1.pdf"),width = 20,height=20)
DynamicPlot(
  CD8T,
  lineages = c("Lineage1"),
  group.by = "Group",
  features = c(
    "TCF7","IL7R","GZMK","EOMES","HAVCR2", "PDCD1", "LAG3", "TIGIT", "CTLA4", "ENTPD1","TOX","ITGAE"
  ),
  point_palcolor = c("Pre_CR" = "#053E7A", 
                    "Pre_NR" = "#A52A2A", 
                    "Post_NR" = "#E46B4F",
                    "Post_CR" = "#8B8B00"),
  compare_lineages = TRUE,
  compare_features = FALSE
)
dev.off()

pdf(file.path(output_dir,'Figures',"Fig3D_CD8T_slingshot_DynamicPlot_Group_Lineage3.pdf"),width = 20,height=20)
DynamicPlot(
  CD8T,
  lineages = c("Lineage3"),
  group.by = "Group",
  features = c(
    "TCF7","IL7R","GZMK","EOMES","HAVCR2", "PDCD1", "LAG3", "TIGIT", "CTLA4", "ENTPD1","TOX","ITGAE"
  ),
  point_palcolor = c("Pre_CR" = "#053E7A", 
                    "Pre_NR" = "#A52A2A", 
                    "Post_NR" = "#E46B4F",
                    "Post_CR" = "#8B8B00"),
  compare_lineages = TRUE,
  compare_features = FALSE
)
dev.off()

df_lineage1 <- data.frame(
  pseudotime = CD8T$Lineage1,
  Group = CD8T$Group
) %>% 
  dplyr::filter(!is.na(pseudotime))
df_lineage1_Pre = df_lineage1 %>% dplyr::filter(Group == "Pre_CR" | Group == "Pre_NR")
df_lineage1_Post = df_lineage1 %>% dplyr::filter(Group == "Post_CR" | Group == "Post_NR")
p=ggplot(df_lineage1_Pre, aes(x = pseudotime, fill = Group)) +
  geom_density(alpha = 0.5) + # fill transparency
  theme_classic() +
  # Key fix: use scale_fill_manual to match the fill aesthetic
  scale_fill_manual(values = c("Pre_CR" = "#053E7A", 
                              "Pre_NR" = "#A52A2A")) +
  # If you also want the line colors to match, add color mapping together with scale_color_manual
  aes(color = Group) +
  scale_color_manual(values = c("Pre_CR" = "#053E7A", 
                                "Pre_NR" = "#A52A2A")) +
  labs(title = "Lineage 1 Pre",
      x = "Pseudotime",
      y = "Density")
ggsave(filename = file.path(output_dir,'Figures',"Fig3D_CD8T_slingshot_RidgePlot_lineage1_Pre.pdf"),plot=p,width=4,height=2)
p=ggplot(df_lineage1_Post, aes(x = pseudotime, fill = Group)) +
  geom_density(alpha = 0.5) + # fill transparency
  theme_classic() +
  # Key fix: use scale_fill_manual to match the fill aesthetic
  scale_fill_manual(values = c("Post_CR" = "#8B8B00", 
                              "Post_NR" = "#E46B4F")) +
  # If you also want the line colors to match, add color mapping together with scale_color_manual
  aes(color = Group) +
  scale_color_manual(values = c("Post_CR" = "#8B8B00", 
                                "Post_NR" = "#E46B4F")) +
  labs(title = "Lineage 1 Post",
      x = "Pseudotime",
      y = "Density")
ggsave(filename = file.path(output_dir,'Figures',"Fig3Ds_CD8T_slingshot_RidgePlot_lineage1_Post.pdf"),plot=p,width=4,height=2)


df_lineage3 <- data.frame(
  pseudotime = CD8T$Lineage3,
  Group = CD8T$Group
) %>% 
  dplyr::filter(!is.na(pseudotime))
df_lineage3_Pre = df_lineage3 %>% dplyr::filter(Group == "Pre_CR" | Group == "Pre_NR")
df_lineage3_Post = df_lineage3 %>% dplyr::filter(Group == "Post_CR" | Group == "Post_NR")
p=ggplot(df_lineage3_Pre, aes(x = pseudotime, fill = Group)) +
  geom_density(alpha = 0.5) + # fill transparency
  theme_classic() +
  # Key fix: use scale_fill_manual to match the fill aesthetic
  scale_fill_manual(values = c("Pre_CR" = "#053E7A", 
                              "Pre_NR" = "#A52A2A")) +
  # If you also want the line colors to match, add color mapping together with scale_color_manual
  aes(color = Group) +
  scale_color_manual(values = c("Pre_CR" = "#053E7A", 
                                "Pre_NR" = "#A52A2A")) +
  labs(title = "Lineage 3 Pre",
      x = "Pseudotime",
      y = "Density")
ggsave(filename = file.path(output_dir,'Figures',"Fig3D_CD8T_slingshot_RidgePlot_lineage3_Pre.pdf"),plot=p,width=4,height=2)
p=ggplot(df_lineage3_Post, aes(x = pseudotime, fill = Group)) +
  geom_density(alpha = 0.5) + # fill transparency
  theme_classic() +
  # Key fix: use scale_fill_manual to match the fill aesthetic
  scale_fill_manual(values = c("Post_CR" = "#8B8B00", 
                              "Post_NR" = "#E46B4F")) +
  # If you also want the line colors to match, add color mapping together with scale_color_manual
  aes(color = Group) +
  scale_color_manual(values = c("Post_CR" = "#8B8B00", 
                                "Post_NR" = "#E46B4F")) +
  labs(title = "Lineage 3 Post",
      x = "Pseudotime",
      y = "Density")
ggsave(filename = file.path(output_dir,'Figures',"Fig3Ds_CD8T_slingshot_RidgePlot_lineage3_Post.pdf"),plot=p,width=4,height=2)


#Cell–cell communication among positively correlated subclusters
subcelltypes = c("m01_Macro_B2M","m03_Macro_SPP1","s06_Fib_CCL19","m04_Macro_CXCL10","t05_CD8_HAVCR2")
sub.obj = subset(seu.combined.IT,idents=subcelltypes)
sub.obj$samples = sub.obj$sample
obj.PreCR = subset(sub.obj,subset=Group=="Pre_CR")
obj.PreNR = subset(sub.obj,subset=Group=="Pre_NR")

cellchat.analysis=function(obj,output,smooth=TRUE,trim=NULL) {
  cellchat <- createCellChat(object = obj, group.by ="subcelltype", assay = "RNA")
  CellChatDB <- CellChatDB.human
  CellChatDB.use <- subsetDB(CellChatDB)
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  #future::plan("multisession", workers = 16)
  future::plan(sequential)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  if (smooth) {
    cellchat <- smoothData(cellchat, adj = PPI.human)
    if (is.null(trim)) {
      cellchat <- computeCommunProb(cellchat, type = "triMean", trim = NULL,raw.use = FALSE) 
    } else {
      cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = trim,raw.use = FALSE)
    }
  } else {
    if (is.null(trim)) {
      cellchat <- computeCommunProb(cellchat, type = "triMean", trim = NULL,raw.use = TRUE)      
    } else {
      cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = trim,raw.use = TRUE)
    }
  }

  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  future::plan(sequential)
  qs_save(cellchat,output)
  gc()
  return(cellchat)
}

cellchat.Pre_CR = cellchat.analysis(obj.PreCR,'cellchat.Pre_CR_NULL.qs',smooth=FALSE)
gc()
cellchat.Pre_NR = cellchat.analysis(obj.PreNR,'cellchat.Pre_NR_NULL.qs',smooth=FALSE)
gc()
sources = c("m01_Macro_B2M","m03_Macro_SPP1","s06_Fib_CCL19","m04_Macro_CXCL10")
target = "t05_CD8_HAVCR2"
#netVisual_bubble(cellchat.Pre_NR, sources.use = sources, targets.use = target, remove.isolate = TRUE,angle.x = 45)
pdf(file.path(output_dir,'Figures',"Fig3E_CD8T_cellchat_Pre_NR.pdf"),width = 10,height=10)
netVisual_chord_gene(cellchat.Pre_NR, sources.use = sources, targets.use = target,signaling=c('CXCL',"ADGRG","COLLAGEN","NECTIN","GALECTIN"),lab.cex = 0.8,legend.pos.x = 0)
dev.off()
pdf(file.path(output_dir,'Figures',"Fig3E_CD8T_cellchat_Pre_CR.pdf"),width = 10,height=10)
netVisual_chord_gene(cellchat.Pre_CR, sources.use = sources, targets.use = target,signaling=c('CXCL',"ADGRG","COLLAGEN","NECTIN","GALECTIN"),lab.cex = 0.8,legend.pos.x = 0)
dev.off()
object.list = list(Pre_NR = cellchat.Pre_NR, Pre_CR = cellchat.Pre_CR)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
netVisual_bubble(cellchat, sources.use = sources, targets.use = target,  comparison = c(1, 2), angle.x = 45)

#Fig3F Violin plots showing differential expression of myeloid CXCL16, LGALS9, and NECTIN2, fibroblast COL3A1, and T-cell CXCR6 and ADGRG1
myeloid.sub = subset(seu.combined.IT,idents=subcelltypes[c(1,2,4)],subset=Time_point=="Pre")
pdf(file.path(output_dir,"Figures","Fig",'Fig3F_Myeloid_CXCL16_LGALS9_NECTIN2.pdf'),width=8,height=3)
#VlnPlot(myeloid.sub,features = c('CXCL16','LGALS9','NECTIN2'),fill.by='ident',split.by='Group',cols = colors_15,raster = F,pt.size = 0)
VlnPlot2(myeloid.sub,features = c('CXCL16','LGALS9','NECTIN2'),fill.by='ident',split.by='Group',cols = colors_15,raster = F,pt=F,stat.method="wilcox.test")
dev.off()
fib.sub = subset(seu.combined.IT,idents=subcelltypes[3],subset=Time_point=="Pre")
pdf(file.path(output_dir,"Figures","Fig",'Fig3F_Fib_COL3A1.pdf'),width=3,height=3)
VlnPlot(fib.sub,features = c('COL3A1'),fill.by='ident',split.by='Group',cols = colors_15,raster = F,pt.size = 0)
dev.off()
T.sub = subset(seu.combined.IT,idents=subcelltypes[5],subset=Time_point=="Pre")
pdf(file.path(output_dir,"Figures","Fig",'Fig3F_T_CXCR6_ADGRG1.pdf'),width=5,height=3)
VlnPlot(T.sub,features = c('CXCR6','ADGRG1'),fill.by='ident',split.by='Group',cols = colors_15,raster = F,pt.size = 0)
dev.off()