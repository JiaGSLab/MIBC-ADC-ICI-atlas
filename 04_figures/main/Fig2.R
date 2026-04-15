# ==================================================
# Script: Figure 2 plotting
# Author: Wenhui Chen
# Purpose: Generates T cell subtype figures, marker plots and related statistics for Figure 2.
# Input: Main Seurat or summary objects already generated in the project, plus any file paths referenced in the script.
# Output: Intermediate result files and/or figure files written to output_dir.
# Note: Review file paths before running. Keep original analysis logic unchanged.
# ==================================================
#Fig. 2A T cell UMAP
output_dir = "~/scRNAseq/GPU_IT/QC"
Tcell.obj = qs_read("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/Tcell_obj_intergration_new.qs")
table(Idents(Tcell.obj))
umap_df = data.frame(Embeddings(Tcell.obj, "umap"),
                    Tcell.obj$subcelltype,
                    Tcell.obj$Group)
colnames(umap_df)=c('UMAP_1','UMAP_2','subcelltype',"Group")
subcl=as.character(sort(unique(umap_df$subcelltype)))
umap_df$subcelltype=factor(umap_df$subcelltype,levels = rev(subcl))
p1<- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = subcelltype)) +
      geom_point(size = 1,shape = 16, stroke = 0) +
      scale_color_manual(values = c(colors_91[1:16])) +
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
ggsave(filename ="./Figures/Fig2A_T_subcelltype_Umap.png",height = 4,width = 4,plot = p1,dpi = 600,bg="white")


p2 <- ggplot(data = umap_df, aes(x = UMAP_1, y = UMAP_2, color = subcelltype)) +
      geom_point(size = 0, shape = 16,stroke = 0) +
      theme_void()+
      scale_color_manual(values = c(colors_91[1:16]), name ='') +
      theme( aspect.ratio = 1,
            legend.position ='bottom',
            plot.margin = margin(0,0,0,0) ) +
      guides(color = guide_legend( ncol = 1, override.aes = list(size = 2, alpha = 1) )) +
      theme(legend.text = element_text(size = 5),
            legend.spacing.y = unit(0,'cm'),
            legend.key.height = unit(0,"cm"),
            legend.box.spacing = unit(0,'cm'))

ggsave("./Figures/Fig2A_T_subcelltype_legend.pdf", plot = p2,width = 4, height = 6, dpi = 500 )


#Fig. 2B T cell markers
CD4T = subset(Tcell.obj,idents=subcl[3:9])
CD8T = subset(Tcell.obj,idents=subcl[10:16])
gdT = subset(Tcell.obj,idents=subcl[1:2])
markers.CD4T = cosg(CD4T,
                groups='all',
                assay = 'RNA',
                mu=100,
                remove_lowly_expressed = TRUE,
                expressed_pct = 0.3,
                n_genes_user=30)
View(markers.CD4T$names)
markers.CD8T = cosg(CD8T,
                groups='all',
                assay = 'RNA',
                mu=100,
                remove_lowly_expressed = TRUE,
                expressed_pct = 0.3,
                n_genes_user=30)
View(markers.CD8T$names)
markers.gdT = cosg(gdT,
                groups='all',
                assay = 'RNA',
                mu=100,
                remove_lowly_expressed = TRUE,
                expressed_pct = 0.3,
                n_genes_user=30)
View(markers.gdT$names)

markers.Tcell = c("IL7R","TCF7","B2M","AKNA","GZMK","EOMES","ZNF683","ITGA1","TNFRSF9","HAVCR2","IFIT1","IFIT3","MKI67","TOP2A",
                  "IKZF1","AHNAK","CCR7","S1PR1","CXCL13","PDCD1","KLRB1","RORA","FOXP3","ETS1","CCR8","BATF","MX1","MX2",
                  "CX3CR1","FGFBP2","KLRC1","FCER1G")
p=jjDotPlot(object = Tcell.obj,
          gene = markers.Tcell,
          id = 'subcelltype',
          #gene.order = rev(markers.Tcell),
          #cluster.order = colnames(ma4rkers.Tcell$names),
          ytree = F,dot.min=0,dot.max=3,legend.position="bottom") + 
          coord_flip() +
          theme(axis.text.x = element_text(angle = 45,vjust = 0.2))
p
ggsave( './Figures/Fig/Fig2B_Tcell_markers.pdf',plot=p,width = 4.4,height =7.5)

#Fig. 2C Signature scoring using the gene sets from Wang Linghua et al.
Tcell.signiture = read.csv("~/scRNAseq/T_cell_gene_list.csv",header = T)
#convert dataframe to list
Tcell.signiture = as.list(Tcell.signiture)
Tcell.signiture = lapply(Tcell.signiture,function(x) {
  x = gsub("\\.","-",x)
  x = intersect(x,rownames(Tcell.obj))
})

features = paste0(names(Tcell.signiture),'_score')
do.SignatureHeatmap(Tcell.obj,features = features,
                    fig_output = file.path(output_dir, 'Figures','Fig2C_Tcell_signiture_heatmap.pdf'),
                    plot_margin = c(1.5,1.5,1.5,1.5))

#Fig. 2D Tex/Trm ratio
plot_tex_teff_ratio_final <- function(seurat_obj, 
                                      subcelltype_col = "subcelltype", 
                                      sample_col = "orig.ident", 
                                      group_col = "Group",
                                      tex_name = "t05_CD8_HAVCR2",
                                      teff_name = "t04_CD8_ZNF683",
                                      remove_outliers = FALSE) {
  
  # 1. Calculate the ratio
  ratio_df <- seurat_obj@meta.data %>%
    dplyr::group_by(!!sym(sample_col), !!sym(group_col), !!sym(subcelltype_col)) %>%
    dplyr::summarise(n = n(), .groups = 'drop') %>%
    dplyr::mutate(n = as.numeric(n)) %>% 
    tidyr::pivot_wider(names_from = !!sym(subcelltype_col), values_from = n, values_fill = 0.1) %>%
    dplyr::mutate(ratio = !!sym(tex_name) / !!sym(teff_name))

  # 2. Identify and remove outliers (1.5 × IQR criterion)
  if (remove_outliers) {
    ratio_df <- ratio_df %>%
      dplyr::group_by(!!sym(group_col)) %>%
      dplyr::mutate(
        Q1 = quantile(ratio, 0.25),
        Q3 = quantile(ratio, 0.75),
        IQR = Q3 - Q1,
        is_outlier = ratio < (Q1 - 1.5 * IQR) | ratio > (Q3 + 1.5 * IQR)
      ) %>%
      dplyr::filter(!is_outlier) %>%
      dplyr::ungroup()
  }
  ratio_df[[group_col]] <- factor(ratio_df[[group_col]], levels = c('Pre_NR', 'Pre_CR', 'Post_NR', 'Post_CR'))
  # 3. Define comparison groups
  my_comparisons <- list(c('Pre_NR', 'Pre_CR'), c('Pre_NR', 'Post_NR'))
  #ratio_df[,group_col] = factor(ratio_df[,group_col], levels = c('Pre_NR', 'Pre_CR', 'Post_NR', 'Post_CR'))
  
  # 4. Plot
  p <- ggplot2::ggplot(ratio_df, aes(x = !!sym(group_col), y = ratio, fill = !!sym(group_col))) +
    # Use boxplots with reduced color distraction outside the outlines
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.5, color = "black", size = 0.4) +
    # Modify scatter points: shape = 16 gives solid points without borders
    ggplot2::geom_jitter(width = 0.15, 
                        size = 2, 
                        alpha = 0.9, 
                        shape = 16) + # Remove borders entirely
    ggplot2::scale_fill_manual(values = c("Pre_CR" = "#053E7A", 
                                          "Pre_NR" = "#A52A2A", 
                                          "Post_NR" = "#E46B4F",
                                          "Post_CR"="#8B8B00")) +
    ggpubr::stat_compare_means(comparisons = my_comparisons, 
                              method = "wilcox.test", 
                              label = "p.signif",
                              step.increase = 0.12) +
    ggplot2::theme_classic() +
    ggplot2::labs(y = "Ratio (Tex / Trm)", x = NULL, 
                  title = "") +
    ggplot2::theme(
      text = element_text(size = 12, family = "sans"),
      axis.text.x = element_text(color = "black", size = 10, angle = 45, vjust = 1, hjust = 1),
      axis.text.y = element_text(color = "black", size = 10),
      legend.position = "none"
    )
  
  return(p)
}

# --- Run the analysis ---
p_ratio <- plot_tex_teff_ratio_final(Tcell.obj,remove_outliers = TRUE)
ggsave(file.path(output_dir,"Figures","Fig2D_CD8T_Tex_Trm_Ratio.pdf"), p_ratio, width = 3, height = 4)


plot_subcelltype_pct <- function(seurat_obj,target_subtype = "t05_CD8_HAVCR2", 
                                subcelltype_col = "subcelltype",
                                sample_col = "sample", 
                                group_col = "Group",
                                target_levels = c('Pre_NR', 'Pre_CR', 'Post_NR', 'Post_CR'),
                                remove_outliers = FALSE) {
  
  # 1. Extract metadata
  meta <- seurat_obj@meta.data
  
  # 2. Calculate the total number of cells in each sample (denominator)
  sample_totals <- meta %>%
    dplyr::group_by(!!sym(sample_col)) %>%
    dplyr::summarise(total_cells = dplyr::n(), .groups = 'drop')
  
  # 3. Calculate the number of target subclusters in each sample (numerator)
  subtype_counts <- meta %>%
    dplyr::filter(!!sym(subcelltype_col) == target_subtype) %>%
    dplyr::group_by(!!sym(sample_col), !!sym(group_col)) %>%
    dplyr::summarise(subtype_n = dplyr::n(), .groups = 'drop')
  
  # 4. Merge and calculate percentages
  plot_df <- sample_totals %>%
    dplyr::left_join(subtype_counts, by = sample_col) %>%
    # Fill missing values with 0 when the subcluster was not detected in that sample
    dplyr::mutate(subtype_n = tidyr::replace_na(subtype_n, 0)) %>%
    dplyr::mutate(pct = (subtype_n / total_cells) * 100)

  plot_df <- plot_df %>% 
    dplyr::filter(!base::is.na(!!sym(group_col))) %>%
    dplyr::filter(!!sym(group_col) != "NA")
  

  plot_df[[group_col]] <- factor(plot_df[[group_col]], levels = target_levels)
  
  # Define default comparison groups, keeping only those present in the data
  my_comparisons <- list(c('Pre_NR', 'Pre_CR'), c('Pre_NR', 'Post_NR'))
  #my_comparisons <- base::Filter(function(x) base::all(x %in% final_levels), potential_comparisons)
  
  # 6. Plot
  p <- ggplot2::ggplot(plot_df, aes(x = !!sym(group_col), y = pct, fill = !!sym(group_col))) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.5, color = "black", size = 0.4) +
    ggplot2::geom_jitter(width = 0.15, size = 2, alpha = 0.9, shape = 16) +
    ggplot2::scale_fill_manual(values = c("Pre_CR" = "#053E7A", 
                                          "Pre_NR" = "#A52A2A", 
                                          "Post_NR" = "#E46B4F",
                                          "Post_CR" = "#8B8B00")) +
    ggplot2::theme_classic() +
    ggplot2::labs(y = paste0("Percentage of ", target_subtype, " (%)"), 
                  x = NULL, 
                  title = paste("Distribution of", target_subtype)) +
    ggplot2::theme(
      text = element_text(size = 12, family = "sans"),
      axis.text.x = element_text(color = "black", size = 10, angle = 45, vjust = 1, hjust = 1),
      axis.text.y = element_text(color = "black", size = 10),
      legend.position = "none"
    )
  
  # Add significance annotations
  if (base::length(my_comparisons) > 0) {
    p <- p + ggpubr::stat_compare_means(comparisons = my_comparisons, 
                                        method = "wilcox.test", 
                                        label = "p.signif",
                                        step.increase = 0.12)
  }
  
  return(p)
}
p_Tex_ratio <- plot_subcelltype_pct(Tcell.obj,target_subtype = "t05_CD8_HAVCR2",remove_outliers = FALSE)
ggsave(file.path(output_dir,"Figures","Fig2D_CD8T_Tex_ratio.pdf"), p_Tex_ratio, width = 3, height = 4)
p_Trm_ratio <- plot_subcelltype_pct(Tcell.obj,target_subtype = "t04_CD8_ZNF683",remove_outliers = FALSE)
ggsave(file.path(output_dir,"Figures","Fig2D_CD8T_Trm_ratio.pdf"), p_Trm_ratio, width = 3, height = 4)
p_T03_ratio <- plot_subcelltype_pct(Tcell.obj,target_subtype = "t03_CD8_EOMES",remove_outliers = FALSE)
ggsave(file.path(output_dir,"Figures","Fig2D_CD8T_T03_ratio.pdf"), p_T03_ratio, width = 3, height = 4)




#Fig. 2E comparison of grouped gene-set scores
plot_exhaustion_violin_full <- function(seurat_obj, 
                                        score_col = "Exhaustion_score", 
                                        group_col = "Group") {
  

  # 1. Extract data and ensure Group order is logical (CR first, NR second)
  plot_data <- seurat_obj@meta.data %>%
    dplyr::select(!!sym(group_col), !!sym(score_col))
  
  # Explicitly set factor levels for a cleaner plot
  target_levels <- c('Pre_NR', 'Pre_CR', 'Post_NR', 'Post_CR')
  plot_data[[group_col]] <- factor(plot_data[[group_col]], levels = target_levels)

  # 2. Define the comparison pairs for displaying P values
  # Although Post_CR is shown in the figure, significance annotations are applied only to the following two comparisons
  my_comparisons <- list(c("Pre_CR", "Pre_NR"), c("Post_CR", "Post_NR"),c("Pre_NR", "Post_NR"))

  # 3. Plot
  p <- ggplot2::ggplot(plot_data, aes(x = !!sym(group_col), y = !!sym(score_col), fill = !!sym(group_col))) +
    # Draw violin plots
    ggplot2::geom_violin(trim = FALSE, alpha = 0.7, color = "black", size = 0.3) +
    # Draw internal median and boxplot layers
    ggplot2::geom_boxplot(width = 0.1, color = "black", outlier.shape = NA, fill = "white", size = 0.3) +
    # Use a custom four-group palette (NC style)
    ggplot2::scale_fill_manual(values = c("Pre_CR" = "#053E7A", 
                                          "Pre_NR" = "#A52A2A", 
                                          "Post_NR" = "#E46B4F",
                                          "Post_CR" = "#8B8B00")) +
    # Differential analysis: calculate and annotate only the comparisons listed in my_comparisons
    ggpubr::stat_compare_means(comparisons = my_comparisons, 
                              method = "wilcox.test", 
                              label = "p.format",
                              step.increase = 0.1) +
    ggplot2::theme_classic() +
    ggplot2::labs(y = score_col, x = NULL, title = "") +
    ggplot2::theme(
      text = element_text(size = 12, family = "sans"),
      axis.text.x = element_text(color = "black", size = 10, angle = 45, vjust = 1, hjust = 1),
      legend.position = "none"
    )
  
  return(p)
}

# --- Execute ---
Tex = subset(Tcell.obj,subset=subcelltype=="t05_CD8_HAVCR2")
p_ex <- plot_exhaustion_violin_full(Tex)
ggplot2::ggsave(file.path(output_dir,"Figures","Fig2E_Tex_Exhaustion_score.pdf"), p_ex, width = 3, height = 4)
p_tex_eff <- plot_exhaustion_violin_full(Tex,score_col = "Activation.Effector.function_score")
ggplot2::ggsave(file.path(output_dir,"Figures","Fig2E_Tex_Effector_score.pdf"), p_tex_eff, width = 3, height = 4)
p_tex_cyto <- plot_exhaustion_violin_full(Tex,score_col = "Cytotoxicity_score")
ggplot2::ggsave(file.path(output_dir,"Figures","Fig2E_Tex_Cytotoxicity_score.pdf"), p_tex_cyto, width = 3, height = 4)
Trm = subset(Tcell.obj,subset=subcelltype=="t04_CD8_ZNF683")
Tcm = subset(Tcell.obj,subset=subcelltype=="t03_CD8_EOMES")
p_trm_eff <- plot_exhaustion_violin_full(Trm,score_col = "Activation.Effector.function_score")
ggplot2::ggsave(file.path(output_dir,"Figures","Fig2E_Trm_Effector_score.pdf"), p_trm_eff, width = 3, height = 4)
p_trm_cyto <- plot_exhaustion_violin_full(Trm,score_col = "Cytotoxicity_score")
ggplot2::ggsave(file.path(output_dir,"Figures","Fig2E_Trm_Cytotoxicity_score.pdf"), p_trm_cyto, width = 3, height = 4)
p_trm_exhau <- plot_exhaustion_violin_full(Trm,score_col = "Exhaustion_score")
ggplot2::ggsave(file.path(output_dir,"Figures","Fig2E_Trm_Exhaustion_score.pdf"), p_trm_exhau, width = 3, height = 4)

p_tcm_eff <- plot_exhaustion_violin_full(Tcm,score_col = "Activation.Effector.function_score")
ggplot2::ggsave(file.path(output_dir,"Figures","Fig2E_Tcm_Effector_score.pdf"), p_tcm_eff, width = 3, height = 4)
p_tcm_cyto <- plot_exhaustion_violin_full(Tcm,score_col = "Cytotoxicity_score")
ggplot2::ggsave(file.path(output_dir,"Figures","Fig2E_Tcm_Cytotoxicity_score.pdf"), p_tcm_cyto, width = 3, height = 4)
p_tcm_exhau <- plot_exhaustion_violin_full(Tcm,score_col = "Exhaustion_score")
ggplot2::ggsave(file.path(output_dir,"Figures","Fig2E_Tcm_Exhaustion_score.pdf"), p_tcm_exhau, width = 3, height = 4)

features = c("PDCD1","LAG3","HAVCR2","TIGIT","CTLA4","ENTPD1","BTLA","VSIR",
            "CD28","ICOS","TNFRSF4","TNFRSF9","CD40LG",
            "GZMB","GZMA","PRF1","GNLY","IFNG","TNF")

#Fig2F
features = c("PDCD1","LAG3","HAVCR2","TIGIT","CTLA4","ENTPD1","BTLA","VSIR",
            "CD28","ICOS","TNFRSF4","TNFRSF9","CD40LG",
            "GZMB","GZMA","PRF1","GNLY","IFNG","TNF")
CD8T$newCluster = paste0(CD8T$Group,"_",CD8T$subcelltype)
CD8T.sub = subset(CD8T,idents=c("t03_CD8_EOMES",'t04_CD8_ZNF683','t05_CD8_HAVCR2'))
p=jjDotPlot(object = CD8T.sub,
          gene = features,
          id = 'newCluster',
          #gene.order = rev(markers.Tcell),
          #cluster.order = colnames(markers.Tcell$names),
          ytree = F) + theme(axis.text.x = element_text(angle = 45,vjust = 0.1))
p
ggsave(file.path(output_dir, 'Figures',"Fig",'Fig2F_CD8T_exhaustive_markers.pdf'),plot=p,width = 12,height =5)



