###############################################################################
# Figure 2 and Supplementary Figure 2
#
# Consolidated from the original Fig2.R and FigS2.R scripts.
# Plotting logic is retained; only panel organization and output numbering
# are corrected below.
###############################################################################

###############################################################################
# Figure 2
###############################################################################

# 如确实需要复现旧版，可在运行前设置：Sys.setenv(RUN_FIG2_LEGACY = "1")。
if (identical(Sys.getenv("RUN_FIG2_LEGACY", unset = "0"), "1")) {

# Fig2A: T细胞UMAP
output_dir = "~/scRNAseq/GPU_IT/QC"
Tcell.obj = qs_read("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/Tcell_obj_intergration_new_new_new.qs")
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
      theme_classic() + # 使用简洁主题
      theme( plot.background = element_blank(), # 移除背景
            panel.grid.major = element_blank(), # 移除主要网格线
            panel.grid.minor = element_blank(), # 移除次要网格线
            axis.title.x = element_blank(), # 移除x轴标题
            axis.title.y = element_blank(), # 移除y轴标题
            axis.text = element_blank(), # 移除坐标轴刻度标签
            axis.ticks = element_blank(), # 移除坐标轴刻度线
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


# Fig2B: T细胞marker
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

# Fig2C: 打分，参考王凌华文章的基因集
Tcell.signiture = read.csv("~/scRNAseq/T_cell_gene_list.csv",header = T)
#dataframe to list
Tcell.signiture = as.list(Tcell.signiture)
Tcell.signiture = lapply(Tcell.signiture,function(x) {
  x = gsub("\\.","-",x)
  x = intersect(x,rownames(Tcell.obj))
})

features = paste0(names(Tcell.signiture),'_score')
do.SignatureHeatmap(Tcell.obj,features = features,
                    fig_output = file.path(output_dir, 'Figures','Fig2C_Tcell_signiture_heatmap.pdf'),
                    plot_margin = c(1.5,1.5,1.5,1.5))

# Fig2D: CD8 T细胞亚群比例（同时保留辅助Tex/Trm ratio分析）
plot_tex_teff_ratio_final <- function(seurat_obj, 
                                      subcelltype_col = "subcelltype", 
                                      sample_col = "orig.ident", 
                                      group_col = "Group",
                                      tex_name = "t05_CD8_HAVCR2",
                                      teff_name = "t04_CD8_ZNF683",
                                      remove_outliers = FALSE) {
  
  # 1. 计算比值
  ratio_df <- seurat_obj@meta.data %>%
    dplyr::group_by(!!sym(sample_col), !!sym(group_col), !!sym(subcelltype_col)) %>%
    dplyr::summarise(n = n(), .groups = 'drop') %>%
    dplyr::mutate(n = as.numeric(n)) %>% 
    tidyr::pivot_wider(names_from = !!sym(subcelltype_col), values_from = n, values_fill = 0.1) %>%
    dplyr::mutate(ratio = !!sym(tex_name) / !!sym(teff_name))

  # 2. 识别并剔除离群值 (1.5*IQR 准则)
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
  # 3. 统计比较组
  my_comparisons <- list(c('Pre_NR', 'Pre_CR'), c('Pre_NR', 'Post_NR'))
  #ratio_df[,group_col] = factor(ratio_df[,group_col], levels = c('Pre_NR', 'Pre_CR', 'Post_NR', 'Post_CR'))
  
  # 4. 绘图
  p <- ggplot2::ggplot(ratio_df, aes(x = !!sym(group_col), y = ratio, fill = !!sym(group_col))) +
    # 箱线图移除线条外的颜色干扰
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.5, color = "black", size = 0.4) +
    # 散点修改：shape = 16 是纯实心点，没有外环(edge)
    ggplot2::geom_jitter(width = 0.15, 
                        size = 2, 
                        alpha = 0.9, 
                        shape = 16) + # 彻底去除 edge
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

} # end of the first preserved legacy section

###############################################################################
#
# Only panels affected by Outcome / Group are redrawn here:
#   Fig2D-F; FigS2A (Outcome UMAP); FigS2D-F.
# The export dimensions below were measured directly from the current
# Figure 2.ai and Figure S2.ai A4 layouts. All text is Arial 8 pt.
# UMAP points alone are rasterized at 600 dpi; text and legend stay vector.
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

tcell_object_path <- "/home/data/chenwenhui/scRNAseq/GPU_IT/QC/Tcell_obj_intergration_new_new_new.qs"
final_output_root <- "/home/data/chenwenhui/scRNAseq/GPU_IT/QC/Figures/Fig2_S2"
fig2_final_dir <- file.path(final_output_root, "Fig2")
figs2_final_dir <- file.path(final_output_root, "FigS2")
source_data_dir <- file.path(final_output_root, "source_data")

base::dir.create(fig2_final_dir, recursive = TRUE, showWarnings = FALSE)
base::dir.create(figs2_final_dir, recursive = TRUE, showWarnings = FALSE)
base::dir.create(source_data_dir, recursive = TRUE, showWarnings = FALSE)

# AI 中可见内容的实测尺寸，以及本次留约 0.5 mm 安全边距后的导出尺寸。
ai_size_spec <- data.frame(
  panel = c("Fig2D", "Fig2E", "Fig2F", "FigS2A_Outcome", "FigS2D", "FigS2E", "FigS2F"),
  ai_content_width_mm = c(83.29, 85.25, 196.25, 40.64, 113.40, 203.16, 181.95),
  ai_content_height_mm = c(47.38, 47.00, 64.71, 48.08, 79.18, 97.11, 49.82),
  export_width_mm = c(84, 86, 197, 41, 114, 204, 183),
  export_height_mm = c(48, 48, 66, 49, 80, 98, 50),
  font_family = "Arial",
  font_size_pt = 8,
  stringsAsFactors = FALSE
)
utils::write.csv(
  ai_size_spec,
  file.path(source_data_dir, "AI_layout_dimensions_and_export_sizes.csv"),
  row.names = FALSE
)

save_panel <- function(plot_object, file_stem, width_mm, height_mm,
                       preview_dpi = 300) {
  ggplot2::ggsave(
    filename = paste0(file_stem, ".pdf"),
    plot = plot_object,
    device = grDevices::cairo_pdf,
    width = width_mm,
    height = height_mm,
    units = "mm",
    bg = "white",
    limitsize = FALSE
  )
  ggplot2::ggsave(
    filename = paste0(file_stem, ".png"),
    plot = plot_object,
    device = ragg::agg_png,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = preview_dpi,
    bg = "white",
    limitsize = FALSE
  )
}

theme_article_8 <- function() {
  ggplot2::theme_classic(base_size = 8, base_family = "Arial") +
    ggplot2::theme(
      text = ggplot2::element_text(family = "Arial", size = 8, colour = "black"),
      plot.title = ggplot2::element_text(
        family = "Arial", size = 8, face = "plain", hjust = 0.5,
        margin = ggplot2::margin(b = 0.8, unit = "mm")
      ),
      axis.title = ggplot2::element_text(family = "Arial", size = 8, colour = "black"),
      axis.text = ggplot2::element_text(family = "Arial", size = 8, colour = "black"),
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.35),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.35),
      axis.ticks.length = grid::unit(1.1, "mm"),
      legend.title = ggplot2::element_text(family = "Arial", size = 8),
      legend.text = ggplot2::element_text(family = "Arial", size = 8),
      plot.margin = ggplot2::margin(0.6, 0.6, 0.6, 0.6, unit = "mm")
    )
}

display_subtype <- function(x) {
  base::gsub("_gdT_", "_γδT_", x, fixed = TRUE)
}

format_p_plotmath <- function(p_value) {
  if (base::is.na(p_value)) {
    return("p = NA")
  }
  if (p_value < .Machine$double.eps) {
    return("p < 2.22e-16")
  }
  if (p_value < 0.001) {
    p_text <- base::formatC(p_value, digits = 1, format = "e")
    p_text <- base::sub("\\.0e", "e", p_text)
    p_text <- base::sub("e-0", "e-", p_text)
    p_text <- base::sub("e\\+0", "e", p_text)
    p_text <- base::sub("e\\+", "e", p_text)
    return(base::paste0("p = ", p_text))
  }
  if (p_value < 0.01) {
    p_text <- base::sprintf("%.4f", p_value)
  } else {
    p_text <- base::sprintf("%.2f", p_value)
  }
  p_text <- base::sub("0+$", "", p_text)
  p_text <- base::sub("\\.$", "", p_text)
  base::paste0("p = ", p_text)
}

bracket_layers <- function(xmin, xmax, y, p_value, span,
                           text_colour = "black", text_x = NULL) {
  tip <- 0.035 * span
  text_offset <- 0.035 * span
  if (base::is.null(text_x)) {
    text_x <- (xmin + xmax) / 2
  }
  bracket_only <- list(
    ggplot2::annotate(
      "segment", x = xmin, xend = xmax, y = y, yend = y,
      colour = "black", linewidth = 0.4
    ),
    ggplot2::annotate(
      "segment", x = xmin, xend = xmin, y = y, yend = y - tip,
      colour = "black", linewidth = 0.4
    ),
    ggplot2::annotate(
      "segment", x = xmax, xend = xmax, y = y, yend = y - tip,
      colour = "black", linewidth = 0.4
    )
  )
  if (base::is.na(p_value)) {
    return(bracket_only)
  }
  c(
    bracket_only,
    list(
    ggplot2::annotate(
      "text",
      x = text_x,
      y = y + text_offset,
      label = format_p_plotmath(p_value),
      parse = FALSE,
      family = "Arial",
      size = 8 / 2.845276,
      colour = text_colour,
      vjust = 0
    )
    )
  )
}

p_text_colour <- function(p_value) {
  if (!base::is.na(p_value) && p_value < 0.05) "#D55E00" else "black"
}

Tcell.obj <- qs2::qs_read(tcell_object_path)
metadata <- Tcell.obj@meta.data
metadata$cell_barcode <- base::rownames(metadata)



# -----------------------------------------------------------------------------
# Fig2D and FigS2E: sample-level T-cell subtype proportions
# One-sided exact permutation test; direction follows the observed median shift,
# matching the test used in the current AI figure. Two-sided p is also exported.
# -----------------------------------------------------------------------------

prop_data <- metadata |>
  dplyr::transmute(
    sample = base::as.character(.data$sample),
    Time_point = base::as.character(.data$Time_point),
    Outcome = base::as.character(.data$Outcome),
    subcelltype = base::as.character(.data$subcelltype)
  ) |>
  dplyr::filter(
    !base::is.na(.data$sample),
    !base::is.na(.data$Time_point),
    !base::is.na(.data$Outcome),
    !base::is.na(.data$subcelltype)
  ) |>
  dplyr::count(
    .data$sample, .data$Time_point, .data$Outcome, .data$subcelltype,
    name = "count"
  ) |>
  dplyr::group_by(.data$sample) |>
  dplyr::mutate(
    total_T_cells = base::sum(.data$count),
    percentage = .data$count / .data$total_T_cells * 100
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    Time_point = base::factor(.data$Time_point, levels = c("Pre", "Post")),
    Outcome = base::factor(.data$Outcome, levels = c("CR", "NR"))
  )

run_one_permutation <- function(test_data, group_column, level1, level2,
                                contrast_name) {
  local_data <- data.frame(
    percentage = base::as.numeric(test_data$percentage),
    test_group = base::factor(
      base::as.character(test_data[[group_column]]),
      levels = c(level1, level2)
    )
  )
  local_data <- local_data[!base::is.na(local_data$percentage) &
                             !base::is.na(local_data$test_group), , drop = FALSE]
  value1 <- local_data$percentage[local_data$test_group == level1]
  value2 <- local_data$percentage[local_data$test_group == level2]

  if (base::length(value1) < 2L || base::length(value2) < 2L) {
    return(data.frame(
      contrast = contrast_name,
      level1 = level1,
      level2 = level2,
      n_level1 = base::length(value1),
      n_level2 = base::length(value2),
      median_level1 = if (base::length(value1)) stats::median(value1) else NA_real_,
      median_level2 = if (base::length(value2)) stats::median(value2) else NA_real_,
      mean_level1 = if (base::length(value1)) base::mean(value1) else NA_real_,
      mean_level2 = if (base::length(value2)) base::mean(value2) else NA_real_,
      difference_level2_minus_level1 = NA_real_,
      alternative_one_sided = NA_character_,
      p_one_sided = NA_real_,
      p_two_sided = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  difference <- stats::median(value2) - stats::median(value1)
  alternative <- if (difference >= 0) "less" else "greater"
  one_sided_test <- coin::independence_test(
    percentage ~ test_group,
    data = local_data,
    distribution = "exact",
    alternative = alternative
  )
  two_sided_test <- coin::independence_test(
    percentage ~ test_group,
    data = local_data,
    distribution = "exact",
    alternative = "two.sided"
  )

  data.frame(
    contrast = contrast_name,
    level1 = level1,
    level2 = level2,
    n_level1 = base::length(value1),
    n_level2 = base::length(value2),
    median_level1 = stats::median(value1),
    median_level2 = stats::median(value2),
    mean_level1 = base::mean(value1),
    mean_level2 = base::mean(value2),
    difference_level2_minus_level1 = difference,
    alternative_one_sided = alternative,
    p_one_sided = base::as.numeric(coin::pvalue(one_sided_test)),
    p_two_sided = base::as.numeric(coin::pvalue(two_sided_test)),
    stringsAsFactors = FALSE
  )
}

build_proportion_stats <- function(one_subtype) {
  subtype_data <- dplyr::filter(prop_data, .data$subcelltype == one_subtype)
  result <- dplyr::bind_rows(
    run_one_permutation(
      dplyr::filter(subtype_data, .data$Outcome == "CR"),
      "Time_point", "Pre", "Post", "CR_Pre_vs_Post"
    ),
    run_one_permutation(
      dplyr::filter(subtype_data, .data$Outcome == "NR"),
      "Time_point", "Pre", "Post", "NR_Pre_vs_Post"
    ),
    run_one_permutation(
      dplyr::filter(subtype_data, .data$Time_point == "Pre"),
      "Outcome", "CR", "NR", "Pre_CR_vs_NR"
    ),
    run_one_permutation(
      dplyr::filter(subtype_data, .data$Time_point == "Post"),
      "Outcome", "CR", "NR", "Post_CR_vs_NR"
    )
  )
  result$subcelltype <- one_subtype
  result[, c("subcelltype", base::setdiff(base::colnames(result), "subcelltype"))]
}

all_subtypes <- base::levels(metadata$subcelltype)
if (base::is.null(all_subtypes)) {
  all_subtypes <- base::sort(base::unique(base::as.character(metadata$subcelltype)))
}
proportion_stats <- dplyr::bind_rows(base::lapply(all_subtypes, build_proportion_stats))

fig2d_subtypes <- c("t03_CD8_EOMES", "t04_CD8_ZNF683", "t05_CD8_HAVCR2")
figs2e_subtypes <- c(
  "t01_CD8_IL7R", "t02_CD8_AKNA", "t06_CD8_ISG", "t07_CD8_MKI67",
  "t08_CD4_IL7R", "t09_CD4_TCF7", "t10_CD4_CXCL13",
  "t11_CD4_RORA", "t12_CD4_ETS1", "t13_CD4_TNFRSF9", "t14_CD4_ISG",
  "t15_gdT_CX3CR1", "t16_gdT_FCER1G"
)

utils::write.csv(
  dplyr::filter(prop_data, .data$subcelltype %in% fig2d_subtypes),
  file.path(source_data_dir, "Fig2D_Tcell_subtype_proportions_plot_data.csv"),
  row.names = FALSE
)
utils::write.csv(
  dplyr::filter(proportion_stats, .data$subcelltype %in% fig2d_subtypes),
  file.path(source_data_dir, "Fig2D_Tcell_subtype_proportions_permutation_stats.csv"),
  row.names = FALSE
)
utils::write.csv(
  dplyr::filter(prop_data, .data$subcelltype %in% figs2e_subtypes),
  file.path(source_data_dir, "FigS2E_Tcell_subtype_proportions_plot_data.csv"),
  row.names = FALSE
)
utils::write.csv(
  dplyr::filter(proportion_stats, .data$subcelltype %in% figs2e_subtypes),
  file.path(source_data_dir, "FigS2E_Tcell_subtype_proportions_permutation_stats.csv"),
  row.names = FALSE
)

timepoint_colours <- c("Pre" = "#E69F00", "Post" = "#56B4E9")

build_proportion_plot <- function(one_subtype, show_y_title = FALSE) {
  plot_data <- dplyr::filter(prop_data, .data$subcelltype == one_subtype)
  plot_summary <- plot_data |>
    dplyr::group_by(.data$Outcome, .data$Time_point) |>
    dplyr::summarise(
      mean_percentage = base::mean(.data$percentage),
      se_percentage = stats::sd(.data$percentage) / base::sqrt(dplyr::n()),
      .groups = "drop"
    )
  plot_stats <- dplyr::filter(proportion_stats, .data$subcelltype == one_subtype)
  get_p <- function(contrast_name) {
    plot_stats$p_one_sided[base::match(contrast_name, plot_stats$contrast)]
  }

  data_max <- base::max(
    plot_data$percentage,
    plot_summary$mean_percentage + plot_summary$se_percentage,
    na.rm = TRUE
  )
  span <- base::max(data_max, 1)
  y_cr <- data_max + 0.08 * span
  y_nr <- data_max + 0.22 * span
  y_pre <- data_max + 0.42 * span
  y_upper <- data_max + 0.62 * span

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$Outcome, y = .data$percentage)
  ) +
    ggplot2::geom_col(
      data = plot_summary,
      ggplot2::aes(
        x = .data$Outcome,
        y = .data$mean_percentage,
        fill = .data$Time_point,
        group = .data$Time_point
      ),
      position = ggplot2::position_dodge(width = 0.72),
      width = 0.64,
      colour = "black",
      linewidth = 0.35,
      alpha = 0.7,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_errorbar(
      data = plot_summary,
      ggplot2::aes(
        x = .data$Outcome,
        ymin = .data$mean_percentage - .data$se_percentage,
        ymax = .data$mean_percentage + .data$se_percentage,
        group = .data$Time_point
      ),
      position = ggplot2::position_dodge(width = 0.72),
      width = 0.16,
      linewidth = 0.35,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = .data$Time_point, group = .data$Time_point),
      position = ggplot2::position_jitterdodge(
        jitter.width = 0.035,
        jitter.height = 0,
        dodge.width = 0.72,
        seed = 20260803
      ),
      shape = 16,
      size = 0.9,
      alpha = 0.95,
      stroke = 0,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = timepoint_colours,
      breaks = c("Pre", "Post"),
      name = "Time point",
      drop = FALSE
    ) +
    ggplot2::scale_colour_manual(values = timepoint_colours, guide = "none") +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::coord_cartesian(ylim = c(0, y_upper), clip = "off") +
    ggplot2::labs(
      x = NULL,
      y = if (show_y_title) "% in T cells" else NULL,
      title = display_subtype(one_subtype)
    ) +
    theme_article_8() +
    ggplot2::theme(
      legend.key.height = grid::unit(3.2, "mm"),
      legend.key.width = grid::unit(3.2, "mm"),
      legend.spacing.x = grid::unit(1.1, "mm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5, unit = "mm")
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        title.position = "top", nrow = 1,
        override.aes = list(alpha = 0.7, colour = "black")
      )
    )

  p <- p + bracket_layers(
    0.82, 1.18, y_cr, get_p("CR_Pre_vs_Post"), span,
    p_text_colour(get_p("CR_Pre_vs_Post")), text_x = 1.12
  )
  p <- p + bracket_layers(
    1.82, 2.18, y_nr, get_p("NR_Pre_vs_Post"), span,
    p_text_colour(get_p("NR_Pre_vs_Post")), text_x = 2.05
  )
  p <- p + bracket_layers(
    0.82, 1.82, y_pre, get_p("Pre_CR_vs_NR"), span,
    p_text_colour(get_p("Pre_CR_vs_NR")), text_x = 1.40
  )
  p
}

fig2d_plot_list <- base::Map(
  function(x, show_y) build_proportion_plot(x, show_y_title = show_y),
  fig2d_subtypes,
  c(TRUE, FALSE, FALSE)
)
fig2d_plot <- patchwork::wrap_plots(
  fig2d_plot_list,
  nrow = 1,
  guides = "collect"
) &
  ggplot2::theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = ggplot2::element_blank()
  )

save_panel(
  fig2d_plot,
  file.path(fig2_final_dir, "Fig2D_Tcell_subtype_proportions"),
  width_mm = 84,
  height_mm = 48
)

figs2e_plot_list <- base::Map(
  function(x, show_y) build_proportion_plot(x, show_y_title = show_y),
  figs2e_subtypes,
  base::rep(FALSE, base::length(figs2e_subtypes))
)
figs2e_top <- patchwork::wrap_plots(figs2e_plot_list[1:7], ncol = 7)
figs2e_bottom <- patchwork::wrap_plots(
  c(figs2e_plot_list[8:13], list(patchwork::guide_area())),
  ncol = 7
)
figs2e_plot <- (figs2e_top / figs2e_bottom) +
  patchwork::plot_layout(guides = "collect", heights = c(1, 1)) &
  ggplot2::theme(
    legend.position = "right",
    legend.direction = "vertical"
  )

save_panel(
  figs2e_plot,
  file.path(figs2_final_dir, "FigS2E_Tcell_subtype_proportions"),
  width_mm = 204,
  height_mm = 98
)


# -----------------------------------------------------------------------------
# Fig2E and FigS2F: signature score violin plots and Wilcoxon rank-sum tests
# -----------------------------------------------------------------------------

score_columns <- c(
  "Activation.Effector.function_score",
  "Cytotoxicity_score",
  "Exhaustion_score"
)
score_titles <- c(
  "Activation.Effector.function_score" = "Activation Effector function",
  "Cytotoxicity_score" = "Cytotoxicity",
  "Exhaustion_score" = "Exhaustion"
)
score_subtypes <- c("t03_CD8_EOMES", "t04_CD8_ZNF683", "t05_CD8_HAVCR2")

score_data <- metadata |>
  dplyr::filter(.data$subcelltype %in% score_subtypes) |>
  dplyr::select(
    dplyr::all_of(c(
      "cell_barcode", "orig.ident", "sample", "Sample_ID", "Time_point",
      "Outcome", "Group", "subcelltype"
    )),
    dplyr::all_of(score_columns)
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(score_columns),
    names_to = "score_name",
    values_to = "score_value"
  ) |>
  dplyr::mutate(
    Group = base::factor(
      base::as.character(.data$Group),
      levels = c("Pre_NR", "Pre_CR", "Post_NR", "Post_CR")
    ),
    score_title = base::unname(score_titles[.data$score_name])
  )

score_comparisons <- list(
  c("Pre_NR", "Pre_CR"),
  c("Post_NR", "Post_CR"),
  c("Pre_NR", "Post_NR")
)

build_score_stats <- function(one_subtype, one_score) {
  local_data <- dplyr::filter(
    score_data,
    .data$subcelltype == one_subtype,
    .data$score_name == one_score
  )
  dplyr::bind_rows(base::lapply(score_comparisons, function(one_comparison) {
    value1 <- local_data$score_value[local_data$Group == one_comparison[1]]
    value2 <- local_data$score_value[local_data$Group == one_comparison[2]]
    test_result <- stats::wilcox.test(value1, value2, exact = FALSE)
    data.frame(
      subcelltype = one_subtype,
      score_name = one_score,
      score_title = base::unname(score_titles[one_score]),
      contrast = base::paste(one_comparison, collapse = "_vs_"),
      level1 = one_comparison[1],
      level2 = one_comparison[2],
      n_level1 = base::length(value1),
      n_level2 = base::length(value2),
      median_level1 = stats::median(value1),
      median_level2 = stats::median(value2),
      mean_level1 = base::mean(value1),
      mean_level2 = base::mean(value2),
      p_value = test_result$p.value,
      method = "Wilcoxon rank-sum test (unpaired, two-sided)",
      stringsAsFactors = FALSE
    )
  }))
}

score_stats <- dplyr::bind_rows(base::lapply(score_subtypes, function(one_subtype) {
  dplyr::bind_rows(base::lapply(score_columns, function(one_score) {
    build_score_stats(one_subtype, one_score)
  }))
}))

utils::write.csv(
  dplyr::filter(score_data, .data$subcelltype == "t05_CD8_HAVCR2"),
  file.path(source_data_dir, "Fig2E_t05_signature_scores_plot_data.csv"),
  row.names = FALSE
)
utils::write.csv(
  dplyr::filter(score_stats, .data$subcelltype == "t05_CD8_HAVCR2"),
  file.path(source_data_dir, "Fig2E_t05_signature_scores_Wilcoxon_stats.csv"),
  row.names = FALSE
)
utils::write.csv(
  dplyr::filter(score_data, .data$subcelltype %in% c("t03_CD8_EOMES", "t04_CD8_ZNF683")),
  file.path(source_data_dir, "FigS2F_t03_t04_signature_scores_plot_data.csv"),
  row.names = FALSE
)
utils::write.csv(
  dplyr::filter(score_stats, .data$subcelltype %in% c("t03_CD8_EOMES", "t04_CD8_ZNF683")),
  file.path(source_data_dir, "FigS2F_t03_t04_signature_scores_Wilcoxon_stats.csv"),
  row.names = FALSE
)

group_colours <- c(
  "Pre_NR" = "#A52A2A",
  "Pre_CR" = "#053E7A",
  "Post_NR" = "#E46B4F",
  "Post_CR" = "#8B8B00"
)

build_score_violin <- function(one_subtype, one_score, show_y_title = FALSE) {
  local_data <- dplyr::filter(
    score_data,
    .data$subcelltype == one_subtype,
    .data$score_name == one_score
  )
  local_stats <- dplyr::filter(
    score_stats,
    .data$subcelltype == one_subtype,
    .data$score_name == one_score
  )
  get_p <- function(contrast_name) {
    local_stats$p_value[base::match(contrast_name, local_stats$contrast)]
  }

  data_min <- base::min(local_data$score_value, na.rm = TRUE)
  data_max <- base::max(local_data$score_value, na.rm = TRUE)
  span <- base::max(data_max - data_min, 0.5)
  y_pre <- data_max + 0.08 * span
  y_post <- data_max + 0.36 * span
  y_time <- data_max + 0.64 * span
  y_upper <- data_max + 0.90 * span
  y_lower <- data_min - 0.04 * span

  p <- ggplot2::ggplot(
    local_data,
    ggplot2::aes(x = .data$Group, y = .data$score_value, fill = .data$Group)
  ) +
    ggplot2::geom_violin(
      trim = FALSE,
      scale = "width",
      alpha = 0.7,
      colour = "black",
      linewidth = 0.3
    ) +
    ggplot2::geom_boxplot(
      width = 0.11,
      outlier.shape = NA,
      fill = "white",
      alpha = 0.8,
      colour = "black",
      linewidth = 0.3
    ) +
    ggplot2::scale_fill_manual(values = group_colours, drop = FALSE) +
    ggplot2::coord_cartesian(ylim = c(y_lower, y_upper), clip = "off") +
    ggplot2::labs(
      x = NULL,
      y = if (show_y_title) "Signature Score" else NULL,
      title = base::unname(score_titles[one_score])
    ) +
    theme_article_8() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        family = "Arial", size = 8, angle = 45, hjust = 1, vjust = 1
      ),
      legend.position = "none",
      plot.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5, unit = "mm")
    )

  p <- p + bracket_layers(
    1, 2, y_pre, get_p("Pre_NR_vs_Pre_CR"), span, "black", text_x = 1.65
  )
  p <- p + bracket_layers(
    3, 4, y_post, get_p("Post_NR_vs_Post_CR"), span, "black", text_x = 3.45
  )
  p <- p + bracket_layers(
    1, 3, y_time, get_p("Pre_NR_vs_Post_NR"), span, "black", text_x = 2.05
  )
  p
}

fig2e_plot_list <- base::Map(
  function(x, show_y) build_score_violin("t05_CD8_HAVCR2", x, show_y),
  score_columns,
  c(TRUE, FALSE, FALSE)
)
fig2e_plot <- patchwork::wrap_plots(fig2e_plot_list, nrow = 1) +
  patchwork::plot_annotation(
    caption = "t05_CD8_HAVCR2",
    theme = ggplot2::theme(
      plot.caption = ggplot2::element_text(
        family = "Arial", size = 8, colour = "black", hjust = 0.5,
        margin = ggplot2::margin(t = 0.8, unit = "mm")
      )
    )
  )

save_panel(
  fig2e_plot,
  file.path(fig2_final_dir, "Fig2E_t05_signature_scores"),
  width_mm = 86,
  height_mm = 48
)

figs2f_t03_list <- base::Map(
  function(x, show_y) build_score_violin("t03_CD8_EOMES", x, show_y),
  score_columns,
  c(TRUE, FALSE, FALSE)
)
figs2f_t04_list <- base::Map(
  function(x, show_y) build_score_violin("t04_CD8_ZNF683", x, show_y),
  score_columns,
  c(FALSE, FALSE, FALSE)
)
figs2f_t03 <- patchwork::wrap_plots(figs2f_t03_list, nrow = 1) +
  patchwork::plot_annotation(
    caption = "t03_CD8_EOMES",
    theme = ggplot2::theme(
      plot.caption = ggplot2::element_text(
        family = "Arial", size = 8, colour = "black", hjust = 0.5,
        margin = ggplot2::margin(t = 0.8, unit = "mm")
      )
    )
  )
figs2f_t04 <- patchwork::wrap_plots(figs2f_t04_list, nrow = 1) +
  patchwork::plot_annotation(
    caption = "t04_CD8_ZNF683",
    theme = ggplot2::theme(
      plot.caption = ggplot2::element_text(
        family = "Arial", size = 8, colour = "black", hjust = 0.5,
        margin = ggplot2::margin(t = 0.8, unit = "mm")
      )
    )
  )
figs2f_plot <- patchwork::wrap_plots(list(figs2f_t03, figs2f_t04), nrow = 1)

save_panel(
  figs2f_plot,
  file.path(figs2_final_dir, "FigS2F_t03_t04_signature_scores"),
  width_mm = 183,
  height_mm = 50
)


# -----------------------------------------------------------------------------
# FigS2A Outcome UMAP: raster points (600 dpi) + vector title/legend
# UMAP coordinates are intentionally not written to CSV.
# -----------------------------------------------------------------------------

umap_embedding <- base::as.data.frame(Seurat::Embeddings(Tcell.obj, reduction = "umap"))
base::colnames(umap_embedding)[1:2] <- c("UMAP_1", "UMAP_2")
umap_embedding$Outcome <- base::factor(
  base::as.character(metadata$Outcome),
  levels = c("CR", "NR")
)
outcome_colours <- c("CR" = "#D55E00", "NR" = "#0072B2")

figs2a_outcome_plot <- ggplot2::ggplot(
  umap_embedding,
  ggplot2::aes(x = .data$UMAP_1, y = .data$UMAP_2, colour = .data$Outcome)
) +
  ggrastr::geom_point_rast(
    shape = 16,
    size = 0.18,
    alpha = 0.8,
    stroke = 0,
    raster.dpi = 600
  ) +
  ggplot2::scale_colour_manual(
    values = outcome_colours,
    breaks = c("CR", "NR"),
    drop = FALSE
  ) +
  ggplot2::coord_fixed() +
  ggplot2::labs(title = "Outcome", colour = NULL) +
  ggplot2::theme_void(base_size = 8, base_family = "Arial") +
  ggplot2::theme(
    text = ggplot2::element_text(family = "Arial", size = 8, colour = "black"),
    plot.title = ggplot2::element_text(
      family = "Arial", size = 8, face = "plain", hjust = 0.5,
      margin = ggplot2::margin(b = 0.5, unit = "mm")
    ),
    legend.position = "inside",
    legend.position.inside = c(0.5, 0.03),
    legend.direction = "horizontal",
    legend.justification = c(0.5, 0),
    legend.background = ggplot2::element_rect(
      fill = scales::alpha("white", 0.78), colour = NA
    ),
    legend.key = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(family = "Arial", size = 8),
    legend.spacing.x = grid::unit(1.2, "mm"),
    legend.margin = ggplot2::margin(0, 0, 0, 0),
    plot.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm")
  ) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(
      nrow = 1,
      override.aes = list(size = 1.6, alpha = 1)
    )
  )

save_panel(
  figs2a_outcome_plot,
  file.path(figs2_final_dir, "FigS2A_Outcome_UMAP"),
  width_mm = 41,
  height_mm = 49,
  preview_dpi = 600
)


# -----------------------------------------------------------------------------
# FigS2D: Ro/e heatmap with the top stacked-fraction bars
# -----------------------------------------------------------------------------

roe_cluster_levels <- base::levels(metadata$subcelltype)
roe_group_levels <- base::levels(base::factor(base::as.character(metadata$Group)))
roe_count_matrix <- base::table(
  base::factor(base::as.character(metadata$subcelltype), levels = roe_cluster_levels),
  base::factor(base::as.character(metadata$Group), levels = roe_group_levels)
)

roe_fisher_rows <- list()
roe_index <- 0L
for (cluster_index in base::seq_along(roe_cluster_levels)) {
  for (group_index in base::seq_along(roe_group_levels)) {
    this_count <- roe_count_matrix[cluster_index, group_index]
    other_in_cluster <- base::sum(roe_count_matrix[cluster_index, ]) - this_count
    other_in_group <- base::sum(roe_count_matrix[, group_index]) - this_count
    other_cells <- base::sum(roe_count_matrix) - this_count - other_in_cluster - other_in_group
    fisher_result <- stats::fisher.test(base::matrix(
      c(this_count, other_in_cluster, other_in_group, other_cells),
      ncol = 2
    ))
    roe_index <- roe_index + 1L
    roe_fisher_rows[[roe_index]] <- data.frame(
      subcelltype = roe_cluster_levels[cluster_index],
      Group = roe_group_levels[group_index],
      count = base::as.numeric(this_count),
      odds_ratio = base::as.numeric(fisher_result$estimate),
      p_value = fisher_result$p.value,
      stringsAsFactors = FALSE
    )
  }
}
roe_stats <- dplyr::bind_rows(roe_fisher_rows)
roe_stats$FDR_BH <- stats::p.adjust(roe_stats$p_value, method = "BH")

roe_count_data <- base::as.data.frame(roe_count_matrix, stringsAsFactors = FALSE)
base::colnames(roe_count_data) <- c("subcelltype", "Group", "count")
roe_fraction_data <- roe_count_data |>
  dplyr::group_by(.data$subcelltype) |>
  dplyr::mutate(fraction = .data$count / base::sum(.data$count)) |>
  dplyr::ungroup()

utils::write.csv(
  roe_count_data,
  file.path(source_data_dir, "FigS2D_Roe_cell_counts.csv"),
  row.names = FALSE
)
utils::write.csv(
  roe_fraction_data,
  file.path(source_data_dir, "FigS2D_Roe_topbar_fraction_data.csv"),
  row.names = FALSE
)
utils::write.csv(
  roe_stats,
  file.path(source_data_dir, "FigS2D_Roe_Fisher_exact_stats.csv"),
  row.names = FALSE
)

roe_group_colours <- c(
  "Pre_NR" = "#A52A2A",
  "Pre_CR" = "#053E7A",
  "Post_NR" = "#E46B4F",
  "Post_CR" = "#8B8B00"
)
roe_fraction_data$subcelltype <- base::factor(
  roe_fraction_data$subcelltype,
  levels = roe_cluster_levels
)
roe_fraction_data$Group <- base::factor(
  roe_fraction_data$Group,
  levels = base::rev(roe_group_levels)
)
roe_stats$subcelltype <- base::factor(roe_stats$subcelltype, levels = roe_cluster_levels)
roe_stats$Group <- base::factor(roe_stats$Group, levels = base::rev(roe_group_levels))
roe_stats$odds_ratio_capped <- base::pmin(base::pmax(roe_stats$odds_ratio, 0), 3)
roe_stats$significant_label <- base::ifelse(
  roe_stats$FDR_BH < 0.05,
  base::sprintf("%.2f", roe_stats$odds_ratio),
  ""
)
roe_stats$label_colour <- base::ifelse(
  roe_stats$odds_ratio > 2 | roe_stats$odds_ratio < 0.3,
  "white",
  "black"
)

roe_top_plot <- ggplot2::ggplot(
  roe_fraction_data,
  ggplot2::aes(
    x = .data$subcelltype,
    y = .data$fraction,
    fill = .data$Group
  )
) +
  ggplot2::geom_col(
    width = 0.98,
    colour = "white",
    linewidth = 0.2,
    position = ggplot2::position_stack(reverse = TRUE)
  ) +
  ggplot2::scale_fill_manual(values = roe_group_colours, guide = "none") +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    expand = c(0, 0)
  ) +
  ggplot2::scale_x_discrete(expand = c(0, 0)) +
  ggplot2::labs(x = NULL, y = "Fraction") +
  theme_article_8() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0.4, 0.4, 0, 0.4, unit = "mm")
  )

roe_heat_plot <- ggplot2::ggplot(
  roe_stats,
  ggplot2::aes(
    x = .data$subcelltype,
    y = .data$Group,
    fill = .data$odds_ratio_capped
  )
) +
  ggplot2::geom_tile(colour = "grey90", linewidth = 0.25) +
  ggplot2::geom_text(
    ggplot2::aes(label = .data$significant_label, colour = .data$label_colour),
    family = "Arial",
    size = 8 / 2.845276,
    show.legend = FALSE
  ) +
  ggplot2::scale_colour_identity() +
  ggplot2::scale_fill_gradientn(
    colours = c("#2166AC", "white", "#B2182B"),
    values = scales::rescale(c(0, 1, 3)),
    limits = c(0, 3),
    breaks = c(0, 1, 2, 3),
    labels = c("0", "1", "2", "≥3"),
    oob = scales::squish,
    name = "Ro/e",
    guide = ggplot2::guide_colourbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(28, "mm"),
      barheight = grid::unit(3, "mm"),
      ticks = FALSE
    )
  ) +
  ggplot2::scale_x_discrete(
    labels = function(x) display_subtype(x),
    expand = c(0, 0)
  ) +
  ggplot2::scale_y_discrete(expand = c(0, 0)) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::theme_minimal(base_size = 8, base_family = "Arial") +
  ggplot2::theme(
    text = ggplot2::element_text(family = "Arial", size = 8, colour = "black"),
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      family = "Arial", size = 8, colour = "black",
      angle = 45, hjust = 1, vjust = 1
    ),
    axis.text.y = base::suppressWarnings(ggplot2::element_text(
      family = "Arial", size = 8,
      colour = roe_group_colours[base::rev(roe_group_levels)]
    )),
    axis.ticks = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(family = "Arial", size = 8),
    legend.text = ggplot2::element_text(family = "Arial", size = 8),
    legend.margin = ggplot2::margin(0, 0, 0, 0),
    plot.margin = ggplot2::margin(0, 0.4, 0.4, 0.4, unit = "mm")
  )

figs2d_plot <- roe_top_plot / roe_heat_plot +
  patchwork::plot_layout(heights = c(0.74, 2.26))

save_panel(
  figs2d_plot,
  file.path(figs2_final_dir, "FigS2D_Tcell_subtype_Roe"),
  width_mm = 114,
  height_mm = 80
)


# -----------------------------------------------------------------------------
# Fig2F: CD8 inhibitory-receptor / effector-gene dot plot
# -----------------------------------------------------------------------------

fig2f_features <- c(
  "PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4", "ENTPD1", "BTLA", "VSIR",
  "CD28", "ICOS", "TNFRSF4", "TNFRSF9", "CD40LG",
  "GZMB", "GZMA", "PRF1", "GNLY", "IFNG", "TNF"
)
fig2f_group_order <- c("Pre_NR", "Pre_CR", "Post_NR", "Post_CR")
fig2f_subtype_order <- c("t05_CD8_HAVCR2", "t04_CD8_ZNF683", "t03_CD8_EOMES")
fig2f_id_order <- base::unlist(base::lapply(fig2f_group_order, function(one_group) {
  base::paste(one_group, fig2f_subtype_order, sep = "___")
}))

fig2f_cells <- metadata$cell_barcode[
  base::as.character(metadata$subcelltype) %in% fig2f_subtype_order
]
fig2f_object <- Tcell.obj[, fig2f_cells]
fig2f_object$newCluster <- base::factor(
  base::paste(
    base::as.character(fig2f_object$Group),
    base::as.character(fig2f_object$subcelltype),
    sep = "___"
  ),
  levels = fig2f_id_order
)

# Keep the original Fig2.R calculation exactly: jjDotPlot() and its scaled
# expression values.  Only the final ggplot styling/export below is customized.
fig2f_dot_data <- scRNAtoolVis::jjDotPlot(
  object = fig2f_object,
  assay = "RNA",
  slot = "data",
  gene = fig2f_features,
  id = "newCluster",
  ytree = FALSE,
  scale = TRUE,
  col.min = -2.5,
  col.max = 2.5
)$data
base::colnames(fig2f_dot_data)[
  base::colnames(fig2f_dot_data) == "gene"
] <- "features.plot"
fig2f_dot_data$avg.exp <- base::as.numeric(fig2f_dot_data$avg.exp)
fig2f_dot_data$avg.exp.scaled <- base::as.numeric(fig2f_dot_data$avg.exp.scaled)
fig2f_dot_data$pct.exp <- base::as.numeric(fig2f_dot_data$pct.exp)
fig2f_dot_data$pct.exp.plot <- base::pmin(fig2f_dot_data$pct.exp, 80)
fig2f_dot_data$id <- base::factor(
  base::as.character(fig2f_dot_data$id),
  levels = base::rev(fig2f_id_order)
)
fig2f_dot_data$features.plot <- base::factor(
  base::as.character(fig2f_dot_data$features.plot),
  levels = fig2f_features
)
fig2f_dot_data$Group <- base::factor(
  base::sub("___.*$", "", base::as.character(fig2f_dot_data$id)),
  levels = fig2f_group_order
)
fig2f_dot_data$subcelltype <- base::factor(
  base::sub("^.*___", "", base::as.character(fig2f_dot_data$id)),
  levels = base::rev(fig2f_subtype_order)
)

utils::write.csv(
  fig2f_dot_data,
  file.path(source_data_dir, "Fig2F_CD8_marker_dotplot_source_data.csv"),
  row.names = FALSE
)

fig2f_highlight_data <- data.frame(
  Group = base::factor(c("Pre_NR", "Post_NR"), levels = fig2f_group_order),
  xmin = 0.5,
  xmax = 6.5,
  ymin = 2.5,
  ymax = 3.5
)

fig2f_plot <- ggplot2::ggplot(
  fig2f_dot_data,
  ggplot2::aes(
    x = .data$features.plot,
    y = .data$subcelltype,
    size = .data$pct.exp.plot,
    fill = .data$avg.exp.scaled
  )
) +
  ggplot2::geom_point(shape = 21, colour = "black", stroke = 0.35) +
  ggplot2::geom_rect(
    data = fig2f_highlight_data,
    ggplot2::aes(
      xmin = .data$xmin, xmax = .data$xmax,
      ymin = .data$ymin, ymax = .data$ymax
    ),
    inherit.aes = FALSE,
    fill = NA,
    colour = "red",
    linewidth = 0.5
  ) +
  ggplot2::scale_size_continuous(
    name = "Fraction of cells in clusters (%)",
    range = c(0.15, 3.4),
    limits = c(0, 80),
    breaks = c(0, 20, 40, 60, 80)
  ) +
  ggplot2::scale_fill_gradientn(
    name = "Mean expression in clusters",
    colours = c("#FFF7F3", "#E5A08D", "#B2182B"),
    values = scales::rescale(c(-1.5, 0, 2.2)),
    limits = c(-1.5, 2.2),
    breaks = c(-1, 0, 1, 2),
    oob = scales::squish,
    guide = ggplot2::guide_colourbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(31, "mm"),
      barheight = grid::unit(3.2, "mm"),
      ticks = TRUE,
      order = 2
    )
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(.data$Group),
    scales = "free_y",
    space = "free_y",
    switch = "y",
    drop = FALSE
  ) +
  ggplot2::scale_x_discrete(
    drop = FALSE,
    expand = ggplot2::expansion(add = 0.5)
  ) +
  ggplot2::scale_y_discrete(
    drop = FALSE,
    labels = function(x) display_subtype(x),
    expand = ggplot2::expansion(add = 0.5)
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::theme_classic(base_size = 8, base_family = "Arial") +
  ggplot2::theme(
    text = ggplot2::element_text(family = "Arial", size = 8, colour = "black"),
    panel.grid.major = ggplot2::element_line(colour = "grey85", linewidth = 0.25),
    panel.border = ggplot2::element_rect(
      fill = NA, colour = "black", linewidth = 0.5
    ),
    axis.line = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      family = "Arial", size = 8, angle = 45, hjust = 1, vjust = 1
    ),
    axis.text.y = ggplot2::element_text(family = "Arial", size = 8),
    strip.background = ggplot2::element_blank(),
    strip.placement = "outside",
    strip.text.y.left = ggplot2::element_text(
      family = "Arial", size = 8, angle = 0, hjust = 1,
      margin = ggplot2::margin(r = 2, unit = "mm")
    ),
    panel.spacing.y = grid::unit(0, "mm"),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = ggplot2::element_text(family = "Arial", size = 8),
    legend.text = ggplot2::element_text(family = "Arial", size = 8),
    legend.margin = ggplot2::margin(0, 0, 0, 1, unit = "mm"),
    plot.margin = ggplot2::margin(1, 1, 2, 1, unit = "mm")
  ) +
  ggplot2::guides(
    size = ggplot2::guide_legend(
      title.position = "top",
      nrow = 1,
      byrow = TRUE,
      order = 1,
      override.aes = list(fill = "grey60")
    )
  )

save_panel(
  fig2f_plot,
  file.path(fig2_final_dir, "Fig2F_CD8_marker_dotplot"),
  width_mm = 197,
  height_mm = 66
)


# Continue preserving the remainder of the historical merged script.
if (identical(Sys.getenv("RUN_FIG2_LEGACY", unset = "0"), "1")) {

# --- 运行分析 ---
p_ratio <- plot_tex_teff_ratio_final(Tcell.obj,remove_outliers = TRUE)
ggsave(file.path(output_dir,"Figures","Fig2D_CD8T_Tex_Trm_Ratio.pdf"), p_ratio, width = 3, height = 4)


plot_subcelltype_pct <- function(seurat_obj,target_subtype = "t05_CD8_HAVCR2", 
                                subcelltype_col = "subcelltype",
                                sample_col = "sample", 
                                group_col = "Group",
                                target_levels = c('Pre_NR', 'Pre_CR', 'Post_NR', 'Post_CR'),
                                remove_outliers = FALSE) {
  
  # 1. 提取元数据
  meta <- seurat_obj@meta.data
  
  # 2. 计算每个样本的总细胞数 (分母)
  sample_totals <- meta %>%
    dplyr::group_by(!!sym(sample_col)) %>%
    dplyr::summarise(total_cells = dplyr::n(), .groups = 'drop')
  
  # 3. 计算目标亚群在每个样本中的数量 (分子)
  subtype_counts <- meta %>%
    dplyr::filter(!!sym(subcelltype_col) == target_subtype) %>%
    dplyr::group_by(!!sym(sample_col), !!sym(group_col)) %>%
    dplyr::summarise(subtype_n = dplyr::n(), .groups = 'drop')
  
  # 4. 合并并计算百分比
  plot_df <- sample_totals %>%
    dplyr::left_join(subtype_counts, by = sample_col) %>%
    # 将缺失值（该样本没检测到该亚群）填补为 0
    dplyr::mutate(subtype_n = tidyr::replace_na(subtype_n, 0)) %>%
    dplyr::mutate(pct = (subtype_n / total_cells) * 100)

  plot_df <- plot_df %>% 
    dplyr::filter(!base::is.na(!!sym(group_col))) %>%
    dplyr::filter(!!sym(group_col) != "NA")
  

  plot_df[[group_col]] <- factor(plot_df[[group_col]], levels = target_levels)
  
  # 定义默认比较组（仅保留数据中存在的）
  my_comparisons <- list(c('Pre_NR', 'Pre_CR'), c('Pre_NR', 'Post_NR'))
  #my_comparisons <- base::Filter(function(x) base::all(x %in% final_levels), potential_comparisons)
  
  # 6. 绘图
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
  
  # 添加显著性
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




# Fig2E: t05_CD8_HAVCR2分组基因集评分比较
plot_exhaustion_violin_full <- function(seurat_obj, 
                                        score_col = "Exhaustion_score", 
                                        group_col = "Group") {
  

  # 1. 提取并确保 Group 的顺序 (按照逻辑排列：CR在前, NR在后)
  plot_data <- seurat_obj@meta.data %>%
    dplyr::select(!!sym(group_col), !!sym(score_col))
  
  # 显式设置 Level 顺序，确保出图美观
  target_levels <- c('Pre_NR', 'Pre_CR', 'Post_NR', 'Post_CR')
  plot_data[[group_col]] <- factor(plot_data[[group_col]], levels = target_levels)

  # 2. 定义需要显示 P 值的比较对
  # 虽然图中显示了 Post_CR，但我们只对以下两组做显著性标注
  my_comparisons <- list(c("Pre_CR", "Pre_NR"), c("Post_CR", "Post_NR"),c("Pre_NR", "Post_NR"))

  # 3. 绘图
  p <- ggplot2::ggplot(plot_data, aes(x = !!sym(group_col), y = !!sym(score_col), fill = !!sym(group_col))) +
    # 绘制小提琴图
    ggplot2::geom_violin(trim = FALSE, alpha = 0.7, color = "black", size = 0.3) +
    # 绘制内部中位数线和箱线
    ggplot2::geom_boxplot(width = 0.1, color = "black", outlier.shape = NA, fill = "white", size = 0.3) +
    # 自定义四组配色 (NC 风格)
    ggplot2::scale_fill_manual(values = c("Pre_CR" = "#053E7A", 
                                          "Pre_NR" = "#A52A2A", 
                                          "Post_NR" = "#E46B4F",
                                          "Post_CR" = "#8B8B00")) +
    # 差异分析：仅对 my_comparisons 中的组合进行计算和标注
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

# --- 执行 ---
Tex = subset(Tcell.obj,subset=subcelltype=="t05_CD8_HAVCR2")
p_ex <- plot_exhaustion_violin_full(Tex)
ggplot2::ggsave(file.path(output_dir,"Figures","Fig2E_Tex_Exhaustion_score.pdf"), p_ex, width = 3, height = 4)
p_tex_eff <- plot_exhaustion_violin_full(Tex,score_col = "Activation.Effector.function_score")
ggplot2::ggsave(file.path(output_dir,"Figures","Fig2E_Tex_Effector_score.pdf"), p_tex_eff, width = 3, height = 4)
p_tex_cyto <- plot_exhaustion_violin_full(Tex,score_col = "Cytotoxicity_score")
ggplot2::ggsave(file.path(output_dir,"Figures","Fig2E_Tex_Cytotoxicity_score.pdf"), p_tex_cyto, width = 3, height = 4)
features = c("PDCD1","LAG3","HAVCR2","TIGIT","CTLA4","ENTPD1","BTLA","VSIR",
            "CD28","ICOS","TNFRSF4","TNFRSF9","CD40LG",
            "GZMB","GZMA","PRF1","GNLY","IFNG","TNF")

# Fig2F: CD8 T细胞耗竭/激活/细胞毒marker
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




###############################################################################
# Supplementary Figure 2
###############################################################################

Tcell.obj = qs_read("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/Tcell_obj_intergration_new_new_new.qs")
table(Idents(Tcell.obj))

# FigS2A: sample、Time point和Outcome UMAP

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
       theme_classic() + # 使用简洁主题
       theme( plot.background = element_blank(), # 移除背景
              panel.grid.major = element_blank(), # 移除主要网格线
              panel.grid.minor = element_blank(), # 移除次要网格线
              axis.title.x = element_blank(), # 移除x轴标题
              axis.title.y = element_blank(), # 移除y轴标题
              axis.text = element_blank(), # 移除坐标轴刻度标签
              axis.ticks = element_blank(), # 移除坐标轴刻度线
              axis.line = element_blank(),
              legend.position ="none",
              aspect.ratio = 1)
p1
ggsave(filename ="./Figures/FigS/FigS02A_sample_Umap.png",height = 4,width = 4,plot = p1,dpi = 600,bg="white")


p1<- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = Time_point)) +
       geom_point(size = 1,shape = 16, stroke = 0) +
       scale_color_manual(values = colors_15) +
       theme_classic() + # 使用简洁主题
       theme( plot.background = element_blank(), # 移除背景
              panel.grid.major = element_blank(), # 移除主要网格线
              panel.grid.minor = element_blank(), # 移除次要网格线
              axis.title.x = element_blank(), # 移除x轴标题
              axis.title.y = element_blank(), # 移除y轴标题
              axis.text = element_blank(), # 移除坐标轴刻度标签
              axis.ticks = element_blank(), # 移除坐标轴刻度线
              axis.line = element_blank(),
              legend.position ="none",
              aspect.ratio = 1)
p1
ggsave(filename ="./Figures/FigS/FigS02A_timepoint_Umap.png",height = 4,width = 4,plot = p1,dpi = 600,bg="white")

p1<- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = Outcome)) +
       geom_point(size = 1,shape = 16, stroke = 0) +
       scale_color_manual(values = colors_15[3:4]) +
       theme_classic() + # 使用简洁主题
       theme( plot.background = element_blank(), # 移除背景
              panel.grid.major = element_blank(), # 移除主要网格线
              panel.grid.minor = element_blank(), # 移除次要网格线
              axis.title.x = element_blank(), # 移除x轴标题
              axis.title.y = element_blank(), # 移除y轴标题
              axis.text = element_blank(), # 移除坐标轴刻度标签
              axis.ticks = element_blank(), # 移除坐标轴刻度线
              axis.line = element_blank(),
              legend.position ="none",
              aspect.ratio = 1)
p1
ggsave(filename ="./Figures/FigS/FigS02A_Outcome_Umap.png",height = 4,width = 4,plot = p1,dpi = 600,bg="white")

##############
# FigS2B: CD4、FOXP3、CD8A和TRDC FeaturePlot
pdf(file = file.path(output_dir,"Figures","FigS","FigS02B_Tcell_marker_FeaturePlot.pdf"),width = 4,height = 4)
FeaturePlot(Tcell.obj,features = c("CD4","FOXP3","CD8A","TRDC"),
            cols=c('lightgrey',"red"),ncol = 2,raster = T,pt.size=3) & NoLegend() & NoAxes()
dev.off()

pdf(file = file.path(output_dir,"Figures","FigS","FigS02B_Tcell_marker_colorbar.pdf"),width = 4,height = 4)
FeaturePlot(Tcell.obj,features = c("CD4"),
            cols=c('lightgrey',"red"),raster = T,pt.size=3)
dev.off()

####################
# FigS2C: 和https://spica.unil.ch/refs的CD8T细胞进行比较
CD8T.ref = readRDS("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/CD8T_human_ref_v1.rds")
CD8T.subcelltype = c("t01_CD8_IL7R","t02_CD8_AKNA","t03_CD8_EOMES","t04_CD8_ZNF683","t05_CD8_HAVCR2","t06_CD8_ISG","t07_CD8_MKI67")
CD8T.obj = subset(Tcell.obj,idents=CD8T.subcelltype)
# 1. 獲取兩個對象的共有高變基因 (Variable Features)
common_genes <- intersect(
  Seurat::VariableFeatures(CD8T.obj), 
  Seurat::VariableFeatures(CD8T.ref)
)

# 2. 計算用戶數據的亞群平均表達量
# 使用 Idents() 作為分組依據
user_avg <- Seurat::AverageExpression(
  CD8T.obj, 
  features = common_genes, 
  group.by = "ident", 
  slot = "data"
)$RNA

# 3. 計算參考數據的亞群平均表達量
ref_avg <- Seurat::AverageExpression(
  CD8T.ref, 
  features = common_genes, 
  group.by = "ident", 
  slot = "data"
)$RNA
# 4. 計算 Spearman 相關性矩陣
# 行為用戶亞群，列為參考亞群
user_matrix <- as.matrix(user_avg[common_genes, ])
ref_matrix <- as.matrix(ref_avg[common_genes, ])
cor_mat <- stats::cor(user_matrix, ref_matrix, method = "spearman")

z_mat <- t(scale(t(cor_mat)))

# 定義色盲友好色系：深藍 (低) -> 白色 (中) -> 深紅 (高)
# 我們將顏色中心點設為中位數，而非 0
col_fun <- circlize::colorRamp2(
  c(min(z_mat), 0, max(z_mat)), 
  c("#4575B4", "#FFFFBF", "#D73027")
)

# 7. 繪製熱圖
hp <- ComplexHeatmap::Heatmap(
  z_mat,
  name = "z-score scaled\nsimilarity",
  col = col_fun,
  cluster_rows = TRUE,          # 對你的亞群進行聚類
  cluster_columns = TRUE,       # 對參考亞群進行聚類
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  rect_gp = grid::gpar(col = "white", lwd = 0.5), # 方塊白邊
  column_names_rot = 45,        # 底部標籤旋轉 45 度
  row_names_side = "left",      # 你的亞群名稱放左邊
  column_title = "Reference Subtypes (CD8T.ref)",
  row_title = "Your Clusters (CD8T.obj)",
  cell_fun = function(j, i, x, y, width, height, fill) {
    # 如果你想在格子裡顯示原始 P 值或相關係數，可以在這裡添加
  }
)

# 渲染熱圖
ComplexHeatmap::draw(hp)
pdf(file=file.path(output_dir, 'Figures',"FigS", 'FigS02C_Tcell_similarity_with_ref.pdf'),width=5,height=3)
ComplexHeatmap::draw(hp)
dev.off()

###########################
# FigS2D: T细胞亚群Ro/e
do.tissueDist.full(
    Tcell.obj@meta.data,
    meta.cluster = 'subcelltype',
    Group.info = 'Group',
    bar_colors = colors_15, # 传入自定义颜色
    pdf_width = 8,         # 增加宽度以容纳 18 个亚群
    fig.output = file.path(output_dir, 'Figures',"FigS", 'FigS02D_Tcell_subcelltype_Roe_with_bar.pdf')
)

###############
# FigS2E: 其余T细胞亚群比例
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
  
  # 確保因子順序
  df$Time_point <- factor(df$Time_point, levels = c("Pre", "Post"))
  df$Outcome <- factor(df$Outcome, levels = c("CR", "NR"))
  return(df)
}

# --- 函数 2：自動方向的置換檢驗 ---
auto_perm_p <- function(data, formula, grp_col) {
  grps <- levels(factor(data[[grp_col]]))
  # 計算兩組的中位數
  val1 <- data$percentage[data[[grp_col]] == grps[1]]
  val2 <- data$percentage[data[[grp_col]] == grps[2]]
  
  # 樣本量檢查，n<2 則不進行檢驗
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
    # 提取当前亚群数据
    sub_df <- prop_df[prop_df[[cluster_col]] == cl, ]
    
    # 定义计算逻辑，确保返回 numeric
    calc_p <- function(data_subset, formula, grp_col) {
      p_raw <- auto_perm_p(data_subset, formula, grp_col)
      return(as.numeric(p_raw)) # 核心修复：强制转为数值
    }
    
    # 构造单行数据框
    data.frame(
      cluster = as.character(cl),
      p_CR_Time = calc_p(dplyr::filter(sub_df, Outcome == "CR"), percentage ~ Time_point, "Time_point"),
      p_NR_Time = calc_p(dplyr::filter(sub_df, Outcome == "NR"), percentage ~ Time_point, "Time_point"),
      p_Pre_Out = calc_p(dplyr::filter(sub_df, Time_point == "Pre"), percentage ~ Outcome, "Outcome"),
      p_Post_Out = calc_p(dplyr::filter(sub_df, Time_point == "Post"), percentage ~ Outcome, "Outcome"),
      stringsAsFactors = FALSE
    )
  })
  
  # 使用 bind_rows 合并
  final_stats <- dplyr::bind_rows(results_list)
  return(final_stats)
}

plot_clusters_with_p <- function(prop_df, stats_df, cluster_col, title_prefix = "") {
  cb_colors <- c("Pre" = "#E69F00", "Post" = "#56B4E9") # 色盲友好色
  
  # 計算每個子圖的最大 Y 值以確定標注高度
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

  # 注意：批量標注 P 值在分面圖中較複雜，建議手動對顯著的類型單獨繪圖或使用此循環
  # 這裡展示標注 CR 組 Pre vs Post 的邏輯示例
  p <- p + ggsignif::geom_signif(
    data = prop_df,
    ggplot2::aes(xmin = 0.8, xmax = 1.2, 
                annotations = sprintf("p=%.3f", stats_df$p_CR_Time[match(get(cluster_col), stats_df$cluster)])),
    manual = TRUE, tip_length = 0.01, vjust = -0.2, size = 0.3, textsize = 2.5
  )
  
  return(p)
}

# 1. 提取數據 (假設你想用目前的 Idents)
prop_t <- get_prop_df(Tcell.obj)
# 2. 計算統計
stats_t <- run_full_stats(prop_t, "current_clusters")
# 3. 繪圖
plot_t <- plot_clusters_with_p(prop_t, stats_t, "current_clusters", "All T Cell Subsets")
print(plot_t)
ggsave(filename ="./Figures/FigS/FigS02E_tcell_proportion_histogram.pdf",height = 8,width = 8,plot = plot_t,dpi = 600,bg="white")

# 1. 提取數據
prop_cd8 <- get_prop_df(CD8T.obj)
# 2. 計算統計
stats_cd8 <- run_full_stats(prop_cd8, "current_clusters")
# 3. 繪圖
plot_cd8 <- plot_clusters_with_p(prop_cd8, stats_cd8, "current_clusters", "CD8+ T Cell Sub-clusters")
print(plot_cd8)
ggsave(filename ="./Figures/FigS/FigS02E_CD8T_proportion_histogram.pdf",height = 8,width = 8,plot = plot_cd8,dpi = 600,bg="white")

###############
# FigS2F: t03_CD8_EOMES和t04_CD8_ZNF683分组基因集评分比较
Tcm = subset(Tcell.obj,subset=subcelltype=="t03_CD8_EOMES")
Trm = subset(Tcell.obj,subset=subcelltype=="t04_CD8_ZNF683")

p_tcm_eff <- plot_exhaustion_violin_full(Tcm,score_col = "Activation.Effector.function_score")
ggplot2::ggsave(file.path(output_dir,"Figures","FigS","FigS02F_Tcm_Effector_score.pdf"), p_tcm_eff, width = 3, height = 4)
p_tcm_cyto <- plot_exhaustion_violin_full(Tcm,score_col = "Cytotoxicity_score")
ggplot2::ggsave(file.path(output_dir,"Figures","FigS","FigS02F_Tcm_Cytotoxicity_score.pdf"), p_tcm_cyto, width = 3, height = 4)
p_tcm_exhau <- plot_exhaustion_violin_full(Tcm,score_col = "Exhaustion_score")
ggplot2::ggsave(file.path(output_dir,"Figures","FigS","FigS02F_Tcm_Exhaustion_score.pdf"), p_tcm_exhau, width = 3, height = 4)

p_trm_eff <- plot_exhaustion_violin_full(Trm,score_col = "Activation.Effector.function_score")
ggplot2::ggsave(file.path(output_dir,"Figures","FigS","FigS02F_Trm_Effector_score.pdf"), p_trm_eff, width = 3, height = 4)
p_trm_cyto <- plot_exhaustion_violin_full(Trm,score_col = "Cytotoxicity_score")
ggplot2::ggsave(file.path(output_dir,"Figures","FigS","FigS02F_Trm_Cytotoxicity_score.pdf"), p_trm_cyto, width = 3, height = 4)
p_trm_exhau <- plot_exhaustion_violin_full(Trm,score_col = "Exhaustion_score")
ggplot2::ggsave(file.path(output_dir,"Figures","FigS","FigS02F_Trm_Exhaustion_score.pdf"), p_trm_exhau, width = 3, height = 4)

}
