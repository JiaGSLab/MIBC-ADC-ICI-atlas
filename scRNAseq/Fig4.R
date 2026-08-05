###############################################################################
# Figure 4 and Supplementary Figures 4-6 — ordered review draft
#
# This file only reorganizes the user's original plotting code by the panel
# letters in the current manuscript:
#   /home/data/chenwenhui/codex project/manuscript 20260803/
#
# Original sources:
#   scRNAseq.R, Figure4.R, FigS3.R, FigS4.R, plot.R, kegg_plot.R
#
# IMPORTANT
# 1. This is a review draft. No panel runs by default.
# 2. Change panels_to_run only after the static comparison is confirmed.
# 3. No integration, PCA, UMAP, clustering, CellChat, or enrichment query is
#    recalculated here. Existing QC results are read directly.
# 4. Outcome-independent spatial maps (Fig4g/h and FigS5a/b) do not need
#    recomputation. Outcome is used there only in labels/file names.
###############################################################################

base::options(stringsAsFactors = FALSE)
base::set.seed(123)


panels_to_run <- c(
  "Fig4a", "Fig4b", "Fig4c", "Fig4i",
  "FigS4c", "FigS4d", "FigS4e"
)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(qs2)
  library(data.table)
  library(ggpubr)
  library(scRNAtoolVis)
})

qc_dir <- "/home/data/chenwenhui/scRNAseq/GPU_IT/QC"
script_dir <- "/home/data/chenwenhui/scRNAseq/GPU_IT/scripts"
output_dir <- qc_dir

malignant_path <- base::file.path(qc_dir, "malignant_obj_intergration_new.qs")
malignant_figs3_path <- base::file.path(qc_dir, "malignant_obj_intergration.qs")
spatial_path <- "/home/data/chenwenhui/visiumHD/ST.combined.obj.qs"
mp_signature_path <- base::file.path(qc_dir, "Final_MP_Top30_Signatures.qs")
robust_program_path <- base::file.path(qc_dir, "Robust_Programs_Ordered_with_ID.csv")
enrichment_path <- base::file.path(qc_dir, "All_MP_enrichR_Results.qs")
sample_info_path <- "/home/data/chenwenhui/sample_info.csv"
fig3_mp_mean_path <- base::file.path(qc_dir, "malignant_mp_mean_matrix.csv")
fig3_proportion_path <- base::file.path(qc_dir, "subcelltype_proportion_matrix.csv")

output_root <- base::file.path(
  qc_dir, "Figures", "Fig4_S4_S5_S6_ordered_review_20260803"
)
fig4_dir <- base::file.path(output_root, "Fig4")
figs4_dir <- base::file.path(output_root, "FigS4")
figs5_dir <- base::file.path(output_root, "FigS5")
figs6_dir <- base::file.path(output_root, "FigS6")
source_data_dir <- base::file.path(output_root, "source_data")
stats_dir <- base::file.path(output_root, "stats")

# Directories are created only when at least one reviewed panel is selected.
if (base::length(panels_to_run) > 0L) {
  for (one_dir in c(
    fig4_dir, figs4_dir, figs5_dir, figs6_dir, source_data_dir, stats_dir
  )) {
    base::dir.create(one_dir, recursive = TRUE, showWarnings = FALSE)
  }
}

font_family <- "Arial"
font_size <- 8
outcome_colors <- c("CR" = "#D55E00", "NR" = "#0072B2")
group_levels <- c("Pre_NR", "Pre_CR", "Post_NR", "Post_CR")
group_colors <- c(
  "Pre_NR" = "#A52A2A",
  "Pre_CR" = "#053E7A",
  "Post_NR" = "#E46B4F",
  "Post_CR" = "#8B8B00"
)
mp_names <- c(
  "Cycle_G2M", "Surface_Trafficking", "RNA_Splicing", "Basal_pEMT",
  "Stress_IEG", "IFN_Response", "Hypoxia", "Metabolic",
  "Uro_Differentiation", "Ep_Inflammation"
)

panel_review <- data.frame(
  panel = c(
    paste0("Fig4", letters[1:14]),
    paste0("FigS4", letters[1:5]),
    paste0("FigS5", letters[1:2]),
    paste0("FigS6", letters[1:6])
  ),
  source = c(
    "scRNAseq.R:1641-1910 + QC/Robust_Programs_Ordered_with_ID.csv",
    "scRNAseq.R:2082-2101,2169-2196",
    "scRNAseq.R:2329-2388",
    "plot.R:21-25 + kegg_plot.R",
    "external mIHC images",
    "Figure4.R:1-64",
    "Figure4.R:213-232 (Slide2 only in current Fig4)",
    "Figure4.R:234-269 (pt23 and pt42 selected in AI)",
    "Figure4.R:67-212",
    "Figure4.R:272-329",
    "Figure4.R:474-592 response block",
    "external experimental schematic",
    "Figure4.R:922-970",
    "Figure4.R:1230-1478",
    "scRNAseq.R:1463-1478 + FigS4.R:1-23",
    "FigS4.R:25-32",
    "scRNAseq.R:2367-2388",
    "FigS3.R:1-139; reads Fig3 result matrices in QC",
    "FigS3.R:142-218; reads malignant_obj_intergration.qs",
    "Figure4.R:213-232 (Slide1 and Slide3)",
    "Figure4.R:234-269",
    "Figure4.R:333-592",
    "Figure4.R:596-730",
    "survival source not found in the supplied R scripts",
    "survival source not found in the supplied R scripts",
    "Figure4.R:871-921",
    "Figure4.R:971-1229"
  ),
  stringsAsFactors = FALSE
)

if (base::length(panels_to_run) > 0L) {
  utils::write.csv(
    panel_review,
    base::file.path(source_data_dir, "Fig4_S4_S5_S6_panel_source.csv"),
    row.names = FALSE
  )
}


###############################################################################
# Small plotting/statistics helpers
###############################################################################

theme_article_8 <- function() {
  ggplot2::theme_classic(base_size = font_size, base_family = font_family) +
    ggplot2::theme(
      text = ggplot2::element_text(
        family = font_family, size = font_size, colour = "black"
      ),
      axis.title = ggplot2::element_text(
        family = font_family, size = font_size, colour = "black"
      ),
      axis.text = ggplot2::element_text(
        family = font_family, size = font_size, colour = "black"
      ),
      legend.title = ggplot2::element_text(
        family = font_family, size = font_size, colour = "black"
      ),
      legend.text = ggplot2::element_text(
        family = font_family, size = font_size, colour = "black"
      ),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        family = font_family, size = font_size, face = "plain"
      ),
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      plot.margin = ggplot2::margin(1.2, 1.2, 1.2, 1.2, unit = "mm")
    )
}

save_panel <- function(plot_object, file_stem, width_mm, height_mm, dpi = 300) {
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
    dpi = dpi,
    bg = "white",
    limitsize = FALSE
  )
}

format_p <- function(p_value) {
  if (base::is.na(p_value)) return("p = NA")
  if (p_value < 0.001) {
    return(paste0("p = ", base::formatC(p_value, format = "e", digits = 1)))
  }
  paste0("p = ", base::trimws(base::formatC(p_value, format = "fg", digits = 2)))
}

format_cor_label <- function(r_value, p_value) {
  p_text <- if (p_value < 2.2e-16) {
    "p < 2.2e-16"
  } else if (p_value < 0.001) {
    paste0("p = ", base::formatC(p_value, format = "e", digits = 1))
  } else {
    paste0("p = ", base::trimws(base::formatC(p_value, format = "fg", digits = 2)))
  }
  paste0("R = ", base::formatC(r_value, format = "f", digits = 2), ", ", p_text)
}

read_qs_object <- function(path) {
  if (!base::file.exists(path)) base::stop("Missing qs file: ", path)
  qs2::qs_read(path)
}



add_named_mp_scores <- function(object) {
  # malignant_obj_intergration_new.qs 是 Fig4 的正确恶性细胞对象。
  # 该对象保存于重新整合完成、MP 打分之前；这里严格复用 scRNAseq.R
  # 的 Final_MP_Top30_Signatures.qs、参数与 MP1-MP10 顺序计算 MP_Final。
  columns_to_remove <- base::grep(
    paste0(
      "^(MP_Final_[0-9]+|",
      base::paste(mp_names, collapse = "|"),
      ")$"
    ),
    base::colnames(object@meta.data),
    value = TRUE
  )
  if (base::length(columns_to_remove) > 0L) {
    object@meta.data <- object@meta.data[
      , base::setdiff(base::colnames(object@meta.data), columns_to_remove),
      drop = FALSE
    ]
  }

  mp_gene_lists <- read_qs_object(mp_signature_path)
  if (base::length(mp_gene_lists) != base::length(mp_names)) {
    base::stop("Final_MP_Top30_Signatures.qs must contain 10 signatures.")
  }

  object <- Seurat::AddModuleScore(
    object,
    features = mp_gene_lists,
    name = "MP_Final_",
    assay = "RNA",
    nbin = 24,
    seed = 123
  )

  for (i in base::seq_along(mp_names)) {
    old_name <- paste0("MP_Final_", i)
    if (!old_name %in% base::colnames(object@meta.data)) {
      base::stop("Missing recomputed MP column: ", old_name)
    }
    base::colnames(object@meta.data)[
      base::colnames(object@meta.data) == old_name
    ] <- mp_names[i]
  }
  object
}

calculate_pairwise_mp_stats <- function(patient_long) {
  comparisons <- list(
    c("Pre_CR", "Pre_NR"),
    c("Pre_NR", "Post_NR")
  )
  result <- list()
  result_index <- 1L
  for (one_mp in base::unique(patient_long$MP)) {
    one_df <- patient_long[patient_long$MP == one_mp, , drop = FALSE]
    score_range <- base::diff(base::range(one_df$Average_Score, na.rm = TRUE))
    if (!base::is.finite(score_range) || score_range == 0) score_range <- 1
    score_max <- base::max(one_df$Average_Score, na.rm = TRUE)
    for (comparison_index in base::seq_along(comparisons)) {
      groups <- comparisons[[comparison_index]]
      x_value <- one_df$Average_Score[one_df$Group == groups[1]]
      y_value <- one_df$Average_Score[one_df$Group == groups[2]]
      test_result <- stats::wilcox.test(
        x = x_value,
        y = y_value,
        alternative = "two.sided"
      )
      result[[result_index]] <- data.frame(
        MP = one_mp,
        group1 = groups[1],
        group2 = groups[2],
        n_group1 = base::length(x_value),
        n_group2 = base::length(y_value),
        method = "two-sided Wilcoxon rank-sum test (stats::wilcox.test defaults)",
        statistic = base::unname(test_result$statistic),
        p_value = test_result$p.value,
        p_label = format_p(test_result$p.value),
        y.position = score_max + score_range * (0.10 + 0.17 * (comparison_index - 1)),
        stringsAsFactors = FALSE
      )
      result_index <- result_index + 1L
    }
  }
  dplyr::bind_rows(result)
}

calculate_correlation_stats <- function(data, x_column, y_columns) {
  result <- lapply(y_columns, function(y_column) {
    complete <- stats::complete.cases(data[, c(x_column, y_column), drop = FALSE])
    x_value <- data[[x_column]][complete]
    y_value <- data[[y_column]][complete]
    cor_result <- stats::cor.test(x_value, y_value, method = "pearson")
    lm_result <- stats::lm(y_value ~ x_value)
    data.frame(
      x = x_column,
      y = y_column,
      n = base::length(x_value),
      pearson_r = base::unname(cor_result$estimate),
      p_value = cor_result$p.value,
      slope = base::unname(stats::coef(lm_result)[2]),
      intercept = base::unname(stats::coef(lm_result)[1]),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(result)
}

make_correlation_plot <- function(
    data, x_column, y_column, x_label, y_label, stats_row,
    show_legend = FALSE) {
  correlation_label <- format_cor_label(
    stats_row$pearson_r[[1]],
    stats_row$p_value[[1]]
  )
  plot_object <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = .data[[x_column]], y = .data[[y_column]])
  ) +
    ggrastr::geom_point_rast(
      ggplot2::aes(colour = .data$Group),
      size = 4,
      alpha = 0.8,
      raster.dpi = 600
    ) +
    ggplot2::geom_smooth(
      ggplot2::aes(group = 1),
      method = "lm",
      formula = y ~ x,
      colour = "black",
      fill = "grey80",
      linewidth = 0.4,
      linetype = "dashed",
      alpha = 0.2,
      se = TRUE
    ) +
    ggplot2::annotate(
      geom = "text",
      x = -Inf,
      y = Inf,
      label = correlation_label,
      hjust = -0.05,
      vjust = 1.15,
      family = font_family,
      size = font_size / ggplot2::.pt
    ) +
    ggplot2::scale_colour_manual(values = group_colors, drop = FALSE) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.15))
    ) +
    ggplot2::labs(x = x_label, y = y_label, colour = NULL) +
    theme_article_8()

  if (!show_legend) {
    plot_object <- plot_object + ggplot2::theme(legend.position = "none")
  } else {
    plot_object <- plot_object + ggplot2::theme(
      legend.position = "right",
      legend.key.height = grid::unit(3.5, "mm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0)
    )
  }
  plot_object
}


###############################################################################
# Fig4a — consensus MP Jaccard heatmap
# Original source: scRNAseq.R:1641-1910; reconstructed from the saved ordered robust-program CSV
###############################################################################
if ("Fig4a" %in% panels_to_run) {
  # Fig4a: consensus meta-program Jaccard heatmap
  ###############################################################################
  
  message("Fig4a: rebuilding the corrected MP annotation heatmap")
  
  robust_raw <- utils::read.csv(
    robust_program_path,
    check.names = FALSE
  )
  # 该文件的第一列是原 write.csv 产生的空行名列；第一条数据行是程序 ID。
  if (base::as.character(robust_raw[[1]][1]) != "ID") {
    base::stop("The first data row of Robust_Programs_Ordered_with_ID.csv must be ID.")
  }
  robust_with_id <- robust_raw[, -1, drop = FALSE]
  robust_programs <- robust_with_id[-1, , drop = FALSE]
  
  calculate_jaccard_matrix <- function(program_matrix) {
    program_count <- base::ncol(program_matrix)
    output <- base::matrix(0, nrow = program_count, ncol = program_count)
    base::colnames(output) <- base::rownames(output) <- base::colnames(program_matrix)
    for (i in base::seq_len(program_count)) {
      genes_i <- base::unique(program_matrix[, i])
      for (j in i:program_count) {
        genes_j <- base::unique(program_matrix[, j])
        similarity <- base::length(base::intersect(genes_i, genes_j)) /
          base::length(base::union(genes_i, genes_j))
        output[i, j] <- similarity
        output[j, i] <- similarity
      }
    }
    output
  }
  
  jaccard_matrix <- calculate_jaccard_matrix(robust_programs)
  
  sample_info <- utils::read.csv(sample_info_path, check.names = FALSE)
  sample_info <- sample_info |>
    dplyr::mutate(
      Outcome_short = dplyr::case_when(
        .data$Outcome == "response" ~ "CR",
        .data$Outcome == "non-response" ~ "NR",
        TRUE ~ base::as.character(.data$Outcome)
      ),
      Time_short = dplyr::case_when(
        .data$Time_point == "pre-treatment" ~ "Pre",
        .data$Time_point == "post-treatment" ~ "Post",
        TRUE ~ base::as.character(.data$Time_point)
      ),
      Group = base::paste(.data$Time_short, .data$Outcome_short, sep = "_")
    )
  
  program_annotation <- data.frame(
    Program = base::colnames(robust_programs),
    Patient = base::sub("_R[0-9].*", "", base::colnames(robust_programs)),
    stringsAsFactors = FALSE
  )
  program_annotation$Group <- sample_info$Group[
    base::match(program_annotation$Patient, sample_info$orig.ident)
  ]
  if (base::anyNA(program_annotation$Group)) {
    missing_patients <- base::unique(program_annotation$Patient[base::is.na(program_annotation$Group)])
    base::stop("Missing clinical Group for: ", base::paste(missing_patients, collapse = ", "))
  }
  if (!base::all(program_annotation$Group[program_annotation$Patient == "LSY-1_1114"] == "Pre_NR")) {
    base::stop("Fig4a annotation still labels LSY-1_1114 incorrectly.")
  }
  
  mp_ranges <- list(
    "Cycle_G2M" = 1:23,
    "Surface_Trafficking" = 24:34,
    "RNA_Splicing" = 37:64,
    "Basal_pEMT" = 69:92,
    "Stress_IEG" = 96:101,
    "IFN_Response" = 111:120,
    "Hypoxia" = 121:134,
    "Metabolic" = 135:155,
    "Uro_Differentiation" = 169:187,
    "Ep_Inflammation" = 188:194
  )
  program_annotation$MP <- "None"
  for (one_mp in base::names(mp_ranges)) {
    program_annotation$MP[mp_ranges[[one_mp]]] <- one_mp
  }
  program_annotation$MP <- base::factor(
    program_annotation$MP,
    levels = c(mp_names, "None")
  )
  program_annotation$Group <- base::factor(
    program_annotation$Group,
    levels = group_levels
  )
  
  palette_environment <- base::new.env(parent = base::emptyenv())
  base::load("/home/data/chenwenhui/scRNAseq/colors_90.RData", envir = palette_environment)
  patient_palette <- palette_environment$colors_91
  unique_patients <- base::unique(program_annotation$Patient)
  patient_colors <- stats::setNames(
    patient_palette[base::seq_along(unique_patients)],
    unique_patients
  )
  mp_colors <- stats::setNames(
    c(RColorBrewer::brewer.pal(10, "Paired"), "white"),
    c(mp_names, "None")
  )
  
  top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    MP = program_annotation$MP,
    Group = program_annotation$Group,
    Patient = program_annotation$Patient,
    col = list(MP = mp_colors, Group = group_colors, Patient = patient_colors),
    show_legend = c(MP = TRUE, Group = TRUE, Patient = FALSE),
    show_annotation_name = TRUE,
    annotation_name_side = "left",
    annotation_name_rot = 0,
    annotation_name_offset = grid::unit(0.8, "mm"),
    annotation_name_gp = grid::gpar(fontfamily = font_family, fontsize = font_size),
    annotation_legend_param = list(
      MP = list(
        title_gp = grid::gpar(fontfamily = font_family, fontsize = font_size),
        labels_gp = grid::gpar(fontfamily = font_family, fontsize = font_size),
        grid_width = grid::unit(3, "mm"),
        grid_height = grid::unit(3, "mm"),
        gap = grid::unit(0.8, "mm"),
        ncol = 1
      ),
      Group = list(
        title_gp = grid::gpar(fontfamily = font_family, fontsize = font_size),
        labels_gp = grid::gpar(fontfamily = font_family, fontsize = font_size),
        grid_width = grid::unit(3, "mm"),
        grid_height = grid::unit(3, "mm")
      )
    ),
    border = FALSE,
    simple_anno_size = grid::unit(4, "mm"),
    gap = grid::unit(0.7, "mm")
  )
  
  fig4a_heatmap <- ComplexHeatmap::Heatmap(
    jaccard_matrix,
    name = "Jaccard Index",
    column_title = "Consensus Meta-Programs (MP1-MP10)",
    column_title_gp = grid::gpar(
      fontfamily = font_family, fontsize = font_size, fontface = "plain"
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_title = NULL,
    top_annotation = top_annotation,
    column_split = program_annotation$MP,
    row_split = program_annotation$MP,
    cluster_column_slices = FALSE,
    cluster_row_slices = FALSE,
    column_gap = grid::unit(0.8, "mm"),
    row_gap = grid::unit(0.8, "mm"),
    col = circlize::colorRamp2(
      c(0.05, 0.10, 0.30),
      c("white", "#fddbc7", "#b2182b")
    ),
    # PDF 中仅将 200 x 200 热图主体栅格化；标题、注释和图例仍为矢量。
    use_raster = TRUE,
    raster_device = "png",
    raster_quality = 4,
    raster_by_magick = FALSE,
    rect_gp = grid::gpar(col = NA),
    heatmap_legend_param = list(
      title_gp = grid::gpar(fontfamily = font_family, fontsize = font_size),
      labels_gp = grid::gpar(fontfamily = font_family, fontsize = font_size),
      legend_height = grid::unit(28, "mm"),
      grid_width = grid::unit(3, "mm")
    )
  )
  
  grDevices::cairo_pdf(
    filename = base::file.path(fig4_dir, "Fig4a_consensus_MP_Jaccard_heatmap.pdf"),
    width = 135 / 25.4,
    height = 84 / 25.4,
    family = font_family
  )
  base::suppressWarnings(
    ComplexHeatmap::draw(
      fig4a_heatmap,
      merge_legend = TRUE,
      annotation_legend_side = "right",
      heatmap_legend_side = "right"
    )
  )
  grDevices::dev.off()
  
  ragg::agg_png(
    filename = base::file.path(fig4_dir, "Fig4a_consensus_MP_Jaccard_heatmap.png"),
    width = 135,
    height = 84,
    units = "mm",
    res = 300,
    background = "white"
  )
  base::suppressWarnings(
    ComplexHeatmap::draw(
      fig4a_heatmap,
      merge_legend = TRUE,
      annotation_legend_side = "right",
      heatmap_legend_side = "right"
    )
  )
  grDevices::dev.off()
  
  jaccard_long <- base::as.data.frame(base::as.table(jaccard_matrix))
  base::colnames(jaccard_long) <- c("Program_row", "Program_column", "Jaccard_index")
  utils::write.csv(
    jaccard_long,
    base::file.path(source_data_dir, "Fig4a_Jaccard_matrix_long.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    program_annotation,
    base::file.path(source_data_dir, "Fig4a_program_annotation.csv"),
    row.names = FALSE
  )
  
  base::rm(
    fig4a_heatmap, top_annotation, robust_raw, robust_with_id, robust_programs,
    jaccard_matrix, jaccard_long
  )
  base::invisible(base::gc())
  
}

###############################################################################
# Shared malignant-cell preparation for Fig4b/c/j and FigS4c
# Original method retained: AddModuleScore(nbin = 24, seed = 123).
###############################################################################

malignant_score_panels <- c("Fig4b", "Fig4c", "Fig4j", "FigS4c")
if (base::any(malignant_score_panels %in% panels_to_run)) {
  malignant_obj <- read_qs_object(malignant_path)
  malignant_obj <- add_named_mp_scores(malignant_obj)
  malignant_meta <- malignant_obj@meta.data
}

###############################################################################
# Fig4b — cell-level Basal_pEMT / Uro_Differentiation ridgelines
# Original source: scRNAseq.R:2169-2196
###############################################################################
if ("Fig4b" %in% panels_to_run) {
  ridge_data <- data.frame(
    cell_id = base::rownames(malignant_meta),
    sample = malignant_meta$sample,
    Group = malignant_meta$Group,
    Basal_pEMT = malignant_meta$Basal_pEMT,
    Uro_Differentiation = malignant_meta$Uro_Differentiation,
    check.names = FALSE
  ) |>
    tidyr::pivot_longer(
      cols = c("Basal_pEMT", "Uro_Differentiation"),
      names_to = "MP",
      values_to = "Score"
    )

  fig4b_plot <- ggplot2::ggplot(
    ridge_data,
    ggplot2::aes(x = .data$Score, y = .data$Group, fill = .data$Group)
  ) +
    ggridges::geom_density_ridges(alpha = 0.7, quantile_lines = TRUE) +
    ggplot2::scale_fill_manual(values = group_colors, drop = FALSE) +
    ggplot2::facet_wrap(~MP, scales = "free") +
    ggridges::theme_ridges(
      font_size = font_size,
      font_family = font_family,
      grid = FALSE
    ) +
    ggplot2::labs(x = "Score", y = NULL) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        family = font_family,
        size = font_size,
        face = "plain"
      ),
      legend.position = "none"
    )

  save_panel(
    fig4b_plot,
    base::file.path(fig4_dir, "Fig4b_MP4_MP9_ridgelines"),
    width_mm = 72,
    height_mm = 35
  )
  utils::write.csv(
    ridge_data,
    base::file.path(source_data_dir, "Fig4b_cell_level_MP_scores.csv"),
    row.names = FALSE
  )
}

###############################################################################
# Fig4c and FigS4c — patient-level MP scores and original Wilcoxon tests
# Original source: scRNAseq.R:2329-2388
###############################################################################
if (base::any(c("Fig4c", "FigS4c") %in% panels_to_run)) {
  patient_average <- malignant_meta |>
    dplyr::group_by(.data$sample, .data$Group) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(mp_names),
        function(x) base::mean(x, na.rm = TRUE)
      ),
      .groups = "drop"
    )
  patient_long <- patient_average |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(mp_names),
      names_to = "MP",
      values_to = "Average_Score"
    )
  mp_stats <- calculate_pairwise_mp_stats(patient_long)

  utils::write.csv(
    patient_long,
    base::file.path(source_data_dir, "Fig4c_FigS4c_patient_level_MP_scores.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    mp_stats,
    base::file.path(stats_dir, "Fig4c_FigS4c_patient_level_Wilcoxon.csv"),
    row.names = FALSE
  )

  original_comparisons <- list(
    c("Pre_CR", "Pre_NR"),
    c("Pre_NR", "Post_NR")
  )

  if ("Fig4c" %in% panels_to_run) {
    fig4c_data <- patient_long[
      patient_long$MP %in% c("Basal_pEMT", "Uro_Differentiation"),
      , drop = FALSE
    ]
    fig4c_plot <- ggpubr::ggboxplot(
      fig4c_data,
      x = "Group",
      y = "Average_Score",
      fill = "Group",
      palette = group_colors,
      add = "jitter",
      facet.by = "MP",
      scales = "free_y",
      ncol = 2
    ) +
      ggpubr::stat_compare_means(
        comparisons = original_comparisons,
        method = "wilcox.test",
        label = "p.format",
        size = font_size / ggplot2::.pt
      ) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0.05, 0.30))
      ) +
      theme_article_8() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "none"
      ) +
      ggplot2::labs(x = NULL, y = "Mean Module Score per Patient")

    save_panel(
      fig4c_plot,
      base::file.path(fig4_dir, "Fig4c_patient_level_Basal_Uro_boxplots"),
      width_mm = 72,
      height_mm = 48
    )
  }

  if ("FigS4c" %in% panels_to_run) {
    figs4c_order <- c(
      "Cycle_G2M", "Ep_Inflammation", "Metabolic", "RNA_Splicing",
      "Hypoxia", "IFN_Response", "Stress_IEG", "Surface_Trafficking"
    )
    figs4c_data <- patient_long[
      patient_long$MP %in% figs4c_order,
      , drop = FALSE
    ]
    figs4c_data$MP <- base::factor(figs4c_data$MP, levels = figs4c_order)
    figs4c_plot <- ggpubr::ggboxplot(
      figs4c_data,
      x = "Group",
      y = "Average_Score",
      fill = "Group",
      palette = group_colors,
      add = "jitter",
      facet.by = "MP",
      scales = "free_y",
      ncol = 2
    ) +
      ggpubr::stat_compare_means(
        comparisons = original_comparisons,
        method = "wilcox.test",
        label = "p.format",
        size = font_size / ggplot2::.pt
      ) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0.05, 0.30))
      ) +
      theme_article_8() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "none"
      ) +
      ggplot2::labs(x = NULL, y = "Mean Module Score per Patient")

    save_panel(
      figs4c_plot,
      base::file.path(figs4_dir, "FigS4c_patient_level_supplementary_MP_boxplots"),
      width_mm = 70,
      height_mm = 184
    )
  }
}


###############################################################################
# Fig4d — enrichment bar plot
# Original source: plot.R:21-25 and kegg_plot.R
# Existing All_MP_enrichR_Results.qs is used; no Enrichr query is rerun.
###############################################################################
if ("Fig4d" %in% panels_to_run) {
  all_mp_enrichr <- read_qs_object(enrichment_path)
  MP4.enrich = all_mp_enrichr[["MP4"]]
  MP9.enrich = all_mp_enrichr[["MP9"]]
  MP9.enrich.GO = MP9.enrich[['GO_Biological_Process_2025']]
  MP4.enrich.HallMark = MP4.enrich$MSigDB_Hallmark_2020
  MP4.enrich.GO = MP4.enrich$GO_Biological_Process_2025
  # 1. 确保已加载所需的库
  library(tidyverse)
  library(forcats)
  
  # 2. 读取您的数据文件 (假设数据结构与之前相同)
  kegg.enriched = rbind(MP4.enrich.HallMark[1:3,],MP4.enrich.GO[1:3,])
  kegg.enriched.negative = MP9.enrich.GO[1:2,]
  
  # 3. 数据预处理 (与之前保持一致)
  kegg_pos <- kegg.enriched %>%
    mutate(
      Regulation = "Positive",
      log_p = -log10(Adjusted.P.value),
      Plot_Value = log_p
    )
  
  kegg_neg <- kegg.enriched.negative %>%
    mutate(
      Regulation = "Negative",
      log_p = -log10(Adjusted.P.value),
      Plot_Value = -log_p
    )
  
  plot_data <- bind_rows(kegg_pos, kegg_neg)
  
  # 4. 排序条目 (与之前保持一致)
  pos_terms_ordered <- kegg_pos %>% arrange(desc(log_p)) %>% pull(Term)
  neg_terms_ordered <- kegg_neg %>% arrange(log_p) %>% pull(Term)
  ordered_terms <- c(pos_terms_ordered, neg_terms_ordered)
  plot_data$Term <- factor(plot_data$Term, levels = ordered_terms)
  
  # 5. 绘图：应用新的色盲友好颜色
  # 定义色盲友好颜色
  COLOR_POSITIVE <- "#D95F02" # 柔和的紫红色/品红色
  COLOR_NEGATIVE <- "#1B9E77" # 柔和的深蓝色
  
  # 定义标签偏移量
  x_label_offset <- max(abs(plot_data$Plot_Value)) * 0.1
  
  p <- ggplot(plot_data, aes(x = Plot_Value, y = Term, fill = Regulation)) +
    
    # 绘制条形图
    geom_bar(stat = "identity") +
    
    # **关键修改：使用色盲友好颜色**
    scale_fill_manual(values = c("Negative" = COLOR_NEGATIVE, "Positive" = COLOR_POSITIVE)) +
    
    # X 轴设置
    scale_x_continuous(
      labels = abs,
      name = expression(-log[10](Adjusted.P.value))
    ) +
    
    # 添加左侧（Positive）标签
    geom_text(
      data = plot_data %>% filter(Regulation == "Positive"),
      aes(
        x = -x_label_offset,
        label = Term
      ),
      hjust = 1,
      size = 3
    ) +
    
    # 添加右侧（Negative）标签
    geom_text(
      data = plot_data %>% filter(Regulation == "Negative"),
      aes(
        x = x_label_offset,
        label = Term
      ),
      hjust = 0,
      size = 3
    ) +
    
    # 移除 Y 轴的默认标签、刻度和标题
    labs(y = NULL, title = "KEGG Enrichment Diverging Bar Plot") +
    
    # 主题美化 (保持一致)
    # 主题美化
    theme_minimal() +
    theme(
      # 移除 Y 轴上的条目文字和刻度线
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      # 调整图例位置
      legend.position = "top",
      # **关键修改：移除所有网格线**
      panel.grid.major = element_blank(), # 移除主网格线
      panel.grid.minor = element_blank()  # 移除次网格线
      # 如果您想保留 Y 轴 (垂直) 的网格线，只移除 X 轴的网格线，则使用：
      # panel.grid.major.x = element_blank(),
      # panel.grid.minor.x = element_blank(),
    )
  
  # 显示图表
  print(p)
  ggsave(file.path(output_dir,"Figures","MP4_MP9_enrich.png"), plot = p, width = 10, height = 7, units = "in", dpi = 300)
  ggsave(file.path(output_dir,"Figures","MP4_MP9_enrich.pdf"), plot = p, width = 8, height = 6, units = "in", dpi = 300)
}

###############################################################################
# Fig4e — representative mIHC images
###############################################################################

###############################################################################
# Fig4f — KRT5+ / UPK2+ mIHC cell proportions
# Original source: Figure4.R:1-64
###############################################################################
if ("Fig4f" %in% panels_to_run) {
  ########################
  #Fig4F KRT5+ UPK2+阳性细胞比例统计
  df = read.csv('/home/data/chenwenhui/scRNAseq/GPU_IT/QC/Figures/Fig/Fig4F.csv')
  # 1. 计算比例并选择需要的列
  plot_df <- df %>%
    dplyr::mutate(
      KRT5_Proportion = (KRT5p / total.cell) * 100,
      UPK2_Proportion = (UPK2p / total.cell) * 100
    ) %>%
    dplyr::select(Group, KRT5_Proportion, UPK2_Proportion)
  
  # 2. 转换为长格式以便分面绘图
  plot_df_long <- plot_df %>%
    tidyr::pivot_longer(
      cols = c("KRT5_Proportion", "UPK2_Proportion"),
      names_to = "CellType",
      values_to = "Percentage"
    )
  
  # 3. 确保 Group 是因子，并设置 CR 在前
  plot_df_long$Group <- factor(plot_df_long$Group, levels = c("CR", "NR"))
  # 定义色盲友好配色 (CR: 朱红, NR: 天蓝)
  cb_palette <- c("CR" = "#E41A1C", "NR" = "#377EB8")
  
  p <- ggpubr::ggbarplot(
    plot_df_long, 
    x = "Group", 
    y = "Percentage",
    fill = "Group", 
    palette = cb_palette,
    add = "mean_se",          # 添加均值和标准误
    add.params = list(width = 0.5),
    position = position_dodge(0.8),
    facet.by = "CellType",     # 按细胞类型分面
    scales = "free_y",
    error.plot = "upper_errorbar"
  ) +
    # 叠加原始数据点，展示 FOV 间的异质性
    ggplot2::geom_jitter(
      width = 0.2, 
      size = 1.5,
      shape = 16, 
      color = "black"
    ) +
    # 添加统计检验 (Wilcoxon 秩和检验)
    ggpubr::stat_compare_means(
      method = "wilcox.test", 
      label = "p.format", 
      label.x = 1.5
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = "Comparison of Cell Proportions within FOVs",
      y = "Positive Cell Proportion (%)",
      x = ""
    ) +
    ggplot2::theme(
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "none" # 因为 X 轴已经标注了组名
    )
  
  print(p)
  ggsave(filename = file.path(output_dir,"Figures","Fig","Fig4F_mIHC_cell_proportion.pdf"),width=3,height = 2.5,plot=p)
}

###############################################################################
# Shared Visium HD preparation for DISPLAY panels only:
# Fig4g, Fig4h, FigS5a and FigS5b.
#
###############################################################################
spatial_map_panels <- c("Fig4g", "Fig4h", "FigS5a", "FigS5b")
if (base::any(spatial_map_panels %in% panels_to_run)) {
  HD.combined.obj <- read_qs_object(spatial_path)
  HD.combined.obj$Outcome <- base::ifelse(
    HD.combined.obj$sample == "ZZJ_8916",
    "CR",
    HD.combined.obj$Outcome
  )
  HD.combined.obj <- base::subset(HD.combined.obj, subset = Therapy == "IT")
  HD.combined.obj <- Seurat::NormalizeData(HD.combined.obj)
  HD.combined.obj$pct_counts_mt <- Seurat::PercentageFeatureSet(
    HD.combined.obj,
    pattern = "^MT-"
  )
  HD.combined.obj <- base::subset(
    HD.combined.obj,
    subset = nCount_Spatial > 50 & nCount_Spatial < 4000 &
      nFeature_Spatial > 50 & nFeature_Spatial < 2000 &
      pct_counts_mt < 25
  )

  MP.genelist <- read_qs_object(mp_signature_path)
  HD.combined.obj <- Seurat::AddModuleScore(
    HD.combined.obj,
    features = MP.genelist,
    name = "MPscore"
  )
  HD.combined.obj@meta.data <- HD.combined.obj@meta.data[
    , base::setdiff(base::colnames(HD.combined.obj@meta.data), mp_names),
    drop = FALSE
  ]
  for (i in base::seq_along(mp_names)) {
    base::colnames(HD.combined.obj@meta.data)[
      base::colnames(HD.combined.obj@meta.data) == paste0("MPscore", i)
    ] <- mp_names[i]
  }

  epi.clusters <- c(
    "e01_Epi_PSCA", "e02_Epi_MKI67", "e03_Epi_KRT5", "e04_Epi_NDRG1",
    "e05_Epi_FOS", "e06_Epi_S100A9", "e07_Epi_KRT13", "e08_Epi_CXCL10",
    "e09_Epi_PIGR", "e10_Epi_MMP7", "e11_Epi_KRT16"
  )
}

###############################################################################
# Fig4g and FigS5a — RCTD cell-type maps
# Original source: Figure4.R:213-232
###############################################################################
if (base::any(c("Fig4g", "FigS5a") %in% panels_to_run)) {
  RCTD.obj <- base::subset(
    HD.combined.obj,
    subset = first_type == "unknown",
    invert = TRUE
  )
  RCTD.obj$celltype <- base::ifelse(
    RCTD.obj$first_type %in% epi.clusters,
    "Epithelial",
    RCTD.obj$first_type
  )
  RCTD.obj$celltype <- base::factor(
    RCTD.obj$celltype,
    levels = c(
      "Epithelial", "Endothelial", "Pericyte", "SMC", "Fibroblast",
      "Macrophage", "Monocyte", "Mast cell", "DC", "Neutrophil", "B",
      "Plasma", "CD8 T", "CD4 T", "Treg", "γδT"
    )
  )

  palette_environment <- base::new.env(parent = base::emptyenv())
  base::load(
    "/home/data/chenwenhui/scRNAseq/colors_90.RData",
    envir = palette_environment
  )
  spatial_celltype_colors <- palette_environment$colors_91

  slides_to_draw <- c(
    if ("Fig4g" %in% panels_to_run) "slide02",
    if ("FigS5a" %in% panels_to_run) c("slide01", "slide03")
  )
  for (slide in slides_to_draw) {
    RCTD.obj.slide <- base::subset(RCTD.obj, subset = orig.ident == slide)
    spatial_celltype_plot <- Seurat::SpatialPlot(
      RCTD.obj.slide,
      group.by = "celltype",
      pt.size.factor = 4,
      image.alpha = 0.2
    ) +
      ggplot2::scale_fill_manual(values = spatial_celltype_colors) +
      ggplot2::theme(legend.position = "none", aspect.ratio = 1)

    ggplot2::ggsave(
      base::file.path(
        if (slide == "slide02") fig4_dir else figs5_dir,
        paste0(if (slide == "slide02") "Fig4g_" else "FigS5a_", slide, "_celltype.png")
      ),
      plot = spatial_celltype_plot,
      device = ragg::agg_png,
      width = 4,
      height = 4,
      units = "in",
      dpi = 600,
      bg = "transparent"
    )
  }
}

###############################################################################
# Fig4h and FigS5b — spatial MP maps
# Original source: Figure4.R:234-269
###############################################################################
if (base::any(c("Fig4h", "FigS5b") %in% panels_to_run)) {
  map_samples <- if ("FigS5b" %in% panels_to_run) {
    base::unique(base::as.character(HD.combined.obj$sample_id))
  } else {
    c("pt23", "pt42")
  }

  for (samp_id in map_samples) {
    sub_obj <- base::subset(HD.combined.obj, subset = sample_id == samp_id)
    group_info <- base::unique(base::as.character(sub_obj$Outcome))
    for (feature in c("Basal_pEMT", "Uro_Differentiation")) {
      score_limits <- if (feature == "Basal_pEMT") c(-0.1, 1.5) else c(-0.1, 0.2)
      spatial_feature_plot <- Seurat::SpatialFeaturePlot(
        sub_obj,
        features = feature,
        pt.size.factor = 5,
        image.alpha = 0.2
      ) +
        ggplot2::scale_fill_viridis_c(option = "magma", limits = score_limits) +
        ggplot2::theme(legend.position = "none")

      ggplot2::ggsave(
        base::file.path(
          if (samp_id %in% c("pt23", "pt42")) fig4_dir else figs5_dir,
          paste0(
            if (samp_id %in% c("pt23", "pt42")) "Fig4h_" else "FigS5b_",
            samp_id, "_", group_info, "_", feature, ".png"
          )
        ),
        plot = spatial_feature_plot,
        device = ragg::agg_png,
        width = 2,
        height = 2,
        units = "in",
        dpi = 600,
        bg = "transparent"
      )
    }
  }
}

###############################################################################
# Fig4i — sample-level spatial MP comparison
# Original source: Figure4.R:67-212
#
# Scope correction found during static comparison:
#   the original line 82 slide02/slide03 filter was upstream of this test and
#   reduced 11 samples to 8. Here the QC object keeps all three slides.
#
# REVIEW before running:
#   the original script calculates Basal with "less" and Uro with "greater",
#   but its plotted stat_compare_means then applies "greater" to both facets.
#   The vector below follows the two separate calculations/comments and is
#   intentionally visible for confirmation.
###############################################################################
if ("Fig4i" %in% panels_to_run) {
  HD.stats <- read_qs_object(spatial_path)
  HD.stats$Outcome <- base::ifelse(
    HD.stats$sample == "ZZJ_8916",
    "CR",
    HD.stats$Outcome
  )
  HD.stats <- base::subset(HD.stats, subset = Therapy == "IT")
  HD.stats <- Seurat::NormalizeData(HD.stats)
  HD.stats$pct_counts_mt <- Seurat::PercentageFeatureSet(
    HD.stats,
    pattern = "^MT-"
  )
  HD.stats <- base::subset(
    HD.stats,
    subset = nCount_Spatial > 50 & nCount_Spatial < 4000 &
      nFeature_Spatial > 50 & nFeature_Spatial < 2000 &
      pct_counts_mt < 25
  )

  epi.clusters <- c(
    "e01_Epi_PSCA", "e02_Epi_MKI67", "e03_Epi_KRT5", "e04_Epi_NDRG1",
    "e05_Epi_FOS", "e06_Epi_S100A9", "e07_Epi_KRT13", "e08_Epi_CXCL10",
    "e09_Epi_PIGR", "e10_Epi_MMP7", "e11_Epi_KRT16"
  )
  HD.epi <- base::subset(
    HD.stats,
    subset = first_type %in% epi.clusters | second_type %in% epi.clusters
  )
  MP.genelist <- read_qs_object(mp_signature_path)
  HD.epi <- Seurat::AddModuleScore(
    HD.epi,
    features = MP.genelist,
    name = "MPscore"
  )
  HD.epi@meta.data <- HD.epi@meta.data[
    , base::setdiff(base::colnames(HD.epi@meta.data), mp_names),
    drop = FALSE
  ]
  for (i in base::seq_along(mp_names)) {
    base::colnames(HD.epi@meta.data)[
      base::colnames(HD.epi@meta.data) == paste0("MPscore", i)
    ] <- mp_names[i]
  }

  hd_sample_avg <- HD.epi@meta.data |>
    dplyr::group_by(.data$sample_id, .data$orig.ident, .data$Outcome) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(mp_names), base::mean, na.rm = TRUE),
      epithelial_spots = dplyr::n(),
      .groups = "drop"
    )
  base::stopifnot(dplyr::n_distinct(hd_sample_avg$sample_id) == 11L)

  hd_long <- hd_sample_avg |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(mp_names),
      names_to = "MP",
      values_to = "Score"
    )
  hd_long$Outcome <- base::factor(
    base::as.character(hd_long$Outcome),
    levels = c("CR", "NR")
  )
  fig4i_data <- hd_long[
    hd_long$MP %in% c("Basal_pEMT", "Uro_Differentiation"),
    , drop = FALSE
  ]

  # 与 Fig1 的原始函数保持完全相同的单侧方向设置：
  # 第二组（NR）的中位数 >= 第一组（CR）时，coin 使用 "less"；
  # 第二组（NR）的中位数 <  第一组（CR）时，coin 使用 "greater"。
  # 这里不再按 Basal/Uro 的预期生物学方向硬编码。
  fig4i_stats <- lapply(
    c("Basal_pEMT", "Uro_Differentiation"),
    function(one_mp) {
    one_data <- fig4i_data[fig4i_data$MP == one_mp, , drop = FALSE]
    median_CR <- stats::median(
      one_data$Score[one_data$Outcome == "CR"],
      na.rm = TRUE
    )
    median_NR <- stats::median(
      one_data$Score[one_data$Outcome == "NR"],
      na.rm = TRUE
    )
    alternative_one_sided <- if (median_NR >= median_CR) {
      "less"
    } else {
      "greater"
    }

    one_sided_test <- coin::independence_test(
      Score ~ Outcome,
      data = one_data,
      distribution = "exact",
      alternative = alternative_one_sided
    )
    two_sided_test <- coin::independence_test(
      Score ~ Outcome,
      data = one_data,
      distribution = "exact",
      alternative = "two.sided"
    )

    data.frame(
      MP = one_mp,
      first_level = "CR",
      second_level = "NR",
      n_CR = base::sum(one_data$Outcome == "CR"),
      n_NR = base::sum(one_data$Outcome == "NR"),
      mean_CR = base::mean(
        one_data$Score[one_data$Outcome == "CR"],
        na.rm = TRUE
      ),
      mean_NR = base::mean(
        one_data$Score[one_data$Outcome == "NR"],
        na.rm = TRUE
      ),
      median_CR = median_CR,
      median_NR = median_NR,
      median_difference_NR_minus_CR = median_NR - median_CR,
      observed_median_direction = if (median_NR >= median_CR) {
        "NR >= CR"
      } else {
        "NR < CR"
      },
      alternative_one_sided = alternative_one_sided,
      p_one_sided = base::as.numeric(coin::pvalue(one_sided_test)),
      p_two_sided = base::as.numeric(coin::pvalue(two_sided_test)),
      stringsAsFactors = FALSE
    )
  }) |>
    dplyr::bind_rows()

  utils::write.csv(
    hd_sample_avg,
    base::file.path(source_data_dir, "Fig4i_spatial_sample_audit.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    fig4i_data,
    base::file.path(source_data_dir, "Fig4i_spatial_MP_scores_by_sample.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    fig4i_stats,
    base::file.path(stats_dir, "Fig4i_spatial_exact_test.csv"),
    row.names = FALSE
  )

  # 图上使用与 Fig1 相同规则得到的单侧精确置换检验 P 值；
  # CSV 同时保留双侧 P 值，便于复核。
  fig4i_plot_stats <- fig4i_stats
  fig4i_plot_stats$label <- base::paste0(
    "p = ",
    base::formatC(fig4i_plot_stats$p_one_sided, format = "fg", digits = 2)
  )
  fig4i_plot <- ggplot2::ggplot(
    fig4i_data,
    ggplot2::aes(x = .data$Outcome, y = .data$Score, fill = .data$Outcome)
  ) +
    ggplot2::geom_boxplot(
      width = 0.62,
      alpha = 0.7,
      outlier.shape = NA,
      linewidth = 0.35
    ) +
    ggplot2::geom_jitter(
      width = 0.12,
      shape = 16,
      size = 1.15,
      colour = "black",
      alpha = 0.8
    ) +
    ggplot2::geom_text(
      data = fig4i_plot_stats,
      ggplot2::aes(x = 1.5, y = Inf, label = .data$label),
      inherit.aes = FALSE,
      vjust = 1.2,
      family = font_family,
      size = font_size / ggplot2::.pt
    ) +
    ggplot2::scale_fill_manual(values = outcome_colors) +
    ggplot2::facet_wrap(~MP, ncol = 1, scales = "free_y") +
    theme_article_8() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(x = NULL, y = "Mean score")

  save_panel(
    fig4i_plot,
    base::file.path(fig4_dir, "Fig4i_spatial_Basal_Uro_by_Outcome"),
    width_mm = 40,
    height_mm = 62
  )
}


###############################################################################
# Fig4j — malignant-cell MP scores versus ERBB2
# Original source: Figure4.R:272-329
###############################################################################
if ("Fig4j" %in% panels_to_run) {
  #######################################################
  # 1. 提取绘图数据 (保持不变)
  plot_df <- data.frame(
    MP4 = malignant.obj$Basal_pEMT,
    MP9 = malignant.obj$Uro_Differentiation,
    ERBB2 = Seurat::FetchData(malignant.obj, vars = "ERBB2")[,1]
  )
  
  # 2. 计算全局斜率 (Slope)
  fit4 <- stats::lm(ERBB2 ~ MP4, data = plot_df)
  slope4 <- stats::coef(fit4)[2]
  
  fit9 <- stats::lm(ERBB2 ~ MP9, data = plot_df)
  slope9 <- stats::coef(fit9)[2]
  p_mp4_global <- ggplot2::ggplot(plot_df, ggplot2::aes(x = MP4, y = ERBB2)) +
    # 栅格化密度层
    ggplot2::geom_bin2d(bins = 100) + 
    ggplot2::scale_fill_viridis_c(option = "magma", trans = "log10") + 
    # 拟合线：青色虚线
    ggplot2::geom_smooth(method = "lm", color = "cyan", linewidth = 1, linetype = "dashed") +
    # 标注 R 和 P 值 [cite: 2025-12-25]
    ggpubr::stat_cor(method = "pearson", color = "black", size = 4, 
                     label.x.npc = "left", label.y.npc = "top") + 
    # 标注斜率 (手动下移 Y 轴位置以防重叠)
    ggplot2::annotate("text", x = min(plot_df$MP4), y = max(plot_df$ERBB2) * 0.85, 
                      label = paste("Slope =", round(slope4, 3)), 
                      color = "black", hjust = 0, fontface = "italic", size = 4) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = "Global Correlation: MP4 vs ERBB2",
                  x = "MP4 (Basal_pEMT) Score", y = "ERBB2 Expression") +
    ggplot2::theme(panel.grid = element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"),aspect.ratio = 1)
  
  # 4. 绘制 MP9 vs ERBB2 (正相关)
  p_mp9_global <- ggplot2::ggplot(plot_df, ggplot2::aes(x = MP9, y = ERBB2)) +
    ggplot2::geom_bin2d(bins = 100) + 
    ggplot2::scale_fill_viridis_c(option = "magma", trans = "log10") + 
    ggplot2::geom_smooth(method = "lm", color = "cyan", linewidth = 1, linetype = "dashed") +
    ggpubr::stat_cor(method = "pearson", color = "black", size = 4, 
                     label.x.npc = "left", label.y.npc = "top") + 
    ggplot2::annotate("text", x = min(plot_df$MP9), y = max(plot_df$ERBB2) * 0.85, 
                      label = paste("Slope =", round(slope9, 3)), 
                      color = "black", hjust = 0, fontface = "italic", size = 4) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = "Global Correlation: MP9 vs ERBB2",
                  x = "MP9 (Uro_Differentiation) Score", y = "ERBB2 Expression") +
    ggplot2::theme(panel.grid = element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"),aspect.ratio = 1)
  
  # 组合展示
  combined_plot <- p_mp4_global / p_mp9_global
  ggplot2::ggsave(
    filename = file.path(output_dir, "Figures","Fig", "Fig4I_MP_vs_ERBB2_Global_Density.pdf"),
    plot = combined_plot, 
    width = 4, 
    height = 8,
    device = grDevices::cairo_pdf
  )
}

###############################################################################
# Fig4k and FigS6a — IMvigor210 response and ERBB2 correlations
# Original source: Figure4.R:333-592
###############################################################################
if (base::any(c("Fig4k", "FigS6a") %in% panels_to_run)) {
  #Fig4H TCGA bulk看pEMT_Basal打分与HER2表达是否负相关
  library(GSVA)
  library(SeuratExtend)
  bulk.TCGA = readRDS("/home/data/chenwenhui/TCGA/data/BLCA/01_BLCA_tpm_filt.rds")
  MP.genelist = qs_read("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/Final_MP_Top30_Signatures.qs")
  library(robustbase)
  library(org.Hs.eg.db)
  
  MP4.gene = MP.genelist$MP4
  MP9.gene = MP.genelist$MP9
  ensembleid = rownames(bulk.TCGA)
  ensembleid = gsub("\\..*","",ensembleid)
  symbol = EnsemblToGenesymbol(ensembleid,spe="human",keep.seq=T)
  bulk.log <- log2(bulk.TCGA + 1)
  
  # 2. 提取 ERBB2 的 Log 表达量
  # 处理 symbol 重复：取所有对应行的平均值
  erbb2_idx <- which(symbol == "ERBB2")
  if(length(erbb2_idx) == 0) stop("ERBB2 not found in symbol vector!")
  erbb2_expr <- colMeans(bulk.log[erbb2_idx, , drop = FALSE])
  
  # 3. 计算 MP4 和 MP9 的 Signature 分数 (基于 Log 表达量的均值)
  # 只计算在矩阵中存在的基因
  mp4_score <- colMeans(bulk.log[symbol %in% MP4.gene, , drop = FALSE])
  mp9_score <- colMeans(bulk.log[symbol %in% MP9.gene, , drop = FALSE])
  
  # 4. 构建绘图 Dataframe
  plot_df <- data.frame(
    Sample = colnames(bulk.log),
    ERBB2 = as.numeric(erbb2_expr),
    MP4_Score = as.numeric(mp4_score),
    MP9_Score = as.numeric(mp9_score)
  )
  threshold_25 <- stats::quantile(plot_df$ERBB2, 0.25)
  
  # 2. 过滤掉最低的 25% 样本，保留“HER2 活跃”群体
  plot_df <- plot_df[plot_df$ERBB2 > threshold_25, ]
  
  # 5. 绘图：展示 Log 转换后的相关性
  # MP4 相关性图
  p1 <- ggplot(plot_df, aes(x = MP4_Score, y = ERBB2)) +
    geom_point(color = "#A52A2A", alpha = 0.5, size = 1.5,shape=16) +
    geom_smooth(method = "glm", color = "black", fill = "gray80", alpha = 0.2) +
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
    theme_classic() +
    labs(title = "",
         x = "MP4 Signature Score",
         y = "ERBB2 Expression") +
    theme(aspect.ratio=1)
  
  # MP9 相关性图
  p2 <- ggplot(plot_df, aes(x = MP9_Score, y = ERBB2)) +
    geom_point(color = "#053E7A", alpha = 0.5, size = 1.5,shape=16) +
    geom_smooth(method = "glm", color = "black", fill = "gray80", alpha = 0.2) +
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
    theme_classic() +
    labs(title = "",
         x = "MP9 Signature Score",
         y = " ERBB2 Expression") + 
    theme(aspect.ratio=1)
  
  
  # 合并展示
  cowplot::plot_grid(p1, p2, ncol = 2)
  
  
  gene_sets <- list(
    MP4_Basal_pEMT = MP4.gene,
    MP9_Differentiation = MP9.gene
  )
  clean_ids <- gsub("\\..*", "", rownames(bulk.TCGA))
  
  # 2. 建立 Ensembl 到 Symbol 的映射表
  gene_map <- AnnotationDbi::select(org.Hs.eg.db, 
                                    keys = clean_ids, 
                                    columns = "SYMBOL", 
                                    keytype = "ENSEMBL") %>%
    dplyr::filter(!is.na(SYMBOL)) # 剔除没匹配上的
  
  # 3. 转换矩阵并处理重复 Symbol (取均值)
  # 这一步将 Ensembl 矩阵转换为 Symbol 矩阵，并合并重复项
  bulk_symbol <- bulk.TCGA %>%
    as.data.frame() %>%
    mutate(ENSEMBL = clean_ids) %>%
    inner_join(gene_map, by = "ENSEMBL") %>%
    dplyr::select(-ENSEMBL) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(across(everything(), mean)) %>% # 重复 Symbol 取平均值
    column_to_rownames("SYMBOL") %>%
    as.matrix()
  ssgsea_params <- ssgseaParam(bulk_symbol, gene_sets) 
  
  # 2. 计算分数
  ssgsea_scores <- gsva(ssgsea_params)
  
  # 3. 提取并过滤 (Top 75%)
  bulk_log <- log2(bulk_symbol + 1)
  
  # 2. 提取 HER2 (ERBB2) 的 Log 表达量
  # 此时提取出的是一个向量
  erbb2_log <- bulk_log["ERBB2", ]
  
  # 3. 重新构造分析 Dataframe
  plot_df_log <- data.frame(
    Sample = colnames(bulk_log),
    ERBB2 = as.numeric(erbb2_log),
    MP4_ssGSEA = as.numeric(ssgsea_scores["MP4_Basal_pEMT", ]),
    MP9_ssGSEA = as.numeric(ssgsea_scores["MP9_Differentiation", ])
  )
  
  # 4. 再次执行 25% 过滤（基于 Log 后的值）
  threshold_25 <- stats::quantile(plot_df_log$ERBB2, 0.25)
  plot_df_final <- plot_df_log[plot_df_log$ERBB2 > threshold_25, ]
  #plot_df_final = plot_df_log
  
  p1 <- ggplot(plot_df_final, aes(x = MP4_ssGSEA, y = ERBB2)) +
    geom_point(color = "#A52A2A", alpha = 0.5, size = 1.5,shape=16) +
    geom_smooth(method = "glm", color = "black", fill = "gray80", alpha = 0.2) +
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
    theme_classic() +
    labs(title = "",
         x = "MP4 Signature Score",
         y = "ERBB2 Expression") +
    theme(aspect.ratio=1)
  
  # MP9 相关性图
  p2 <- ggplot(plot_df_final, aes(x = MP9_ssGSEA, y = ERBB2)) +
    geom_point(color = "#053E7A", alpha = 0.5, size = 1.5,shape=16) +
    geom_smooth(method = "glm", color = "black", fill = "gray80", alpha = 0.2) +
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
    theme_classic() +
    labs(title = "",
         x = "MP9 Signature Score",
         y = " ERBB2 Expression") + 
    theme(aspect.ratio=1)
  
  cowplot::plot_grid(p1, p2, ncol = 2)
  
  ggsave(file.path(output_dir,"Figures","Fig","Fig4J_TCGA_BLCA_bulk.pdf"),width=6,height=3)
  
  #Imvigor210shuj
  patient.data = fread('/home/data/chenwenhui/scRNAseq/public/blca_iatlas_imvigor210_2017/data_clinical_patient.txt',
                            header=T,skip=4)
  sample.data = fread('/home/data/chenwenhui/scRNAseq/public/blca_iatlas_imvigor210_2017/data_clinical_sample.txt',
                            header=T,skip=4)
  # View(sample.data)
  exp_tpm <- fread("/home/data/chenwenhui/scRNAseq/public/blca_iatlas_imvigor210_2017/data_mrna_seq_tpm.txt", header=T,data.table = F)
  rownames(exp_tpm) = exp_tpm$Hugo_Symbol
  exp_tpm = exp_tpm[,2:ncol(exp_tpm)]
  # View(exp_tpm)
  exp_tpm =as.matrix(exp_tpm)
  all(colnames(exp_tpm)==sample.data$PATIENT_ID)
  
  MP.genelist = qs_read("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/Final_MP_Top30_Signatures.qs")
  MP4.gene = MP.genelist$MP4
  MP9.gene = MP.genelist$MP9
  mp_list <- list(
    MP4_Basal_pEMT = MP4.gene,
    MP9_Differentiation = MP9.gene
  )
  ssgsea_param <- ssgseaParam(
    exprData = exp_tpm, 
    geneSets = mp_list, 
    minSize = 5,           # 过滤掉成员过少的基因集
    maxSize = 500, 
    normalize = TRUE       # ssGSEA 通常建议进行内部归一化
  )
  mp_scores_res <- gsva(ssgsea_param)
  mp_scores_df <- as.data.frame(t(mp_scores_res))
  colnames(mp_scores_df) <- c("MP4_Basal_pEMT", "MP9_Uro_Diff")
  mp_scores_df$Sample_ID <- rownames(mp_scores_df)
  # View(mp_scores_df)
  her2_exp <- data.frame(Sample_ID = colnames(exp_tpm), 
                         HER2 = as.numeric(exp_tpm["ERBB2", ]))
  plot_df <- left_join(mp_scores_df, her2_exp, by = "Sample_ID")
  p1 <- ggplot(plot_df, aes(x = MP4_Basal_pEMT, y = HER2)) +
    geom_point(color = "#A52A2A", alpha = 0.5, size = 1.5,shape=16) +
    geom_smooth(method = "glm", color = "black", fill = "gray80", alpha = 0.2) +
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
    theme_classic() +
    labs(title = "",
         x = "MP4 Signature Score",
         y = "ERBB2 Expression") +
    theme(aspect.ratio=1)
  
  # MP9 相关性图
  p2 <- ggplot(plot_df, aes(x = MP9_Uro_Diff, y = HER2)) +
    geom_point(color = "#053E7A", alpha = 0.5, size = 1.5,shape=16) +
    geom_smooth(method = "glm", color = "black", fill = "gray80", alpha = 0.2) +
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
    theme_classic() +
    labs(title = "",
         x = "MP9 Signature Score",
         y = " ERBB2 Expression") + 
    theme(aspect.ratio=1)
  
  cowplot::plot_grid(p1, p2, ncol = 2)
  
  ggsave(file.path(output_dir,"Figures","FigS","Fig4J_IMvigor210_MPvsHER2.pdf"),width=6,height=3)
  
  #################################################
  #MP4评分，MP9评分与疗效的关系
  sample.data$SAMPLE_ID = as.character(sample.data$SAMPLE_ID)
  plot_data <- mp_scores_df %>%
    dplyr::inner_join(sample.data, by = c("Sample_ID" = "SAMPLE_ID")) %>%
    # 剔除没有响应信息的样本
    dplyr::filter(!is.na(RESPONSE) & RESPONSE != "") %>%
    # 将 RESPONSE 转为因子并排序：CR -> PR -> SD -> PD
    dplyr::mutate(RESPONSE = factor(RESPONSE, levels = c("Complete Response", 
                                                        "Partial Response", 
                                                        "Stable Disease", 
                                                        "Progressive Disease")))
  
  # 检查清洗后的数据量
  table(plot_data$RESPONSE)
  # 定义颜色：根据 Outcome 逻辑，使用色盲友好配色 [cite: 2025-12-25]
  # 这里也可以直接使用 ggpubr 的内置色板如 "npg" 或 "jco"
  res_cols <- c("#207122", "#62A765", "#E1C239", "#D55E00") # 绿到红渐变
  
  # 2. 绘制 MP4 (Basal_pEMT)
  p_mp4 <- ggpubr::ggboxplot(
    plot_data, x = "RESPONSE", y = "MP4_Basal_pEMT",
    fill = "RESPONSE", palette = res_cols,
    add = "jitter", # 叠加散点展示分布 [cite: 2026-03-04]
    add.params = list(size = 1.5, alpha = 0.4, shape = 16) # AI 编辑友好型实心圆
  ) +
    # 关键：以 Complete Response 为基准进行比较
    ggpubr::stat_compare_means(
      method = "wilcox.test", 
      label = "p.format", # 显示显著性星号，也可选 "p.format"
      ref.group = "Complete Response"
    ) +
    ggplot2::labs(title = "MP4 (Basal_pEMT) vs Atezolizumab Response",
                  y = "ssGSEA Score", x = "") +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1),
                   legend.position = "none")
  
  # 3. 绘制 MP9 (Uro_Diff)
  p_mp9 <- ggpubr::ggboxplot(
    plot_data, x = "RESPONSE", y = "MP9_Uro_Diff",
    fill = "RESPONSE", palette = res_cols,
    add = "jitter", 
    add.params = list(size = 1.5, alpha = 0.4, shape = 16)
  ) +
    ggpubr::stat_compare_means(
      method = "wilcox.test", 
      label = "p.format",
      ref.group = "Complete Response"
    ) +
    ggplot2::labs(title = "MP9 (Uro_Diff) vs Atezolizumab Response",
                  y = "ssGSEA Score", x = "") +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1),
                   legend.position = "none")
  
  # 合并展示
  library(patchwork)
  p_final <- p_mp4 | p_mp9
  print(p_final)
  ggsave(file.path(output_dir,"Figures","FigS","Fig4J_IMvigor210_MPwithRESPONSE.pdf"),width=6,height=3,plot=p_final)
}

###############################################################################
# Fig4l — experimental schematic
###############################################################################

###############################################################################
# Fig4m — T24 qPCR
# Original source: Figure4.R:922-970
###############################################################################
if ("Fig4m" %in% panels_to_run) {
  #T24
  qpcr.data = read.csv('/home/data/chenwenhui/scRNAseq/GPU_IT/QC/cell_experiment/t24-qPCR.csv')
  df_norm <- qpcr.data %>%
    dplyr::group_by(Marker) %>%
    dplyr::mutate(Ctrl_Mean = mean(Expression[Group == "Control"])) %>%
    dplyr::mutate(Fold_Change = Expression / Ctrl_Mean) %>%
    dplyr::ungroup()
  df_norm$Marker <- factor(df_norm$Marker, levels=c("E-cadherin", "N-cadherin", "Vimentin", "HER2"))
  df_norm$Group <- factor(df_norm$Group, levels=c("Control", "EMT"))
  dodge_w <- 0.5
  p_final <- ggplot2::ggplot(df_norm, aes(x = Marker, y = Fold_Change, fill = Group)) +
    # 图层 1：柱状图展示均值
    # 显式设置 position_dodge(dodge_w)
    ggplot2::geom_bar(stat = "summary", fun = mean, position = ggplot2::position_dodge(width=dodge_w), 
                    color = "black", width = 0.5, linewidth = 0.5) +
    
    # 图层 2：误差棒 (标准误 SE)
    # 必须使用相同的 dodge_w 偏移
    ggplot2::geom_errorbar(stat = "summary", fun.data = mean_se, 
                          position = ggplot2::position_dodge(width=dodge_w), 
                          width = 0.25, color = "black", linewidth = 0.5) +
    ggplot2::geom_point(aes(group = Group), 
                        position = ggplot2::position_dodge(width = dodge_w), 
                        size = 2, shape = 16, fill = "black", color = "black", 
                        stroke = 0.8, alpha = 1) +
  
    # 添加组间对比 P 值 (T 检验标注)
    ggpubr::stat_compare_means(aes(group = Group), method = "t.test", 
                              label = "p.format", hide.ns = FALSE,
                              label.y = max(df_norm$Fold_Change) * 1.2) +
    
    # 色盲友好配色 (Control: 橙色, EMT: 蓝色) [cite: 2025-12-25]
    ggplot2::scale_fill_manual(values = c("Control" = "#A52A2A", "EMT" = "#053E7A")) +
    
    ggplot2::theme_classic() +
    ggplot2::labs(title = "T24", 
                  y = "Relative mRNA Fold Change", x = "") +
    ggplot2::theme(
      axis.text.x = element_text(face = "bold", size = 12, color = "black"),
      axis.title.y = element_text(face = "bold"),
      legend.title = ggplot2::element_blank(),
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  p_final
  # 5. 保存为 PDF
  ggplot2::ggsave(file.path(output_dir,"Figures","Fig","Fig4M_qPCR_T24_EMT_Markers.pdf"), plot = p_final, width = 4, height = 3)
  
  
}

###############################################################################
# Fig4n — T24 viability
# Original source: Figure4.R:1230-1478
###############################################################################
if ("Fig4n" %in% panels_to_run) {
  #T24
  OD.data =  read.csv('/home/data/chenwenhui/scRNAseq/GPU_IT/QC/cell_experiment/T24_CCK8_data2.csv')
  baseline.data = read.csv('/home/data/chenwenhui/scRNAseq/GPU_IT/QC/cell_experiment/T24_CCK8_data2_baseline.csv')
  # View(OD.data)
  # View(baseline.data)
  mean.baseline = mean(baseline.data$baseline_OD)
  control_0 <- c(OD.data$Control1[1], OD.data$Control2[1], OD.data$Control3[1])
  treated_0 <- c(OD.data$Treated1[1], OD.data$Treated2[1], OD.data$Treated3[1])
  # 计算 Cell Viability (%)
  df_viability <- OD.data %>%
    mutate(
      # 各重复单独归一化：(OD - 空白) / (0浓度OD - 空白) * 100
      Control1_viability = (Control1 - mean.baseline) / (control_0[1] - mean.baseline) * 100,
      Control2_viability = (Control2 - mean.baseline) / (control_0[2] - mean.baseline) * 100,
      Control3_viability = (Control3 - mean.baseline) / (control_0[3] - mean.baseline) * 100,
      Treated1_viability = (Treated1 - mean.baseline) / (treated_0[1] - mean.baseline) * 100,
      Treated2_viability = (Treated2 - mean.baseline) / (treated_0[2] - mean.baseline) * 100,
      Treated3_viability = (Treated3 - mean.baseline) / (treated_0[3] - mean.baseline) * 100
    ) %>%
    dplyr::select(inducer_con, tidyselect::ends_with("_viability"))
  #print(df_viability)
  
  df_long <- df_viability %>%
    tidyr::pivot_longer(
      cols = -inducer_con,
      names_to = "Sample",
      values_to = "Viability"
    ) %>%
    dplyr::mutate(
      Group = ifelse(stringr::str_detect(Sample, "Control"), "Control", "Treated"),
      Group = factor(Group, levels = c("Control", "Treated"))
    )
  
  df_summary <- df_long %>%
    dplyr::group_by(inducer_con, Group) %>%
    dplyr::summarise(
      Mean = mean(Viability),
      SD = sd(Viability),
      SE = sd(Viability) / sqrt(dplyr::n()),
      n = dplyr::n(),
      .groups = "drop"
    )
  
  min_conc <- min(df_long$inducer_con[df_long$inducer_con > 0])
  df_long_fit <- df_long %>%
    dplyr::mutate(conc_fit = ifelse(inducer_con == 0, min_conc / 100, inducer_con))
  
  # 分组拟合 4PL 模型
  # LL.4() 是 4 参数 log-logistic 模型: Y = Bottom + (Top - Bottom) / (1 + (X/IC50)^HillSlope)
  fit_4pl <- drc::drm(Viability ~ conc_fit, 
                      curveid = Group,
                      data = df_long_fit, 
                      fct = drc::LL.4(names = c("HillSlope", "Bottom", "Top", "IC50")))
  
  ic50_values <- drc::ED(fit_4pl, 50, interval = "delta")
  fit_combined <- drc::drm(Viability ~ conc_fit, 
                           data = df_long_fit, 
                           fct = drc::LL.4())
  comparison <- anova(fit_combined, fit_4pl)
  print(comparison)
  
  ic50_comparison <- drc::compParm(fit_4pl, "IC50", "-")
  print(ic50_comparison)
  
  conc_seq <- exp(seq(log(min_conc/100), log(max(df_long$inducer_con)), length.out = 100))
  
  pred_control <- predict(fit_4pl, newdata = data.frame(conc_fit = conc_seq, Group = "Control"), 
                          interval = "confidence")
  pred_treated <- predict(fit_4pl, newdata = data.frame(conc_fit = conc_seq, Group = "Treated"), 
                          interval = "confidence")
   
  df_pred <- tibble::tibble(
    conc = rep(conc_seq, 2),
    Group = factor(rep(c("Control", "Treated"), each = length(conc_seq)), 
                   levels = c("Control", "Treated")),
    Predicted = c(pred_control[, "Prediction"], pred_treated[, "Prediction"])
  )
   
  # 如果有置信区间
  if (ncol(pred_control) > 1) {
    df_pred$Lower <- c(pred_control[, "Lower"], pred_treated[, "Lower"])
    df_pred$Upper <- c(pred_control[, "Upper"], pred_treated[, "Upper"])
  }
   
  # ============================================
  # 7. 绘制带 4PL 拟合曲线的图
  # ============================================
   
  # 处理 0 浓度用于绑图
  df_summary_plot <- df_summary %>%
    dplyr::mutate(conc_plot = ifelse(inducer_con == 0, min_conc / 100, inducer_con))
   
  df_long_plot <- df_long %>%
    dplyr::mutate(conc_plot = ifelse(inducer_con == 0, min_conc / 100, inducer_con))
   
  # 提取 IC50 数值用于标注
  ic50_control <- ic50_values["e:Control:50", "Estimate"]
  ic50_treated <- ic50_values["e:Treated:50", "Estimate"]
  
  transform_conc_logbreak <- function(x) {
    log_x <- log10(x)
    
    # 断轴参数
    log_break_start <- -4      # log10(0.0001)
    log_break_end <- -2        # log10(0.01)
    log_break_scale <- 0.20    # 压缩段占 25%
    log_max <- log10(2)               # log10(1)
    
    result <- rep(NA, length(x))
    
    # 左侧：0.0001 到 0.01（被压缩到 0-0.25）
    left_mask <- log_x >= log_break_start & log_x < log_break_end
    if (any(left_mask)) {
      result[left_mask] <- log_break_scale * (log_x[left_mask] - log_break_start) / 
                          (log_break_end - log_break_start)
    }
    
    # 右侧：0.01 到 2（保持对数，放在 0.25-2）
    right_mask <- log_x >= log_break_end
    if (any(right_mask)) {
      result[right_mask] <- log_break_scale + (2 - log_break_scale) * 
                            (log_x[right_mask] - log_break_end) / (log_max - log_break_end)
    }
    
    return(result)
  }
   
  df_summary_plot_transformed <- df_summary_plot %>%
    dplyr::mutate(
      conc_plot_transformed = transform_conc_logbreak(conc_plot)
    )
  
  # 应用变换到散点数据
  df_long_plot_transformed <- df_long_plot %>%
    dplyr::mutate(
      conc_plot_transformed = transform_conc_logbreak(conc_plot)
    )
  
  # 应用变换到拟合曲线
  df_pred_transformed <- df_pred %>%
    dplyr::mutate(
      conc_transformed = transform_conc_logbreak(conc)
    )
  
  # IC50 值变换
  ic50_control_transformed <- transform_conc_logbreak(ic50_control)
  ic50_treated_transformed <- transform_conc_logbreak(ic50_treated)
  
  # ============================================
  # 绘图
  # ============================================
  
  p_4pl <- ggplot2::ggplot() +
    
    # 图层 1：4PL 拟合曲线
    ggplot2::geom_line(data = df_pred_transformed, 
                       ggplot2::aes(x = conc_transformed, y = Predicted, color = Group), 
                       linewidth = 1.2) +
    
    # 图层 2：误差棒（用 SD 代替 SE）
    ggplot2::geom_errorbar(data = df_summary_plot_transformed, 
                           ggplot2::aes(x = conc_plot_transformed, 
                                        ymin = Mean - SD, 
                                        ymax = Mean + SD, 
                                        color = Group),
                           width = 0.05, 
                           linewidth = 0.9,
                           alpha = 0.9) +
    
    # 图层 3：均值点
    ggplot2::geom_point(data = df_summary_plot_transformed, 
                        ggplot2::aes(x = conc_plot_transformed, y = Mean, 
                                     color = Group, shape = Group),
                        size = 2, fill = "white", stroke = 1.5) +
    
    # 图层 4：各重复的原始点（透明，显示变异）
    ggplot2::geom_point(data = df_long_plot_transformed, 
                        ggplot2::aes(x = conc_plot_transformed, y = Viability, 
                                     color = Group),
                        alpha = 0.25, size = 1.5,shape=16) +
    
    # 图层 5：IC50 垂直线
    ggplot2::geom_vline(xintercept = ic50_control_transformed, 
                        color = "#A52A2A", linetype = "dashed", alpha = 0.7, linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = ic50_treated_transformed, 
                        color = "#053E7A", linetype = "dashed", alpha = 0.7, linewidth = 0.8) +
    
    # 配色
    ggplot2::scale_color_manual(values = c("Control" = "#A52A2A", "Treated" = "#053E7A")) +
    ggplot2::scale_fill_manual(values = c("Control" = "#A52A2A", "Treated" = "#053E7A")) +
    ggplot2::scale_shape_manual(values = c("Control" = 16, "Treated" = 22)) +
    
    # X 轴：变换后的坐标（修复：5个breaks对应5个labels）
    ggplot2::scale_x_continuous(
      breaks = c(transform_conc_logbreak(0.0001),
                 transform_conc_logbreak(0.01),
                 transform_conc_logbreak(0.1),
                 transform_conc_logbreak(1),
                 transform_conc_logbreak(2)),
      labels = c("0", "0.01", "0.1", "1","2"),
      limits = c(-0.02, 2.02)
    ) +
    
    # Y 轴范围
    ggplot2::ylim(0, 120) +
    
    # 主题：去掉网格线，简洁风格
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(
      title = "T24",
      x = "Drug Concentration (mg/mL)",
      y = "Cell Viability (%)",
      color = "Group",
      fill = "Group",
      shape = "Group"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black"),
      legend.position = "top",
      # 关键：去掉所有网格线
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    
    # 添加 IC50 和 p 值标注
    ggplot2::annotate("text", 
                      x = 0.02, 
                      y = 110, 
                      label = paste0("IC50 = ", round(ic50_control, 3), " mg/mL"),
                      color = "#A52A2A", size = 3.5, hjust = 0) +
    ggplot2::annotate("text", 
                      x = 0.02, 
                      y = 100, 
                      label = paste0("IC50 = ", round(ic50_treated, 3), " mg/mL"),
                      color = "#053E7A", size = 3.5, hjust = 0) +
    ggplot2::annotate("text", 
                      x = 0.02, 
                      y = 90, 
                      label = paste0("p = ", sprintf("%.2e", comparison$`p value`[2])),
                      color = "black", size = 3.5, hjust = 0)
  
  print(p_4pl)
  
  # 保存图片
  ggplot2::ggsave(file.path(output_dir, "Figures", "Fig", "Fig4N_T24_cell_viability_curve.pdf"), 
                  p_4pl, width = 4, height = 3, dpi = 300)
  
}

###############################################################################
# Supplementary Figure 4
###############################################################################

###############################################################################
# FigS4a — malignant UMAP before/after Harmony
# Original sources: scRNAseq.R:1463-1478 and FigS4.R:1-23
###############################################################################
if ("FigS4a" %in% panels_to_run) {
  malignant_umap_obj <- read_qs_object(malignant_figs3_path)
  available_reductions <- SeuratObject::Reductions(malignant_umap_obj)
  base::stopifnot(
    "umap.raw" %in% available_reductions,
    "umap.harmony" %in% available_reductions
  )
  figs4a_raw <- Seurat::DimPlot(
    malignant_umap_obj,
    reduction = "umap.raw",
    group.by = "sample",
    raster = TRUE,
    raster.dpi = c(600, 600)
  ) +
    ggplot2::ggtitle("unintegration") +
    Seurat::NoLegend() +
    ggplot2::theme_void(base_size = font_size, base_family = font_family)
  figs4a_harmony <- Seurat::DimPlot(
    malignant_umap_obj,
    reduction = "umap.harmony",
    group.by = "sample",
    raster = TRUE,
    raster.dpi = c(600, 600)
  ) +
    ggplot2::ggtitle("harmony integration") +
    Seurat::NoLegend() +
    ggplot2::theme_void(base_size = font_size, base_family = font_family)
  save_panel(
    figs4a_raw / figs4a_harmony,
    base::file.path(figs4_dir, "FigS4a_malignant_UMAP"),
    width_mm = 70,
    height_mm = 105,
    dpi = 600
  )
}

###############################################################################
# FigS4b — MP scores on UMAP
# Original source: FigS4.R:25-32
###############################################################################
if ("FigS4b" %in% panels_to_run) {
  if (!base::exists("malignant_umap_obj")) {
    malignant_umap_obj <- read_qs_object(malignant_figs3_path)
  }
  grDevices::cairo_pdf(
    base::file.path(figs4_dir, "FigS4b_MP_score_UMAP.pdf"),
    width = 137 / 25.4,
    height = 153 / 25.4
  )
  print(
    Seurat::FeaturePlot(
      malignant_umap_obj,
      features = mp_names,
      reduction = "umap.raw",
      cols = c("lightgrey", "red"),
      ncol = 3,
      raster = TRUE,
      pt.size = 0.5,
      order = TRUE
    ) &
      Seurat::NoLegend() &
      Seurat::NoAxes()
  )
  grDevices::dev.off()
}

###############################################################################
# FigS4c — remaining eight patient-level MP boxplots
###############################################################################
# Uses exactly the same patient aggregation and Wilcoxon code as Fig4c and is
# kept with that shared block above to avoid scoring the 5-GB object twice.

###############################################################################
# FigS4d — t05_CD8_HAVCR2 proportion vs Cycle_G2M / IFN_Response
# Original source: FigS3.R:1-139; reads existing Fig3 result matrices.
###############################################################################
if ("FigS4d" %in% panels_to_run) {
  fig3_mp_mean <- utils::read.csv(
    fig3_mp_mean_path,
    row.names = 1,
    check.names = FALSE
  )
  fig3_proportion <- utils::read.csv(
    fig3_proportion_path,
    check.names = FALSE
  )
  base::stopifnot(
    base::nrow(fig3_mp_mean) == base::nrow(fig3_proportion),
    "Cycle_G2M" %in% base::colnames(fig3_mp_mean),
    "IFN_Response" %in% base::colnames(fig3_mp_mean),
    "t05_CD8_HAVCR2" %in% base::colnames(fig3_proportion)
  )

  figs4d_data <- data.frame(
    sample = base::rownames(fig3_mp_mean),
    t05_proportion = fig3_proportion$t05_CD8_HAVCR2,
    Cycle_G2M = fig3_mp_mean$Cycle_G2M,
    IFN_Response = fig3_mp_mean$IFN_Response,
    stringsAsFactors = FALSE
  )
  sample_info_clean <- utils::read.csv(sample_info_path, check.names = FALSE) |>
    dplyr::mutate(
      Outcome_short = dplyr::case_when(
        .data$Outcome == "response" ~ "CR",
        .data$Outcome == "non-response" ~ "NR",
        TRUE ~ base::as.character(.data$Outcome)
      ),
      Time_short = dplyr::case_when(
        .data$Time_point == "pre-treatment" ~ "Pre",
        .data$Time_point == "post-treatment" ~ "Post",
        TRUE ~ base::as.character(.data$Time_point)
      ),
      Group = base::paste(.data$Time_short, .data$Outcome_short, sep = "_")
    )
  sample_group_map <- sample_info_clean |>
    dplyr::select(sample = "orig.ident", "Group") |>
    dplyr::distinct()
  figs4d_data <- figs4d_data |>
    dplyr::left_join(sample_group_map, by = "sample")
  figs4d_data$Group <- base::factor(figs4d_data$Group, levels = group_levels)
  figs4d_stats <- calculate_correlation_stats(
    figs4d_data,
    "t05_proportion",
    c("Cycle_G2M", "IFN_Response")
  )

  utils::write.csv(
    figs4d_data,
    base::file.path(source_data_dir, "FigS4d_sample_level.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    figs4d_stats,
    base::file.path(stats_dir, "FigS4d_Pearson.csv"),
    row.names = FALSE
  )

  figs4d_g2m <- make_correlation_plot(
    figs4d_data,
    "t05_proportion",
    "Cycle_G2M",
    "Proportion of t05_CD8_HAVCR2 (%)",
    "Cycle_G2M MP Median Score",
    figs4d_stats[figs4d_stats$y == "Cycle_G2M", , drop = FALSE],
    FALSE
  )
  figs4d_ifn <- make_correlation_plot(
    figs4d_data,
    "t05_proportion",
    "IFN_Response",
    "Proportion of t05_CD8_HAVCR2 (%)",
    "IFN_Response MP Median Score",
    figs4d_stats[figs4d_stats$y == "IFN_Response", , drop = FALSE],
    FALSE
  )
  save_panel(
    figs4d_g2m | figs4d_ifn,
    base::file.path(figs4_dir, "FigS4d_t05_proportion_MP_correlations"),
    width_mm = 137,
    height_mm = 63
  )
}

###############################################################################
# FigS4e — IFN_Response vs CD274 / CXCL10
# Original source: FigS3.R:142-218.
###############################################################################
if ("FigS4e" %in% panels_to_run) {
  malignant_figs3_obj <- read_qs_object(malignant_figs3_path)
  ifn_expression <- Seurat::FetchData(
    malignant_figs3_obj,
    vars = c("sample", "Group", "IFN_Response", "CD274", "CXCL10"),
    layer = "data"
  )
  figs4e_data <- ifn_expression |>
    dplyr::group_by(.data$sample, .data$Group) |>
    dplyr::summarise(
      Median_IFN = stats::median(.data$IFN_Response, na.rm = TRUE),
      Avg_CD274 = base::mean(.data$CD274, na.rm = TRUE),
      Avg_CXCL10 = base::mean(.data$CXCL10, na.rm = TRUE),
      .groups = "drop"
    )
  figs4e_data$Group <- base::factor(figs4e_data$Group, levels = group_levels)
  figs4e_stats <- calculate_correlation_stats(
    figs4e_data,
    "Median_IFN",
    c("Avg_CD274", "Avg_CXCL10")
  )

  utils::write.csv(
    figs4e_data,
    base::file.path(source_data_dir, "FigS4e_sample_level.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    figs4e_stats,
    base::file.path(stats_dir, "FigS4e_Pearson.csv"),
    row.names = FALSE
  )

  figs4e_cd274 <- make_correlation_plot(
    figs4e_data,
    "Median_IFN",
    "Avg_CD274",
    "Median IFN Response Score",
    "Average CD274 Expression",
    figs4e_stats[figs4e_stats$y == "Avg_CD274", , drop = FALSE],
    FALSE
  )
  figs4e_cxcl10 <- make_correlation_plot(
    figs4e_data,
    "Median_IFN",
    "Avg_CXCL10",
    "Median IFN Response Score",
    "Average CXCL10 Expression",
    figs4e_stats[figs4e_stats$y == "Avg_CXCL10", , drop = FALSE],
    TRUE
  )
  save_panel(
    figs4e_cd274 | figs4e_cxcl10,
    base::file.path(figs4_dir, "FigS4e_IFN_CD274_CXCL10_correlations"),
    width_mm = 137,
    height_mm = 76
  )
}

###############################################################################
# Supplementary Figure 5
###############################################################################

###############################################################################
# FigS5a — Slide1 and Slide3 RCTD cell-type maps
###############################################################################
# Executable code is shared with Fig4g because the input and plotting method
# are identical; only the selected slides differ.

###############################################################################
# FigS5b — all-sample Basal_pEMT / Uro_Differentiation spatial maps
###############################################################################

###############################################################################
# Supplementary Figure 6
###############################################################################

###############################################################################
# FigS6a — TCGA BLCA and IMvigor210 ERBB2 correlations
###############################################################################
# Retained in the shared Fig4k/FigS6a block because the same IMvigor210 ssGSEA
# matrix is the direct input of Fig4k.


###############################################################################
# FigS6b — TCGA ESCA correlations
# Original source: Figure4.R:596-730
###############################################################################
if ("FigS6b" %in% panels_to_run) {
  #####################################################
  bulk.TCGA = readRDS("/home/data/chenwenhui/TCGA/data/ESCA/01_ESCA_tpm_filt.rds")
  MP.genelist = qs_read("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/Final_MP_Top30_Signatures.qs")
  
  MP4.gene = MP.genelist$MP4
  MP9.gene = MP.genelist$MP9
  ensembleid = rownames(bulk.TCGA)
  ensembleid = gsub("\\..*","",ensembleid)
  symbol = EnsemblToGenesymbol(ensembleid,spe="human",keep.seq=T)
  bulk.log <- log2(bulk.TCGA + 1)
  
  # 2. 提取 ERBB2 的 Log 表达量
  # 处理 symbol 重复：取所有对应行的平均值
  erbb2_idx <- which(symbol == "ERBB2")
  if(length(erbb2_idx) == 0) stop("ERBB2 not found in symbol vector!")
  erbb2_expr <- colMeans(bulk.log[erbb2_idx, , drop = FALSE])
  
  # 3. 计算 MP4 和 MP9 的 Signature 分数 (基于 Log 表达量的均值)
  # 只计算在矩阵中存在的基因
  mp4_score <- colMeans(bulk.log[symbol %in% MP4.gene, , drop = FALSE])
  mp9_score <- colMeans(bulk.log[symbol %in% MP9.gene, , drop = FALSE])
  
  # 4. 构建绘图 Dataframe
  plot_df <- data.frame(
    Sample = colnames(bulk.log),
    ERBB2 = as.numeric(erbb2_expr),
    MP4_Score = as.numeric(mp4_score),
    MP9_Score = as.numeric(mp9_score)
  )
  threshold_25 <- stats::quantile(plot_df$ERBB2, 0.25)
  
  # 2. 过滤掉最低的 25% 样本，保留“HER2 活跃”群体
  plot_df <- plot_df[plot_df$ERBB2 > threshold_25, ]
  
  # 5. 绘图：展示 Log 转换后的相关性
  # MP4 相关性图
  p1 <- ggplot(plot_df, aes(x = MP4_Score, y = ERBB2)) +
    geom_point(color = "#A52A2A", alpha = 0.5, size = 1.5,shape=16) +
    geom_smooth(method = "glm", color = "black", fill = "gray80", alpha = 0.2) +
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
    theme_classic() +
    labs(title = "",
         x = "MP4 Signature Score",
         y = "ERBB2 Expression") +
    theme(aspect.ratio=1)
  
  # MP9 相关性图
  p2 <- ggplot(plot_df, aes(x = MP9_Score, y = ERBB2)) +
    geom_point(color = "#053E7A", alpha = 0.5, size = 1.5,shape=16) +
    geom_smooth(method = "glm", color = "black", fill = "gray80", alpha = 0.2) +
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
    theme_classic() +
    labs(title = "",
         x = "MP9 Signature Score",
         y = " ERBB2 Expression") + 
    theme(aspect.ratio=1)
  
  
  # 合并展示
  cowplot::plot_grid(p1, p2, ncol = 2)
  
  
  gene_sets <- list(
    MP4_Basal_pEMT = MP4.gene,
    MP9_Differentiation = MP9.gene
  )
  clean_ids <- gsub("\\..*", "", rownames(bulk.TCGA))
  
  # 2. 建立 Ensembl 到 Symbol 的映射表
  gene_map <- AnnotationDbi::select(org.Hs.eg.db, 
                                    keys = clean_ids, 
                                    columns = "SYMBOL", 
                                    keytype = "ENSEMBL") %>%
    dplyr::filter(!is.na(SYMBOL)) # 剔除没匹配上的
  
  # 3. 转换矩阵并处理重复 Symbol (取均值)
  # 这一步将 Ensembl 矩阵转换为 Symbol 矩阵，并合并重复项
  bulk_symbol <- bulk.TCGA %>%
    as.data.frame() %>%
    mutate(ENSEMBL = clean_ids) %>%
    inner_join(gene_map, by = "ENSEMBL") %>%
    dplyr::select(-ENSEMBL) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(across(everything(), mean)) %>% # 重复 Symbol 取平均值
    column_to_rownames("SYMBOL") %>%
    as.matrix()
  ssgsea_params <- ssgseaParam(bulk_symbol, gene_sets) 
  
  # 2. 计算分数
  ssgsea_scores <- gsva(ssgsea_params)
  
  # 3. 提取并过滤 (Top 75%)
  bulk_log <- log2(bulk_symbol + 1)
  
  # 2. 提取 HER2 (ERBB2) 的 Log 表达量
  # 此时提取出的是一个向量
  erbb2_log <- bulk_log["ERBB2", ]
  
  # 3. 重新构造分析 Dataframe
  plot_df_log <- data.frame(
    Sample = colnames(bulk_log),
    ERBB2 = as.numeric(erbb2_log),
    MP4_ssGSEA = as.numeric(ssgsea_scores["MP4_Basal_pEMT", ]),
    MP9_ssGSEA = as.numeric(ssgsea_scores["MP9_Differentiation", ])
  )
  
  # 4. 再次执行 25% 过滤（基于 Log 后的值）
  # threshold_25 <- stats::quantile(plot_df_log$ERBB2, 0.25)
  # plot_df_final <- plot_df_log[plot_df_log$ERBB2 > threshold_25, ]
  plot_df_final = plot_df_log
  
  p1 <- ggplot(plot_df_final, aes(x = MP4_ssGSEA, y = ERBB2)) +
    geom_point(color = "#A52A2A", alpha = 0.5, size = 1.5,shape=16) +
    geom_smooth(method = "glm", color = "black", fill = "gray80", alpha = 0.2) +
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
    theme_classic() +
    labs(title = "",
         x = "MP4 Signature Score",
         y = "ERBB2 Expression") +
    theme(aspect.ratio=1)
  
  # MP9 相关性图
  p2 <- ggplot(plot_df_final, aes(x = MP9_ssGSEA, y = ERBB2)) +
    geom_point(color = "#053E7A", alpha = 0.5, size = 1.5,shape=16) +
    geom_smooth(method = "glm", color = "black", fill = "gray80", alpha = 0.2) +
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
    theme_classic() +
    labs(title = "",
         x = "MP9 Signature Score",
         y = " ERBB2 Expression") + 
    theme(aspect.ratio=1)
  
  cowplot::plot_grid(p1, p2, ncol = 2)
  
  ggsave(file.path(output_dir,"Figures","FigS","Fig4J_TCGA_ESCA_bulk.pdf"),width=6,height=3)
}

###############################################################################
# FigS6c — Basal_pEMT survival panels
###############################################################################
# Exact R source was not found in the supplied GPU_IT/scripts directory.
# Current AI panel is retained; it is Outcome-independent.

###############################################################################
# FigS6d — Uro_Differentiation survival panels
###############################################################################
# Exact R source was not found in the supplied GPU_IT/scripts directory.
# Current AI panel is retained; it is Outcome-independent.

###############################################################################
# FigS6e — RT4 qPCR
# Original source: Figure4.R:871-921
###############################################################################
if ("FigS6e" %in% panels_to_run) {
  #######################################
  #细胞实验
  #Fig4M EMT qPCR RT4
  qpcr.data = read.csv('/home/data/chenwenhui/scRNAseq/GPU_IT/QC/cell_experiment/rt4-qPCR.csv')
  # View(qpcr.data)
  df_norm <- qpcr.data %>%
    dplyr::group_by(Marker) %>%
    dplyr::mutate(Ctrl_Mean = mean(Expression[Group == "Control"])) %>%
    dplyr::mutate(Fold_Change = Expression / Ctrl_Mean) %>%
    dplyr::ungroup()
  df_norm$Marker <- factor(df_norm$Marker, levels=c("E-cadherin", "N-cadherin", "Vimentin", "HER2"))
  df_norm$Group <- factor(df_norm$Group, levels=c("Control", "EMT"))
  dodge_w <- 0.5
  p_final <- ggplot2::ggplot(df_norm, aes(x = Marker, y = Fold_Change, fill = Group)) +
    # 图层 1：柱状图展示均值
    # 显式设置 position_dodge(dodge_w)
    ggplot2::geom_bar(stat = "summary", fun = mean, position = ggplot2::position_dodge(width=dodge_w), 
                    color = "black", width = 0.5, linewidth = 0.5) +
    
    # 图层 2：误差棒 (标准误 SE)
    # 必须使用相同的 dodge_w 偏移
    ggplot2::geom_errorbar(stat = "summary", fun.data = mean_se, 
                          position = ggplot2::position_dodge(width=dodge_w), 
                          width = 0.25, color = "black", linewidth = 0.5) +
    ggplot2::geom_point(aes(group = Group), 
                        position = ggplot2::position_dodge(width = dodge_w), 
                        size = 2, shape = 16, fill = "black", color = "black", 
                        stroke = 0.8, alpha = 1) +
  
    # 添加组间对比 P 值 (T 检验标注)
    ggpubr::stat_compare_means(aes(group = Group), method = "t.test", 
                              label = "p.format", hide.ns = FALSE,
                              label.y = max(df_norm$Fold_Change) * 1.2) +
    
    # 色盲友好配色 (Control: 橙色, EMT: 蓝色) [cite: 2025-12-25]
    ggplot2::scale_fill_manual(values = c("Control" = "#A52A2A", "EMT" = "#053E7A")) +
    
    ggplot2::theme_classic() +
    ggplot2::labs(title = "RT4", 
                  y = "Relative mRNA Fold Change", x = "") +
    ggplot2::theme(
      axis.text.x = element_text(face = "bold", size = 12, color = "black"),
      axis.title.y = element_text(face = "bold"),
      legend.title = ggplot2::element_blank(),
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  p_final
  # 5. 保存为 PDF
  ggplot2::ggsave(file.path(output_dir,"Figures","Fig","Fig4M_qPCR_RT4_EMT_Markers.pdf"), plot = p_final, width = 4, height = 3)
  
}

###############################################################################
# FigS6f — RT4 viability
# Original source: Figure4.R:971-1229
###############################################################################
if ("FigS6f" %in% panels_to_run) {
  ################################
  # Fig4N CCK-8 Cell viablility 曲线
  #RT4
  library(drc)
  library(ggbreak)
  OD.data =  read.csv('/home/data/chenwenhui/scRNAseq/GPU_IT/QC/cell_experiment/RT4_CCK8_data2.csv')
  baseline.data = read.csv('/home/data/chenwenhui/scRNAseq/GPU_IT/QC/cell_experiment/RT4_CCK8_data2_baseline.csv')
  # View(OD.data)
  # View(baseline.data)
  mean.baseline = mean(baseline.data$baseline_OD)
  control_0 <- c(OD.data$Control1[1], OD.data$Control2[1], OD.data$Control3[1])
  treated_0 <- c(OD.data$Treated1[1], OD.data$Treated2[1], OD.data$Treated3[1])
  # 计算 Cell Viability (%)
  df_viability <- OD.data %>%
    mutate(
      # 各重复单独归一化：(OD - 空白) / (0浓度OD - 空白) * 100
      Control1_viability = (Control1 - mean.baseline) / (control_0[1] - mean.baseline) * 100,
      Control2_viability = (Control2 - mean.baseline) / (control_0[2] - mean.baseline) * 100,
      Control3_viability = (Control3 - mean.baseline) / (control_0[3] - mean.baseline) * 100,
      Treated1_viability = (Treated1 - mean.baseline) / (treated_0[1] - mean.baseline) * 100,
      Treated2_viability = (Treated2 - mean.baseline) / (treated_0[2] - mean.baseline) * 100,
      Treated3_viability = (Treated3 - mean.baseline) / (treated_0[3] - mean.baseline) * 100
    ) %>%
    dplyr::select(inducer_con, tidyselect::ends_with("_viability"))
  #print(df_viability)
  
  df_long <- df_viability %>%
    tidyr::pivot_longer(
      cols = -inducer_con,
      names_to = "Sample",
      values_to = "Viability"
    ) %>%
    dplyr::mutate(
      Group = ifelse(stringr::str_detect(Sample, "Control"), "Control", "Treated"),
      Group = factor(Group, levels = c("Control", "Treated"))
    )
  
  df_summary <- df_long %>%
    dplyr::group_by(inducer_con, Group) %>%
    dplyr::summarise(
      Mean = mean(Viability),
      SD = sd(Viability),
      SE = sd(Viability) / sqrt(dplyr::n()),
      n = dplyr::n(),
      .groups = "drop"
    )
  
  min_conc <- min(df_long$inducer_con[df_long$inducer_con > 0])
  df_long_fit <- df_long %>%
    dplyr::mutate(conc_fit = ifelse(inducer_con == 0, min_conc / 100, inducer_con))
  
  # 分组拟合 4PL 模型
  # LL.4() 是 4 参数 log-logistic 模型: Y = Bottom + (Top - Bottom) / (1 + (X/IC50)^HillSlope)
  fit_4pl <- drc::drm(Viability ~ conc_fit, 
                      curveid = Group,
                      data = df_long_fit, 
                      fct = drc::LL.4(names = c("HillSlope", "Bottom", "Top", "IC50")))
  
  ic50_values <- drc::ED(fit_4pl, 50, interval = "delta")
  fit_combined <- drc::drm(Viability ~ conc_fit, 
                           data = df_long_fit, 
                           fct = drc::LL.4())
  comparison <- anova(fit_combined, fit_4pl)
  print(comparison)
  
  ic50_comparison <- drc::compParm(fit_4pl, "IC50", "-")
  print(ic50_comparison)
  
  conc_seq <- exp(seq(log(min_conc/100), log(max(df_long$inducer_con)), length.out = 100))
  
  pred_control <- predict(fit_4pl, newdata = data.frame(conc_fit = conc_seq, Group = "Control"), 
                          interval = "confidence")
  pred_treated <- predict(fit_4pl, newdata = data.frame(conc_fit = conc_seq, Group = "Treated"), 
                          interval = "confidence")
   
  df_pred <- tibble::tibble(
    conc = rep(conc_seq, 2),
    Group = factor(rep(c("Control", "Treated"), each = length(conc_seq)), 
                   levels = c("Control", "Treated")),
    Predicted = c(pred_control[, "Prediction"], pred_treated[, "Prediction"])
  )
   
  # 如果有置信区间
  if (ncol(pred_control) > 1) {
    df_pred$Lower <- c(pred_control[, "Lower"], pred_treated[, "Lower"])
    df_pred$Upper <- c(pred_control[, "Upper"], pred_treated[, "Upper"])
  }
   
  # ============================================
  # 7. 绘制带 4PL 拟合曲线的图
  # ============================================
   
  # 处理 0 浓度用于绑图
  df_summary_plot <- df_summary %>%
    dplyr::mutate(conc_plot = ifelse(inducer_con == 0, min_conc / 100, inducer_con))
   
  df_long_plot <- df_long %>%
    dplyr::mutate(conc_plot = ifelse(inducer_con == 0, min_conc / 100, inducer_con))
   
  # 提取 IC50 数值用于标注
  ic50_control <- ic50_values["e:Control:50", "Estimate"]
  ic50_treated <- ic50_values["e:Treated:50", "Estimate"]
  
  
  transform_conc_logbreak <- function(x) {
    log_x <- log10(x)
    
    # 断轴参数
    log_break_start <- -4      # log10(0.0001)
    log_break_end <- -2        # log10(0.01)
    log_break_scale <- 0.20    # 压缩段占 25%
    log_max <- 0               # log10(1)
    
    result <- rep(NA, length(x))
    
    # 左侧：0.0001 到 0.01（被压缩到 0-0.25）
    left_mask <- log_x >= log_break_start & log_x < log_break_end
    if (any(left_mask)) {
      result[left_mask] <- log_break_scale * (log_x[left_mask] - log_break_start) / 
                           (log_break_end - log_break_start)
    }
    
    # 右侧：0.01 到 1（保持对数，放在 0.25-1）
    right_mask <- log_x >= log_break_end
    if (any(right_mask)) {
      result[right_mask] <- log_break_scale + (1 - log_break_scale) * 
                            (log_x[right_mask] - log_break_end) / (log_max - log_break_end)
    }
    
    return(result)
  }
  
  # ============================================
  # 数据变换
  # ============================================
  
  # 应用变换到摘要数据（用 SD 作为误差）
  df_summary_plot_transformed <- df_summary_plot %>%
    dplyr::mutate(
      conc_plot_transformed = transform_conc_logbreak(conc_plot)
    )
  
  # 应用变换到散点数据
  df_long_plot_transformed <- df_long_plot %>%
    dplyr::mutate(
      conc_plot_transformed = transform_conc_logbreak(conc_plot)
    )
  
  # 应用变换到拟合曲线
  df_pred_transformed <- df_pred %>%
    dplyr::mutate(
      conc_transformed = transform_conc_logbreak(conc)
    )
  
  # IC50 值变换
  ic50_control_transformed <- transform_conc_logbreak(ic50_control)
  ic50_treated_transformed <- transform_conc_logbreak(ic50_treated)
  
  # ============================================
  # 绘图
  # ============================================
  
  p_4pl <- ggplot2::ggplot() +
    
    # 图层 1：4PL 拟合曲线
    ggplot2::geom_line(data = df_pred_transformed, 
                       ggplot2::aes(x = conc_transformed, y = Predicted, color = Group), 
                       linewidth = 1.2) +
    
    # 图层 2：误差棒（用 SD 代替 SE）
    ggplot2::geom_errorbar(data = df_summary_plot_transformed, 
                           ggplot2::aes(x = conc_plot_transformed, 
                                        ymin = Mean - SD, 
                                        ymax = Mean + SD, 
                                        color = Group),
                           width = 0.05, 
                           linewidth = 0.9,
                           alpha = 0.9) +
    
    # 图层 3：均值点
    ggplot2::geom_point(data = df_summary_plot_transformed, 
                        ggplot2::aes(x = conc_plot_transformed, y = Mean, 
                                     color = Group, shape = Group),
                        size = 2, fill = "white", stroke = 1.5) +
    
    # 图层 4：各重复的原始点（透明，显示变异）
    ggplot2::geom_point(data = df_long_plot_transformed, 
                        ggplot2::aes(x = conc_plot_transformed, y = Viability, 
                                     color = Group),
                        alpha = 0.5, size = 1.5,shape=16) +
    
    # 图层 5：IC50 垂直线
    ggplot2::geom_vline(xintercept = ic50_control_transformed, 
                        color = "#A52A2A", linetype = "dashed", alpha = 0.7, linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = ic50_treated_transformed, 
                        color = "#053E7A", linetype = "dashed", alpha = 0.7, linewidth = 0.8) +
    
    # 配色
    ggplot2::scale_color_manual(values = c("Control" = "#A52A2A", "Treated" = "#053E7A")) +
    ggplot2::scale_fill_manual(values = c("Control" = "#A52A2A", "Treated" = "#053E7A")) +
    ggplot2::scale_shape_manual(values = c("Control" = 16, "Treated" = 22)) +
    
    # X 轴：变换后的坐标（修复：5个breaks对应5个labels）
    ggplot2::scale_x_continuous(
      breaks = c(transform_conc_logbreak(0.0001),
                 transform_conc_logbreak(0.01),
                 transform_conc_logbreak(0.1),
                 transform_conc_logbreak(1)),
      labels = c("0", "0.01", "0.1", "1"),
      limits = c(-0.02, 1.02)
    ) +
    
    # Y 轴范围
    ggplot2::ylim(0, 120) +
    
    # 主题：去掉网格线，简洁风格
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(
      title = "RT4",
      x = "Drug Concentration (mg/mL)",
      y = "Cell Viability (%)",
      color = "Group",
      fill = "Group",
      shape = "Group"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black"),
      legend.position = "top",
      # 关键：去掉所有网格线
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    
    # 添加 IC50 和 p 值标注
    ggplot2::annotate("text", 
                      x = 0.02, 
                      y = 110, 
                      label = paste0("IC50 = ", round(ic50_control, 3), " mg/mL"),
                      color = "#A52A2A", size = 3.5, hjust = 0) +
    ggplot2::annotate("text", 
                      x = 0.02, 
                      y = 100, 
                      label = paste0("IC50 = ", round(ic50_treated, 3), " mg/mL"),
                      color = "#053E7A", size = 3.5, hjust = 0) +
    ggplot2::annotate("text", 
                      x = 0.02, 
                      y = 90, 
                      label = paste0("p = ", sprintf("%.2e", comparison$`p value`[2])),
                      color = "black", size = 3.5, hjust = 0)
  
  print(p_4pl)
  
  # 保存图片
  ggplot2::ggsave(file.path(output_dir, "Figures", "Fig", "Fig4N_RT4_cell_viability_curve.pdf"), 
                  p_4pl, width = 4, height = 3, dpi = 300)
  
  
}

###############################################################################
# End of review draft
###############################################################################
if (base::length(panels_to_run) == 0L) {
  message(
    "Fig4.R review mode: no panels selected; no data were read and no plots were run."
  )
}
