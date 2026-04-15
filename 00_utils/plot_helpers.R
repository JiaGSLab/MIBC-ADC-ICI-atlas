# ==================================================
# Script: Utility plotting functions
# Author: Wenhui Chen
# Purpose: Custom plotting utilities reused across multiple analysis and figure scripts.
# Input: Main Seurat or summary objects already generated in the project, plus any file paths referenced in the script.
# Output: Intermediate result files and/or figure files written to output_dir.
# ==================================================
if (TRUE) {
  suppressPackageStartupMessages({library(Seurat)
    library(patchwork)
    library(reticulate)
    library(ggrepel)
    library(SoupX)
    library(Matrix)
    library(dplyr)
    library(tidyr)
    library(rlang)
    library(future)
    library(scCustomize)
    library(BiocParallel)
    library(SeuratWrappers)
    library(ggplot2)
    library(circlize)
    library(RColorBrewer)
    library(ComplexHeatmap)
    library(DescTools)
    library(ggsci)
    library(COSG)
    #library(ROGUE)
    library(foreach)
    library(doParallel)
    library(furrr)
    library(clustree)
    library(future.apply)
    library(gplots)
    library(biomaRt)
    library(ggalluvial)
    library(msigdbr)
    library(ggsignif)
    library(ggthemr)
    library(ggunchull)
    library(corrplot)
    library(Hmisc)
    library(sscVis)
    library(data.table)
    library(data.table)
    library(circlize)
    library(grid)
    library(Startrac)
    library(readr)
    library(sscClust)
    library(pheatmap)
    library(scRNAtoolVis)
    library(DESeq2)
    library(miloR)
    library(ggpubr)
    library(enrichR)
    library(qs2)
    library(celda)
    library(sceasy)
    library(progeny)
    library(scMetabolism)
    library(rsvd)
    library(spacexr)
    library(SCNT)
    library(scStarPlot)
    library(scop)
    library(scproportion)
    library(CellChat)}
  )
  colors_15 <- c("#A52A2A", "#053E7A", "#E46B4F", "#8B8B00", "#B8860B",  
                "#CD661D", "#8B5A2B", "#CDAD00", "#2F4F4F", "#483D8B", 
                "#7A378B", "#4A7088", "#556B2F", "#008B8B", "#CD7054")
  load('~/scRNAseq/colors_90.RData') #this is a color list and can be commented out if not needed
  colors_90 = colors_91
}


do.tissueDist.full <- function(
  cellInfo.tb,
  meta.cluster = "subcelltype",
  Group.info = "Group",
  fig.output,
  show_star = TRUE,
  p_threshold = 0.05,

  # Output canvas size (compact by default)
  pdf_width = 8.5,
  pdf_height = 5.2,

  # Top bar-plot appearance parameters
  top_bar_height_cm = 3.2,
  top_bar_width = 0.98,

  # Heatmap gap parameters
  column_gap_mm = 0.8,
  row_gap_mm = 0.8,

  # Stacking order: TRUE means the bar-plot stacking order is the reverse of the heatmap row order (recommended)
  bar_stack_reverse = TRUE,

  # Colors: may be unnamed as long as length >= number of groups, or named to match Group labels
  bar_colors = NULL
) {
  # ---- Dependency check ----
  req_pkgs <- c(
    "data.table",
    "tibble",
    "ComplexHeatmap",
    "circlize",
    "grid",
    "ggsci"
  )
  miss <- req_pkgs[
    !vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(miss) > 0) {
    stop("Missing required packages: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  # ---- Basic checks ----
  cellInfo.dt <- data.table::as.data.table(cellInfo.tb)
  if (!meta.cluster %in% colnames(cellInfo.dt)) {
    stop("meta.cluster does not exist: ", meta.cluster, call. = FALSE)
  }
  if (!Group.info %in% colnames(cellInfo.dt)) {
    stop("Group.info does not exist: ", Group.info, call. = FALSE)
  }

  # Convert to factors while preserving existing levels; if not originally factors, alphabetical order will be used
  cellInfo.dt[[meta.cluster]] <- factor(cellInfo.dt[[meta.cluster]])
  cellInfo.dt[[Group.info]] <- factor(cellInfo.dt[[Group.info]])

  # Record the desired Group order for the heatmap (top to bottom)
  group_levels <- levels(cellInfo.dt[[Group.info]])

  # ---- 1) Count matrix: rows = cell types, columns = groups ----
  count_dist <- unclass(cellInfo.dt[, table(
    get(meta.cluster),
    get(Group.info)
  )])[,
    group_levels
  ]

  # Force the column order of count_dist to match group_levels
  count_dist <- count_dist[, group_levels, drop = FALSE]

  # ---- 2) Statistics (depends on your existing test.dist.table function) ----
  dist_stats <- test.dist.table(count_dist)

  OR_mtx <- data.table::dcast(dist_stats, rid ~ cid, value.var = "OR") |>
    tibble::column_to_rownames("rid") |>
    as.matrix()

  FDR_mtx <- data.table::dcast(
    dist_stats,
    rid ~ cid,
    value.var = "adj.p.value"
  ) |>
    tibble::column_to_rownames("rid") |>
    as.matrix()

  # Align matrices to the row and column order of count_dist
  OR_mtx <- OR_mtx[rownames(count_dist), colnames(count_dist), drop = FALSE]
  FDR_mtx <- FDR_mtx[rownames(count_dist), colnames(count_dist), drop = FALSE]

  # ---- 3) Fraction matrix: rows = cell types, columns = groups ----
  frac_mtx <- sweep(count_dist, 1, rowSums(count_dist), "/")

  # ---- 4) Heatmap orientation: rows = groups, columns = cell types ----
  plot_OR <- t(OR_mtx)
  plot_FDR <- t(FDR_mtx)

  # Force the heatmap row order (top to bottom) to match group_levels
  plot_OR <- plot_OR[group_levels, , drop = FALSE]
  plot_FDR <- plot_FDR[group_levels, , drop = FALSE]

  # ---- 5) Top annotation matrix: cell type × group, with row order matching the heatmap column order ----
  anno_mat <- frac_mtx[colnames(plot_OR), , drop = FALSE]
  # Force the column order to match group_levels
  anno_mat <- anno_mat[, group_levels, drop = FALSE]

  if (nrow(anno_mat) != ncol(plot_OR)) {
    stop(
      "Top annotation dimensions do not match: nrow(anno_mat)=",
      nrow(anno_mat),
      " but ncol(plot_OR)=",
      ncol(plot_OR),
      call. = FALSE
    )
  }
  if (!identical(rownames(anno_mat), colnames(plot_OR))) {
    stop("Row names of anno_mat are not aligned with column names of plot_OR.", call. = FALSE)
  }

  # ---- 6) Colors: bar_colors may be unnamed as long as length >= number of groups; named colors will be matched to Group labels ----
  n_groups <- length(group_levels)

  if (is.null(bar_colors)) {
    bar_colors_use <- ggsci::pal_npg()(n_groups)
    names(bar_colors_use) <- group_levels
  } else {
    if (!is.null(names(bar_colors)) && any(nzchar(names(bar_colors)))) {
      bar_colors_use <- bar_colors[group_levels]
      if (any(is.na(bar_colors_use))) {
        miss_g <- group_levels[is.na(bar_colors_use)]
        stop(
          "bar_colors is missing named colors for the following Groups: ",
          paste(miss_g, collapse = ", "),
          " Alternatively, provide an unnamed color vector with length >= number of groups.",
          call. = FALSE
        )
      }
      names(bar_colors_use) <- group_levels
    } else {
      if (length(bar_colors) < n_groups) {
        stop(
          "bar_colors is too short: requires >= ",
          n_groups,
          "; current length is ",
          length(bar_colors),
          "。",
          call. = FALSE
        )
      }
      bar_colors_use <- bar_colors[seq_len(n_groups)]
      names(bar_colors_use) <- group_levels
    }
  }

  # ---- 7) Key point: reverse the bar-plot stacking order relative to the heatmap row order for visual alignment ----
  stack_levels <- if (isTRUE(bar_stack_reverse)) {
    rev(group_levels)
  } else {
    group_levels
  }
  anno_mat_bar <- anno_mat[, stack_levels, drop = FALSE]

  top_anno <- ComplexHeatmap::HeatmapAnnotation(
    which = "column",
    Fraction = ComplexHeatmap::anno_barplot(
      anno_mat_bar,
      gp = grid::gpar(
        fill = bar_colors_use[colnames(anno_mat_bar)],
        col = "white",
        lwd = 0.4
      ),
      bar_width = top_bar_width,
      height = grid::unit(top_bar_height_cm, "cm")
    ),
    annotation_name_side = "left"
  )

  # ---- 8) Heatmap ----
  ht <- ComplexHeatmap::Heatmap(
    plot_OR,
    name = "Ro/e",
    col = circlize::colorRamp2(c(0, 1, 3), c("#2166ac", "white", "#b2182b")),
    top_annotation = top_anno,

    cluster_rows = FALSE,
    cluster_columns = FALSE,

    column_gap = grid::unit(column_gap_mm, "mm"),
    row_gap = grid::unit(row_gap_mm, "mm"),

    rect_gp = grid::gpar(col = "grey90", lwd = 0.5),
    row_names_side = "left",
    row_names_gp = grid::gpar(
      col = bar_colors_use[rownames(plot_OR)],
      fontface = "bold"
    ),
    column_names_rot = 45,

    cell_fun = function(j, i, x, y, w, h, fill) {
      if (isTRUE(show_star)) {
        p_val <- plot_FDR[i, j]
        if (!is.na(p_val) && p_val < p_threshold) {
          val_text <- sprintf("%.2f", plot_OR[i, j])
          text_col <- ifelse(
            plot_OR[i, j] > 2 | plot_OR[i, j] < 0.3,
            "white",
            "black"
          )
          grid::grid.text(
            val_text,
            x,
            y,
            gp = grid::gpar(fontsize = 9, col = text_col)
          )
        }
      }
    },

    heatmap_legend_param = list(
      title = "Ro/e",
      direction = "horizontal",
      at = c(0, 1, 2, 3),
      labels = c("0", "1", "2", ">=3")
    )
  )

  # ---- 9) Save ----
  grDevices::pdf(fig.output, width = pdf_width + 1, height = pdf_height + 1)
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "right",
    padding = grid::unit(c(0.5, 0.5, 0.5, 0.5), "in")
  )
  grDevices::dev.off()

  invisible(list(
    count_dist = count_dist,
    OR_mtx = OR_mtx,
    FDR_mtx = FDR_mtx,
    frac_mtx = frac_mtx,
    plot_OR = plot_OR,
    plot_FDR = plot_FDR,
    anno_mat = anno_mat,
    anno_mat_bar = anno_mat_bar,
    bar_colors = bar_colors_use,
    group_levels = group_levels,
    stack_levels = stack_levels
  ))
}


do.tissueDist <- function(
  cellInfo.tb = metadata,
  meta.cluster = cellInfo.tb$subcelltype,
  Group.info = cellInfo.tb$Group,
  fig.output,
  pdf.width = 4,
  pdf.height = 6,
  verbose = 1
) {
  dir.create(dirname(fig.output), showWarnings = FALSE, recursive = TRUE)

  cellInfo.tb <- data.table::as.data.table(cellInfo.tb)
  cellInfo.tb$meta.cluster <- as.character(meta.cluster)

  if (is.factor(Group.info)) {
    cellInfo.tb$Group.info <- Group.info
  } else {
    cellInfo.tb$Group.info <- as.factor(Group.info)
  }

  Group.info.avai.vec <- levels(cellInfo.tb[["Group.info"]])
  count.dist <- unclass(cellInfo.tb[, table(meta.cluster, Group.info)])[,
    Group.info.avai.vec
  ]
  freq.dist <- sweep(count.dist, 1, rowSums(count.dist), "/")
  freq.dist.bin <- floor(freq.dist * 100 / 10)
  print(freq.dist.bin)

  count.dist.melt.ext.tb <- test.dist.table(count.dist)

  p.dist.tb <- data.table::dcast(count.dist.melt.ext.tb, rid ~ cid, value.var = "p.value")
  OR.dist.tb <- data.table::dcast(count.dist.melt.ext.tb, rid ~ cid, value.var = "OR")
  FDR.dist.tb <- data.table::dcast(
    count.dist.melt.ext.tb,
    rid ~ cid,
    value.var = "adj.p.value"
  )
  OR.dist.mtx <- as.matrix(OR.dist.tb[, -1])
  rownames(OR.dist.mtx) <- OR.dist.tb[[1]]
  p.dist.mtx <- as.matrix(p.dist.tb[, -1])
  rownames(p.dist.mtx) <- p.dist.tb[[1]]
  FDR.dist.mtx <- as.matrix(FDR.dist.tb[, -1])
  rownames(FDR.dist.mtx) <- FDR.dist.tb[[1]]

  plot.mtx <- OR.dist.mtx
  plot.mtx <- round(plot.mtx, 2)
  upper <- min(4, max(plot.mtx, na.rm = TRUE))
  OR_for_color <- pmax(pmin(plot.mtx, upper), 0)
  piecewise_map <- function(x, upper = upper) {
    x <- pmax(pmin(x, upper), 0)
    ifelse(x <= 1, 0.5 * (x - 0) / (1 - 0), 0.5 + 0.5 * (x - 1) / (upper - 1))
  }
  OR_col_z <- piecewise_map(OR_for_color, upper)
  pal <- colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100)
  col_fun <- circlize::colorRamp2(seq(0, 1, length.out = length(pal)), pal)

  legend_breaks_raw <- c(0, 1, upper)
  legend_breaks_mapped <- piecewise_map(legend_breaks_raw, upper)

  ht <- ComplexHeatmap::Heatmap(
    OR_col_z,
    col = col_fun,
    name = "OR",
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    rect_gp = grid::gpar(col = "white", lwd = 1),
    heatmap_legend_param = list(
      at = legend_breaks_mapped,
      labels = legend_breaks_raw,
      title = 'Ro/e'
    ),
    column_names_rot = 45,
    cell_fun = function(j, i, x, y, w, h, fill) {
      if (!is.na(FDR.dist.mtx[i, j]) && FDR.dist.mtx[i, j] < 0.05) {
        grid::grid.text(
          sprintf("%.1f", OR.dist.mtx[i, j]),
          x,
          y,
          gp = grid::gpar(fontsize = 8)
        )
      }
    }
  )
  grDevices::pdf(fig.output, width = pdf.width + 1, height = pdf.height + 1)
  ComplexHeatmap::draw(ht, padding = grid::unit(c(0.5, 0.5, 0.5, 0.5), "in"))
  grDevices::dev.off()
  if (verbose == 1) {
    return(list(
      "count.dist.melt.ext.tb" = count.dist.melt.ext.tb,
      "p.dist.tb" = p.dist.tb,
      "OR.dist.tb" = OR.dist.tb,
      'FDR.dist.tb' = FDR.dist.tb,
      "OR.dist.mtx" = OR.dist.mtx,
      "p.dist.mtx" = p.dist.mtx,
      'FDR.dist.mtx' = FDR.dist.mtx
    ))
  } else {
    return(OR.dist.mtx)
  }
}


#' Perform Fisher's exact test for cell-type distribution analysis
test.dist.table <- function(count.dist, min.rowSum = 0) {
  count.dist <- count.dist[rowSums(count.dist) >= min.rowSum, , drop = FALSE]
  sum.col <- colSums(count.dist)
  sum.row <- rowSums(count.dist)
  count.dist.tb <- as.data.frame(count.dist)
  data.table::setDT(count.dist.tb, keep.rownames = TRUE)
  count.dist.melt.tb <- data.table::melt(count.dist.tb, id.vars = "rn")
  colnames(count.dist.melt.tb) <- c("rid", "cid", "count")
  
  count.dist.melt.ext.tb <- data.table::as.data.table(plyr::ldply(
    seq_len(nrow(count.dist.melt.tb)),
    function(i) {
      this.row <- count.dist.melt.tb$rid[i]
      this.col <- count.dist.melt.tb$cid[i]
      this.c <- count.dist.melt.tb$count[i]
      other.col.c <- sum.col[this.col] - this.c
      this.m <- matrix(
        c(
          this.c,
          sum.row[this.row] - this.c,
          other.col.c,
          sum(sum.col) - sum.row[this.row] - other.col.c
        ),
        ncol = 2
      )
      res.test <- fisher.test(this.m)
      data.frame(
        rid = this.row,
        cid = this.col,
        p.value = res.test$p.value,
        OR = res.test$estimate
      )
    }
  ))
  
  count.dist.melt.ext.tb <- merge(
    count.dist.melt.tb,
    count.dist.melt.ext.tb,
    by = c("rid", "cid")
  )
  
  # Ensure the object is in data.table format before using :=
  data.table::setDT(count.dist.melt.ext.tb)
  count.dist.melt.ext.tb[, adj.p.value := p.adjust(p.value, "BH")]
  
  return(count.dist.melt.ext.tb)
}


#' Plot cell-proportion distributions and compare groups
plot_cell_fraction_expert <- function(
  input_data,
  sample_col = "sample",
  group_col = "Group",
  celltype_col = "celltype",
  comparisons = NULL,
  method = "wilcox.test",
  re_normalize = FALSE,
  exclude_types = NULL,
  min_fraction = 0,
  fig_output = "cell_fraction_results.pdf",
  pdf_width = 13,
  pdf_height = 11,
  colors = NULL,
  plot_margin = c(0.3, 0.3, 0.3, 0.3)
) {
  # 1. Extract metadata
  df <- if (inherits(input_data, "Seurat")) {
    as.data.frame(input_data@meta.data)
  } else {
    as.data.frame(input_data)
  }

  # 2. Initial proportion calculation (based on raw totals)
  plot_df <- df %>%
    dplyr::group_by(
      .data[[sample_col]],
      .data[[group_col]],
      .data[[celltype_col]]
    ) %>%
    dplyr::summarise(count = dplyr::n(), .groups = 'drop') %>%
    dplyr::group_by(.data[[sample_col]]) %>%
    dplyr::mutate(fraction = count / sum(count)) %>%
    dplyr::ungroup()

  # 3. Filtering logic
  if (!is.null(exclude_types)) {
    plot_df <- plot_df %>%
      dplyr::filter(!(.data[[celltype_col]] %in% exclude_types))
  }

  if (min_fraction > 0) {
    keep <- plot_df %>%
      dplyr::group_by(.data[[celltype_col]]) %>%
      dplyr::summarise(m = mean(fraction), .groups = 'drop') %>%
      dplyr::filter(m >= min_fraction) %>%
      dplyr::pull(1)
    plot_df <- plot_df %>% dplyr::filter(.data[[celltype_col]] %in% keep)
  }

  # 4. Recalculate proportions if requested by the user
  if (re_normalize) {
    plot_df <- plot_df %>%
      dplyr::group_by(.data[[sample_col]]) %>%
      dplyr::mutate(fraction = count / sum(count)) %>%
      dplyr::ungroup()
    message(
      "--- Notice: Fractions have been re-normalized based on the remaining cell types."
    )
  }

  # 5. Fill in zeros to ensure complete statistics (revised code)
  all_samples <- unique(df[[sample_col]])
  all_celltypes <- unique(plot_df[[celltype_col]])
  sample_map <- unique(df[, c(sample_col, group_col)])

  # Use setNames instead of the !!var := syntax
  expand_df <- expand.grid(
    sample = all_samples,
    celltype = all_celltypes,
    stringsAsFactors = FALSE
  )
  colnames(expand_df) <- c(sample_col, celltype_col)
  
  plot_df <- expand_df %>%
    dplyr::left_join(sample_map, by = sample_col) %>%
    dplyr::left_join(plot_df, by = c(sample_col, celltype_col, group_col)) %>%
    dplyr::mutate(
      fraction = tidyr::replace_na(fraction, 0),
      count = tidyr::replace_na(count, 0)
    )

  # 6. Plot and save
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data[[group_col]],
      y = fraction,
      fill = .data[[group_col]]
    )
  ) +
    ggplot2::geom_boxplot(
      outlier.shape = NA,
      alpha = 0.6,
      width = 0.7,
      color = "black"
    ) +
    ggplot2::geom_jitter(
      width = 0.15,
      size = 1.2,
      shape = 21,
      color = "white",
      alpha = 0.8
    ) +
    ggplot2::facet_wrap(
      stats::as.formula(paste("~", celltype_col)),
      scales = "free_y"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      plot.margin = ggplot2::unit(plot_margin, "cm"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  # Add colors
 if (!is.null(colors)) {
    p <- p + ggplot2::scale_fill_manual(values = colors)
  }

  # Add statistical comparisons
  if (!is.null(comparisons)) {
    p <- p +
      ggpubr::stat_compare_means(
        comparisons = comparisons,
        method = method,
        label = "p.signif",
        method.args = list(exact = FALSE)
      )
  }

  # Automatic export
  utils::write.csv(
    plot_df,
    file = gsub("\\.pdf$", ".csv", fig_output),
    row.names = FALSE
  )
  ggplot2::ggsave(fig_output, p, width = pdf_width, height = pdf_height)

  return(p)
}

#' Generate a feature heatmap
#'
#' Calculate mean feature values for each subcluster based on the given feature set and grouping variable, and generate a heatmap
#'
#' @param obj single-cell object containing metadata; metadata should include the feature columns
#' @param features vector of features to display; these must already exist in metadata
#' @param group.by metadata column used for grouping; default is "subcelltype"
#' @param split_list list used for row splitting, for example list(Differentiation = c("Naive", ...))
#' @param use_split whether to enable row splitting
#' @param cluster_cols whether to cluster columns (subclusters); default is TRUE
#' @param cluster_rows whether to cluster rows (gene sets); default is FALSE
#' @param col_palette heatmap color gradient; default is blue-white-red
#' @param fig_output output PDF filename
#' @param pdf_width output PDF width (inches)
#' @param pdf_height output PDF height (inches)
#' @param plot_margin plot margins (inches)
#'
#' @return normalized feature matrix with features as rows and subclusters as columns
#'
#' @details this function will:
#' 1. calculate mean feature values for each subcluster
#' 2. Z-score normalize the feature values
#' 3. generate the heatmap according to the specified parameters and save it as a PDF
#' 4. support row splitting and clustering options
#'
#' @note make sure the features have been computed in advance, for example with AddModuleScore, and stored in metadata
do.SignatureHeatmap <- function(obj, 
                                features, 
                                group.by = "subcelltype",
                                split_list = NULL,     # Pass a list, for example list(Differentiation = c("Naive", ...))
                                use_split = FALSE,     # Whether to enable row splitting
                                cluster_cols = FALSE,   # Whether to cluster subclusters
                                cluster_rows = FALSE,  # Whether to cluster gene sets
                                col_palette = c("#2166ac", "white", "#b2182b"), # Blue-white-red palette
                                fig_output = "Signature_Heatmap.pdf",
                                pdf_width = 9,
                                pdf_height = 9,
                                plot_margin = c(0.6, 0.6, 0.6, 0.6)) {
  
  # 1. Extract score data and compute means by subcluster
  # Assume features are already stored in metadata
  meta <- obj@meta.data
  if(!all(features %in% colnames(meta))) {
    stop("Some features are not present in metadata; please run AddModuleScore or another mapping step first.")
  }
  
  plot_mtx <- meta %>%
    dplyr::group_by(.data[[group.by]]) %>%
    dplyr::summarise(across(all_of(features), mean), .groups = 'drop') %>%
    tibble::column_to_rownames(group.by) %>%
    as.matrix() %>%
    t() # Transpose so rows are gene sets and columns are subclusters

  # 2. Normalize with Z-scores to enhance contrast
  plot_mtx_scaled <- t(scale(t(plot_mtx)))
  
  # 3. Handle row-splitting logic
  row_split_vec <- NULL
  if (use_split && !is.null(split_list)) {
    # Convert the list to the corresponding vector mapping
    feature_to_group <- stack(split_list)
    row_split_vec <- feature_to_group$ind[match(rownames(plot_mtx_scaled), feature_to_group$values)]
    # Preserve the order of the input list
    row_split_vec <- factor(row_split_vec, levels = names(split_list))
  }

  # 4. Set color mapping
  # Automatically compute symmetric limits for the color scale
  limit <- max(abs(plot_mtx_scaled), na.rm = TRUE)
  col_fun <- colorRamp2(c(-limit, 0, limit), col_palette)

  # 5. Plot
  ht <- Heatmap(
    plot_mtx_scaled,
    name = "Z-score",
    col = col_fun,
    row_split = row_split_vec,
    row_title_rot = 0,               # Display row-split titles horizontally
    cluster_rows = cluster_rows,
    cluster_columns = cluster_cols,
    rect_gp = gpar(col = "grey30", lwd = 0.5), # Thin gray grid lines
    column_names_rot = 45,
    column_names_side = "bottom",
    row_names_side = "right",        # Place names on the right, following the reference figure
    row_title_gp = gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param = list(
      title = "Signature score",
      direction = "horizontal",
      title_position = "lefttop"
    ),
    width = unit(ncol(plot_mtx_scaled) * 0.8, "cm"), # Preserve tile shape
    height = unit(nrow(plot_mtx_scaled) * 0.6, "cm")
  )

  # 6. Save PDF
  pdf(fig_output, width = pdf_width, height = pdf_height)
  draw(ht, 
      heatmap_legend_side = "top",  # Place the legend on top to match the reference figure
      padding = unit(plot_margin, "in"))
  dev.off()
  
  message("--- Heatmap saved to: ", fig_output)
  return(plot_mtx_scaled)
}

plot_density_dimplot_new <- function(seurat_object, reduction = "umap", group.by, split.by = NULL, 
                                ncol = 3, density_alpha = 0.5, contour_alpha = 0.3, contour_size = 0.3, 
                                color_palette = "viridis", point_size = NULL, 
                                # New parameter to control legend point size
                                legend_point_size = 3, 
                                title = "Cell Clusters with Density Distribution", 
                                theme_style = ggplot2::theme_bw(), show_axis_text = FALSE, 
                                show_grid = FALSE, bins = 100, assay = NULL, slot = "data",
                                group_colors = NULL) 
{
    # --- 1. Parameter checks and initialization ---
    if (!inherits(seurat_object, "Seurat")) {
        stop("Input must be a Seurat object")
    }
    if (!reduction %in% names(seurat_object@reductions)) {
        stop(paste0("Reduction '", reduction, "' not found in Seurat object."))
    }

    # --- 2. Identify the type of group.by (gene or metadata) ---
    is_gene <- FALSE
    is_continuous <- FALSE
    group_data_vector <- NULL

    if (group.by %in% colnames(seurat_object@meta.data)) {
        group_data_vector <- seurat_object@meta.data[[group.by]]
        if (is.numeric(group_data_vector)) {
            is_continuous <- TRUE
            message(paste0("Metadata variable '", group.by, "' is numeric (continuous). Using gradient color scale."))
        } else {
            message(paste0("Metadata variable '", group.by, "' is categorical. Using discrete color scale."))
        }
    } else {
        if (is.null(assay)) {
            assay <- Seurat::DefaultAssay(seurat_object)
        }
        tryCatch({
            gene_data <- Seurat::FetchData(seurat_object, vars = group.by, slot = slot)
            is_gene <- TRUE
            is_continuous <- TRUE 
            message(paste0("Using '", group.by, "' as a gene/feature (continuous)."))
        }, error = function(e) {
            stop(paste0("'", group.by, "' not found in metadata or as a gene/feature."))
        })
    }

    # --- 3. Collect plotting data ---
    reduction_upper <- toupper(reduction)
    coord_data <- as.data.frame(Seurat::Embeddings(seurat_object, reduction))
    if (ncol(coord_data) < 2) {
        stop("Reduction must have at least 2 dimensions")
    }
    colnames(coord_data)[1:2] <- c("Dim_1", "Dim_2")

    if (is_gene) {
        gene_expr <- Seurat::FetchData(seurat_object, vars = group.by, slot = slot)
        coord_data[[group.by]] <- gene_expr[[group.by]]
    } else {
        coord_data[[group.by]] <- group_data_vector
        if (!is_continuous) {
            coord_data[[group.by]] <- factor(coord_data[[group.by]])
        }
    }

    if (!is.null(split.by)) {
        if (!split.by %in% colnames(seurat_object@meta.data)) {
            stop(paste0("Column '", split.by, "' not found in metadata"))
        }
        coord_data[[split.by]] <- seurat_object@meta.data[[split.by]]
    }
    
    # --- 4. Check custom color settings (for discrete variables only) ---
    if (!is_continuous && !is.null(group_colors)) {
        unique_groups <- unique(as.character(coord_data[[group.by]]))
        num_groups <- length(unique_groups)
        num_colors <- length(group_colors)
        if (num_colors < num_groups) {
            warning(paste0("Provided 'group_colors' (", num_colors, ") is less than the number of unique groups (", 
                          num_groups, "). Colors will be recycled."))
        }
    }


    # --- 5. Build the base scatter plot ---
    base_plot <- ggplot2::ggplot(coord_data, ggplot2::aes(x = .data$Dim_1, y = .data$Dim_2))
    
    # Add facets
    if (!is.null(split.by)) {
        base_plot <- base_plot + ggplot2::facet_wrap(as.formula(paste0("~", split.by)), ncol = ncol)
    }

    # Add points and colors
    base_plot <- base_plot + ggplot2::geom_point(ggplot2::aes(color = .data[[group.by]]), 
        size = ifelse(is.null(point_size), 0.5, point_size)) + 
        ggplot2::labs(x = paste0(reduction_upper, "_1"), 
                      y = paste0(reduction_upper, "_2"), color = group.by)

    # Color settings: use a gradient for continuous variables and discrete colors for categorical variables
    if (is_continuous) {
        # Continuous variable: use gradient colors
        color_scale_discrete <- ggplot2::scale_color_gradientn(
            colors = c("lightgrey", "blue", "red", "yellow"), name = group.by
        )
    } else {
        # Discrete variable: use custom colors if provided; otherwise use default discrete colors
        if (!is.null(group_colors)) {
            color_scale_discrete <- ggplot2::scale_color_manual(values = group_colors)
        } else {
            color_scale_discrete <- ggplot2::scale_color_discrete(name = group.by)
        }
        
        # ***** Key fix: use guides() to explicitly control point size in the legend *****
        # For discrete variables, explicitly set legend point size to legend_point_size
        base_plot <- base_plot + ggplot2::guides(
            color = ggplot2::guide_legend(override.aes = list(size = legend_point_size))
        )
    }
    
    base_plot <- base_plot + color_scale_discrete


    # --- 6. Apply theme and density overlay ---
    # (same as the previous version; applies the theme and adds density overlays)
    base_plot <- base_plot + theme_style
    theme_elements <- list()
    if (!show_grid) {
        theme_elements$panel.grid <- ggplot2::element_blank()
    }
    if (!show_axis_text) {
        theme_elements$axis.ticks <- ggplot2::element_blank()
        theme_elements$axis.text <- ggplot2::element_blank()
    }
    theme_elements$strip.background <- ggplot2::element_blank()
    theme_elements$strip.text <- ggplot2::element_text(face = "bold.italic", size = ggplot2::rel(1))
    base_plot <- base_plot + do.call(ggplot2::theme, theme_elements)
    
    # Overlay density plots
    final_plot <- base_plot + ggplot2::stat_density_2d(data = coord_data, 
        ggplot2::aes(x = .data$Dim_1, y = .data$Dim_2, fill = ggplot2::after_stat(density)), 
        geom = "raster", contour = FALSE, alpha = density_alpha, 
        bins = bins) + ggplot2::geom_density_2d(data = coord_data, 
        ggplot2::aes(x = .data$Dim_1, y = .data$Dim_2), color = "black", 
        size = contour_size, alpha = contour_alpha, bins = bins)

    # Density-plot color settings
    color_scale_density <- switch(tolower(color_palette), 
        viridis = viridis::scale_fill_viridis(option = "viridis"), 
        magma = viridis::scale_fill_viridis(option = "magma"), 
        plasma = viridis::scale_fill_viridis(option = "plasma"), 
        inferno = viridis::scale_fill_viridis(option = "inferno"), 
        cividis = viridis::scale_fill_viridis(option = "cividis"), 
        viridis::scale_fill_viridis(option = "viridis"))
    
    final_plot <- final_plot + color_scale_density + ggplot2::labs(title = title)
    
    return(final_plot)
}



##########################################
#' Generate a stacked bar plot of major cell types (supports custom colors)
#' @param metadata data frame containing cell information
#' @param celltype_col column name storing cell types
#' @param group_col column name storing group labels (for example, "Response")
#' @param output_name output figure filename
#' @param colors custom color vector (for example, c("T"="#E64B35", "B"="#4DBBD5")); plain color values are also accepted
plot_stacked_bar <- function(metadata, celltype_col, group_col, output_name, colors = NULL,pdf_width = 4,pdf_height = 5) {
  
  # 1. Explicitly load required packages
  require(ggplot2)
  require(dplyr)
  require(scales)
  require(grDevices)

  # 2. Data preprocessing
  plot_data <- metadata %>%
    dplyr::group_by(!!sym(group_col), !!sym(celltype_col)) %>%
    dplyr::summarise(n = n(), .groups = 'drop') %>%
    dplyr::group_by(!!sym(group_col)) %>%
    dplyr::mutate(proportion = n / sum(n))
  
  
  # 3. Handle color logic
  # Get the number of unique cell types
  unique_cells <- unique(as.character(plot_data[[celltype_col]]))
  num_cells <- length(unique_cells)
  plot_data[[celltype_col]] <- factor(plot_data[[celltype_col]], levels = levels(metadata[[celltype_col]]))

  if (is.null(colors)) {
    # Default NC-style palette
    base_palette <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", 
                      "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#B09C85FF")
    # If the number of cell types exceeds the preset palette, interpolate additional colors
    final_colors <- grDevices::colorRampPalette(base_palette)(num_cells)
  } else {
    # Check the number of custom colors; interpolate if insufficient, otherwise use them directly
    if (length(colors) < num_cells) {
      warning("The number of provided colors is smaller than the number of cell types; colors were interpolated automatically.")
      final_colors <- grDevices::colorRampPalette(colors)(num_cells)
    } else {
      final_colors <- colors
    }
  }

  # 4. Generate the plot
  p <- ggplot2::ggplot(plot_data, aes(x = !!sym(group_col), y = proportion, fill = !!sym(celltype_col))) +
    ggplot2::geom_bar(stat = "identity", width = 0.7, color = "black", size = 0.3) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
    ggplot2::scale_fill_manual(values = final_colors) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = NULL, y = "Cell Composition (%)", fill = "Cell Type") +
    ggplot2::theme(
      text = element_text(size = 12, family = "sans"),
      axis.text.x = ggplot2::element_text(color = "black", angle = 45, hjust = 1),,
      axis.text.y = element_text(color = "black", size = 10),
      axis.line = element_line(size = 0.5),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 10)
    )

  # 5. Save PDF
  ggplot2::ggsave(filename = output_name, plot = p, width = pdf_width, height = pdf_height, useDingbats = FALSE)
  
  message(paste0("Success: Stacked bar plot saved to ", output_name))
  return(p)
}



#' Draw violin plots that include Post_CR but apply differential analysis selectively
#' @param seurat_obj Seurat object of CD8 T cells
#' @param score_col column name for the exhaustion score
#' @param group_col group column name, expected to include Pre_CR, Pre_NR, Post_CR, and Post_NR
plot_exhaustion_violin_full <- function(seurat_obj, 
                                        score_col = "Exhaustion_score", 
                                        group_col = "Group") {
  
  require(Seurat)
  require(ggplot2)
  require(ggpubr)
  require(dplyr)

  # 1. Extract data and ensure Group order is logical (CR first, NR second)
  plot_data <- seurat_obj@meta.data %>%
    dplyr::select(!!sym(group_col), !!sym(score_col))
  
  # Explicitly set factor levels for a cleaner plot
  # target_levels <- c("Pre_CR", "Post_CR", "Pre_NR", "Post_NR")
  # plot_data[[group_col]] <- factor(plot_data[[group_col]], levels = target_levels)

  # 2. Define the comparison pairs for displaying P values
  # Although Post_CR is shown in the figure, significance annotations are applied only to the following two comparisons
  my_comparisons <- list(c("Pre_CR", "Pre_NR"), c("Pre_NR", "Post_NR"))

  # 3. Plot
  p <- ggplot2::ggplot(plot_data, aes(x = !!sym(group_col), y = !!sym(score_col), fill = !!sym(group_col))) +
    # Draw violin plots
    ggplot2::geom_violin(trim = FALSE, alpha = 0.7, color = "black", size = 0.3) +
    # Draw internal median and boxplot layers
    ggplot2::geom_boxplot(width = 0.1, color = "black", outlier.shape = NA, fill = "white", size = 0.3) +
    # Use a custom four-group palette (NC style)
    ggplot2::scale_fill_manual(values = c("Pre_CR" = "#00A087FF", 
                                          "Post_CR" = "#91D1C2FF",  # light green
                                          "Pre_NR" = "#E64B35FF", 
                                          "Post_NR" = "#3C5488FF")) +
    # Differential analysis: calculate and annotate only the comparisons listed in my_comparisons
    ggpubr::stat_compare_means(comparisons = my_comparisons, 
                               method = "wilcox.test", 
                               label = "p.signif",
                               step.increase = 0.1) +
    ggplot2::theme_classic() +
    ggplot2::labs(y = "Exhaustion Score", x = NULL, title = "T Cell Exhaustion Landscape") +
    ggplot2::theme(
      text = element_text(size = 12, family = "sans"),
      axis.text.x = element_text(color = "black", size = 10, angle = 45, vjust = 1, hjust = 1),
      legend.position = "none"
    )

  # Export
  ggplot2::ggsave(filename = "Fig1_Exhaustion_Violin_Full.pdf", 
                  plot = p, width = 5, height = 5.5, useDingbats = FALSE)
  
  return(p)
}

# --- Execute ---
# p_ex <- plot_exhaustion_violin_full(CD8T)
#ggplot2::ggsave(file.path(output_dir,"Figures","CD8T_Tex_score.pdf"), p_ex, width = 3, height = 4)
# Tex = subset(CD8T,subset=subcelltype=="t03_CD8_LAG3 (Tex)")
# p_ex <- plot_exhaustion_violin_full(Tex)
# ggplot2::ggsave(file.path(output_dir,"Figures","CD8Tex_Tex_score.pdf"), p_ex, width = 3, height = 4)
# p_gly <- plot_exhaustion_violin_full(Tex,score_col = 'Glycolysis_score')

#' Draw a clean correlation heatmap without numeric annotations
#' @param seurat_obj Seurat object of CD8 T cells
plot_score_correlation_clean <- function(seurat_obj) {
  
  require(dplyr)
  require(ggcorrplot)

  # 1. Extract all columns ending in _score
  score_data <- seurat_obj@meta.data %>%
    dplyr::select(ends_with("_score"))
  
  # 2. Calculate the correlation matrix
  corr_matrix <- cor(score_data, method = "pearson")
  
  # 3. Plot
  p <- ggcorrplot(corr_matrix, 
                  hc.order = TRUE,       # Use hierarchical clustering to group related scores automatically
                  type = "lower",        # Display only the lower triangle to avoid visual redundancy
                  outline.color = "white",
                  lab = FALSE,           # Key change: do not display numeric labels
                  colors = c("#2166AC", "white", "#B2182B"), # Classic RdBu palette
                  title = "Functional Program Correlation Matrix",
                  ggtheme = theme_classic) +
    ggplot2::theme(
      axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  return(p)
}

# Run and save
# p_corr <- plot_score_correlation_clean(CD8T)
# ggplot2::ggsave("Fig2_Score_Correlation_Heatmap.pdf", p_corr, width = 6, height = 6)
