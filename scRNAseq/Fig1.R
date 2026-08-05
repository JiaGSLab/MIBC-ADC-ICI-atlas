# 加载必要的包
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(patchwork)
library(stringr)
library(scales)
library(stats)

output_dir <- "~/scRNAseq/GPU_IT/QC"
source("/home/data/chenwenhui/scRNAseq/GPU_IT/scripts/function_for_plot.R", echo = TRUE)
if (!dir.exists(file.path(output_dir, "Figures", "FigS"))) {
  dir.create(file.path(output_dir, "Figures", "FigS"), recursive = TRUE)
}

#Fig1C UMAP美化
seu.combined.IT = qs_read("/home/data/chenwenhui/scRNAseq/GPU_IT/QC/seu_combined_IT_filtered_obj_intergration_new.qs")
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
p
ggsave(filename ="./Figures/Fig/Fig1D_major_Umap.png",height = 2,width = 2,plot = p,dpi = 600,bg="white")

p1<- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = subcelltype)) +
       geom_point(size = 0.15,shape = 16, stroke = 0) +
       scale_color_manual(values = c(colors_91[1:54],colors_15[5:6])) +
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
ggsave(filename ="./Figures/Fig/Fig1D_subcelltype_Umap.png",height = 4,width = 4,plot = p1,dpi = 600,bg="white")


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

ggsave("./Figures/Fig/Fig1D_subcelltype_legend.pdf", plot = p2,width = 10, height = 6, dpi = 500 )


#################
#figs


# BEGIN FigS01A_PATHOLOGY_OUTCOME
# FigS01A: patient-level pathology features by clinical outcome.
# This section is intentionally self-contained so it can be extracted and run
# without loading the Seurat object used by the earlier Fig1 panels.

first_non_missing <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  value <- x[!is.na(x)]
  if (length(value)) value[[1]] else NA_character_
}

prepare_figs01a_data <- function(sample_info_path) {
  required_columns <- c(
    "Sample_ID", "Therapy", "Outcome", "Sex", "Age", "Grade",
    "Invasive", "Differentiation"
  )
  patient_fields <- c(
    "Therapy", "Outcome", "Sex", "Age", "Grade", "Invasive",
    "Differentiation"
  )

  raw_data <- utils::read.csv(
    sample_info_path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fill = TRUE,
    strip.white = TRUE
  )

  missing_columns <- setdiff(required_columns, colnames(raw_data))
  if (length(missing_columns)) {
    stop(
      "sample_info.csv is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  normalized_data <- raw_data |>
    dplyr::filter(trimws(.data$Therapy) == "IT") |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(patient_fields),
        ~ dplyr::na_if(trimws(as.character(.x)), "")
      )
    )

  conflicts <- normalized_data |>
    dplyr::group_by(.data$Sample_ID) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(patient_fields),
        ~ dplyr::n_distinct(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )
  if (any(as.matrix(conflicts[, patient_fields, drop = FALSE]) > 1L)) {
    stop("Conflicting patient-level values detected in the IT cohort")
  }

  patient_data <- normalized_data |>
    dplyr::group_by(.data$Sample_ID) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(patient_fields), first_non_missing),
      n_specimens = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      Outcome = dplyr::recode(
        .data$Outcome,
        response = "CR",
        `non-response` = "NR",
        .default = NA_character_
      ),
      Outcome = factor(.data$Outcome, levels = c("CR", "NR")),
      Age = suppressWarnings(as.numeric(.data$Age)),
      Sex = factor(.data$Sex, levels = c("female", "male")),
      Grade = factor(.data$Grade, levels = c("high", "low")),
      Invasive = factor(
        .data$Invasive,
        levels = c("invasive", "non-invasive")
      ),
      Differentiation = factor(
        .data$Differentiation,
        levels = c("no", "yes")
      )
    ) |>
    dplyr::arrange(.data$Outcome, .data$Sample_ID)

  outcome_counts <- table(patient_data$Outcome)
  valid_counts <- identical(names(outcome_counts), c("CR", "NR")) &&
    identical(as.integer(outcome_counts), c(27L, 21L))
  if (nrow(patient_data) != 48L || !valid_counts) {
    stop("Unexpected IT cohort size or Outcome counts")
  }
  patient_data
}

figs01a_palette <- c("CR" = "#D55E00", "NR" = "#0072B2")
figs01a_width_mm <- 189
figs01a_height_mm <- 59.3
figs01a_base_size <- 8
figs01a_font_family <- "Arial"

format_figs01a_p <- function(p_value) {
  if (!is.finite(p_value)) {
    return("P = NA")
  }
  if (p_value < 0.001) {
    return("P < 0.001")
  }
  paste0("P = ", format(round(p_value, 3), nsmall = 0, trim = TRUE))
}

theme_figs01a <- function() {
  ggplot2::theme_classic(
    base_size = figs01a_base_size,
    base_family = figs01a_font_family
  ) +
    ggplot2::theme(
      text = ggplot2::element_text(
        family = figs01a_font_family,
        size = figs01a_base_size,
        colour = "black"
      ),
      panel.border = ggplot2::element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.35
      ),
      axis.line = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(
        family = figs01a_font_family,
        size = figs01a_base_size,
        face = "plain"
      ),
      axis.text = ggplot2::element_text(
        family = figs01a_font_family,
        size = figs01a_base_size,
        colour = "black"
      ),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      axis.ticks = ggplot2::element_line(linewidth = 0.3),
      plot.title = ggplot2::element_text(
        family = figs01a_font_family,
        size = figs01a_base_size,
        face = "plain",
        hjust = 0.5,
        margin = ggplot2::margin(b = 2, unit = "pt")
      ),
      plot.margin = ggplot2::margin(2, 2, 1, 1, unit = "pt"),
      legend.title = ggplot2::element_text(
        family = figs01a_font_family,
        size = figs01a_base_size
      ),
      legend.text = ggplot2::element_text(
        family = figs01a_font_family,
        size = figs01a_base_size
      ),
      legend.key.size = grid::unit(8, "pt"),
      legend.spacing.y = grid::unit(1, "pt")
    )
}

plot_figs01a_age <- function(data) {
  plot_data <- data |>
    dplyr::filter(!is.na(.data$Age), !is.na(.data$Outcome))

  p_value <- stats::wilcox.test(
    Age ~ Outcome,
    data = plot_data,
    alternative = "two.sided",
    exact = FALSE
  )$p.value

  age_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$Outcome, y = .data$Age, fill = .data$Outcome)
  ) +
    ggplot2::geom_boxplot(
      width = 0.55,
      outlier.shape = NA,
      colour = "black",
      linewidth = 0.35,
      alpha = 0.7
    ) +
    ggplot2::geom_jitter(
      shape = 21,
      colour = "black",
      stroke = 0.15,
      width = 0.12,
      height = 0,
      size = 1.2,
      alpha = 1
    ) +
    ggplot2::scale_fill_manual(
      values = figs01a_palette,
      limits = names(figs01a_palette),
      drop = FALSE,
      name = "Outcome",
      guide = "none"
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.04, 0.18))
    ) +
    ggplot2::annotate(
      "text",
      x = 1.5,
      y = max(plot_data$Age, na.rm = TRUE),
      vjust = -1.1,
      label = format_figs01a_p(p_value),
      family = figs01a_font_family,
      size = figs01a_base_size / ggplot2::.pt
    ) +
    ggplot2::labs(title = "Age", x = NULL, y = "Years") +
    theme_figs01a()

  attr(age_plot, "p_value") <- unname(p_value)
  age_plot
}

plot_figs01a_categorical <- function(data, variable, title = variable) {
  plot_data <- data |>
    dplyr::filter(!is.na(.data[[variable]]), !is.na(.data$Outcome))

  contingency_table <- table(plot_data[[variable]], plot_data$Outcome)
  chi_result <- suppressWarnings(stats::chisq.test(contingency_table))
  if (any(chi_result$expected < 5)) {
    p_value <- stats::fisher.test(contingency_table)$p.value
  } else {
    p_value <- chi_result$p.value
  }

  proportion_data <- plot_data |>
    dplyr::count(.data[[variable]], .data$Outcome, name = "n") |>
    dplyr::group_by(.data[[variable]]) |>
    dplyr::mutate(proportion = .data$n / sum(.data$n)) |>
    dplyr::ungroup()

  category_count <- nlevels(droplevels(plot_data[[variable]]))
  categorical_plot <- ggplot2::ggplot(
    proportion_data,
    ggplot2::aes(
      x = .data[[variable]],
      y = .data$proportion,
      fill = .data$Outcome
    )
  ) +
    ggplot2::geom_col(
      width = 0.7,
      colour = "black",
      linewidth = 0.35,
      alpha = 0.7
    ) +
    ggplot2::scale_fill_manual(
      values = figs01a_palette,
      limits = names(figs01a_palette),
      drop = FALSE,
      name = "Outcome"
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 1, by = 0.25),
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1.1),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::annotate(
      "text",
      x = mean(seq_len(category_count)),
      y = 1.055,
      label = format_figs01a_p(p_value),
      family = figs01a_font_family,
      size = figs01a_base_size / ggplot2::.pt
    ) +
    ggplot2::labs(title = title, x = NULL, y = "Proportion") +
    theme_figs01a()

  attr(categorical_plot, "p_value") <- unname(p_value)
  categorical_plot
}

build_figs01a <- function(
  sample_info_path = "/home/data/chenwenhui/sample_info.csv",
  output_directory = "/home/data/chenwenhui/scRNAseq/GPU_IT/QC/Figures/FigS"
) {
  required_packages <- c("dplyr", "ggplot2", "patchwork", "scales")
  unavailable_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(unavailable_packages)) {
    stop(
      "Missing required R packages: ",
      paste(unavailable_packages, collapse = ", ")
    )
  }
  if (!capabilities("cairo")) {
    stop("This R installation does not support Cairo PDF output")
  }

  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  patient_data <- prepare_figs01a_data(sample_info_path)

  panels <- list(
    Age = plot_figs01a_age(patient_data),
    Sex = plot_figs01a_categorical(patient_data, "Sex", "Gender"),
    Grade = plot_figs01a_categorical(patient_data, "Grade", "Grade"),
    Invasive = plot_figs01a_categorical(
      patient_data,
      "Invasive",
      "Invasive"
    ),
    Differentiation = plot_figs01a_categorical(
      patient_data,
      "Differentiation",
      "Differentiation"
    )
  )
  statistics <- vapply(
    panels,
    function(panel) attr(panel, "p_value"),
    numeric(1)
  )

  combined_plot <- patchwork::wrap_plots(
    panels,
    nrow = 1,
    widths = c(1.0, 1.0, 1.0, 1.15, 1.15),
    guides = "collect"
  ) &
    ggplot2::theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.justification = "center",
      text = ggplot2::element_text(
        family = figs01a_font_family,
        size = figs01a_base_size
      )
    )

  pdf_path <- file.path(
    output_directory,
    "FigS01A_pathology_outcome.pdf"
  )
  csv_path <- file.path(
    output_directory,
    "FigS01A_pathology_outcome_source_data.csv"
  )

  export_data <- patient_data |>
    dplyr::mutate(
      dplyr::across(
        where(is.factor),
        as.character
      )
    )
  utils::write.csv(
    export_data,
    csv_path,
    row.names = FALSE,
    na = "",
    fileEncoding = "UTF-8"
  )

  set.seed(20260731)
  ggplot2::ggsave(
    filename = pdf_path,
    plot = combined_plot,
    width = figs01a_width_mm,
    height = figs01a_height_mm,
    units = "mm",
    device = grDevices::cairo_pdf,
    bg = "white"
  )

  list(
    plot = combined_plot,
    panels = panels,
    data = patient_data,
    statistics = statistics,
    palette = figs01a_palette,
    pdf_path = pdf_path,
    csv_path = csv_path,
    width_mm = figs01a_width_mm,
    height_mm = figs01a_height_mm
  )
}

if (!identical(Sys.getenv("FIGS01A_TEST_MODE"), "true")) {
  figs01a_result <- build_figs01a()
  message("FigS01A PDF: ", figs01a_result$pdf_path)
  message("FigS01A source data: ", figs01a_result$csv_path)
}

# END FigS01A_PATHOLOGY_OUTCOME

#####################################################

# BEGIN FigS1_SIMPLE_REVISED
# FigS1: keep the original sequential plotting style; only standardize
# colours, font, size and source-data export.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(qs2)
  library(Seurat)
  library(coin)
  library(ggsignif)
  library(ComplexHeatmap)
  library(circlize)
})

figs1_output_dir <- "/home/data/chenwenhui/scRNAseq/GPU_IT/QC/Figures/FigS/FigS01_revised"
figs1_object_path <- "/home/data/chenwenhui/scRNAseq/GPU_IT/QC/seu_combined_IT_filtered_obj_intergration_new.qs"
dir.create(figs1_output_dir, recursive = TRUE, showWarnings = FALSE)

if (!exists("test.dist.table", mode = "function") ||
    !exists("colors_15", inherits = TRUE) ||
    !exists("colors_91", inherits = TRUE)) {
  suppressWarnings(source(
    "/home/data/chenwenhui/scRNAseq/GPU_IT/scripts/function_for_plot.R",
    echo = FALSE
  ))
}

figs1_font <- "Arial"
figs1_size <- 8
figs1_outcome_colors <- c(CR = "#D55E00", NR = "#0072B2")
figs1_time_colors <- c(Pre = "#D55E00", Post = "#0072B2")
figs1_group_levels <- c("Pre_NR", "Pre_CR", "Post_NR", "Post_CR")
figs1_major_levels <- c(
  "T cell", "B cell", "Myeloid", "Fibroblast",
  "Endothelial", "Mural cells", "Epithelial", "malignant"
)
figs1_major_colors <- setNames(colors_15[seq_along(figs1_major_levels)], figs1_major_levels)

theme_figs1 <- function() {
  ggplot2::theme_classic(base_size = figs1_size, base_family = figs1_font) +
    ggplot2::theme(
      text = ggplot2::element_text(family = figs1_font, size = figs1_size, colour = "black"),
      axis.title = ggplot2::element_text(family = figs1_font, size = figs1_size),
      axis.text = ggplot2::element_text(family = figs1_font, size = figs1_size, colour = "black"),
      legend.title = ggplot2::element_text(family = figs1_font, size = figs1_size),
      legend.text = ggplot2::element_text(family = figs1_font, size = figs1_size),
      plot.title = ggplot2::element_text(family = figs1_font, size = figs1_size, hjust = 0.5),
      plot.margin = ggplot2::margin(2, 2, 2, 2, unit = "pt")
    )
}

write_figs1_csv <- function(x, filename) {
  x[] <- lapply(x, function(z) if (is.factor(z)) as.character(z) else z)
  write.csv(
    x,
    file.path(figs1_output_dir, filename),
    row.names = FALSE,
    na = "",
    fileEncoding = "UTF-8"
  )
}

# FigS1A -----------------------------------------------------------------
# Reuse the already verified FigS01A functions above in this same Fig1.R.
figs1_A_result <- build_figs01a(
  sample_info_path = "/home/data/chenwenhui/sample_info.csv",
  output_directory = figs1_output_dir
)
FigS1A <- figs1_A_result$plot
figs1_A_statistics <- data.frame(
  variable = names(figs1_A_result$statistics),
  p_value = unname(figs1_A_result$statistics),
  row.names = NULL
)
write_figs1_csv(figs1_A_statistics, "FigS01A_pathology_outcome_statistics.csv")

# Load the corrected Seurat object once for FigS1B--G --------------------
seu.combined.IT <- qs2::qs_read(figs1_object_path)
required_meta <- c("sample", "Sample_ID", "Outcome", "Time_point", "Group", "major", "subcelltype")
stopifnot(all(required_meta %in% colnames(seu.combined.IT[[]])))
seu.combined.IT$Outcome <- factor(as.character(seu.combined.IT$Outcome), levels = c("CR", "NR"))
seu.combined.IT$Time_point <- factor(as.character(seu.combined.IT$Time_point), levels = c("Pre", "Post"))
seu.combined.IT$Group <- factor(as.character(seu.combined.IT$Group), levels = figs1_group_levels)

# FigS1B: Sample, Time point and Outcome UMAP ----------------------------
seu.combined.IT$sampleid <- paste0(seu.combined.IT$Sample_ID, "-", seu.combined.IT$Time_point)
sample_levels <- sort(unique(seu.combined.IT$sampleid))
sample_palette <- setNames(
  rep(unname(unlist(colors_91)), length.out = length(sample_levels)),
  sample_levels
)

umap_xy <- SeuratObject::Embeddings(seu.combined.IT, "umap")
x_range <- range(umap_xy[, 1], finite = TRUE)
y_range <- range(umap_xy[, 2], finite = TRUE)
x0 <- x_range[1] + 0.04 * diff(x_range)
y0 <- y_range[1] + 0.05 * diff(y_range)
x1 <- x0 + 0.18 * diff(x_range)
y1 <- y0 + 0.18 * diff(y_range)

umap_theme <- ggplot2::theme_void(base_size = figs1_size, base_family = figs1_font) +
  ggplot2::theme(
    text = ggplot2::element_text(family = figs1_font, size = figs1_size, colour = "black"),
    plot.title = ggplot2::element_text(family = figs1_font, size = figs1_size, hjust = 0.5),
    legend.position = "bottom",
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(family = figs1_font, size = figs1_size),
    legend.key.size = grid::unit(5, "pt"),
    plot.margin = ggplot2::margin(1, 1, 1, 1, unit = "pt")
  )

FigS1B_sample <- Seurat::DimPlot(
  seu.combined.IT,
  reduction = "umap",
  group.by = "sampleid",
  cols = sample_palette,
  pt.size = 0.15,
  raster = TRUE,
  raster.dpi = c(600, 600)
) +
  Seurat::NoAxes() + ggplot2::ggtitle("Sample") + umap_theme + Seurat::NoLegend() +
  ggplot2::annotate(
    "segment", x = x0, xend = x1, y = y0, yend = y0,
    linewidth = 0.3, arrow = grid::arrow(length = grid::unit(1.4, "mm"), type = "closed")
  ) +
  ggplot2::annotate(
    "segment", x = x0, xend = x0, y = y0, yend = y1,
    linewidth = 0.3, arrow = grid::arrow(length = grid::unit(1.4, "mm"), type = "closed")
  ) +
  ggplot2::annotate("text", x = x1, y = y0, label = "UMAP_1", hjust = 0.9, vjust = 1.8,
           family = figs1_font, size = figs1_size / .pt) +
  ggplot2::annotate("text", x = x0, y = y1, label = "UMAP_2", angle = 90, hjust = 0.9, vjust = 1.8,
           family = figs1_font, size = figs1_size / .pt)

FigS1B_timepoint <- Seurat::DimPlot(
  seu.combined.IT,
  reduction = "umap",
  group.by = "Time_point",
  cols = figs1_time_colors,
  pt.size = 0.15,
  raster = TRUE,
  raster.dpi = c(600, 600)
) +
  Seurat::NoAxes() + ggplot2::ggtitle("Time point") + umap_theme +
  ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1, override.aes = list(size = 1.5, alpha = 1)))

FigS1B_outcome <- Seurat::DimPlot(
  seu.combined.IT,
  reduction = "umap",
  group.by = "Outcome",
  cols = figs1_outcome_colors,
  pt.size = 0.15,
  raster = TRUE,
  raster.dpi = c(600, 600)
) +
  Seurat::NoAxes() + ggplot2::ggtitle("Outcome") + umap_theme +
  ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1, override.aes = list(size = 1.5, alpha = 1)))

FigS1B <- FigS1B_sample | FigS1B_timepoint | FigS1B_outcome
ggplot2::ggsave(
  file.path(figs1_output_dir, "FigS01B_umap.pdf"),
  FigS1B,
  width = 126,
  height = 50,
  units = "mm",
  device = grDevices::cairo_pdf,
  bg = "white"
)

# FigS1C: marker FeaturePlot ---------------------------------------------
figs1_markers <- c("KRT19", "CD3D", "CD79A", "CD14", "CSF3R", "COL1A1", "VWF", "ACTA2", "RGS5")
stopifnot(all(figs1_markers %in% rownames(seu.combined.IT)))
marker_values <- SeuratObject::FetchData(seu.combined.IT, vars = figs1_markers, layer = "data")
marker_limits <- range(unlist(marker_values, use.names = FALSE), finite = TRUE)
rm(marker_values)

figs1_C_panels <- Seurat::FeaturePlot(
  seu.combined.IT,
  features = figs1_markers,
  cols = c("grey90", "#EF3B2C"),
  keep.scale = "all",
  combine = FALSE,
  raster = TRUE,
  raster.dpi = c(600, 600),
  pt.size = 0.15
)

for (i in seq_along(figs1_C_panels)) {
  figs1_C_panels[[i]] <- figs1_C_panels[[i]] +
    ggplot2::scale_colour_gradientn(
      colours = c("grey90", "#EF3B2C"),
      limits = marker_limits,
      breaks = marker_limits,
      labels = c("Min", "Max"),
      oob = scales::squish
    ) +
    ggplot2::labs(title = figs1_markers[i], colour = NULL) +
    ggplot2::theme_void(base_size = figs1_size, base_family = figs1_font) +
    ggplot2::theme(
      text = ggplot2::element_text(family = figs1_font, size = figs1_size, colour = "black"),
      plot.title = ggplot2::element_text(family = figs1_font, size = figs1_size, hjust = 0.5),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(family = figs1_font, size = figs1_size),
      legend.key.width = grid::unit(7, "mm"),
      legend.key.height = grid::unit(1.4, "mm"),
      plot.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt")
    )
}

figs1_C_panels[[7]] <- figs1_C_panels[[7]] +
  ggplot2::annotate(
    "segment", x = x0, xend = x1, y = y0, yend = y0,
    linewidth = 0.3, arrow = grid::arrow(length = grid::unit(1.2, "mm"), type = "closed")
  ) +
  ggplot2::annotate(
    "segment", x = x0, xend = x0, y = y0, yend = y1,
    linewidth = 0.3, arrow = grid::arrow(length = grid::unit(1.2, "mm"), type = "closed")
  ) +
  ggplot2::annotate("text", x = x1, y = y0, label = "UMAP_1", hjust = 0.9, vjust = 1.8,
           family = figs1_font, size = figs1_size / .pt) +
  ggplot2::annotate("text", x = x0, y = y1, label = "UMAP_2", angle = 90, hjust = 0.9, vjust = 1.8,
           family = figs1_font, size = figs1_size / .pt)

FigS1C <- patchwork::wrap_plots(figs1_C_panels, ncol = 3, guides = "collect") &
  ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(
  file.path(figs1_output_dir, "FigS01C_marker_featureplot.pdf"),
  FigS1C,
  width = 64,
  height = 65,
  units = "mm",
  device = grDevices::cairo_pdf,
  bg = "white"
)
rm(figs1_C_panels)
gc()

# Metadata used by D--G --------------------------------------------------
figs1_meta <- seu.combined.IT[[]] |>
  dplyr::mutate(
    Group = factor(as.character(Group), levels = figs1_group_levels),
    Outcome = factor(as.character(Outcome), levels = c("CR", "NR")),
    Time_point = factor(as.character(Time_point), levels = c("Pre", "Post")),
    major = factor(as.character(major), levels = figs1_major_levels)
  )

# FigS1D: major-cell composition ----------------------------------------
figs1_D_data <- figs1_meta |>
  dplyr::count(Group, major, .drop = FALSE, name = "n") |>
  dplyr::group_by(Group) |>
  dplyr::mutate(proportion = n / sum(n)) |>
  dplyr::ungroup()
write_figs1_csv(figs1_D_data, "FigS01D_major_composition_source_data.csv")

FigS1D <- ggplot2::ggplot(figs1_D_data, ggplot2::aes(Group, proportion, fill = major)) +
  ggplot2::geom_col(width = 0.72, colour = "black", linewidth = 0.25, alpha = 0.7) +
  ggplot2::scale_fill_manual(values = figs1_major_colors, drop = FALSE) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
  ggplot2::labs(x = NULL, y = "Proportion of cell clusters", fill = "Clusters") +
  theme_figs1() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.key.size = grid::unit(7, "pt"))
ggplot2::ggsave(
  file.path(figs1_output_dir, "FigS01D_major_composition.pdf"),
  FigS1D,
  width = 58,
  height = 65,
  units = "mm",
  device = grDevices::cairo_pdf,
  bg = "white"
)

# FigS1E: four selected major-cell proportions --------------------------
sample_info_e <- figs1_meta |>
  dplyr::distinct(sample, Sample_ID, Time_point, Outcome)
figs1_E_all_data <- figs1_meta |>
  dplyr::count(sample, major, .drop = FALSE, name = "count") |>
  dplyr::left_join(sample_info_e, by = "sample") |>
  dplyr::group_by(sample) |>
  dplyr::mutate(sample_total = sum(count), percentage = 100 * count / sample_total) |>
  dplyr::ungroup()
write_figs1_csv(figs1_E_all_data, "FigS01E_major_proportion_source_data.csv")

selected_majors <- c("malignant", "Fibroblast", "Epithelial", "B cell")
major_titles <- c(malignant = "Malignant", Fibroblast = "Fibroblast", Epithelial = "Epithelial", `B cell` = "B cells")
figs1_E_data <- figs1_E_all_data |>
  dplyr::filter(major %in% selected_majors) |>
  dplyr::mutate(major = factor(as.character(major), levels = selected_majors))

run_coin_comparison <- function(data, group_column, first_level, second_level) {
  test_data <- data |>
    dplyr::filter(.data[[group_column]] %in% c(first_level, second_level)) |>
    dplyr::mutate(test_group = factor(as.character(.data[[group_column]]), levels = c(first_level, second_level)))
  medians <- tapply(test_data$percentage, test_data$test_group, median, na.rm = TRUE)
  alternative <- if (medians[2] >= medians[1]) "less" else "greater"
  one_sided <- coin::independence_test(
    percentage ~ test_group,
    data = test_data,
    distribution = "exact",
    alternative = alternative
  )
  two_sided <- coin::independence_test(
    percentage ~ test_group,
    data = test_data,
    distribution = "exact",
    alternative = "two.sided"
  )
  data.frame(
    first_level = first_level,
    second_level = second_level,
    n_first = sum(test_data$test_group == first_level),
    n_second = sum(test_data$test_group == second_level),
    median_first = unname(medians[1]),
    median_second = unname(medians[2]),
    alternative_one_sided = alternative,
    p_one_sided = unname(coin::pvalue(one_sided)),
    p_two_sided = unname(coin::pvalue(two_sided))
  )
}

figs1_E_statistics <- lapply(selected_majors, function(m) {
  dm <- dplyr::filter(figs1_E_data, major == m)
  dplyr::bind_rows(
    cbind(major = m, comparison = "CR_Pre_vs_Post", displayed = TRUE,
          run_coin_comparison(dplyr::filter(dm, Outcome == "CR"), "Time_point", "Pre", "Post")),
    cbind(major = m, comparison = "NR_Pre_vs_Post", displayed = TRUE,
          run_coin_comparison(dplyr::filter(dm, Outcome == "NR"), "Time_point", "Pre", "Post")),
    cbind(major = m, comparison = "Pre_CR_vs_NR", displayed = TRUE,
          run_coin_comparison(dplyr::filter(dm, Time_point == "Pre"), "Outcome", "CR", "NR")),
    cbind(major = m, comparison = "Post_CR_vs_NR", displayed = FALSE,
          run_coin_comparison(dplyr::filter(dm, Time_point == "Post"), "Outcome", "CR", "NR"))
  )
}) |>
  dplyr::bind_rows()
write_figs1_csv(figs1_E_statistics, "FigS01E_major_proportion_statistics.csv")

FigS1E_panels <- lapply(selected_majors, function(m) {
  dm <- dplyr::filter(figs1_E_data, major == m)
  sm <- dplyr::filter(figs1_E_statistics, major == m, displayed)
  p_labels <- paste0("p = ", formatC(sm$p_one_sided, digits = 3, format = "f"))
  y_data_max <- max(dm$percentage, na.rm = TRUE)
  y_step <- max(5, 0.15 * y_data_max)
  y_positions <- y_data_max + y_step * c(0.55, 1.2, 2.15)

  ggplot2::ggplot(dm, ggplot2::aes(Outcome, percentage, fill = Time_point)) +
    ggplot2::stat_summary(
      fun = mean, geom = "bar", position = ggplot2::position_dodge(0.72),
      width = 0.62, colour = "black", linewidth = 0.3, alpha = 0.7
    ) +
    ggplot2::stat_summary(
      fun.data = ggplot2::mean_se, geom = "errorbar", position = ggplot2::position_dodge(0.72),
      width = 0.18, linewidth = 0.3
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = Time_point),
      position = ggplot2::position_jitterdodge(jitter.width = 0.04, dodge.width = 0.72),
      size = 0.9, show.legend = FALSE
    ) +
    ggsignif::geom_signif(
      annotations = p_labels,
      xmin = c(0.82, 1.82, 0.82),
      xmax = c(1.18, 2.18, 1.82),
      y_position = y_positions,
      tip_length = 0.01,
      textsize = figs1_size / .pt,
      size = 0.3,
      vjust = 0.2
    ) +
    ggplot2::scale_fill_manual(values = figs1_time_colors, drop = FALSE) +
    ggplot2::scale_colour_manual(values = figs1_time_colors, drop = FALSE) +
    ggplot2::scale_y_continuous(limits = c(0, max(y_positions) + y_step), expand = c(0, 0)) +
    ggplot2::labs(title = unname(major_titles[m]), x = NULL, y = "Proportion (%)", fill = NULL) +
    theme_figs1()
})
names(FigS1E_panels) <- selected_majors
FigS1E <- patchwork::wrap_plots(FigS1E_panels, nrow = 1, guides = "collect") &
  ggplot2::theme(legend.position = "right")
ggplot2::ggsave(
  file.path(figs1_output_dir, "FigS01E_major_proportion.pdf"),
  FigS1E,
  width = 126,
  height = 42,
  units = "mm",
  device = grDevices::cairo_pdf,
  bg = "white"
)

# Ro/e calculation shared by FigS1F and FigS1G --------------------------
calculate_roe <- function(metadata, cluster_column, cluster_levels = NULL) {
  cluster <- as.character(metadata[[cluster_column]])
  if (is.null(cluster_levels)) cluster_levels <- sort(unique(cluster))
  cluster <- factor(cluster, levels = cluster_levels)
  group <- factor(as.character(metadata$Group), levels = figs1_group_levels)
  count_matrix <- unclass(table(cluster, group))
  test_table <- test.dist.table(count_matrix)

  to_matrix <- function(value_column) {
    wide <- data.table::dcast(data.table::as.data.table(test_table), rid ~ cid, value.var = value_column)
    matrix_value <- as.matrix(wide[, -1])
    rownames(matrix_value) <- wide$rid
    matrix_value[rownames(count_matrix), colnames(count_matrix), drop = FALSE]
  }

  list(
    count = count_matrix,
    OR = to_matrix("OR"),
    p = to_matrix("p.value"),
    FDR = to_matrix("adj.p.value")
  )
}

roe_long_table <- function(result, cluster_name) {
  out <- expand.grid(
    cluster = rownames(result$count),
    Group = colnames(result$count),
    stringsAsFactors = FALSE
  )
  names(out)[1] <- cluster_name
  out$count <- as.vector(result$count)
  out$RoE <- as.vector(result$OR)
  out$p_value <- as.vector(result$p)
  out$FDR <- as.vector(result$FDR)
  out
}

roe_colors <- circlize::colorRamp2(c(0, 1, 3), c("#2166AC", "white", "#B2182B"))

# FigS1F: major Ro/e heatmap --------------------------------------------
figs1_F_result <- calculate_roe(figs1_meta, "major", figs1_major_levels)
write_figs1_csv(roe_long_table(figs1_F_result, "major"), "FigS01F_major_RoE_source_data.csv")
figs1_F_OR <- t(figs1_F_result$OR)
figs1_F_FDR <- t(figs1_F_result$FDR)

FigS1F <- ComplexHeatmap::Heatmap(
  figs1_F_OR,
  name = "Ro/e",
  col = roe_colors,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  rect_gp = grid::gpar(col = "grey90", lwd = 0.4),
  row_names_gp = grid::gpar(
    fontfamily = figs1_font, fontsize = figs1_size,
    col = c(Pre_NR = "#A52A2A", Pre_CR = "#053E7A", Post_NR = "#E46B4F", Post_CR = "#8B8B00")[rownames(figs1_F_OR)]
  ),
  column_names_gp = grid::gpar(fontfamily = figs1_font, fontsize = figs1_size),
  column_names_rot = 45,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(figs1_F_FDR[i, j]) && figs1_F_FDR[i, j] < 0.05) {
      grid::grid.text(
        sprintf("%.2f", figs1_F_OR[i, j]), x, y,
        gp = grid::gpar(
          fontfamily = figs1_font, fontsize = figs1_size,
          col = ifelse(figs1_F_OR[i, j] > 2 || figs1_F_OR[i, j] < 0.3, "white", "black")
        )
      )
    }
  },
  heatmap_legend_param = list(
    title = "Ro/e", direction = "horizontal",
    title_gp = grid::gpar(fontfamily = figs1_font, fontsize = figs1_size),
    labels_gp = grid::gpar(fontfamily = figs1_font, fontsize = figs1_size)
  )
)
grDevices::cairo_pdf(
  file.path(figs1_output_dir, "FigS01F_major_RoE.pdf"),
  width = 100 / 25.4,
  height = 50 / 25.4,
  family = figs1_font
)
ComplexHeatmap::draw(FigS1F, heatmap_legend_side = "bottom", padding = grid::unit(rep(2, 4), "mm"))
dev.off()

# FigS1G: subcelltype Ro/e heatmap --------------------------------------
figs1_G_result <- calculate_roe(figs1_meta, "subcelltype")
write_figs1_csv(
  roe_long_table(figs1_G_result, "subcelltype"),
  "FigS01G_subcelltype_RoE_source_data.csv"
)
figs1_G_OR <- figs1_G_result$OR
figs1_G_FDR <- figs1_G_result$FDR

FigS1G <- ComplexHeatmap::Heatmap(
  figs1_G_OR,
  name = "Ro/e",
  col = roe_colors,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_side = "right",
  rect_gp = grid::gpar(col = "white", lwd = 0.5),
  row_names_gp = grid::gpar(fontfamily = figs1_font, fontsize = figs1_size),
  column_names_gp = grid::gpar(fontfamily = figs1_font, fontsize = figs1_size),
  column_names_rot = 45,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(figs1_G_FDR[i, j]) && figs1_G_FDR[i, j] < 0.05) {
      grid::grid.text(
        sprintf("%.1f", figs1_G_OR[i, j]), x, y,
        gp = grid::gpar(fontfamily = figs1_font, fontsize = figs1_size)
      )
    }
  },
  heatmap_legend_param = list(
    title = "Ro/e",
    title_gp = grid::gpar(fontfamily = figs1_font, fontsize = figs1_size),
    labels_gp = grid::gpar(fontfamily = figs1_font, fontsize = figs1_size)
  )
)
grDevices::cairo_pdf(
  file.path(figs1_output_dir, "FigS01G_subcelltype_RoE.pdf"),
  width = 59 / 25.4,
  height = 207 / 25.4,
  family = figs1_font
)
FigS1G_drawn <- ComplexHeatmap::draw(FigS1G, heatmap_legend_side = "bottom", padding = grid::unit(rep(2, 4), "mm"))
dev.off()

figs1_G_order <- data.frame(
  plot_order = seq_along(ComplexHeatmap::row_order(FigS1G_drawn)),
  subcelltype = rownames(figs1_G_OR)[ComplexHeatmap::row_order(FigS1G_drawn)]
)
write_figs1_csv(figs1_G_order, "FigS01G_subcelltype_RoE_order.csv")

message("FigS1 A--G saved in: ", figs1_output_dir)
# END FigS1_SIMPLE_REVISED
