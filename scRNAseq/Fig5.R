###############################################################################
# Figure 5 and Supplementary Figure 7 revised panels
#
# Final manuscript letters are used below. The old scripts used shifted names:
# current Fig5C = old Fig5B, current Fig5D = old Fig5C,
# current Fig5E = old Fig5D, and current Fig5H = old Fig5G.
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
})

required_packages <- c(
  "dplyr", "tidyr", "tibble", "ggplot2", "ggrepel", "ggrastr", "RANN",
  "GSEABase", "future", "coin", "ComplexHeatmap", "circlize",
  "car", "emmeans", "sandwich", "lmtest", "ragg"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

set.seed(123)
options(future.globals.maxSize = 600 * 1024^3)

# independent of those fields, so they are not rebuilt in the default run.
# Set either option to "1" only when that panel itself needs to be redrawn.
run_spatial_panels <- identical(
  Sys.getenv("RUN_FIG5_SPATIAL", unset = "0"),
  "1"
)
run_mouse_panel <- identical(
  Sys.getenv("RUN_FIG5_MOUSE", unset = "0"),
  "1"
)

fig5_output_dir <- "/home/data/chenwenhui/scRNAseq/GPU_IT/QC/Figures/Fig/Fig5_revised"
figs7_output_dir <- "/home/data/chenwenhui/scRNAseq/GPU_IT/QC/Figures/FigS/FigS7_revised"
dir.create(fig5_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figs7_output_dir, recursive = TRUE, showWarnings = FALSE)

fig5_font <- "Arial"
fig5_size <- 8
outcome_colors <- c(CR = "#D55E00", NR = "#0072B2")
program_colors <- c(Basal_pEMT = "#D55E00", Uro_Differentiation = "#56B4E9")

# ComplexHeatmap measures text through the PDF font database before drawing.
# Register Arial metrics there; cairo_pdf still embeds the installed Arial TTF.
if (!fig5_font %in% names(grDevices::pdfFonts())) {
  arial_pdf_font <- grDevices::pdfFonts()[["Helvetica"]]
  arial_pdf_font$family <- fig5_font
  do.call(grDevices::pdfFonts, stats::setNames(list(arial_pdf_font), fig5_font))
}

theme_fig5 <- ggplot2::theme_classic(base_size = fig5_size, base_family = fig5_font) +
  ggplot2::theme(
    text = ggplot2::element_text(family = fig5_font, size = fig5_size, colour = "black"),
    axis.title = ggplot2::element_text(family = fig5_font, size = fig5_size),
    axis.text = ggplot2::element_text(family = fig5_font, size = fig5_size, colour = "black"),
    legend.title = ggplot2::element_text(family = fig5_font, size = fig5_size),
    legend.text = ggplot2::element_text(family = fig5_font, size = fig5_size),
    plot.margin = ggplot2::margin(2, 2, 2, 2, unit = "pt")
  )

###############################################################################
# Shared spatial preparation for Fig5B, Fig5C and FigS7C
###############################################################################

if (run_spatial_panels) {

# Exact input used by the original Figure4.R/Fig5.R workflow. These panels do
spatial_object_path <- "/home/data/chenwenhui/visiumHD/ST.combined.obj.qs"
mp_signature_path <- "/home/data/chenwenhui/scRNAseq/GPU_IT/QC/Final_MP_Top30_Signatures.qs"

HD.combined.obj <- qs2::qs_read(spatial_object_path)
spatial_meta <- HD.combined.obj[[]]
required_spatial_meta <- c(
  "sample_id", "sample", "Outcome", "Therapy", "orig.ident",
  "nCount_Spatial", "nFeature_Spatial", "first_type", "second_type"
)
stopifnot(all(required_spatial_meta %in% colnames(spatial_meta)))


# Preserve the original Figure4.R metadata correction for ZZJ_8916.
spatial_outcome[as.character(HD.combined.obj$sample) == "ZZJ_8916"] <- "CR"
HD.combined.obj$Outcome <- spatial_outcome

HD.combined.obj <- base::subset(HD.combined.obj, subset = Therapy == "IT")
HD.combined.obj <- Seurat::NormalizeData(HD.combined.obj, verbose = FALSE)
HD.combined.obj$pct_counts_mt <- Seurat::PercentageFeatureSet(
  HD.combined.obj,
  pattern = "^MT-"
)
HD.combined.obj <- base::subset(
  HD.combined.obj,
  subset = nCount_Spatial > 50 & nCount_Spatial < 4000 &
    nFeature_Spatial > 50 & nFeature_Spatial < 2000 & pct_counts_mt < 25
)
stopifnot(length(unique(as.character(HD.combined.obj$sample_id))) == 11L)

MP.genelist <- qs2::qs_read(mp_signature_path)
mp_names <- c(
  "Cycle_G2M", "Surface_Trafficking", "RNA_Splicing", "Basal_pEMT",
  "Stress_IEG", "IFN_Response", "Hypoxia", "Metabolic",
  "Uro_Differentiation", "Ep_Inflammation"
)
stopifnot(length(MP.genelist) == length(mp_names))

HD.combined.obj@meta.data <- HD.combined.obj@meta.data[
  , setdiff(colnames(HD.combined.obj@meta.data), mp_names), drop = FALSE
]
HD.combined.obj <- Seurat::AddModuleScore(
  object = HD.combined.obj,
  features = MP.genelist,
  name = "MPscore",
  seed = 1
)
module_columns <- paste0("MPscore", seq_along(mp_names))
stopifnot(all(module_columns %in% colnames(HD.combined.obj[[]])))
colnames(HD.combined.obj@meta.data)[
  match(module_columns, colnames(HD.combined.obj@meta.data))
] <- mp_names

known_cells <- rownames(HD.combined.obj[[]])[
  !is.na(HD.combined.obj$first_type) & HD.combined.obj$first_type != "unknown"
]
RCTD.obj <- base::subset(HD.combined.obj, cells = known_cells)

epi_clusters <- c(
  "e01_Epi_PSCA", "e02_Epi_MKI67", "e03_Epi_KRT5", "e04_Epi_NDRG1",
  "e05_Epi_FOS", "e06_Epi_S100A9", "e07_Epi_KRT13", "e08_Epi_CXCL10",
  "e09_Epi_PIGR", "e10_Epi_MMP7", "e11_Epi_KRT16"
)
RCTD.obj$first_type <- ifelse(
  RCTD.obj$first_type %in% epi_clusters,
  "Epithelial",
  as.character(RCTD.obj$first_type)
)
RCTD.obj$second_type <- ifelse(
  RCTD.obj$second_type %in% epi_clusters,
  "Epithelial",
  as.character(RCTD.obj$second_type)
)

###############################################################################
# Fig5B: continuous bivariate MP-score spatial scatter
###############################################################################

pt32_cells <- rownames(RCTD.obj[[]])[
  as.character(RCTD.obj$sample_id) == "pt32"
]
stopifnot(length(pt32_cells) > 0L)
pt32_obj <- base::subset(RCTD.obj, cells = pt32_cells)
pt32_coordinates <- Seurat::GetTissueCoordinates(pt32_obj)
stopifnot(
  nrow(pt32_coordinates) > 0L,
  all(c("x", "y") %in% colnames(pt32_coordinates))
)

pt32_meta <- pt32_obj[[]][
  rownames(pt32_coordinates),
  c("Basal_pEMT", "Uro_Differentiation", "first_type", "second_type"),
  drop = FALSE
]
pt32_data <- data.frame(
  spot = rownames(pt32_coordinates),
  new_x = pt32_coordinates$y,
  new_y = -pt32_coordinates$x,
  is_endothelial = as.character(pt32_meta$first_type) == "Endothelial" |
    as.character(pt32_meta$second_type) == "Endothelial",
  Basal_pEMT = as.numeric(pt32_meta$Basal_pEMT),
  Uro_Differentiation = as.numeric(pt32_meta$Uro_Differentiation),
  stringsAsFactors = FALSE
)

valid_score <- is.finite(pt32_data$Basal_pEMT) &
  is.finite(pt32_data$Uro_Differentiation)
pt32_data <- pt32_data[
  pt32_data$is_endothelial | valid_score,
  , drop = FALSE
]
non_endothelial <- !pt32_data$is_endothelial
stopifnot(any(non_endothelial), any(pt32_data$is_endothelial))

score_percentile <- function(values, reference) {
  sorted_reference <- sort(reference[is.finite(reference)])
  if (length(sorted_reference) < 2L) {
    stop("At least two finite reference scores are required.")
  }
  percentile <- (findInterval(values, sorted_reference) - 1) /
    (length(sorted_reference) - 1)
  pmin(1, pmax(0, percentile))
}

bivariate_corner_colors <- c(
  low_low = "#E8E8E8",
  high_basal = "#E69F00",
  high_uro = "#56B4E9",
  high_both = "#7030A0"
)
interpolate_bivariate_colour <- function(basal, uro) {
  stopifnot(
    length(basal) == length(uro),
    all(is.finite(basal)),
    all(is.finite(uro)),
    all(basal >= 0 & basal <= 1),
    all(uro >= 0 & uro <= 1)
  )
  corner_matrix <- do.call(rbind, lapply(
    unname(bivariate_corner_colors),
    function(colour) grDevices::col2rgb(colour)[, 1] / 255
  ))
  weights <- cbind(
    (1 - basal) * (1 - uro),
    basal * (1 - uro),
    (1 - basal) * uro,
    basal * uro
  )
  rgb_value <- weights %*% corner_matrix
  grDevices::rgb(rgb_value[, 1], rgb_value[, 2], rgb_value[, 3])
}

pt32_data$Basal_percentile <- NA_real_
pt32_data$Uro_percentile <- NA_real_
pt32_data$bivariate_color <- NA_character_
pt32_data$Basal_percentile[non_endothelial] <- score_percentile(
  pt32_data$Basal_pEMT[non_endothelial],
  pt32_data$Basal_pEMT[non_endothelial]
)
pt32_data$Uro_percentile[non_endothelial] <- score_percentile(
  pt32_data$Uro_Differentiation[non_endothelial],
  pt32_data$Uro_Differentiation[non_endothelial]
)
pt32_data$bivariate_color[non_endothelial] <- interpolate_bivariate_colour(
  pt32_data$Basal_percentile[non_endothelial],
  pt32_data$Uro_percentile[non_endothelial]
)
stopifnot(
  all(pt32_data$Basal_percentile[non_endothelial] >= 0 &
    pt32_data$Basal_percentile[non_endothelial] <= 1),
  all(pt32_data$Uro_percentile[non_endothelial] >= 0 &
    pt32_data$Uro_percentile[non_endothelial] <= 1),
  length(unique(pt32_data$bivariate_color[non_endothelial])) > 256L
)

utils::write.csv(
  pt32_data[, c(
    "spot", "new_x", "new_y", "is_endothelial",
    "Basal_pEMT", "Uro_Differentiation",
    "Basal_percentile", "Uro_percentile", "bivariate_color"
  )],
  file.path(fig5_output_dir, "Fig5B_bivariate_continuous_data.csv"),
  row.names = FALSE
)

scatter_data_B <- pt32_data[non_endothelial, , drop = FALSE]
scatter_data_B <- scatter_data_B[order(
  scatter_data_B$Basal_percentile + scatter_data_B$Uro_percentile
), , drop = FALSE]
endothelial_data_B <- pt32_data[pt32_data$is_endothelial, , drop = FALSE]

Fig5B <- ggplot2::ggplot() +
  ggplot2::geom_point(
    data = scatter_data_B,
    ggplot2::aes(
      x = .data$new_x,
      y = .data$new_y,
      colour = .data$bivariate_color
    ),
    size = 0.10, alpha = 0.70
  ) +
  ggplot2::geom_point(
    data = endothelial_data_B,
    ggplot2::aes(x = .data$new_x, y = .data$new_y),
    colour = "black", size = 0.14
  ) +
  ggplot2::scale_colour_identity() +
  ggplot2::coord_fixed(
    xlim = range(pt32_data$new_x, na.rm = TRUE),
    ylim = range(pt32_data$new_y, na.rm = TRUE),
    expand = FALSE,
    clip = "on"
  ) +
  ggplot2::theme_void(base_size = fig5_size, base_family = fig5_font) +
  ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt"))

ggplot2::ggsave(
  filename = file.path(
    fig5_output_dir,
    "Fig5B_pt32_bivariate_scatter_continuous.png"
  ),
  plot = Fig5B, width = 32, height = 44, units = "mm", dpi = 600,
  device = ragg::agg_png, bg = "white"
)

legend_data_B <- expand.grid(
  Basal_percentile = (seq_len(100L) - 0.5) / 100,
  Uro_percentile = (seq_len(100L) - 0.5) / 100
)
legend_data_B$bivariate_color <- interpolate_bivariate_colour(
  legend_data_B$Basal_percentile,
  legend_data_B$Uro_percentile
)

Fig5B_legend <- ggplot2::ggplot(
  legend_data_B,
  ggplot2::aes(x = .data$Basal_percentile, y = .data$Uro_percentile)
) +
  ggplot2::geom_raster(ggplot2::aes(fill = .data$bivariate_color)) +
  ggplot2::scale_fill_identity() +
  ggplot2::scale_x_continuous(
    breaks = NULL, limits = c(0, 1),
    expand = ggplot2::expansion(mult = 0)
  ) +
  ggplot2::scale_y_continuous(
    breaks = NULL, limits = c(0, 1),
    expand = ggplot2::expansion(mult = 0)
  ) +
  ggplot2::coord_fixed(clip = "off") +
  ggplot2::labs(x = "Basal →", y = "Uro →") +
  ggplot2::theme_minimal(base_size = fig5_size, base_family = fig5_font) +
  ggplot2::theme(
    text = ggplot2::element_text(
      family = fig5_font, size = fig5_size, colour = "black"
    ),
    axis.title = ggplot2::element_text(
      family = fig5_font, size = fig5_size, colour = "black"
    ),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_text(
      margin = ggplot2::margin(t = 1, unit = "pt")
    ),
    axis.title.y = ggplot2::element_text(
      margin = ggplot2::margin(r = 1, unit = "pt")
    ),
    plot.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt")
  )

ggplot2::ggsave(
  filename = file.path(fig5_output_dir, "Fig5B_bivariate_legend_continuous.pdf"),
  plot = Fig5B_legend, width = 20, height = 20, units = "mm",
  device = grDevices::cairo_pdf, family = fig5_font, bg = "white"
)

rm(pt32_obj, pt32_meta, pt32_coordinates)

###############################################################################
# Fig5C and FigS7C: MP scores along distance to vasculature
###############################################################################

calculate_niche_stats_500 <- function(seurat_obj, sample_name) {
  sample_cells <- rownames(seurat_obj[[]])[
    as.character(seurat_obj$sample_id) == sample_name
  ]
  sub_obj <- base::subset(seurat_obj, cells = sample_cells)
  coordinates <- Seurat::GetTissueCoordinates(sub_obj)
  sub_meta <- sub_obj[[]][rownames(coordinates), , drop = FALSE]

  is_endothelial <- sub_meta$first_type == "Endothelial" |
    sub_meta$second_type == "Endothelial"
  is_epithelial <- sub_meta$first_type == "Epithelial" |
    sub_meta$second_type == "Epithelial"

  endothelial_coordinates <- as.matrix(coordinates[is_endothelial, c("x", "y"), drop = FALSE])
  epithelial_coordinates <- as.matrix(coordinates[is_epithelial, c("x", "y"), drop = FALSE])
  if (nrow(endothelial_coordinates) == 0L || nrow(epithelial_coordinates) == 0L) {
    return(NULL)
  }

  minimum_distance <- RANN::nn2(
    data = endothelial_coordinates,
    query = epithelial_coordinates,
    k = 1
  )$nn.dists[, 1]

  sample_data <- data.frame(
    sample = sample_name,
    distance = minimum_distance,
    Basal_pEMT = sub_meta$Basal_pEMT[is_epithelial],
    Uro_Differentiation = sub_meta$Uro_Differentiation[is_epithelial],
    stringsAsFactors = FALSE
  )

  sample_data |>
    dplyr::filter(.data$distance <= 500) |>
    dplyr::mutate(
      bin = cut(
        .data$distance,
        breaks = seq(0, 500, by = 50),
        labels = FALSE,
        include.lowest = TRUE
      ),
      dist_bin = (.data$bin - 0.5) * 50
    ) |>
    dplyr::filter(!is.na(.data$bin)) |>
    dplyr::group_by(.data$sample, .data$dist_bin) |>
    dplyr::summarise(
      mean_Basal = mean(.data$Basal_pEMT, na.rm = TRUE),
      mean_Uro = mean(.data$Uro_Differentiation, na.rm = TRUE),
      n_epithelial_spots = dplyr::n(),
      .groups = "drop"
    )
}

sample_list <- sort(unique(as.character(RCTD.obj$sample_id)))
all_sample_stats_500 <- dplyr::bind_rows(lapply(
  sample_list,
  function(sample_name) calculate_niche_stats_500(RCTD.obj, sample_name)
))
stopifnot(
  nrow(all_sample_stats_500) > 0L,
  dplyr::n_distinct(all_sample_stats_500$sample) == 11L
)

summary_500 <- all_sample_stats_500 |>
  dplyr::group_by(.data$dist_bin) |>
  dplyr::summarise(
    n_samples = dplyr::n(),
    grand_mean_Basal = mean(.data$mean_Basal, na.rm = TRUE),
    grand_mean_Uro = mean(.data$mean_Uro, na.rm = TRUE),
    se_Basal = stats::sd(.data$mean_Basal, na.rm = TRUE) / sqrt(.data$n_samples),
    se_Uro = stats::sd(.data$mean_Uro, na.rm = TRUE) / sqrt(.data$n_samples),
    .groups = "drop"
  )

ratio <- max(summary_500$grand_mean_Basal, na.rm = TRUE) /
  max(summary_500$grand_mean_Uro, na.rm = TRUE)
stopifnot(is.finite(ratio), ratio > 0)

spatial_sample_audit <- RCTD.obj[[]] |>
  dplyr::distinct(
    .data$orig.ident, .data$sample_id, .data$sample, .data$Outcome, .data$Therapy
  ) |>
  dplyr::arrange(.data$orig.ident, .data$sample_id)

utils::write.csv(
  all_sample_stats_500,
  file.path(fig5_output_dir, "Fig5C_FigS7C_per_sample_bins.csv"),
  row.names = FALSE
)
utils::write.csv(
  summary_500,
  file.path(fig5_output_dir, "Fig5C_FigS7C_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  spatial_sample_audit,
  file.path(fig5_output_dir, "Fig5C_sample_audit.csv"),
  row.names = FALSE
)

Fig5C <- ggplot2::ggplot(summary_500, ggplot2::aes(x = .data$dist_bin)) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = .data$grand_mean_Basal - .data$se_Basal,
      ymax = .data$grand_mean_Basal + .data$se_Basal
    ),
    fill = program_colors[["Basal_pEMT"]], alpha = 0.10
  ) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = (.data$grand_mean_Uro - .data$se_Uro) * ratio,
      ymax = (.data$grand_mean_Uro + .data$se_Uro) * ratio
    ),
    fill = program_colors[["Uro_Differentiation"]], alpha = 0.10
  ) +
  ggplot2::geom_line(
    ggplot2::aes(y = .data$grand_mean_Basal, colour = "Basal_pEMT"),
    linewidth = 0.8
  ) +
  ggplot2::geom_line(
    ggplot2::aes(y = .data$grand_mean_Uro * ratio, colour = "Uro_Differentiation"),
    linewidth = 0.8
  ) +
  ggplot2::geom_point(
    ggplot2::aes(y = .data$grand_mean_Basal),
    colour = program_colors[["Basal_pEMT"]], size = 1.4, shape = 16
  ) +
  ggplot2::geom_point(
    ggplot2::aes(y = .data$grand_mean_Uro * ratio),
    colour = program_colors[["Uro_Differentiation"]], size = 1.4, shape = 16
  ) +
  ggplot2::scale_colour_manual(values = program_colors) +
  ggplot2::scale_y_continuous(
    name = "Basal_pEMT score",
    sec.axis = ggplot2::sec_axis(~ . / ratio, name = "Uro_Differentiation score")
  ) +
  ggplot2::labs(x = "Distance to vasculature (µm)", colour = NULL) +
  theme_fig5 +
  ggplot2::theme(
    axis.title.y.left = ggplot2::element_text(colour = program_colors[["Basal_pEMT"]]),
    axis.title.y.right = ggplot2::element_text(colour = program_colors[["Uro_Differentiation"]]),
    legend.position = "top",
    legend.direction = "horizontal"
  ) +
  ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1))

FigS7C <- ggplot2::ggplot() +
  ggplot2::geom_line(
    data = all_sample_stats_500,
    ggplot2::aes(
      x = .data$dist_bin, y = .data$mean_Basal,
      group = .data$sample
    ),
    colour = program_colors[["Basal_pEMT"]], alpha = 0.20, linewidth = 0.35
  ) +
  ggplot2::geom_line(
    data = all_sample_stats_500,
    ggplot2::aes(
      x = .data$dist_bin, y = .data$mean_Uro * ratio,
      group = .data$sample
    ),
    colour = program_colors[["Uro_Differentiation"]], alpha = 0.20, linewidth = 0.35
  ) +
  ggplot2::geom_line(
    data = summary_500,
    ggplot2::aes(
      x = .data$dist_bin, y = .data$grand_mean_Basal,
      colour = "Basal_pEMT"
    ),
    linewidth = 0.8
  ) +
  ggplot2::geom_line(
    data = summary_500,
    ggplot2::aes(
      x = .data$dist_bin, y = .data$grand_mean_Uro * ratio,
      colour = "Uro_Differentiation"
    ),
    linewidth = 0.8
  ) +
  ggplot2::geom_point(
    data = summary_500,
    ggplot2::aes(x = .data$dist_bin, y = .data$grand_mean_Basal),
    colour = program_colors[["Basal_pEMT"]], size = 1.2, shape = 16
  ) +
  ggplot2::geom_point(
    data = summary_500,
    ggplot2::aes(x = .data$dist_bin, y = .data$grand_mean_Uro * ratio),
    colour = program_colors[["Uro_Differentiation"]], size = 1.2, shape = 16
  ) +
  ggplot2::scale_colour_manual(values = program_colors) +
  ggplot2::scale_y_continuous(
    name = "Basal_pEMT score",
    sec.axis = ggplot2::sec_axis(~ . / ratio, name = "Uro_Differentiation score")
  ) +
  ggplot2::labs(x = "Distance to vasculature (µm)", colour = NULL) +
  theme_fig5 +
  ggplot2::theme(
    axis.title.y.left = ggplot2::element_text(colour = program_colors[["Basal_pEMT"]]),
    axis.title.y.right = ggplot2::element_text(colour = program_colors[["Uro_Differentiation"]]),
    legend.position = "top",
    legend.direction = "horizontal"
  ) +
  ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1))

ggplot2::ggsave(
  filename = file.path(fig5_output_dir, "Fig5C_MPscore_distance_to_vasculature.pdf"),
  plot = Fig5C, width = 80, height = 58, units = "mm",
  device = grDevices::cairo_pdf, family = fig5_font, bg = "white"
)
ggplot2::ggsave(
  filename = file.path(figs7_output_dir, "FigS7C_MPscore_distance_to_vasculature.pdf"),
  plot = FigS7C, width = 80, height = 58, units = "mm",
  device = grDevices::cairo_pdf, family = fig5_font, bg = "white"
)

rm(HD.combined.obj, RCTD.obj, MP.genelist)
invisible(gc())

} else {
  base::message(
    "Skipping outcome-independent Fig5B/Fig5C/FigS7C; ",
    "set RUN_FIG5_SPATIAL=1 only when a spatial redraw is required."
  )
}

###############################################################################
# Corrected malignant-cell input shared by outcome/group-dependent panels
###############################################################################

malignant_object_path <- paste0(
  "/home/data/chenwenhui/scRNAseq/GPU_IT/QC/",
  "malignant_obj_intergration_new_new.qs"
)
malignant.obj <- qs2::qs_read(malignant_object_path)
stopifnot(ncol(malignant.obj[[]]) >= 155L)

###############################################################################
# Fig5A: Basal_pEMT and Uro_Differentiation density by Outcome/Time point
###############################################################################

plot_data_A <- malignant.obj[[]] |>
  dplyr::transmute(
    cell = base::rownames(malignant.obj[[]]),
    Basal_pEMT = base::as.numeric(.data$Basal_pEMT),
    Uro_Differentiation = base::as.numeric(.data$Uro_Differentiation),
    Outcome = base::factor(
      base::as.character(.data$Outcome),
      levels = c("CR", "NR")
    ),
    Time_point = base::factor(
      base::as.character(.data$Time_point),
      levels = c("Pre", "Post")
    )
  ) |>
  dplyr::filter(
    base::is.finite(.data$Basal_pEMT),
    base::is.finite(.data$Uro_Differentiation),
    !base::is.na(.data$Outcome),
    !base::is.na(.data$Time_point)
  )

utils::write.csv(
  plot_data_A,
  file.path(fig5_output_dir, "Fig5A_MP4_MP9_density_source_data.csv"),
  row.names = FALSE
)

Fig5A <- ggplot2::ggplot(
  plot_data_A,
  ggplot2::aes(
    x = .data$Basal_pEMT,
    y = .data$Uro_Differentiation
  )
) +
  ggplot2::stat_density_2d(
    ggplot2::aes(fill = ggplot2::after_stat(level)),
    geom = "polygon"
  ) +
  ggplot2::scale_fill_viridis_c(option = "magma", name = "Cell density") +
  ggplot2::facet_grid(Outcome ~ Time_point, drop = FALSE) +
  ggplot2::coord_cartesian(
    xlim = base::range(plot_data_A$Basal_pEMT, na.rm = TRUE),
    ylim = base::range(plot_data_A$Uro_Differentiation, na.rm = TRUE)
  ) +
  ggplot2::labs(
    x = "Basal_pEMT score",
    y = "Uro_Differentiation score"
  ) +
  ggplot2::theme_bw(base_size = fig5_size, base_family = fig5_font) +
  ggplot2::theme(
    text = ggplot2::element_text(
      family = fig5_font, size = fig5_size, colour = "black"
    ),
    strip.background = ggplot2::element_rect(
      fill = "#FFF5F0", colour = "black", linewidth = 0.3
    ),
    strip.text = ggplot2::element_text(
      family = fig5_font, size = fig5_size, colour = "black"
    ),
    panel.grid = ggplot2::element_blank(),
    legend.title = ggplot2::element_text(
      family = fig5_font, size = fig5_size
    ),
    legend.text = ggplot2::element_text(
      family = fig5_font, size = fig5_size
    ),
    plot.margin = ggplot2::margin(1, 1, 1, 1, unit = "mm")
  )

ggplot2::ggsave(
  filename = file.path(fig5_output_dir, "Fig5A_MP4_MP9_density.pdf"),
  plot = Fig5A,
  width = 62,
  height = 54,
  units = "mm",
  device = grDevices::cairo_pdf,
  family = fig5_font,
  bg = "white"
)

###############################################################################
# Fig5D: pathway rank plot
###############################################################################

meta_sub <- malignant.obj[[]] |>
  dplyr::filter(.data$Time_point == "Pre", .data$MP4_Status_99 == "MP4_High")

kegg_columns <- colnames(meta_sub)[46:130]
progeny_columns <- colnames(meta_sub)[c(131:143, 154)]
mp_columns <- colnames(meta_sub)[144:153]
target_columns <- c(kegg_columns, progeny_columns, mp_columns)
stopifnot(length(target_columns) == 109L, all(target_columns %in% colnames(meta_sub)))

diff_results <- dplyr::bind_rows(lapply(target_columns, function(feature_name) {
  feature_value <- meta_sub[[feature_name]]
  if (all(is.na(feature_value)) || stats::sd(feature_value, na.rm = TRUE) == 0) {
    return(NULL)
  }

  standardized_value <- as.numeric(scale(feature_value))
  cr_value <- standardized_value[meta_sub$Outcome == "CR"]
  nr_value <- standardized_value[meta_sub$Outcome == "NR"]
  test_result <- stats::wilcox.test(nr_value, cr_value, exact = FALSE)

  category <- dplyr::case_when(
    feature_name %in% kegg_columns ~ "KEGG",
    feature_name %in% progeny_columns ~ "PROGENy",
    feature_name %in% mp_columns ~ "MP"
  )

  data.frame(
    Feature = feature_name,
    Category = category,
    EffectSize = mean(nr_value, na.rm = TRUE) - mean(cr_value, na.rm = TRUE),
    p_value = test_result$p.value,
    n_CR = sum(!is.na(cr_value)),
    n_NR = sum(!is.na(nr_value)),
    stringsAsFactors = FALSE
  )
}))

plot_data_D <- diff_results |>
  dplyr::mutate(
    p_adj = stats::p.adjust(.data$p_value, method = "BH"),
    minus_log10_p_adj = -log10(.data$p_adj),
    Clean_Name = gsub(
      " and ", " & ", gsub("_", " ", .data$Feature), fixed = TRUE
    ),
    Bias = ifelse(.data$EffectSize > 0, "NR", "CR")
  ) |>
  dplyr::arrange(dplyr::desc(.data$EffectSize)) |>
  dplyr::mutate(Rank = dplyr::row_number())

label_data_D_top <- plot_data_D |>
  dplyr::filter(.data$Rank <= 10L, .data$p_adj < 0.05) |>
  dplyr::arrange(dplyr::desc(.data$EffectSize)) |>
  dplyr::mutate(
    label_x = 20,
    label_y = seq(0.48, 0.05, length.out = dplyr::n())
  )
label_data_D_bottom <- plot_data_D |>
  dplyr::filter(.data$Rank > max(.data$Rank) - 10L, .data$p_adj < 0.05) |>
  dplyr::arrange(dplyr::desc(.data$EffectSize)) |>
  dplyr::mutate(
    label_x = 0,
    label_y = seq(-0.47, -0.98, length.out = dplyr::n())
  )
label_data_D <- dplyr::bind_rows(label_data_D_top, label_data_D_bottom)

fig5D_group_counts <- meta_sub |>
  dplyr::count(.data$Outcome, name = "n_cells")
utils::write.csv(
  plot_data_D,
  file.path(fig5_output_dir, "Fig5D_pathway_rank_source_data.csv"),
  row.names = FALSE
)
utils::write.csv(
  label_data_D,
  file.path(fig5_output_dir, "Fig5D_pathway_rank_labels.csv"),
  row.names = FALSE
)
utils::write.csv(
  fig5D_group_counts,
  file.path(fig5_output_dir, "Fig5D_group_counts.csv"),
  row.names = FALSE
)

Fig5D <- ggplot2::ggplot(
  plot_data_D,
  ggplot2::aes(x = .data$Rank, y = .data$EffectSize)
) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.35) +
  ggplot2::geom_hline(
    yintercept = c(-0.2, 0.2), linetype = "dashed",
    colour = "grey60", linewidth = 0.35
  ) +
  ggplot2::geom_point(
    ggplot2::aes(shape = .data$Category, colour = .data$Bias),
    size = 1.35, alpha = 0.70
  ) +
  ggplot2::geom_segment(
    data = label_data_D,
    ggplot2::aes(
      x = .data$Rank, y = .data$EffectSize,
      xend = .data$label_x, yend = .data$label_y,
      colour = .data$Bias
    ),
    inherit.aes = FALSE, linewidth = 0.25, alpha = 0.45,
    show.legend = FALSE
  ) +
  ggplot2::geom_text(
    data = label_data_D,
    ggplot2::aes(
      x = .data$label_x, y = .data$label_y,
      label = .data$Clean_Name, colour = .data$Bias
    ),
    inherit.aes = FALSE, family = fig5_font, size = 2.8,
    hjust = 0, vjust = 0.5, show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(
    values = outcome_colors,
    breaks = c("CR", "NR"),
    labels = c(CR = "CR higher", NR = "NR higher"),
    name = NULL
  ) +
  ggplot2::scale_shape_manual(
    values = c(KEGG = 15, MP = 16, PROGENy = 17),
    name = "Category"
  ) +
  ggplot2::labs(x = "Rank of features", y = "Effect size (NR − CR)") +
  theme_fig5 +
  ggplot2::theme(legend.position = "right") +
  ggplot2::guides(
    colour = ggplot2::guide_legend(override.aes = list(shape = 16, alpha = 0.7, size = 1.8)),
    shape = ggplot2::guide_legend(override.aes = list(colour = "black", alpha = 0.7, size = 1.8))
  )

ggplot2::ggsave(
  filename = file.path(fig5_output_dir, "Fig5D_pathway_rank.pdf"),
  plot = Fig5D, width = 100, height = 75, units = "mm",
  device = grDevices::cairo_pdf, family = fig5_font, bg = "white"
)

###############################################################################
# Fig5E: MP4-high Pre-treatment differential-expression scatter
###############################################################################

mp4_pre <- base::subset(
  malignant.obj,
  subset = MP4_Status_99 == "MP4_High" & Time_point == "Pre"
)
SeuratObject::Idents(mp4_pre) <- "Outcome"
future::plan(future::sequential)
set.seed(123)
DEG.MP4.wilcox <- Seurat::FindMarkers(
  object = mp4_pre,
  ident.1 = "NR",
  ident.2 = "CR",
  test.use = "wilcox",
  logfc.threshold = 0.1,
  min.pct = 0.05,
  random.seed = 123,
  verbose = FALSE
)

plot_data_E <- DEG.MP4.wilcox |>
  tibble::rownames_to_column("gene") |>
  dplyr::mutate(
    diff_pct = .data$pct.1 - .data$pct.2,
    group = dplyr::case_when(
      .data$avg_log2FC > 1 & .data$diff_pct > 0.1 ~ "NR_high",
      .data$avg_log2FC < -1 & .data$diff_pct < -0.1 ~ "CR_high",
      TRUE ~ "NS"
    )
  )

read_gmt_genes <- function(gmt_path) {
  gmt_object <- suppressWarnings(GSEABase::getGmt(gmt_path))
  GSEABase::geneIds(gmt_object[[1]])
}

hsa00980_genes <- read_gmt_genes(
  "/home/data/chenwenhui/scRNAseq/GPU_IT/QC/hsa00980.gmt"
)
hsa00982_genes <- read_gmt_genes(
  "/home/data/chenwenhui/scRNAseq/GPU_IT/QC/hsa00982.gmt"
)
hsa00983_genes <- read_gmt_genes(
  "/home/data/chenwenhui/scRNAseq/GPU_IT/QC/hsa00983.gmt"
)
glutathione_genes <- read_gmt_genes(
  "/home/data/chenwenhui/scRNAseq/GPU_IT/QC/GLUTATHIONE_METABOLISM.gmt"
)
hypoxia_hallmark <- read_gmt_genes(
  "/home/data/chenwenhui/scRNAseq/GPU_IT/QC/HALLMARK_HYPOXIA.v2026.1.Hs.gmt"
)
hypoxia_harris <- read_gmt_genes(
  "/home/data/chenwenhui/scRNAseq/GPU_IT/QC/HARRIS_HYPOXIA.v2026.1.Hs.gmt"
)

drug_metabolism_genes <- unique(c(
  hsa00980_genes, hsa00982_genes, hsa00983_genes, glutathione_genes
))
drug_metabolism_genes <- intersect(
  drug_metabolism_genes,
  plot_data_E$gene[plot_data_E$group == "NR_high"]
)
hypoxia_genes <- intersect(
  unique(c(hypoxia_hallmark, hypoxia_harris)),
  plot_data_E$gene[plot_data_E$group == "CR_high"]
)
preferred_hypoxia_labels <- c(
  "FGF3", "CA9", "DDIT4", "ANGPTL4", "PPP1R3C",
  "HK2", "PGK1", "PDK1", "PTGS2", "PGM1"
)
target_genes <- unique(c(
  "TNFSF10", "IGF1R", "IGF1",
  drug_metabolism_genes,
  intersect(preferred_hypoxia_labels, hypoxia_genes)
))
priority_genes_E <- c("IGF1", "IGF1R")
label_data_E <- plot_data_E |>
  dplyr::filter(.data$gene %in% target_genes) |>
  dplyr::mutate(
    group = factor(.data$group, levels = c("NR_high", "CR_high", "NS")),
    label_tier = factor(
      dplyr::if_else(.data$gene %in% priority_genes_E, "priority", "standard"),
      levels = c("standard", "priority")
    ),
    label_fontface = dplyr::if_else(
      .data$gene %in% priority_genes_E, "bold", "plain"
    ),
    label_size = dplyr::if_else(
      .data$gene %in% priority_genes_E, 3.0, 2.8
    ),
    label_background = dplyr::if_else(
      .data$gene %in% priority_genes_E, "white", NA_character_
    )
  ) |>
  dplyr::arrange(.data$label_tier)
standard_label_data_E <- label_data_E |>
  dplyr::filter(.data$label_tier == "standard")
priority_label_data_E <- label_data_E |>
  dplyr::filter(.data$label_tier == "priority")
stopifnot(setequal(priority_label_data_E$gene, priority_genes_E))

utils::write.csv(
  plot_data_E,
  file.path(fig5_output_dir, "Fig5E_MP4_high_DEG_source_data.csv"),
  row.names = FALSE
)
utils::write.csv(
  label_data_E,
  file.path(fig5_output_dir, "Fig5E_MP4_high_DEG_labels.csv"),
  row.names = FALSE
)

direction_colors <- c(
  NR_high = outcome_colors[["NR"]],
  CR_high = outcome_colors[["CR"]],
  NS = "grey85"
)
direction_and_label_colors <- c(
  direction_colors,
  standard = "black",
  priority = "#C2185B"
)

Fig5E <- ggplot2::ggplot(
  plot_data_E,
  ggplot2::aes(x = .data$diff_pct, y = .data$avg_log2FC)
) +
  ggrastr::geom_point_rast(
    ggplot2::aes(colour = .data$group),
    alpha = 0.40, size = 1.0, shape = 16,
    raster.dpi = 600
  ) +
  ggplot2::geom_vline(
    xintercept = c(-0.1, 0.1), linetype = "dashed",
    colour = "grey50", linewidth = 0.35
  ) +
  ggplot2::geom_hline(
    yintercept = c(-1, 1), linetype = "dashed",
    colour = "grey50", linewidth = 0.35
  ) +
  ggplot2::geom_vline(xintercept = 0, colour = "black", linewidth = 0.35) +
  ggplot2::geom_hline(yintercept = 0, colour = "black", linewidth = 0.35) +
  ggrastr::geom_point_rast(
    data = standard_label_data_E,
    ggplot2::aes(fill = .data$group),
    shape = 21, colour = "black", size = 2.2, stroke = 0.4,
    raster.dpi = 600
  ) +
  ggrastr::geom_point_rast(
    data = priority_label_data_E,
    shape = 21, fill = "grey70", colour = "#C2185B",
    size = 3.2, stroke = 0.8, raster.dpi = 600
  ) +
  ggrepel::geom_text_repel(
    data = label_data_E,
    ggplot2::aes(
      label = .data$gene,
      colour = .data$label_tier,
      segment.colour = .data$label_tier,
      size = .data$label_size,
      fontface = .data$label_fontface,
      bg.colour = .data$label_background
    ),
    family = fig5_font, seed = 123,
    box.padding = 0.20, point.padding = 0.20, bg.r = 0.12,
    min.segment.length = 0, segment.size = 0.25,
    max.overlaps = Inf
  ) +
  ggplot2::scale_colour_manual(
    values = direction_and_label_colors,
    aesthetics = c("colour", "segment.colour")
  ) +
  ggplot2::scale_fill_manual(values = direction_colors) +
  ggplot2::scale_size_identity() +
  ggplot2::scale_discrete_identity(aesthetics = "bg.colour") +
  ggplot2::labs(
    x = "Difference in pct (NR − CR)",
    y = "Log2 fold change (NR / CR)"
  ) +
  theme_fig5 +
  ggplot2::theme(legend.position = "none")

ggplot2::ggsave(
  filename = file.path(fig5_output_dir, "Fig5E_MP4_high_DEG.pdf"),
  plot = Fig5E, width = 86, height = 75, units = "mm",
  device = grDevices::cairo_pdf, family = fig5_font, bg = "white"
)

###############################################################################
# FigS7B: sample-level MP4-high tumour-cell proportion
###############################################################################

figs7b_data <- malignant.obj[[]] |>
  dplyr::transmute(
    sample = base::as.character(.data$sample),
    Outcome = base::factor(
      base::as.character(.data$Outcome),
      levels = c("CR", "NR")
    ),
    Time_point = base::factor(
      base::as.character(.data$Time_point),
      levels = c("Pre", "Post")
    ),
    MP4_Status_99 = base::as.character(.data$MP4_Status_99)
  ) |>
  dplyr::group_by(.data$sample, .data$Outcome, .data$Time_point) |>
  dplyr::summarise(
    Total_Cells = dplyr::n(),
    High_Count = base::sum(.data$MP4_Status_99 == "MP4_High"),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    Percentage = .data$High_Count / .data$Total_Cells * 100
  )

run_figs7b_test <- function(data, group_column, level1, level2, contrast) {
  test_data <- data.frame(
    value = base::as.numeric(data$Percentage),
    test_group = base::factor(
      base::as.character(data[[group_column]]),
      levels = c(level1, level2)
    )
  )
  test_data <- test_data[
    base::is.finite(test_data$value) & !base::is.na(test_data$test_group),
    ,
    drop = FALSE
  ]
  value1 <- test_data$value[test_data$test_group == level1]
  value2 <- test_data$value[test_data$test_group == level2]
  test_result <- coin::independence_test(
    value ~ test_group,
    data = test_data,
    distribution = "exact",
    alternative = "two.sided"
  )
  data.frame(
    contrast = contrast,
    level1 = level1,
    level2 = level2,
    n_level1 = base::length(value1),
    n_level2 = base::length(value2),
    median_level1 = stats::median(value1),
    median_level2 = stats::median(value2),
    p_value = base::as.numeric(coin::pvalue(test_result)),
    method = "Exact two-sided permutation test",
    stringsAsFactors = FALSE
  )
}

figs7b_stats <- dplyr::bind_rows(
  run_figs7b_test(
    dplyr::filter(figs7b_data, .data$Outcome == "CR"),
    "Time_point", "Pre", "Post", "CR_Pre_vs_Post"
  ),
  run_figs7b_test(
    dplyr::filter(figs7b_data, .data$Outcome == "NR"),
    "Time_point", "Pre", "Post", "NR_Pre_vs_Post"
  ),
  run_figs7b_test(
    dplyr::filter(figs7b_data, .data$Time_point == "Pre"),
    "Outcome", "CR", "NR", "Pre_CR_vs_NR"
  ),
  run_figs7b_test(
    dplyr::filter(figs7b_data, .data$Time_point == "Post"),
    "Outcome", "CR", "NR", "Post_CR_vs_NR"
  )
)

utils::write.csv(
  figs7b_data,
  file.path(figs7_output_dir, "FigS7B_MP4_high_proportion_source_data.csv"),
  row.names = FALSE
)
utils::write.csv(
  figs7b_stats,
  file.path(figs7_output_dir, "FigS7B_MP4_high_proportion_stats.csv"),
  row.names = FALSE
)

figs7b_summary <- figs7b_data |>
  dplyr::group_by(.data$Outcome, .data$Time_point) |>
  dplyr::summarise(
    mean_percentage = base::mean(.data$Percentage),
    se_percentage = stats::sd(.data$Percentage) / base::sqrt(dplyr::n()),
    .groups = "drop"
  )
figs7b_max <- base::max(figs7b_data$Percentage, na.rm = TRUE)
figs7b_step <- base::max(8, figs7b_max * 0.14)
figs7b_p <- function(contrast) {
  figs7b_stats$p_value[base::match(contrast, figs7b_stats$contrast)]
}
format_figs7b_p <- function(value) {
  if (value < 0.001) {
    base::paste0("p = ", base::formatC(value, format = "e", digits = 1))
  } else {
    base::paste0("p = ", base::sprintf("%.3f", value))
  }
}

FigS7B <- ggplot2::ggplot(
  figs7b_data,
  ggplot2::aes(x = .data$Outcome, y = .data$Percentage)
) +
  ggplot2::geom_col(
    data = figs7b_summary,
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
    data = figs7b_summary,
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
      dodge.width = 0.72,
      seed = 123
    ),
    shape = 16,
    size = 1.0,
    alpha = 0.85,
    stroke = 0
  ) +
  ggplot2::annotate(
    "segment", x = 0.82, xend = 1.18,
    y = figs7b_max + 0.35 * figs7b_step,
    yend = figs7b_max + 0.35 * figs7b_step,
    linewidth = 0.35
  ) +
  ggplot2::annotate(
    "text", x = 1, y = figs7b_max + 0.60 * figs7b_step,
    label = format_figs7b_p(figs7b_p("CR_Pre_vs_Post")),
    family = fig5_font, size = fig5_size / 2.845276
  ) +
  ggplot2::annotate(
    "segment", x = 1.82, xend = 2.18,
    y = figs7b_max + 1.05 * figs7b_step,
    yend = figs7b_max + 1.05 * figs7b_step,
    linewidth = 0.35
  ) +
  ggplot2::annotate(
    "text", x = 2, y = figs7b_max + 1.30 * figs7b_step,
    label = format_figs7b_p(figs7b_p("NR_Pre_vs_Post")),
    family = fig5_font, size = fig5_size / 2.845276
  ) +
  ggplot2::annotate(
    "segment", x = 0.82, xend = 1.82,
    y = figs7b_max + 1.85 * figs7b_step,
    yend = figs7b_max + 1.85 * figs7b_step,
    linewidth = 0.35
  ) +
  ggplot2::annotate(
    "text", x = 1.32, y = figs7b_max + 2.10 * figs7b_step,
    label = format_figs7b_p(figs7b_p("Pre_CR_vs_NR")),
    family = fig5_font, size = fig5_size / 2.845276
  ) +
  ggplot2::scale_fill_manual(
    values = c("Pre" = "#D5A04A", "Post" = "#85C1E9"),
    name = NULL
  ) +
  ggplot2::scale_colour_manual(
    values = c("Pre" = "#D5A04A", "Post" = "#85C1E9"),
    guide = "none"
  ) +
  ggplot2::coord_cartesian(
    ylim = c(0, figs7b_max + 2.45 * figs7b_step),
    clip = "off"
  ) +
  ggplot2::labs(x = NULL, y = "Cell proportion %") +
  theme_fig5 +
  ggplot2::theme(legend.position = "right")

ggplot2::ggsave(
  filename = file.path(figs7_output_dir, "FigS7B_MP4_high_proportion.pdf"),
  plot = FigS7B,
  width = 72,
  height = 56,
  units = "mm",
  device = grDevices::cairo_pdf,
  family = fig5_font,
  bg = "white"
)

###############################################################################
# FigS7D/F: stored MP and PROGENy scores in MP4-high cells
###############################################################################

figs7_group_levels <- c("Pre_NR", "Pre_CR", "Post_NR", "Post_CR")
mp4_high_metadata <- malignant.obj[[]] |>
  dplyr::filter(.data$MP4_Status_99 == "MP4_High") |>
  dplyr::mutate(
    Group = base::factor(
      base::as.character(.data$Group),
      levels = figs7_group_levels
    )
  )

mp_columns_s7 <- base::colnames(mp4_high_metadata)[144:153]
progeny_columns_s7 <- base::colnames(mp4_high_metadata)[c(131:143, 154)]

calculate_group_score_matrix <- function(metadata, features) {
  group_means <- metadata |>
    dplyr::group_by(.data$Group) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(features),
        function(value) base::mean(value, na.rm = TRUE)
      ),
      .groups = "drop"
    )
  mean_matrix <- base::t(base::as.matrix(group_means[, features, drop = FALSE]))
  base::colnames(mean_matrix) <- base::as.character(group_means$Group)
  mean_matrix <- mean_matrix[, figs7_group_levels, drop = FALSE]
  z_matrix <- base::t(base::scale(base::t(mean_matrix)))
  list(mean = mean_matrix, z = z_matrix)
}

figs7d_scores <- calculate_group_score_matrix(
  mp4_high_metadata,
  mp_columns_s7
)
figs7f_scores <- calculate_group_score_matrix(
  mp4_high_metadata,
  progeny_columns_s7
)

write_score_matrix <- function(value, filename) {
  output <- data.frame(
    feature = base::rownames(value),
    value,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  utils::write.csv(output, filename, row.names = FALSE)
}

write_score_matrix(
  figs7d_scores$mean,
  file.path(figs7_output_dir, "FigS7D_MP_program_group_means.csv")
)
write_score_matrix(
  figs7d_scores$z,
  file.path(figs7_output_dir, "FigS7D_MP_program_zscores.csv")
)
write_score_matrix(
  figs7f_scores$mean,
  file.path(figs7_output_dir, "FigS7F_PROGENy_group_means.csv")
)
write_score_matrix(
  figs7f_scores$z,
  file.path(figs7_output_dir, "FigS7F_PROGENy_zscores.csv")
)

heatmap_colours <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("#56B4E9", "white", "#D55E00")
)

grDevices::cairo_pdf(
  filename = file.path(figs7_output_dir, "FigS7D_MP_program_heatmap.pdf"),
  width = 72 / 25.4,
  height = 75 / 25.4,
  family = fig5_font
)
ComplexHeatmap::draw(ComplexHeatmap::Heatmap(
  figs7d_scores$z,
  name = "Z-score",
  col = heatmap_colours,
  cluster_columns = FALSE,
  column_order = figs7_group_levels,
  column_names_side = "bottom",
  column_names_rot = 45,
  row_names_gp = grid::gpar(fontfamily = fig5_font, fontsize = fig5_size),
  column_names_gp = grid::gpar(fontfamily = fig5_font, fontsize = fig5_size),
  heatmap_legend_param = list(
    title_gp = grid::gpar(fontfamily = fig5_font, fontsize = fig5_size),
    labels_gp = grid::gpar(fontfamily = fig5_font, fontsize = fig5_size)
  )
))
grDevices::dev.off()

grDevices::cairo_pdf(
  filename = file.path(figs7_output_dir, "FigS7F_PROGENy_heatmap.pdf"),
  width = 66 / 25.4,
  height = 75 / 25.4,
  family = fig5_font
)
ComplexHeatmap::draw(ComplexHeatmap::Heatmap(
  figs7f_scores$z,
  name = "Z-score",
  col = heatmap_colours,
  cluster_columns = FALSE,
  column_order = figs7_group_levels,
  column_names_side = "bottom",
  column_names_rot = 45,
  row_names_gp = grid::gpar(fontfamily = fig5_font, fontsize = fig5_size),
  column_names_gp = grid::gpar(fontfamily = fig5_font, fontsize = fig5_size),
  heatmap_legend_param = list(
    title_gp = grid::gpar(fontfamily = fig5_font, fontsize = fig5_size),
    labels_gp = grid::gpar(fontfamily = fig5_font, fontsize = fig5_size)
  )
))
grDevices::dev.off()

###############################################################################
# FigS7E: all Pre-treatment malignant cells, NR versus CR
###############################################################################

pre_tumour <- base::subset(malignant.obj, subset = Time_point == "Pre")
SeuratObject::Idents(pre_tumour) <- "Outcome"
future::plan(future::sequential)
set.seed(123)
DEG.pre.wilcox <- Seurat::FindMarkers(
  object = pre_tumour,
  ident.1 = "NR",
  ident.2 = "CR",
  test.use = "wilcox",
  logfc.threshold = 0.1,
  min.pct = 0.05,
  random.seed = 123,
  verbose = FALSE
)

plot_data_S7E <- DEG.pre.wilcox |>
  tibble::rownames_to_column("gene") |>
  dplyr::mutate(
    diff_pct = .data$pct.1 - .data$pct.2,
    group = dplyr::case_when(
      .data$avg_log2FC > 1 & .data$diff_pct > 0.1 ~ "NR_high",
      .data$avg_log2FC < -1 & .data$diff_pct < -0.1 ~ "CR_high",
      TRUE ~ "NS"
    )
  )

s7e_drug_genes <- base::intersect(
  base::unique(c(
    hsa00980_genes,
    hsa00982_genes,
    hsa00983_genes,
    glutathione_genes
  )),
  plot_data_S7E$gene[plot_data_S7E$group == "NR_high"]
)
s7e_hypoxia_genes <- base::intersect(
  base::unique(c(hypoxia_hallmark, hypoxia_harris)),
  plot_data_S7E$gene[plot_data_S7E$group == "CR_high"]
)
s7e_target_genes <- base::unique(c(
  "IGF1R",
  "IGF1",
  s7e_drug_genes,
  base::intersect(preferred_hypoxia_labels, s7e_hypoxia_genes)
))
label_data_S7E <- plot_data_S7E |>
  dplyr::filter(.data$gene %in% s7e_target_genes) |>
  dplyr::mutate(
    group = base::factor(
      .data$group,
      levels = c("NR_high", "CR_high", "NS")
    ),
    priority = .data$gene %in% c("IGF1", "IGF1R")
  )

utils::write.csv(
  plot_data_S7E,
  file.path(figs7_output_dir, "FigS7E_all_Pre_DEG_source_data.csv"),
  row.names = FALSE
)
utils::write.csv(
  label_data_S7E,
  file.path(figs7_output_dir, "FigS7E_all_Pre_DEG_labels.csv"),
  row.names = FALSE
)

FigS7E <- ggplot2::ggplot(
  plot_data_S7E,
  ggplot2::aes(x = .data$diff_pct, y = .data$avg_log2FC)
) +
  ggrastr::geom_point_rast(
    ggplot2::aes(colour = .data$group),
    alpha = 0.40,
    size = 1.0,
    shape = 16,
    raster.dpi = 600
  ) +
  ggplot2::geom_vline(
    xintercept = c(-0.1, 0.1),
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.35
  ) +
  ggplot2::geom_hline(
    yintercept = c(-1, 1),
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.35
  ) +
  ggplot2::geom_vline(xintercept = 0, colour = "black", linewidth = 0.35) +
  ggplot2::geom_hline(yintercept = 0, colour = "black", linewidth = 0.35) +
  ggrastr::geom_point_rast(
    data = label_data_S7E,
    ggplot2::aes(fill = .data$group),
    shape = 21,
    colour = "black",
    size = 2.2,
    stroke = 0.4,
    raster.dpi = 600
  ) +
  ggrepel::geom_text_repel(
    data = label_data_S7E,
    ggplot2::aes(
      label = .data$gene,
      fontface = base::ifelse(.data$priority, "bold", "plain")
    ),
    family = fig5_font,
    size = fig5_size / 2.845276,
    seed = 123,
    box.padding = 0.20,
    point.padding = 0.20,
    min.segment.length = 0,
    segment.size = 0.25,
    max.overlaps = Inf
  ) +
  ggplot2::scale_colour_manual(values = direction_colors) +
  ggplot2::scale_fill_manual(values = direction_colors) +
  ggplot2::labs(
    x = "Difference in pct (NR − CR)",
    y = "Log2 fold change (NR / CR)"
  ) +
  theme_fig5 +
  ggplot2::theme(legend.position = "none")

ggplot2::ggsave(
  filename = file.path(figs7_output_dir, "FigS7E_all_Pre_DEG.pdf"),
  plot = FigS7E,
  width = 102,
  height = 74,
  units = "mm",
  device = grDevices::cairo_pdf,
  family = fig5_font,
  bg = "white"
)

rm(
  malignant.obj,
  mp4_pre,
  pre_tumour,
  mp4_high_metadata
)
invisible(gc())

###############################################################################
# Fig5H: horizontal paired tumour-volume plot
###############################################################################

if (run_mouse_panel) {

tumour_volume_path <- "/home/data/chenwenhui/scRNAseq/GPU_IT/QC/mouse_experiment/Tumor_volume.csv"
tumour_raw <- utils::read.csv(
  tumour_volume_path,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8-BOM",
  check.names = FALSE
)
stopifnot(all(c("ID", "Time", "Genotype", "Group", "Tumor_volume") %in% colnames(tumour_raw)))

group_levels <- c("Control1", "Control2", "Treated1", "Treated2")
group_positions <- c(Control1 = 4, Control2 = 3, Treated1 = 2, Treated2 = 1)
group_labels <- c(
  Control1 = "Group 1\nPBS",
  Control2 = "Group 3\nPBS",
  Treated1 = "Group 2\nRC-48",
  Treated2 = "Group 4\nRC-48"
)
genotype_colors <- c(HER2_OE = "#E64B35", Igf1r_KD = "#4DBBD5")

tumour_raw$Group <- factor(tumour_raw$Group, levels = group_levels)
tumour_raw$Time <- factor(tumour_raw$Time, levels = c("Pre", "OneWeek"))

tumour_wide <- tumour_raw |>
  tidyr::pivot_wider(
    id_cols = c("ID", "Genotype", "Group"),
    names_from = "Time",
    values_from = "Tumor_volume"
  ) |>
  dplyr::mutate(
    Growth_Rate = .data$OneWeek / .data$Pre,
    Delta_Volume = .data$OneWeek - .data$Pre,
    Percent_Change = (.data$OneWeek - .data$Pre) / .data$Pre * 100
  ) |>
  dplyr::group_by(.data$Group) |>
  dplyr::arrange(.data$ID, .by_group = TRUE) |>
  dplyr::mutate(
    y_base = unname(group_positions[as.character(.data$Group)]),
    subject_offset = seq(-0.18, 0.18, length.out = dplyr::n()),
    y_position = .data$y_base + .data$subject_offset
  ) |>
  dplyr::ungroup()

# Baseline-adjusted Figure 5H analysis.
tumour_ancova <- tumour_wide |>
  dplyr::mutate(
    Treatment = factor(
      ifelse(.data$Group %in% c("Treated1", "Treated2"), "RC48", "PBS"),
      levels = c("PBS", "RC48")
    ),
    Igf1r = factor(
      ifelse(.data$Group %in% c("Control2", "Treated2"), "KD", "WT"),
      levels = c("WT", "KD")
    )
  )

ancova_fit <- stats::lm(
  log(OneWeek) ~ log(Pre) + Treatment * Igf1r,
  data = tumour_ancova,
  contrasts = list(
    Treatment = stats::contr.sum(2),
    Igf1r = stats::contr.sum(2)
  )
)

ancova_type3 <- as.data.frame(car::Anova(ancova_fit, type = 3))
ancova_global_tests <- data.frame(
  term = rownames(ancova_type3),
  sum_sq = as.numeric(ancova_type3[["Sum Sq"]]),
  df = as.numeric(ancova_type3[["Df"]]),
  f_value = as.numeric(ancova_type3[["F value"]]),
  p_value = as.numeric(ancova_type3[["Pr(>F)"]]),
  model = "log(OneWeek) ~ log(Pre) + Treatment * Igf1r",
  stringsAsFactors = FALSE
)

ancova_emmeans <- emmeans::emmeans(
  ancova_fit,
  specs = ~ Treatment * Igf1r
)
baseline_reference <- exp(mean(log(tumour_ancova$Pre)))
ancova_adjusted_means <- as.data.frame(base::summary(
  ancova_emmeans,
  infer = c(TRUE, FALSE),
  type = "response"
)) |>
  dplyr::mutate(
    Group = dplyr::case_when(
      .data$Treatment == "PBS" & .data$Igf1r == "WT" ~ "Group1",
      .data$Treatment == "RC48" & .data$Igf1r == "WT" ~ "Group2",
      .data$Treatment == "PBS" & .data$Igf1r == "KD" ~ "Group3",
      .data$Treatment == "RC48" & .data$Igf1r == "KD" ~ "Group4"
    ),
    source_group = dplyr::recode(
      .data$Group,
      Group1 = "Control1", Group2 = "Treated1",
      Group3 = "Control2", Group4 = "Treated2"
    ),
    baseline_reference_mm3 = baseline_reference
  ) |>
  dplyr::transmute(
    Group = .data$Group,
    source_group = .data$source_group,
    Treatment = as.character(.data$Treatment),
    Igf1r = as.character(.data$Igf1r),
    baseline_reference_mm3 = .data$baseline_reference_mm3,
    adjusted_geometric_mean_mm3 = .data$response,
    standard_error = .data$SE,
    df = .data$df,
    lower_95_CI = .data$lower.CL,
    upper_95_CI = .data$upper.CL
  )

primary_contrast <- list(
  Treatment_x_Igf1r = c(1, -1, -1, 1)
)
secondary_contrasts <- list(
  Group2_vs_Group1 = c(-1, 1, 0, 0),
  Group4_vs_Group3 = c(0, 0, -1, 1),
  Group4_vs_Group2 = c(0, -1, 0, 1)
)

ancova_primary <- emmeans::contrast(
  ancova_emmeans,
  method = primary_contrast,
  adjust = "none"
)
ancova_secondary_raw <- emmeans::contrast(
  ancova_emmeans,
  method = secondary_contrasts,
  adjust = "none"
)
ancova_secondary_holm <- emmeans::contrast(
  ancova_emmeans,
  method = secondary_contrasts,
  adjust = "holm"
)

ancova_primary_summary <- as.data.frame(base::summary(
  ancova_primary, infer = c(TRUE, TRUE), type = "response"
))
ancova_secondary_raw_summary <- as.data.frame(base::summary(
  ancova_secondary_raw, infer = c(TRUE, TRUE), type = "response"
))
ancova_secondary_holm_summary <- as.data.frame(base::summary(
  ancova_secondary_holm, infer = c(FALSE, TRUE), type = "response"
))

ancova_hc3_vcov <- sandwich::vcovHC(ancova_fit, type = "HC3")
ancova_emmeans_hc3 <- emmeans::emmeans(
  ancova_fit,
  specs = ~ Treatment * Igf1r,
  vcov. = ancova_hc3_vcov
)
hc3_primary <- emmeans::contrast(
  ancova_emmeans_hc3,
  method = primary_contrast,
  adjust = "none"
)
hc3_secondary_raw <- emmeans::contrast(
  ancova_emmeans_hc3,
  method = secondary_contrasts,
  adjust = "none"
)
hc3_secondary_holm <- emmeans::contrast(
  ancova_emmeans_hc3,
  method = secondary_contrasts,
  adjust = "holm"
)

hc3_primary_summary <- as.data.frame(base::summary(
  hc3_primary, infer = c(TRUE, TRUE), type = "response"
))
hc3_secondary_raw_summary <- as.data.frame(base::summary(
  hc3_secondary_raw, infer = c(TRUE, TRUE), type = "response"
))
hc3_secondary_holm_summary <- as.data.frame(base::summary(
  hc3_secondary_holm, infer = c(FALSE, TRUE), type = "response"
))

exact_comparison_groups <- list(
  Group2_vs_Group1 = c("Treated1", "Control1"),
  Group4_vs_Group3 = c("Treated2", "Control2"),
  Group4_vs_Group2 = c("Treated2", "Treated1")
)
exact_wilcoxon_raw_p <- vapply(
  names(exact_comparison_groups),
  function(contrast_name) {
    comparison_groups <- exact_comparison_groups[[contrast_name]]
    first_value <- log(tumour_ancova$Growth_Rate[
      tumour_ancova$Group == comparison_groups[1]
    ])
    second_value <- log(tumour_ancova$Growth_Rate[
      tumour_ancova$Group == comparison_groups[2]
    ])
    stats::wilcox.test(
      first_value, second_value,
      alternative = "two.sided", exact = TRUE
    )$p.value
  },
  numeric(1)
)
exact_wilcoxon_holm_p <- stats::p.adjust(
  exact_wilcoxon_raw_p,
  method = "holm"
)

ancova_primary_table <- data.frame(
  contrast = as.character(ancova_primary_summary$contrast),
  contrast_type = "primary_interaction",
  effect_definition = "(Group4 / Group3) / (Group2 / Group1)",
  effect_ratio = ancova_primary_summary$ratio,
  ancova_lower_95_CI = ancova_primary_summary$lower.CL,
  ancova_upper_95_CI = ancova_primary_summary$upper.CL,
  ancova_raw_p = ancova_primary_summary$p.value,
  ancova_holm_p = NA_real_,
  hc3_lower_95_CI = hc3_primary_summary$lower.CL,
  hc3_upper_95_CI = hc3_primary_summary$upper.CL,
  hc3_raw_p = hc3_primary_summary$p.value,
  hc3_holm_p = NA_real_,
  exact_wilcoxon_raw_p = NA_real_,
  exact_wilcoxon_holm_p = NA_real_,
  stringsAsFactors = FALSE
)

secondary_effect_definitions <- c(
  Group2_vs_Group1 = "Group2 / Group1",
  Group4_vs_Group3 = "Group4 / Group3",
  Group4_vs_Group2 = "Group4 / Group2"
)
ancova_secondary_table <- data.frame(
  contrast = as.character(ancova_secondary_raw_summary$contrast),
  contrast_type = "planned_pairwise",
  effect_definition = unname(secondary_effect_definitions[
    as.character(ancova_secondary_raw_summary$contrast)
  ]),
  effect_ratio = ancova_secondary_raw_summary$ratio,
  ancova_lower_95_CI = ancova_secondary_raw_summary$lower.CL,
  ancova_upper_95_CI = ancova_secondary_raw_summary$upper.CL,
  ancova_raw_p = ancova_secondary_raw_summary$p.value,
  ancova_holm_p = ancova_secondary_holm_summary$p.value,
  hc3_lower_95_CI = hc3_secondary_raw_summary$lower.CL,
  hc3_upper_95_CI = hc3_secondary_raw_summary$upper.CL,
  hc3_raw_p = hc3_secondary_raw_summary$p.value,
  hc3_holm_p = hc3_secondary_holm_summary$p.value,
  exact_wilcoxon_raw_p = unname(exact_wilcoxon_raw_p[
    as.character(ancova_secondary_raw_summary$contrast)
  ]),
  exact_wilcoxon_holm_p = unname(exact_wilcoxon_holm_p[
    as.character(ancova_secondary_raw_summary$contrast)
  ]),
  stringsAsFactors = FALSE
)
ancova_planned_contrasts <- dplyr::bind_rows(
  ancova_primary_table,
  ancova_secondary_table
)

ancova_shapiro <- stats::shapiro.test(stats::residuals(ancova_fit))
ancova_breusch_pagan <- lmtest::bptest(ancova_fit)
ancova_cooks_distance <- stats::cooks.distance(ancova_fit)
ancova_max_cook_index <- which.max(ancova_cooks_distance)
ancova_diagnostics <- data.frame(
  n_mice = nrow(tumour_ancova),
  n_per_group_min = min(as.integer(table(tumour_ancova$Group))),
  n_per_group_max = max(as.integer(table(tumour_ancova$Group))),
  shapiro_W = unname(ancova_shapiro$statistic),
  shapiro_p = ancova_shapiro$p.value,
  breusch_pagan_statistic = unname(ancova_breusch_pagan$statistic),
  breusch_pagan_df = unname(ancova_breusch_pagan$parameter),
  breusch_pagan_p = ancova_breusch_pagan$p.value,
  max_cooks_distance = ancova_cooks_distance[ancova_max_cook_index],
  cooks_threshold_4_over_n = 4 / nrow(tumour_ancova),
  influential_ID = tumour_ancova$ID[ancova_max_cook_index],
  influential_group = as.character(tumour_ancova$Group[ancova_max_cook_index]),
  influential_Pre = tumour_ancova$Pre[ancova_max_cook_index],
  influential_OneWeek = tumour_ancova$OneWeek[ancova_max_cook_index],
  stringsAsFactors = FALSE
)

tumour_summary <- tumour_raw |>
  dplyr::group_by(.data$Group, .data$Genotype, .data$Time) |>
  dplyr::summarise(
    Mean = mean(.data$Tumor_volume),
    SE = stats::sd(.data$Tumor_volume) / sqrt(dplyr::n()),
    n = dplyr::n(),
    .groups = "drop"
  )

tumour_mean_wide <- tumour_summary |>
  tidyr::pivot_wider(
    id_cols = c("Group", "Genotype"),
    names_from = "Time",
    values_from = c("Mean", "SE", "n")
  ) |>
  dplyr::mutate(
    y_base = unname(group_positions[as.character(.data$Group)]),
    Percent_Change = (.data$Mean_OneWeek - .data$Mean_Pre) / .data$Mean_Pre * 100,
    label_x = ifelse(
      .data$Percent_Change >= 0,
      pmax(.data$Mean_Pre, .data$Mean_OneWeek) + 10,
      pmin(.data$Mean_Pre, .data$Mean_OneWeek) - 10
    ),
    label_y = .data$y_base + 0.22,
    label_hjust = ifelse(.data$Percent_Change >= 0, 0, 1)
  )

tumour_mean_points <- tumour_summary |>
  dplyr::mutate(
    y_position = unname(group_positions[as.character(.data$Group)]),
    Time = factor(.data$Time, levels = c("Pre", "OneWeek"))
  ) |>
  dplyr::select(
    dplyr::all_of(c(
      "Group", "Genotype", "Time", "y_position", "Mean", "SE", "n"
    ))
  )

tumour_points <- tumour_wide |>
  dplyr::select(
    dplyr::all_of(c("ID", "Genotype", "Group", "y_position", "Pre", "OneWeek"))
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(c("Pre", "OneWeek")),
    names_to = "Time",
    values_to = "Tumor_volume"
  ) |>
  dplyr::mutate(Time = factor(.data$Time, levels = c("Pre", "OneWeek")))

treated1_growth <- tumour_wide$Growth_Rate[tumour_wide$Group == "Treated1"]
treated2_growth <- tumour_wide$Growth_Rate[tumour_wide$Group == "Treated2"]
growth_test <- stats::t.test(treated1_growth, treated2_growth)
tumour_statistics <- data.frame(
  comparison = "Treated1 vs Treated2 growth rate",
  method = growth_test$method,
  statistic = unname(growth_test$statistic),
  degrees_of_freedom = unname(growth_test$parameter),
  p_value = growth_test$p.value,
  stringsAsFactors = FALSE
)

utils::write.csv(
  tumour_raw,
  file.path(fig5_output_dir, "Fig5H_tumor_volume_raw.csv"),
  row.names = FALSE
)
utils::write.csv(
  tumour_wide,
  file.path(fig5_output_dir, "Fig5H_tumor_volume_wide.csv"),
  row.names = FALSE
)
utils::write.csv(
  tumour_summary,
  file.path(fig5_output_dir, "Fig5H_tumor_volume_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  tumour_points,
  file.path(fig5_output_dir, "Fig5H_tumor_volume_plotted_points.csv"),
  row.names = FALSE
)
utils::write.csv(
  tumour_mean_points,
  file.path(fig5_output_dir, "Fig5H_tumor_volume_mean_points.csv"),
  row.names = FALSE
)
utils::write.csv(
  tumour_statistics,
  file.path(fig5_output_dir, "Fig5H_statistics.csv"),
  row.names = FALSE
)
utils::write.csv(
  ancova_global_tests,
  file.path(fig5_output_dir, "Fig5H_ANCOVA_global_tests.csv"),
  row.names = FALSE
)
utils::write.csv(
  ancova_adjusted_means,
  file.path(fig5_output_dir, "Fig5H_ANCOVA_adjusted_means.csv"),
  row.names = FALSE
)
utils::write.csv(
  ancova_planned_contrasts,
  file.path(fig5_output_dir, "Fig5H_ANCOVA_planned_contrasts.csv"),
  row.names = FALSE
)
utils::write.csv(
  ancova_diagnostics,
  file.path(fig5_output_dir, "Fig5H_ANCOVA_diagnostics.csv"),
  row.names = FALSE
)

Fig5H <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = tumour_wide,
    ggplot2::aes(
      x = .data$Pre, xend = .data$OneWeek,
      y = .data$y_position, yend = .data$y_position,
      colour = .data$Genotype
    ),
    linewidth = 0.35, alpha = 0.45
  ) +
  ggplot2::geom_point(
    data = tumour_points,
    ggplot2::aes(
      x = .data$Tumor_volume, y = .data$y_position,
      colour = .data$Genotype, shape = .data$Time
    ),
    fill = "white", size = 1.25, stroke = 0.35, alpha = 0.70
  ) +
  ggplot2::geom_segment(
    data = tumour_mean_wide,
    ggplot2::aes(
      x = .data$Mean_Pre, xend = .data$Mean_OneWeek,
      y = .data$y_base, yend = .data$y_base,
      colour = .data$Genotype
    ),
    linewidth = 1.1
  ) +
  ggplot2::geom_errorbar(
    data = tumour_mean_points,
    ggplot2::aes(
      y = .data$y_position,
      xmin = .data$Mean - .data$SE,
      xmax = .data$Mean + .data$SE,
      colour = .data$Genotype
    ),
    orientation = "y", width = 0.08, linewidth = 0.45
  ) +
  ggplot2::geom_point(
    data = tumour_mean_points,
    ggplot2::aes(
      x = .data$Mean, y = .data$y_position,
      colour = .data$Genotype, shape = .data$Time
    ),
    fill = "white", size = 1.8, stroke = 0.5,
    show.legend = FALSE
  ) +
  ggplot2::geom_text(
    data = tumour_mean_wide,
    ggplot2::aes(
      x = .data$label_x, y = .data$label_y,
      label = sprintf("%+.0f%%", .data$Percent_Change),
      colour = .data$Genotype,
      hjust = .data$label_hjust
    ),
    family = fig5_font, size = 2.8, show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(
    values = genotype_colors,
    breaks = c("HER2_OE", "Igf1r_KD"),
    labels = c("HER2 OE", "Igf1r KD"),
    name = NULL
  ) +
  ggplot2::scale_shape_manual(
    values = c(Pre = 21, OneWeek = 16),
    labels = c(Pre = "14 d", OneWeek = "21 d"),
    name = NULL
  ) +
  ggplot2::scale_y_continuous(
    breaks = unname(group_positions[group_levels]),
    labels = unname(group_labels[group_levels]),
    limits = c(0.65, 4.35),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::scale_x_continuous(
    name = expression("Tumor volume (mm"^3 * ")"),
    limits = c(0, 650),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(y = NULL) +
  theme_fig5 +
  ggplot2::theme(
    legend.position = "top",
    legend.box = "vertical",
    legend.direction = "horizontal",
    legend.spacing.y = grid::unit(0, "pt"),
    legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt"),
    axis.text.y = ggplot2::element_text(lineheight = 0.9),
    plot.margin = ggplot2::margin(1, 2, 1, 1, unit = "pt")
  ) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(nrow = 2, order = 1, override.aes = list(shape = 16, alpha = 1)),
    shape = ggplot2::guide_legend(nrow = 1, order = 2)
  )

ggplot2::ggsave(
  filename = file.path(fig5_output_dir, "Fig5H_tumor_volume_horizontal.pdf"),
  plot = Fig5H, width = 42, height = 61.67, units = "mm",
  device = grDevices::cairo_pdf, family = fig5_font, bg = "white"
)

} else {
  base::message(
    "Skipping outcome-independent Fig5H; ",
    "set RUN_FIG5_MOUSE=1 only when the mouse panel needs to be redrawn."
  )
}

message("Figure 5 revised outputs: ", fig5_output_dir)
message("Figure S7 revised outputs: ", figs7_output_dir)
