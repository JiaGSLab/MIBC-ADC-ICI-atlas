#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[[1]] else "/data/home/chenwenhui/mIHC/FigS9_redraw_20260731_v2"
mode <- if (length(args) >= 2) args[[2]] else "all"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_csv <- "/data/home/chenwenhui/mIHC/fig6_pm/deliverables/WP1.3_sample_proportion_long.csv"
dt <- fread(input_csv)

group_order <- c("Pre_NR", "Pre_CR", "Post_NR", "Post_CR")
group_colors <- c(
  Pre_NR = "#A52A2A", Pre_CR = "#053E7A",
  Post_NR = "#E46B4F", Post_CR = "#8B8B00"
)

colors_15 <- c(
  "#A52A2A", "#053E7A", "#E46B4F", "#8B8B00", "#B8860B",
  "#CD661D", "#8B5A2B", "#CDAD00", "#2F4F4F", "#483D8B",
  "#7A378B", "#4A7088", "#556B2F", "#008B8B", "#CD7054"
)

cell_data_order <- c(
  "CD8 T cell", "CD8 exhausted T cell", "CD4 T cell", "Treg",
  "B cell", "CD21+ B cell", "CD163- Macrophage", "CD163+ Macrophage",
  "DC", "Vascular", "Fibroblast / SMA+ stromal", "HER2+ Tumor", "PanCK+ Tumor"
)
cell_display <- c(
  "CD8 T cell", "CD8 Tex cell", "CD4 T cell", "Treg",
  "B cell", "CD21 B cell", "CD163- Macrophage", "CD163+ Macrophage",
  "DC", "Endothelial", expression(alpha*"SMA+"), "HER2 Epithelial", "PanCK Epithelial"
)

make_count_matrix <- function(x, feature_col, feature_order) {
  agg <- x[, .(count = sum(as.numeric(count))), by = c(feature_col, "Group")]
  m <- matrix(0, nrow = length(feature_order), ncol = length(group_order),
              dimnames = list(feature_order, group_order))
  for (i in seq_len(nrow(agg))) {
    r <- match(as.character(agg[[feature_col]][i]), feature_order)
    c <- match(as.character(agg$Group[i]), group_order)
    if (!is.na(r) && !is.na(c)) m[r, c] <- agg$count[i]
  }
  m
}

fisher_matrices <- function(counts) {
  nr <- nrow(counts); nc <- ncol(counts)
  OR <- P <- matrix(NA_real_, nr, nc, dimnames = dimnames(counts))
  rs <- rowSums(counts); cs <- colSums(counts); total <- sum(counts)
  jobs <- expand.grid(i = seq_len(nr), j = seq_len(nc))
  values <- parallel::mclapply(seq_len(nrow(jobs)), function(z) {
    i <- jobs$i[z]; j <- jobs$j[z]
    a <- counts[i, j]
    tab <- matrix(c(a, rs[i] - a, cs[j] - a, total - rs[i] - cs[j] + a), ncol = 2)
    ft <- suppressWarnings(fisher.test(tab))
    c(or = unname(ft$estimate), p = ft$p.value)
  }, mc.cores = min(20L, nrow(jobs)), mc.preschedule = TRUE)
  for (z in seq_len(nrow(jobs))) {
    OR[jobs$i[z], jobs$j[z]] <- values[[z]][["or"]]
    P[jobs$i[z], jobs$j[z]] <- values[[z]][["p"]]
  }
  FDR <- matrix(p.adjust(as.vector(P), method = "BH"), nr, nc, dimnames = dimnames(P))
  list(OR = OR, P = P, FDR = FDR)
}

mix_color <- function(a, b, t) {
  aa <- col2rgb(a); bb <- col2rgb(b)
  rgb((aa[1,] + (bb[1,]-aa[1,])*t)/255,
      (aa[2,] + (bb[2,]-aa[2,])*t)/255,
      (aa[3,] + (bb[3,]-aa[3,])*t)/255)
}

or_color <- function(v) {
  if (!is.finite(v)) return("#BDBDBD")
  if (v <= 1) return(mix_color("#2166AC", "#FFFFFF", max(0, v)))
  mix_color("#FFFFFF", "#B2182B", min(1, (v - 1) / 2))
}

open_pdf <- function(path, width_mm, height_mm) {
  cairo_pdf(path, width = width_mm / 25.4, height = height_mm / 25.4,
            family = "Arial", onefile = TRUE)
}

draw_panel_a <- function(counts, path) {
  prop <- sweep(counts, 2, colSums(counts), "/")
  open_pdf(path, 64, 58)
  grid.newpage()
  plot_vp <- viewport(x = 0.19, y = 0.20, width = 0.28, height = 0.75,
                      just = c("left", "bottom"), xscale = c(0.5, 4.5), yscale = c(0, 1))
  pushViewport(plot_vp)
  bottom <- rep(0, 4)
  for (idx in rev(seq_len(nrow(prop)))) {
    h <- prop[idx,]
    for (j in seq_len(4)) {
      grid.rect(x = unit(j, "native"), y = unit(bottom[j], "native"),
                width = unit(0.68, "native"), height = unit(h[j], "native"),
                just = c("center", "bottom"),
                gp = gpar(fill = colors_15[idx], col = "black", lwd = 0.45))
    }
    bottom <- bottom + h
  }
  grid.lines(x = unit(c(0.5, 4.5), "native"), y = unit(c(0, 0), "native"), gp = gpar(lwd = 0.65))
  grid.lines(x = unit(c(0.5, 0.5), "native"), y = unit(c(0, 1), "native"), gp = gpar(lwd = 0.65))
  for (v in seq(0, 1, 0.25)) {
    grid.lines(x = unit(c(0.47, 0.5), "native"), y = unit(c(v, v), "native"), gp = gpar(lwd = 0.5))
    grid.text(sprintf("%d%%", round(100*v)), x = unit(-0.04, "npc"), y = unit(v, "native"),
              just = "right", gp = gpar(fontfamily = "Arial", fontsize = 8))
  }
  for (j in seq_len(4)) {
    grid.text(group_order[j], x = unit(j, "native"), y = unit(-0.055, "npc"),
              rot = 45, just = "right", gp = gpar(fontfamily = "Arial", fontsize = 8))
  }
  grid.text("Cell Composition (%)", x = unit(-0.48, "npc"), y = unit(0.5, "npc"),
            rot = 90, gp = gpar(fontfamily = "Arial", fontsize = 8))
  popViewport()

  grid.text("Cell Type", x = unit(0.545, "npc"), y = unit(0.94, "npc"),
            just = "left", gp = gpar(fontfamily = "Arial", fontsize = 8))
  y0 <- 0.885; dy <- 0.065
  for (i in seq_along(cell_display)) {
    yy <- y0 - (i-1)*dy
    grid.rect(x = unit(0.55, "npc"), y = unit(yy, "npc"), width = unit(0.032, "npc"),
              height = unit(0.043, "npc"), gp = gpar(fill = colors_15[i], col = "black", lwd = 0.35))
    grid.text(cell_display[[i]], x = unit(0.575, "npc"), y = unit(yy, "npc"),
              just = "left", gp = gpar(fontfamily = "Arial", fontsize = 8))
  }
  dev.off()
}

draw_distribution <- function(counts, labels, path, width_mm, height_mm) {
  stats <- fisher_matrices(counts)
  frac <- sweep(counts, 1, rowSums(counts), "/")
  nfeat <- nrow(counts)
  open_pdf(path, width_mm, height_mm)
  grid.newpage()

  x0 <- if (nfeat > 10) 0.12 else 0.18
  x1 <- if (nfeat > 10) 0.88 else 0.82
  top0 <- 0.68; top1 <- 0.96
  hm0 <- 0.30; hm1 <- 0.66

  top_vp <- viewport(x = x0, y = top0, width = x1-x0, height = top1-top0,
                     just = c("left", "bottom"), xscale = c(0.5, nfeat+0.5), yscale = c(0,1))
  pushViewport(top_vp)
  bottom <- rep(0, nfeat)
  for (g in rev(group_order)) {
    vals <- frac[, g]
    for (j in seq_len(nfeat)) {
      grid.rect(x = unit(j, "native"), y = unit(bottom[j], "native"),
                width = unit(0.98, "native"), height = unit(vals[j], "native"),
                just = c("center", "bottom"), gp = gpar(fill = group_colors[g], col = "white", lwd = 0.35))
    }
    bottom <- bottom + vals
  }
  grid.lines(x = unit(c(0.5,0.5), "native"), y = unit(c(0,1), "native"), gp = gpar(lwd=0.55))
  for (v in c(0,0.5,1)) {
    grid.lines(x = unit(c(0.47,0.5), "native"), y = unit(c(v,v), "native"), gp = gpar(lwd=0.45))
    grid.text(format(v, trim=TRUE), x=unit(-0.035,"npc"), y=unit(v,"native"), just="right",
              gp=gpar(fontfamily="Arial", fontsize=8))
  }
  grid.text("Fraction", x=unit(-0.12,"npc"), y=unit(0.5,"npc"), rot=90,
            gp=gpar(fontfamily="Arial", fontsize=8))
  popViewport()

  hm_vp <- viewport(x = x0, y = hm0, width = x1-x0, height = hm1-hm0,
                    just = c("left", "bottom"), xscale = c(0.5,nfeat+0.5), yscale = c(0.5,4.5))
  pushViewport(hm_vp)
  for (i in seq_len(4)) for (j in seq_len(nfeat)) {
    yy <- 5-i
    val <- stats$OR[j,i]
    grid.rect(x=unit(j,"native"), y=unit(yy,"native"), width=unit(1,"native"), height=unit(1,"native"),
              gp=gpar(fill=or_color(val), col="white", lwd=0.45))
    if (is.finite(stats$FDR[j,i]) && stats$FDR[j,i] < 0.05) {
      txt_col <- if (val > 2 || val < 0.3) "white" else "black"
      grid.text(sprintf("%.2f",val), x=unit(j,"native"), y=unit(yy,"native"),
                gp=gpar(fontfamily="Arial", fontsize=7.2, col=txt_col))
    }
  }
  for (i in seq_len(4)) {
    grid.text(group_order[i], x=unit(-0.018,"npc"), y=unit(5-i,"native"), just="right",
              gp=gpar(fontfamily="Arial", fontsize=8, col=group_colors[group_order[i]]))
  }
  for (j in seq_len(nfeat)) {
    grid.text(labels[[j]], x=unit(j,"native"), y=unit(-0.04,"npc"), rot=45, just="right",
              gp=gpar(fontfamily="Arial", fontsize=8))
  }
  popViewport()

  lx <- if (nfeat > 10) 0.925 else 0.90
  grid.text("Ro/e", x=unit(lx,"npc"), y=unit(0.63,"npc"),
            gp=gpar(fontfamily="Arial", fontsize=8))
  for (k in 0:59) {
    val <- 3*k/59
    grid.rect(x=unit(lx,"npc"), y=unit(0.43 + 0.17*k/59,"npc"),
              width=unit(0.035,"npc"), height=unit(0.17/59+0.001,"npc"),
              gp=gpar(fill=or_color(val), col=NA))
  }
  for (val in 0:3) {
    yy <- 0.43 + 0.17*val/3
    grid.text(if (val==3) ">=3" else as.character(val), x=unit(lx+0.03,"npc"), y=unit(yy,"npc"),
              just="left", gp=gpar(fontfamily="Arial", fontsize=7.5))
  }
  dev.off()
}

cell_dt <- dt[feature_domain == "celltype"]
cell_counts <- make_count_matrix(cell_dt, "feature_name", cell_data_order)
cn_order <- paste0("CN", 0:9)
cn_dt <- dt[feature_domain == "CN" & as.character(k) == "10"]
cn_counts <- make_count_matrix(cn_dt, "feature_id", cn_order)

if (mode %in% c("all", "a"))
  draw_panel_a(cell_counts, file.path(out_dir, "FigS9a_B0_server_redraw_v2.pdf"))
if (mode %in% c("all", "b"))
  draw_distribution(cell_counts, cell_display,
                    file.path(out_dir, "FigS9b_B0_server_redraw_v2.pdf"), 115, 58)
if (mode %in% c("all", "f"))
  draw_distribution(cn_counts, as.list(cn_order),
                    file.path(out_dir, "FigS9f_B0_server_redraw_width120.pdf"), 95, 62)

cat("Created:\n",
    file.path(out_dir, "FigS9a_B0_server_redraw_v2.pdf"), "\n",
    file.path(out_dir, "FigS9b_B0_server_redraw_v2.pdf"), "\n",
    file.path(out_dir, "FigS9f_B0_server_redraw_v2.pdf"), "\n", sep="")
