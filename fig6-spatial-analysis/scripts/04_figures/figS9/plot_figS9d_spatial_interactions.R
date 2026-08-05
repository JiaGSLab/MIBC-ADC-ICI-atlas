#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
input_csv <- args[[1]]
output_pdf <- args[[2]]
d <- fread(input_csv)

data_order <- c(
  "B cell", "CD163+ Macrophage", "CD163- Macrophage", "CD21+ B cell",
  "CD4 T cell", "CD8 T cell", "CD8 exhausted T cell", "DC", "Vascular",
  "Fibroblast / SMA+ stromal", "HER2+ Tumor", "PanCK+ Tumor", "Treg"
)
display_order <- c(
  "B cell", "CD163+ Macrophage", "CD163- Macrophage", "CD21 B cell",
  "CD4 T cell", "CD8 T cell", "CD8 Tex cell", "DC", "Endothelial",
  "αSMA+", "HER2 Epithelial", "PanCK Epithelial", "Treg"
)
group_order <- c("Pre_CR", "Pre_NR")
group_colors <- c(Pre_CR = "#053E7A", Pre_NR = "#A52A2A")

mix_color <- function(a, b, t) {
  aa <- col2rgb(a); bb <- col2rgb(b)
  rgb((aa[1,] + (bb[1,]-aa[1,])*t)/255,
      (aa[2,] + (bb[2,]-aa[2,])*t)/255,
      (aa[3,] + (bb[3,]-aa[3,])*t)/255)
}

interaction_color <- function(v) {
  if (!is.finite(v)) return("#BDBDBD")
  if (v <= 0) return(mix_color("#2166AC", "#FFFFFF", max(0, v + 1)))
  mix_color("#FFFFFF", "#B2182B", min(1, v))
}

make_matrix <- function(group) {
  z <- d[Group == group]
  m <- matrix(NA_real_, length(data_order), length(data_order),
              dimnames = list(data_order, data_order))
  for (i in seq_len(nrow(z))) {
    r <- match(z$phenotype[i], data_order)
    c <- match(z$neighbour_phenotype[i], data_order)
    if (!is.na(r) && !is.na(c)) m[r,c] <- z$masked_binary_mean[i]
  }
  m
}

cairo_pdf(output_pdf, width = 100/25.4, height = 62/25.4,
          family = "Arial", onefile = TRUE)
grid.newpage()

x_starts <- c(0.02, 0.34)
x_width <- 0.29
y0 <- 0.47
y1 <- 0.91
n <- length(data_order)

for (gidx in seq_along(group_order)) {
  g <- group_order[gidx]
  m <- make_matrix(g)
  vp <- viewport(x = x_starts[gidx], y = y0, width = x_width, height = y1-y0,
                 just = c("left", "bottom"), xscale = c(0.5,n+0.5), yscale = c(0.5,n+0.5))
  pushViewport(vp)
  for (i in seq_len(n)) for (j in seq_len(n)) {
    yy <- n - i + 1
    grid.rect(x=unit(j,"native"), y=unit(yy,"native"),
              width=unit(1,"native"), height=unit(1,"native"),
              gp=gpar(fill=interaction_color(m[i,j]), col="white", lwd=0.35))
  }
  for (j in seq_len(n)) {
    grid.text(display_order[j], x=unit(j,"native"), y=unit(-0.025,"npc"),
              rot=90, just="right", gp=gpar(fontfamily="Arial", fontsize=8))
  }
  if (gidx == 2) {
    for (i in seq_len(n)) {
      yy <- n-i+1
      grid.text(display_order[i], x=unit(1.018,"npc"), y=unit(yy,"native"),
                just="left", gp=gpar(fontfamily="Arial", fontsize=8))
    }
  }
  popViewport()
  grid.text(g, x=unit(x_starts[gidx] + x_width/2,"npc"), y=unit(0.955,"npc"),
            gp=gpar(fontfamily="Arial", fontsize=8, col=group_colors[g]))
}

grid.text("neighbour phenotype", x=unit(0.325,"npc"), y=unit(0.015,"npc"),
          gp=gpar(fontfamily="Arial", fontsize=8))
grid.text("phenotype", x=unit(0.90,"npc"), y=unit(0.66,"npc"), rot=90,
          gp=gpar(fontfamily="Arial", fontsize=8))

for (k in 0:59) {
  val <- -1 + 2*k/59
  grid.rect(x=unit(0.945,"npc"), y=unit(0.72 + 0.19*k/59,"npc"),
            width=unit(0.025,"npc"), height=unit(0.19/59+0.001,"npc"),
            gp=gpar(fill=interaction_color(val), col=NA))
}
for (val in c(-1,0,1)) {
  yy <- 0.72 + 0.19*(val+1)/2
  grid.text(as.character(val), x=unit(0.962,"npc"), y=unit(yy,"npc"), just="left",
            gp=gpar(fontfamily="Arial", fontsize=8))
}

dev.off()
cat(output_pdf, "\n")
