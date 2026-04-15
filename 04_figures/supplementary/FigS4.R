# ==================================================
# Script: Supplementary Figure S4 plotting
# Author: Wenhui Chen
# Purpose: Generates plots used in Supplementary Figure S4.
# Input: Main Seurat or summary objects already generated in the project, plus any file paths referenced in the script.
# Output: Intermediate result files and/or figure files written to output_dir.
# Note: Review file paths before running. Keep original analysis logic unchanged.
# ==================================================
umap_df = data.frame(Embeddings(malignant.obj, "umap.harmony"),
                    malignant.obj$sample)
colnames(umap_df)=c('UMAP_1','UMAP_2','sample')
umap_df$sample = factor(umap_df$sample,levels = sort(unique(umap_df$sample)))

p1<- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = sample)) +
       geom_point(size = 0.5,shape = 16, stroke = 0) +
       scale_color_manual(values = c(colors_91)) +
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
ggsave(filename ="./Figures/FigS/FigS04B_sample_Umap.png",height = 4,width = 4,plot = p1,dpi = 600,bg="white")

#####################################
#Project MP scores onto the UMAP
new_names = c(
  "Cycle_G2M", "Surface_Trafficking", "RNA_Splicing", "Basal_pEMT", 
  "Stress_IEG", "IFN_Response", "Hypoxia", "Metabolic", 
  "Uro_Differentiation", "Ep_Inflammation"
)
pdf(file = file.path(output_dir,"Figures","FigS","FigS04B_MP_score_UMAP.pdf"),width = 8,height = 5)
FeaturePlot(malignant.obj,features = new_names,
            cols=c('lightgrey',"red"),ncol = 4,raster = T,pt.size=0.5,order = T) & NoLegend() & NoAxes()
dev.off()