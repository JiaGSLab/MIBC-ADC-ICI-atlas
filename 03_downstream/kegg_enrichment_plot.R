# ==================================================
# Script: KEGG enrichment plotting
# Author: Wenhui Chen
# Purpose: Builds a diverging bar plot for selected positively and negatively enriched pathways.
# Input: Main Seurat or summary objects already generated in the project, plus any file paths referenced in the script.
# Output: Intermediate result files and/or figure files written to output_dir.
# Note: Review file paths before running. Keep original analysis logic unchanged.
# ==================================================
# 1. Ensure required libraries are loaded
library(tidyverse)
library(forcats)

# 2. Read the input data file (assuming the same data structure as before)
kegg.enriched = rbind(MP4.enrich.HallMark[1:3,],MP4.enrich.GO[1:3,])
kegg.enriched.negative = MP9.enrich.GO[1:2,]

# 3. Data preprocessing (kept consistent with the previous version)
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

# 4. Order the terms (kept consistent with the previous version)
pos_terms_ordered <- kegg_pos %>% arrange(desc(log_p)) %>% pull(Term)
neg_terms_ordered <- kegg_neg %>% arrange(log_p) %>% pull(Term)
ordered_terms <- c(pos_terms_ordered, neg_terms_ordered)
plot_data$Term <- factor(plot_data$Term, levels = ordered_terms)

# 5. Plot using the updated color-blind-friendly palette
# Define color-blind-friendly colors
COLOR_POSITIVE <- "#D95F02" # soft magenta / reddish-purple
COLOR_NEGATIVE <- "#1B9E77" # soft dark blue-green

# Define label offset
x_label_offset <- max(abs(plot_data$Plot_Value)) * 0.1

p <- ggplot(plot_data, aes(x = Plot_Value, y = Term, fill = Regulation)) +
  
  # Draw bar plot
  geom_bar(stat = "identity") +
  
  # **Key change: use a color-blind-friendly palette**
  scale_fill_manual(values = c("Negative" = COLOR_NEGATIVE, "Positive" = COLOR_POSITIVE)) +
  
  # X-axis settings
  scale_x_continuous(
    labels = abs,
    name = expression(-log[10](Adjusted.P.value))
  ) +
  
  # Add left-side (Positive) labels
  geom_text(
    data = plot_data %>% filter(Regulation == "Positive"),
    aes(
      x = -x_label_offset,
      label = Term
    ),
    hjust = 1,
    size = 3
  ) +
  
  # Add right-side (Negative) labels
  geom_text(
    data = plot_data %>% filter(Regulation == "Negative"),
    aes(
      x = x_label_offset,
      label = Term
    ),
    hjust = 0,
    size = 3
  ) +
  
  # Remove default y-axis labels, ticks, and title
  labs(y = NULL, title = "KEGG Enrichment Diverging Bar Plot") +
  
  # Theme styling (kept consistent)
  # Theme styling
  theme_minimal() +
  theme(
    # Remove y-axis text and tick marks
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    # Adjust legend position
    legend.position = "top",
    # **Key change: remove all grid lines**
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
    # If you want to keep vertical grid lines and only remove the x-axis grid lines, use:
    # panel.grid.major.x = element_blank(),
    # panel.grid.minor.x = element_blank(),
  )

# Display plot
print(p)
ggsave(file.path(output_dir,"Figures","MP4_MP9_enrich.png"), plot = p, width = 10, height = 7, units = "in", dpi = 300)
ggsave(file.path(output_dir,"Figures","MP4_MP9_enrich.pdf"), plot = p, width = 8, height = 6, units = "in", dpi = 300)
