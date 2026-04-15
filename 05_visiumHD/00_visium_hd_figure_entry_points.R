
# ==================================================
# Script: Visium HD figure entry points
# Author: Wenhui Chen
# Purpose: Lists figure scripts that use Visium HD derived objects or RCTD annotations so readers can trace preprocessing dependencies.
# Input: Annotated Visium HD objects and figure-specific intermediate files referenced in the project.
# Output: This note for repository navigation.
# Note: Review file paths before running. Keep original analysis logic unchanged.
# ==================================================
# Main figure scripts with Visium HD related content:
# - scripts/04_figures/main/figure4.R
# - scripts/04_figures/main/figure5.R
# Supplementary figure scripts with Visium HD related content:
# - scripts/04_figures/supplementary/figureS4.R
# - scripts/04_figures/supplementary/figureS5.R
# Recommended order:
# 1. scripts/05_visiumHD/01_build_rctd_reference_from_scrna.R
# 2. scripts/05_visiumHD/02_run_rctd_slide01.R
# 3. scripts/05_visiumHD/03_run_rctd_slide02.R
# 4. scripts/05_visiumHD/04_run_rctd_slide03.R
# 5. scripts/05_visiumHD/05_attach_rctd_to_visium_objects_and_merge.R
# 6. scripts/05_visiumHD/06_visium_hd_spatial_plots_and_celltype_highlights.R
# 7. figure scripts as needed
