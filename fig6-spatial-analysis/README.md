# Fig. 6 spatial cellular-neighbourhood analysis

This folder contains the final data-processing and plotting code used for Fig. 6 and Supplementary Figs. 8-9 in manuscript version `2026-08-03`.

The code was collected from the authoritative GPU-server analysis and excludes superseded CN analyses, the leaky/obsolete prediction workflow, debugging scripts, validation-only utilities and unused figure candidates. No patient-level data or images are included.

## Final analysis scope

- Registered three-panel mIHC sections were converted to cell-level marker measurements and phenotypes.
- The locked B0 dataset contained 72 specimen-level registered stacks from 45 patients and 23,087,820 classified cell segmentations.
- Cellular neighbourhoods (CNs) were defined from phenotype composition within 200 pixels (approximately 100 micrometres), with at least 10 neighbours per cell; `k=10` was primary and `k=11` was retained as a sensitivity analysis.
- Primary response analyses used 36 baseline specimens with at least 30,000 CN-valid cells (19 CR and 17 NR).
- CN coupling used same-CN Delaunay-connected blocks, a 200-pixel edge limit, a minimum block size of 5 cells and within-specimen CN-label permutations.
- Response models used fold-internal feature ranking, direction assignment, median imputation and scaling within LOOCV.

## Repository layout

```text
scripts/
  01_image_to_cell/                         image-to-cell measurements, thresholds and phenotypes
  02_locked_cohort_and_spatial_analysis/    B0 cohort, CNs, proportions, coupling and coherence
  03_response_prediction/                   leakage-free fold-internal response modelling
  04_figures/
    fig6/                                    final Fig. 6 panel sources
    figS8/                                   Supplementary Fig. 8 source
    figS9/                                   final Supplementary Fig. 9 quantitative panels
data/README.md                               expected inputs; data are not distributed
SOURCE_MANIFEST.tsv                         original server path and final use of every script
environment_analysis.yml                    downstream conda environment
```

## Recommended run order

If starting from registered panel-level cell data, run the scripts in `01_image_to_cell` in numerical order. Scripts 01 and 02 must remain together because the panel-merge/crosstalk script imports the registered-cell feature-extraction module.

For the locked downstream analysis:

1. `01_build_locked_b0_cohort.py`
2. `02_compute_cellular_neighborhoods.py`
3. `03_qc_cellular_neighborhoods.py`
4. `04_annotate_cellular_neighborhoods.py`
5. `05_test_celltype_and_cn_proportions.py`
6. `06_compute_cn_coupling.py`
7. `07_compute_celltype_coherence.py`
8. `01_fit_fold_internal_response_models.py`

The exact CN command used for the locked B0 analysis was:

```bash
python scripts/02_locked_cohort_and_spatial_analysis/02_compute_cellular_neighborhoods.py \
  --input_table /path/to/all_samples_cell_phenotypes_B0.parquet \
  --output_dir /path/to/neighborhood_radius200_B0 \
  --radius 200 --min_neighbors 10 --k_values 10 11 --primary_k 10 \
  --sample_col sample --cellid_col cellid --x_col spatial_x --y_col spatial_y \
  --phenotype_col celltype --workers 20 --random_state 0 \
  --max_cells_for_k_selection 100000 --prefer_output parquet
```

## Final figure-to-script map

| Panel | Script(s) retained here | Note |
|---|---|---|
| Fig. 6a | `plot_fig6a_workflow_schematic.py` | Workflow schematic source |
| Fig. 6b | `build_fig6b_c_e_panel_sources.py` | CN composition heatmap |
| Fig. 6c and 6e | `select_fig6c_e_spatial_regions.ipynb`, `export_fig6c_e_spatial_regions.py`, `build_fig6b_c_e_panel_sources.py` | Final ROIs and rasterized-scatter/vector-text exports used in Illustrator |
| Fig. 6d and 6f | `assemble_fig6_panels.py` | TLS-like presence and CN7-coupling effect-size panels; final layout was split/rearranged in Illustrator |
| Fig. 6g | `plot_fig6g_auc_ladder.py`, `plot_fig6_revised_panels.py` | Fair fold-internal AUC ladder; statistics come from `01_fit_fold_internal_response_models.py` |
| Fig. 6h | `plot_fig6h_response_roc.py` | Composite and three individual coupling ROC curves; composite alone has the bootstrap confidence band |
| Fig. S8 | `plot_figS8_processing_qc.py` | Quantitative/assembled source; microscopy images and final labels were arranged in Illustrator |
| Fig. S9a,b,f | `plot_figS9a_b_f_composition.R` | Pooled cell-type and CN composition summaries; panel f uses the final 20% wider export |
| Fig. S9d | `compute_figS9d_spatial_interactions.py`, `plot_figS9d_spatial_interactions.R` | 100-pixel radius, samples processed independently and all registered panels pooled within sample |
| Fig. S9g | `plot_figS9g_spatial_metrics.py` | Final four specimen-level panels: CD21+ B-cell coherence plus CN1-CN7, CN3-CN7 and CN4-CN7 coupling |
| Fig. S9c,e | none | Microscopy panels assembled directly from image exports in Illustrator |

The final `.ai` files are publication-layout assets and are intentionally not duplicated in this code repository.

## Paths and reproducibility

The archived analysis was run under `/data/home/chenwenhui/` and several scripts retain those absolute source/output constants for exact provenance. Before running on another system, replace the path constants at the top of each script or reproduce the same directory layout. Figure helpers in this folder were minimally rewired to locate renamed sibling modules; statistical calculations were not changed.

The crop-selection notebook has had outputs removed for GitHub. Its code and final ROI definitions are retained, including the seven final samples and coordinate windows used in Fig. 6c/e.

Use Arial at 8 pt for regenerated manuscript panels. The plotting code requests Arial and writes editable PDF/SVG text where applicable.

## Environment

Downstream analyses were run with the `scimap` conda environment (Python 3.10.18). Important recorded versions are in `environment_analysis.yml`. Image-to-cell processing used the existing `valis` environment (Python 3.10.13) because it also depends on Bio-Formats/JVM, VALIS, libvips and other system-level imaging libraries.
