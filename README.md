Python and R code for reproducing the statistical analyses, figures and tables for the manuscript "Sound-source composition explains spatially clustered noise at citywide metro entrances".

## Requirements

- Python 3.11+ (developed and verified on 3.13) and R 4.5+ (verified on R 4.5.3)

- Python packages: `numpy`, `pandas`, `scipy`, `statsmodels`, `scikit-learn`, `matplotlib`, `seaborn`, `openpyxl`, `libpysal`, `esda`, `spreg`

- R packages: `tidyverse`, `patchwork`, `scales`, `ggrepel`, `ragg`

- Install:

  ```bash
  pip install numpy pandas scipy statsmodels scikit-learn matplotlib seaborn openpyxl libpysal esda spreg
  ```

  ```r
  install.packages(c("tidyverse", "patchwork", "scales", "ggrepel", "ragg"))
  ```

- No non-standard hardware is required. A complete run takes about two minutes on a normal desktop; the bootstrap correlations of the built-environment screen and the repeated cross-validation account for most of it.

## Data

The analysis reads the Mendeley Data workbook `Shenyang_metro_entrance_acoustic_environment_data.xlsx` — the station-level table, the station-hour panel, the derived acoustic indicators and the land-use name map for the 111 metro-station entrance sites in Shenyang. Download it from Mendeley Data, **DOI [10.17632/7rwt35tchz](https://doi.org/10.17632/7rwt35tchz)** (CC BY 4.0), and place it in `data/` (see `data/README.md`).

One further input ships with the repository: `data/osm_roads_shenyang.csv`, the OpenStreetMap major-road layer drawn behind the point maps of Figures 1 and 5 (© OpenStreetMap contributors, ODbL).

Two small measurement-validation tables are not part of the deposit and are not distributed here (`data/validation/README.md`). Without them `run_all.py` skips Supplementary Table S2 and the recorder-collocation figures quoted in Methods, with a message, and reproduces everything else.

## File structure

- `run_all.py` — master script: extracts the workbook tables, runs every analysis module, then builds every manuscript display item. `python run_all.py --list` prints the steps and what each needs.
- `code/load_data.py` — single data entry: writes the four workbook sheets the analysis reads to `intermediate/clean/` as CSV, cell for cell as the workbook stores them, so every module re-infers column types exactly as the working pipeline did.
- `code/analysis/` — the analysis modules, each writing its result tables to `intermediate/analysis/`: descriptive statistics and the hourly curve (03), line and station-type group differences (04), the multiscale built-environment screen with bootstrap and partial correlations (06), sound-source composition and its correlations with level (08), global and local Moran's I and the spatial models on built-environment predictors (12), the acoustic-metric correlation and PCA structure (18), within-station heat–noise synchrony (19), the repeated five-fold cross-validation of the domain models (21), measurement-timing robustness (22), the station-level source metrics (24), the NDVI greenness gradient (25), influence and method robustness of the source associations (27), the source-composition signal under spatial adjustment (28), the partial correlations of natural-source share with greenness (29), the recorder collocation check (32, optional) and the classifier validation against the human-labelled clips (33, optional).
- `code/display/` — `style.R` (shared plot style and paths), `fig1_study_map.R` … `fig5_spatial_greenness.R` (Figures 1–5), `si1_outcome_structure.R` … `si5_ndvi_gradient.R` (Supplementary Figures S1–S5) and `si_tables.R` (Supplementary Tables S1–S10 as LaTeX fragments).
- `code/fetch_osm.py` — refreshes the OpenStreetMap road cache in `data/`; not part of the run.

## Usage

From the repository root:

```bash
python run_all.py
```

Outputs are written to:

- `output/figures/` — Figures 1–5 (PNG, 600 dpi) and, in `output/figures/si/`, Supplementary Figures S1–S5
- `output/tables/si/` — Supplementary Tables S1–S10 (LaTeX fragments)
- `output/figure_data/` — the per-figure CSVs each figure is drawn from
- `intermediate/clean/` — the four workbook tables; `intermediate/analysis/` — every module's full result tables

`Rscript` must be on the PATH. `intermediate/` and `output/` are produced at run time and are safe to delete.

## Verification

Run against the deposited workbook, this pipeline reproduces all ten Supplementary Table fragments byte-identically, and seven of the ten figures byte-identically. Figures 4 and 5 and Supplementary Figure S5 differ from the manuscript renders in at most 735 of about 13 million pixels, by at most 2 of 255 in intensity — anti-aliasing at point edges, with no mark, label or digit differing.

Workbook cells are capped at 15 significant digits, so a few inputs differ from the working pipeline's values in their sixteenth digit. Every quantity that reaches a figure, a table or the text reproduces exactly at its reported precision. The residual effect is confined to result-table columns the manuscript does not report: the Shannon diversity index of the sound-event labels carries rank ties that flip under the perturbation, moving its correlations by up to 0.002 and one bootstrap bound by 0.001, and the maximum-likelihood spatial models' p-values move by up to 1e-5 relative through the optimiser's tolerance.

Every random element is seeded: the bootstrap and cross-validation draws in the modules, the permutation tests of the global and local Moran's I, and the render-time jitter and label placement of Figures 3 and 1. Repeat runs are byte-identical.

## Notes

- Every step runs from the repository root. The Python steps use the interpreter that runs `run_all.py`; the R steps use `Rscript`. Each module reads `intermediate/clean/` (and, for 29, the output of 24) and writes only its own result tables, so modules can be rerun individually after `code/load_data.py` has run.
- The built-environment screen (06) runs with 1,000 bootstrap resamples for both the zero-order and the partial correlations, as reported in the manuscript.
- The domain models (21) are compared by five-fold cross-validation repeated 100 times on fixed folds; the paired win counts in Table 2 and Supplementary Table S5 are the number of repeats in which the source-only model had the lower error.
- Spatial weights are k = 8 nearest neighbours throughout (12, 28), with 999 permutations for the Moran tests.

## License

Code in this repository is released under the MIT License (see `LICENSE`). The input data are archived separately under CC BY 4.0 at Mendeley Data (DOI as above).
