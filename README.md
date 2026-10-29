Code for data analysis for the paper "Information needs for noise priority assessment at urban metro entrances".

## Requirements

R 4.5 or later, with tidyverse, patchwork, scales, digest and ragg. Python 3.11 or later, with openpyxl, numpy, pandas, scipy, scikit-learn, statsmodels, libpysal, esda and spreg. The figures use Helvetica; install that font or supply a compatible substitute if unavailable.

## Data

Place `Shenyang_metro_entrance_integrated_data.xlsx` from the Mendeley dataset [10.17632/7rwt35tchz](https://doi.org/10.17632/7rwt35tchz) in `data/`. This pipeline requires the integrated workbook, including the validation and collocation sheets. The earlier workbook does not contain the complete input set.

The workbook contains analysis-ready levels, typical-day activity indices, hourly event labels, contextual indicators, recording dates, anonymous validation labels and collocation pairs. Raw audio and original location-service exports are not distributed. This pipeline reproduces analyses from the supplied hourly winners; it does not rerun the audio classifier or reconstruct a historical inference batch.

The background road cache in `data/osm_roads_shenyang.csv` is from OpenStreetMap (© OpenStreetMap contributors, ODbL). It is used only to draw map backgrounds. Code is MIT licensed; deposited data are CC BY 4.0; the road cache retains its separate ODbL terms.

## File structure

- `run_all.R`: one-command reproduction.
- `code/load_data.py`: extracts numeric workbook cells into runtime CSV inputs.
- `code/analysis/`: selection, source-profile, matching, sensitivity and supporting model calculations.
- `code/shared/`: shared definitions and plot style.
- `code/build_displays.R` and `code/support_tables.R`: manuscript display generators.
- `data/`: downloaded workbook and the supplied background road cache.
- `intermediate/`: generated analysis tables; safe to remove and rebuild.
- `output/figures/`: seven main figures and two supplementary figures.
- `output/tables/`: two main and thirteen supplementary table bodies in LaTeX.
- `verification/`: input-read trace and reproduction evidence generated during checking.

## Usage

From the repository directory:

```sh
Rscript run_all.R
```

The command uses `python3` on PATH. To select another Python environment, set `P11_PYTHON` to its interpreter path. No external drive, network request or private project file is needed once the workbook has been downloaded. Seeds are specified in the analysis modules. Intermediate results include quantities used in the manuscript prose as well as those used in the figures and tables.

Hourly records remain linked by entrance identifier. Event categories are grouped only after selecting the winning event. The conventional day/night boundaries and the alternative Chinese boundaries are specified in the shared functions. Persistent-error scenarios are conditional assumptions, not estimates of the probability of future selection.

## Verification

A clean run from the integrated workbook reproduces all nine figure PNGs and all fifteen table bodies byte-for-byte against the manuscript build. Thirteen independent numerical checks pass. Clip-level accuracy and kappa match the retained calculations, and collocation bias and RMSE reproduce 0.41 and 0.45 dB. Workbook values are stored to 15 significant figures; no displayed result changes after the round trip. Very small numerical zeros are formatted as zero.
