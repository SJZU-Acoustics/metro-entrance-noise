#!/usr/bin/env python3
"""Reproduce the manuscript's display items from the Mendeley Data workbook.

    python run_all.py            # everything the shipped inputs allow
    python run_all.py --list     # show the steps and what each needs

Reads   data/Shenyang_metro_entrance_acoustic_environment_data.xlsx   (Mendeley deposit)
        data/osm_roads_shenyang.csv                                  (shipped here)
        data/validation/*.csv                                        (optional, see data/validation/README.md)
Writes  intermediate/clean/       the four workbook tables the analysis reads
        intermediate/analysis/    every analysis module's result tables
        output/figures/           Figures 1-5 and, in output/figures/si/, Supplementary Figures S1-S5
        output/tables/si/         Supplementary Tables S1-S10 (LaTeX fragments)
        output/figure_data/       the per-figure CSVs each figure is drawn from

Every step runs from the repository root. Python steps use the interpreter
running this script; R steps use Rscript on the PATH.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
VALIDATION = ROOT / "data" / "validation"
CLIPS = VALIDATION / "sound_event_validation_clips.csv"
PAIRS = VALIDATION / "collocation_pairs.csv"

# (step, command, needs) -- needs: None, "clips" or "pairs"
STEPS = [
    ("load_data", ["py", "code/load_data.py"], None),
    # analysis modules (results feed the figures, tables and the numbers quoted in the text)
    ("03_descriptive_plots", ["py", "code/analysis/03_descriptive_plots.py"], None),
    ("04_group_differences", ["py", "code/analysis/04_group_differences.py"], None),
    ("06_built_env_multiscale", ["py", "code/analysis/06_built_env_multiscale.py", "--boot", "1000", "--boot-partial", "1000"], None),
    ("08_sound_source_composition", ["py", "code/analysis/08_sound_source_composition.py"], None),
    ("12_spatial_robustness_models", ["py", "code/analysis/12_spatial_robustness_models.py"], None),
    ("28_source_spatial_confounding", ["py", "code/analysis/28_source_spatial_confounding.py"], None),
    ("33_sound_event_validation", ["py", "code/analysis/33_sound_event_validation.py"], "clips"),
    ("18_outcome_structure_pca", ["R", "code/analysis/18_outcome_structure_pca.R"], None),
    ("19_temporal_synchrony", ["R", "code/analysis/19_temporal_synchrony.R"], None),
    ("21_domain_cv_models", ["R", "code/analysis/21_domain_cv_models.R"], None),
    ("22_measurement_timing_robustness", ["R", "code/analysis/22_measurement_timing_robustness.R"], None),
    ("24_acoustic_diversity", ["R", "code/analysis/24_acoustic_diversity.R"], None),
    ("25_ndvi_greenness_gradient", ["R", "code/analysis/25_ndvi_greenness_gradient.R"], None),
    ("27_source_assoc_influence_robustness", ["R", "code/analysis/27_source_assoc_influence_robustness.R"], None),
    ("29_beyond_loudness_partial", ["R", "code/analysis/29_beyond_loudness_partial.R"], None),
    ("32_collocation_check", ["R", "code/analysis/32_collocation_check.R"], "pairs"),
    # display items
    ("fig1_study_map", ["R", "code/display/fig1_study_map.R"], None),
    ("fig2_temporal_burden", ["R", "code/display/fig2_temporal_burden.R"], None),
    ("fig3_synchrony", ["R", "code/display/fig3_synchrony.R"], None),
    ("fig4_source_mechanism", ["R", "code/display/fig4_source_mechanism.R"], None),
    ("fig5_spatial_greenness", ["R", "code/display/fig5_spatial_greenness.R"], None),
    ("si1_outcome_structure", ["R", "code/display/si1_outcome_structure.R"], None),
    ("si2_timing_robustness", ["R", "code/display/si2_timing_robustness.R"], None),
    ("si3_builtenv_multiscale", ["R", "code/display/si3_builtenv_multiscale.R"], None),
    ("si4_lisa_maps", ["R", "code/display/si4_lisa_maps.R"], None),
    ("si5_ndvi_gradient", ["R", "code/display/si5_ndvi_gradient.R"], None),
    ("si_tables", ["R", "code/display/si_tables.R"], None),
]

NEEDS = {
    "clips": (CLIPS, "Supplementary Table S2 (classifier validation)"),
    "pairs": (PAIRS, "the recorder collocation figures quoted in Methods"),
}


def resolve(cmd: list[str]) -> list[str]:
    kind, *rest = cmd
    if kind == "py":
        return [sys.executable, *rest]
    rscript = shutil.which("Rscript")
    if rscript is None:
        sys.exit("Rscript was not found on the PATH; install R 4.5+ and retry.")
    return [rscript, *rest]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true", help="list the steps and exit")
    args = ap.parse_args()
    if args.list:
        for name, cmd, needs in STEPS:
            note = "" if needs is None else f"   [needs {NEEDS[needs][0].relative_to(ROOT)}]"
            print(f"{name:38s} {' '.join(cmd[1:])}{note}")
        return
    start = time.time()
    for name, cmd, needs in STEPS:
        if needs is not None and not NEEDS[needs][0].exists():
            print(f"\n=== {name} === skipped: {NEEDS[needs][0].relative_to(ROOT)} not supplied, so {NEEDS[needs][1]} is not reproduced (see data/validation/README.md)")
            continue
        print(f"\n=== {name} ===", flush=True)
        t0 = time.time()
        subprocess.run(resolve(cmd), cwd=ROOT, check=True)
        print(f"    done in {time.time() - t0:.1f} s", flush=True)
    print(f"\nComplete in {(time.time() - start) / 60:.1f} min. Figures are in output/figures, table fragments in output/tables/si.")


if __name__ == "__main__":
    main()
