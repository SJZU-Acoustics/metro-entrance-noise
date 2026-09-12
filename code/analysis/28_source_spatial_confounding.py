"""28_source_spatial_confounding.py

Robustness analysis: does the source-composition signal survive spatial
adjustment?

The headline finding is that traffic share raises and natural share lowers
station noise. Acoustic outcomes are spatially clustered (global Moran's I
~0.20-0.25 in analysis 12), and the spatial models in analysis 12 only used
built-environment predictors. This script re-estimates the source -> noise
association while accounting for spatial structure, using:

  1. OLS with source predictors (statsmodels), with residual Moran's I.
  2. A spatial trend-surface OLS that adds a quadratic lng/lat surface.
  3. Spatial-lag (SAR / ML_Lag) and spatial-error (SEM / ML_Error) models
     with the same source predictors, using KNN(k=8) weights to match
     analysis 12.

If the traffic and natural coefficients stay in the same direction and remain
significant across these specifications, the source signal is not merely a
by-product of spatial clustering.
"""

from __future__ import annotations

import os
from pathlib import Path

import tempfile
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "mplconfig"))

import numpy as np
import pandas as pd
import statsmodels.api as sm

TOPIC_ROOT = Path(__file__).resolve().parents[2]  # repository root
DATA_CLEAN_DIR = TOPIC_ROOT / "intermediate" / "clean"
OUTPUT_DIR = TOPIC_ROOT / "intermediate" / "analysis"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

K_NEIGHBORS = 8
PERM_SEED = 20260206  # fixes the permutation draws of the residual Moran tests
SOURCE_COLS = ["sound_traffic_share", "sound_natural_share"]
OUTCOMES = {
    "laeq_24h_leq_db_a": "LAeq_24h",
    "lday_leq_db_a": "Lday",
    "lnight_leq_db_a": "Lnight",
}


def build_knn_weights(coords: np.ndarray, k: int):
    from libpysal.weights import KNN

    w = KNN.from_array(coords, k=k)
    w.transform = "R"
    return w


def residual_moran(resid: np.ndarray, w, permutations: int = 999):
    from esda.moran import Moran

    np.random.seed(PERM_SEED)
    m = Moran(resid, w, permutations=permutations)
    return float(m.I), float(m.p_sim)


def extract_spreg(model, model_name: str, outcome_label: str, spatial_term: str) -> pd.DataFrame:
    """Defensive coefficient extraction from a spreg model."""
    betas = np.asarray(model.betas).flatten()
    try:
        se = np.asarray(model.std_err).flatten()
    except Exception:
        se = np.full(betas.shape, np.nan)
    try:
        z = [(float(t[0]), float(t[1])) for t in model.z_stat]
    except Exception:
        z = [(np.nan, np.nan)] * len(betas)
    # spreg's name_x already includes the leading 'CONSTANT'; the spatial
    # parameter (rho for lag, lambda for error) is the last beta.
    names = list(model.name_x)
    while len(names) < len(betas):
        names.append(spatial_term)
    rows = []
    for i, nm in enumerate(names):
        b = betas[i] if i < len(betas) else np.nan
        s = se[i] if i < len(se) else np.nan
        p = z[i][1] if i < len(z) else np.nan
        rows.append({
            "outcome": outcome_label, "model": model_name, "term": nm,
            "coef": b, "se": s, "p_value": p,
        })
    return pd.DataFrame(rows)


def main() -> None:
    master = pd.read_csv(DATA_CLEAN_DIR / "station_master.csv")
    noise = pd.read_csv(DATA_CLEAN_DIR / "noise_metrics_station.csv")[
        ["station_key", "laeq_24h_leq_db_a", "lday_leq_db_a", "lnight_leq_db_a"]
    ]
    df = master.drop(columns=[c for c in ["laeq_24h_leq_db_a"] if c in master.columns]).merge(
        noise, on="station_key", how="inner"
    )

    from spreg import ML_Error, ML_Lag

    coef_frames = []
    moran_rows = []

    for ycol, ylabel in OUTCOMES.items():
        cols = ["station_key", "lng_round", "lat_round", ycol] + SOURCE_COLS
        d = df[cols].copy()
        for c in [ycol, "lng_round", "lat_round"] + SOURCE_COLS:
            d[c] = pd.to_numeric(d[c], errors="coerce")
        d = d.dropna(subset=[ycol, "lng_round", "lat_round"] + SOURCE_COLS).reset_index(drop=True)

        y = d[ycol].to_numpy(dtype=float)
        ym = y.reshape(-1, 1)
        coords = d[["lng_round", "lat_round"]].to_numpy(dtype=float)
        w = build_knn_weights(coords, K_NEIGHBORS)

        # --- (0) Intercept-only residual Moran (raw outcome clustering) ---
        resid_raw = y - y.mean()
        mi_raw, p_raw = residual_moran(resid_raw, w)
        moran_rows.append({"outcome": ylabel, "spec": "intercept_only",
                           "resid_moran_I": mi_raw, "p_sim": p_raw})

        # Spatial trend-surface basis (centred quadratic lng/lat).
        lng = d["lng_round"].to_numpy(dtype=float)
        lat = d["lat_round"].to_numpy(dtype=float)
        lng_c = lng - lng.mean()
        lat_c = lat - lat.mean()
        trend = np.column_stack([lng_c, lat_c, lng_c**2, lat_c**2, lng_c * lat_c])
        trend_names = ["lng", "lat", "lng2", "lat2", "lng_lat"]

        # --- Univariate spatial test for each headline source term ---
        for src in SOURCE_COLS:
            xs = d[[src]].to_numpy(dtype=float)

            # OLS (single source predictor) + residual Moran.
            ols = sm.OLS(y, sm.add_constant(xs, has_constant="add")).fit()
            coef_frames.append(pd.DataFrame([
                {"outcome": ylabel, "model": "OLS", "term": "CONSTANT",
                 "coef": ols.params[0], "se": ols.bse[0], "p_value": ols.pvalues[0]},
                {"outcome": ylabel, "model": "OLS", "term": src,
                 "coef": ols.params[1], "se": ols.bse[1], "p_value": ols.pvalues[1]},
            ]))
            mi_ols, p_ols = residual_moran(ols.resid, w)
            moran_rows.append({"outcome": ylabel, "spec": f"OLS[{src}]",
                               "resid_moran_I": mi_ols, "p_sim": p_ols})

            # OLS + quadratic spatial trend surface.
            Xt = sm.add_constant(np.column_stack([xs, trend]), has_constant="add")
            ols_t = sm.OLS(y, Xt).fit()
            t_names = ["CONSTANT", src] + trend_names
            for i, nm in enumerate(t_names):
                coef_frames.append(pd.DataFrame([{
                    "outcome": ylabel, "model": "OLS_trend", "term": nm,
                    "coef": ols_t.params[i], "se": ols_t.bse[i], "p_value": ols_t.pvalues[i],
                }]))
            mi_t, p_t = residual_moran(ols_t.resid, w)
            moran_rows.append({"outcome": ylabel, "spec": f"OLS_trend[{src}]",
                               "resid_moran_I": mi_t, "p_sim": p_t})

            # SAR (ML_Lag) and SEM (ML_Error) with the single source predictor.
            sar = ML_Lag(ym, xs, w=w, name_y=ylabel, name_x=[src])
            coef_frames.append(extract_spreg(sar, "SAR", ylabel, "W_rho"))
            sem = ML_Error(ym, xs, w=w, name_y=ylabel, name_x=[src])
            coef_frames.append(extract_spreg(sem, "SEM", ylabel, "lambda"))

    coef_tab = pd.concat(coef_frames, ignore_index=True)
    coef_tab.to_csv(OUTPUT_DIR / "28_source_spatial__coefficients.csv", index=False)

    moran_tab = pd.DataFrame(moran_rows)
    moran_tab.to_csv(OUTPUT_DIR / "28_source_spatial__residual_moran.csv", index=False)

    # Build a clean source-term summary across model specs.
    src_tab = coef_tab[coef_tab["term"].isin(SOURCE_COLS)].copy()
    src_tab["sig_05"] = src_tab["p_value"] < 0.05
    src_tab["direction"] = np.sign(src_tab["coef"]).map({1.0: "+", -1.0: "-", 0.0: "0"})
    pivot = src_tab.pivot_table(
        index=["outcome", "term"], columns="model",
        values=["coef", "p_value"], aggfunc="first"
    )
    pivot.columns = [f"{a}__{b}" for a, b in pivot.columns]
    pivot = pivot.reset_index()
    pivot.to_csv(OUTPUT_DIR / "28_source_spatial__source_term_summary.csv", index=False)

    # Console digest.
    print("Residual Moran's I by spec:")
    print(moran_tab.round(3).to_string(index=False))
    print("\nSource-term coefficients across models:")
    show = src_tab[["outcome", "term", "model", "coef", "p_value", "sig_05"]]
    print(show.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
