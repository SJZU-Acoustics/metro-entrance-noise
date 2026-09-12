from __future__ import annotations

import argparse
import math
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "mplconfig"))

TOPIC_ROOT = Path(__file__).resolve().parents[2]  # repository root
DATA_CLEAN_DIR = TOPIC_ROOT / "intermediate" / "clean"
OUTPUT_DIR = TOPIC_ROOT / "intermediate" / "analysis"
PERM_SEED = 20260206  # fixes the permutation draws of the global and local Moran tests


@dataclass(frozen=True)
class PeriodSpec:
    period: str
    y_col: str
    label_en: str


DEFAULT_PERIODS = [
    PeriodSpec(period="day", y_col="lday_leq_db_a", label_en="Lday"),
    PeriodSpec(period="night", y_col="lnight_leq_db_a", label_en="Lnight"),
    PeriodSpec(period="peak", y_col="lpeak_leq_db_a", label_en="Peak (07–09 + 17–19)"),
    PeriodSpec(period="offpeak", y_col="linter_peak_leq_db_a", label_en="Inter-peak (10–16)"),
]


DEFAULT_BASELINE_PREDICTORS = [
    "100m_urban_arterial_road_density_m_m",
    "500m_parks_green_space_land_pct",
    "500m_residential_land_pct",
    "500m_commercial_services_land_pct",
    "300m_parks_green_space_land_pct",
    "100m_shannon_entropy_index",
]


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Step 11: spatial robustness checks and spatial models (kNN weights, Moran's I, LISA, spatial CV, SAR/SEM)."
    )
    p.add_argument(
        "--station-master",
        type=str,
        default=str(DATA_CLEAN_DIR / "station_master.csv"),
        help="Station master CSV (default: data/clean/station_master.csv).",
    )
    p.add_argument(
        "--noise-metrics",
        type=str,
        default=str(DATA_CLEAN_DIR / "noise_metrics_station.csv"),
        help="Noise metrics station CSV (default: data/clean/noise_metrics_station.csv).",
    )
    p.add_argument(
        "--k",
        type=int,
        default=8,
        help="k for kNN spatial weights (default: 8).",
    )
    p.add_argument(
        "--permutations",
        type=int,
        default=999,
        help="Permutations for Moran's I and Local Moran p-values (default: 999).",
    )
    p.add_argument(
        "--periods",
        type=str,
        default=",".join([f"{s.period}:{s.y_col}:{s.label_en}" for s in DEFAULT_PERIODS]),
        help="Comma-separated period specs: period:y_col:label_en.",
    )
    p.add_argument(
        "--baseline-predictors",
        type=str,
        default=",".join(DEFAULT_BASELINE_PREDICTORS),
        help="Comma-separated baseline predictors (fixed structure variable set).",
    )
    p.add_argument(
        "--cv-folds",
        type=int,
        default=5,
        help="Folds for CV comparisons (default: 5).",
    )
    p.add_argument(
        "--spatial-blocks",
        type=int,
        default=8,
        help="Number of spatial blocks for block CV (k-means on coords, default: 8).",
    )
    p.add_argument(
        "--seed",
        type=int,
        default=20260206,
        help="Random seed (default: 20260206).",
    )
    p.add_argument(
        "--out-dir",
        type=str,
        default=str(OUTPUT_DIR),
        help="Output directory (default: output/).",
    )
    return p


def parse_periods(s: str) -> list[PeriodSpec]:
    out: list[PeriodSpec] = []
    for part in s.split(","):
        part = part.strip()
        if not part:
            continue
        bits = part.split(":", 2)
        if len(bits) != 3:
            raise ValueError(f"Invalid period spec: {part} (expected period:y_col:label_en)")
        out.append(PeriodSpec(period=bits[0], y_col=bits[1], label_en=bits[2]))
    if not out:
        raise ValueError("No periods specified.")
    return out


def zscore(series: pd.Series) -> pd.Series:
    s = pd.to_numeric(series, errors="coerce")
    sd = float(s.std(ddof=0))
    if not np.isfinite(sd) or sd == 0:
        return s * 0.0
    return (s - float(s.mean())) / sd


def build_knn_weights(coords: np.ndarray, k: int):
    from libpysal.weights import KNN

    w = KNN.from_array(coords, k=k)
    w.transform = "R"
    return w


def global_moran(values: np.ndarray, w, permutations: int):
    from esda.moran import Moran

    np.random.seed(PERM_SEED)
    moran = Moran(values, w, permutations=permutations)
    return moran


def local_moran(values: np.ndarray, w, permutations: int):
    from esda.moran import Moran_Local

    loc = Moran_Local(values, w, permutations=permutations, seed=PERM_SEED)
    return loc


def lisa_cluster_labels(values_z: np.ndarray, w, loc) -> list[str]:
    # Standard local moran cluster types based on quadrant and significant p.
    # loc.q: 1 HH, 2 LH, 3 LL, 4 HL
    labels = []
    for q, p in zip(loc.q, loc.p_sim, strict=False):
        if not np.isfinite(p) or p > 0.05:
            labels.append("NS")
        else:
            labels.append({1: "HH", 2: "LH", 3: "LL", 4: "HL"}.get(int(q), "NS"))
    return labels


def plot_lisa_map(df: pd.DataFrame, x_col: str, y_col: str, title: str, out_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme(style="whitegrid", palette="colorblind")
    palette = {"HH": "#DD8452", "LL": "#4C72B0", "HL": "#55A868", "LH": "#C44E52", "NS": "#BBBBBB"}
    fig, ax = plt.subplots(figsize=(7.2, 5.2), dpi=300)
    sns.scatterplot(
        data=df,
        x=x_col,
        y=y_col,
        hue="lisa_cluster",
        palette=palette,
        s=45,
        alpha=0.9,
        ax=ax,
        edgecolor="white",
        linewidth=0.3,
    )
    ax.set_title(title)
    ax.set_xlabel("Longitude (rounded)")
    ax.set_ylabel("Latitude (rounded)")
    ax.legend(title="LISA cluster", loc="best", frameon=True, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def fit_baseline(df: pd.DataFrame, y_col: str, x_cols: list[str]) -> tuple[pd.Series, pd.Series, pd.DataFrame]:
    import statsmodels.api as sm

    work = df[["station_key", y_col] + x_cols].copy()
    work[y_col] = pd.to_numeric(work[y_col], errors="coerce")
    for c in x_cols:
        work[c] = zscore(work[c])
    work = work.dropna(subset=[y_col] + x_cols)
    X = sm.add_constant(work[x_cols], has_constant="add").apply(pd.to_numeric, errors="coerce").astype(float)
    y = work[y_col].astype(float)
    model = sm.OLS(y, X).fit(cov_type="HC3")
    pred = pd.Series(model.predict(X), index=work.index, name="baseline_pred")
    resid = y - pred
    coef = model.params.rename("coef").to_frame()
    coef["se_hc3"] = model.bse
    coef["p_value"] = model.pvalues
    coef = coef.reset_index(names="term")
    return pred, resid.rename("baseline_resid"), coef


def cv_compare(df: pd.DataFrame, y_col: str, x_cols: list[str], cv_folds: int, blocks: int, seed: int) -> pd.DataFrame:
    from sklearn.cluster import KMeans
    from sklearn.model_selection import KFold
    from sklearn.metrics import mean_squared_error

    import statsmodels.api as sm

    work = df[["lng_round", "lat_round", y_col] + x_cols].copy()
    work[y_col] = pd.to_numeric(work[y_col], errors="coerce")
    for c in x_cols:
        work[c] = zscore(work[c])
    work = work.dropna(subset=["lng_round", "lat_round", y_col] + x_cols)

    X = sm.add_constant(work[x_cols], has_constant="add").astype(float).to_numpy()
    y = work[y_col].to_numpy(dtype=float)

    kf = KFold(n_splits=cv_folds, shuffle=True, random_state=seed)
    preds = np.full_like(y, np.nan, dtype=float)
    for train, test in kf.split(X):
        model = sm.OLS(y[train], X[train]).fit()
        preds[test] = model.predict(X[test])
    rmse_kfold = float(np.sqrt(mean_squared_error(y, preds)))

    # Spatial block CV via k-means clusters; use clusters as folds (up to cv_folds by merging).
    coords = work[["lng_round", "lat_round"]].to_numpy(dtype=float)
    # Avoid joblib/loky warnings and keep deterministic behavior.
    os.environ.setdefault("LOKY_MAX_CPU_COUNT", "1")
    km = KMeans(n_clusters=int(blocks), random_state=seed, n_init=10)
    block = km.fit_predict(coords)
    # Convert blocks to cv_folds by hashing block id to fold (keeps separation).
    fold = np.mod(block, int(cv_folds))
    preds2 = np.full_like(y, np.nan, dtype=float)
    for f in range(int(cv_folds)):
        test = np.where(fold == f)[0]
        train = np.where(fold != f)[0]
        if train.size < X.shape[1] + 5 or test.size == 0:
            continue
        model = sm.OLS(y[train], X[train]).fit()
        preds2[test] = model.predict(X[test])
    ok = np.isfinite(preds2)
    rmse_block = float(np.sqrt(mean_squared_error(y[ok], preds2[ok]))) if ok.sum() > 0 else math.nan

    return pd.DataFrame(
        [
            {"outcome": y_col, "cv": "kfold", "rmse": rmse_kfold, "n": int(len(work))},
            {"outcome": y_col, "cv": "spatial_block_kmeans", "rmse": rmse_block, "n": int(len(work))},
        ]
    )


def fit_sar_sem(y: np.ndarray, X: np.ndarray, w):
    from spreg import ML_Error, ML_Lag

    # spreg expects 2D arrays.
    y2 = y.reshape((-1, 1))
    # Already includes constant.
    sar = ML_Lag(y2, X, w=w, name_y="y", name_x=[f"x{i}" for i in range(X.shape[1])])
    sem = ML_Error(y2, X, w=w, name_y="y", name_x=[f"x{i}" for i in range(X.shape[1])])
    return sar, sem


def extract_spreg_results(model, model_name: str, outcome: str, terms: list[str]) -> pd.DataFrame:
    # betas includes constant and x; sar includes rho, sem includes lambda at end.
    betas = np.array(model.betas).reshape((-1,))
    std = np.array(model.std_err).reshape((-1,))
    z = np.array(model.z_stat)  # list of tuples
    pvals = np.array([t[1] for t in z], dtype=float)

    # Determine extra parameter name for SAR/SEM.
    extra = []
    if model_name == "SAR":
        extra = ["rho"]
    elif model_name == "SEM":
        extra = ["lambda"]

    names = terms + extra
    rows = []
    for i, name in enumerate(names):
        coef = float(betas[i]) if i < betas.size else math.nan
        se = float(std[i]) if i < std.size else math.nan
        p = float(pvals[i]) if i < pvals.size else math.nan
        rows.append({"outcome": outcome, "model": model_name, "term": name, "coef": coef, "se": se, "p_value": p})
    return pd.DataFrame(rows)


def main() -> None:
    args = build_parser().parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    periods = parse_periods(args.periods)
    x_cols = [c.strip() for c in args.baseline_predictors.split(",") if c.strip()]

    station = pd.read_csv(args.station_master)
    noise = pd.read_csv(args.noise_metrics)
    df0 = station.merge(noise, on="station_key", how="left", suffixes=("", "__noise"))

    needed = ["station_key", "lng_round", "lat_round"] + x_cols
    missing = [c for c in needed if c not in df0.columns]
    if missing:
        raise SystemExit(f"Missing required columns: {missing}")

    # Coordinates and basic checks.
    coords_all = df0[["lng_round", "lat_round"]].apply(pd.to_numeric, errors="coerce")
    mask_coords = coords_all.notna().all(axis=1)
    coords_ok = coords_all.loc[mask_coords].to_numpy(dtype=float)

    # Summary on full coordinate set (before outcome-specific filtering).
    w0 = build_knn_weights(coords_ok, k=int(args.k))
    w_summary = pd.DataFrame(
        [
            {
                "k": int(args.k),
                "n": int(len(coords_ok)),
                "islands": int(len(w0.islands)),
                "components": int(w0.n_components),
            }
        ]
    )
    w_summary.to_csv(out_dir / "12_spatial__weights_summary.csv", index=False)

    global_rows: list[dict] = []
    lisa_outputs: list[pd.DataFrame] = []
    coef_outputs: list[pd.DataFrame] = []
    cv_outputs: list[pd.DataFrame] = []
    spreg_outputs: list[pd.DataFrame] = []

    for spec in periods:
        if spec.y_col not in df0.columns:
            continue
        # Baseline residuals for this outcome (baseline uses only rows with complete X and Y).
        pred, resid, coef = fit_baseline(df0.loc[mask_coords].copy(), y_col=spec.y_col, x_cols=x_cols)
        coef["period"] = spec.period
        coef["label_en"] = spec.label_en
        coef_outputs.append(coef)

        # Build modeling frame for spatial tests.
        work = df0.loc[mask_coords, ["station_key", "lng_round", "lat_round", spec.y_col]].copy()
        work[spec.y_col] = pd.to_numeric(work[spec.y_col], errors="coerce")
        work = work.dropna(subset=[spec.y_col]).copy()
        work["baseline_resid"] = resid.values
        work = work.dropna(subset=["baseline_resid"]).copy()

        if len(work) < 30:
            continue
        coords_y = work[["lng_round", "lat_round"]].to_numpy(dtype=float)
        w_y = build_knn_weights(coords_y, k=int(args.k))
        y_valid = work[spec.y_col].to_numpy(dtype=float)
        y_sd = float(np.std(y_valid, ddof=0))
        y_z = (y_valid - float(np.mean(y_valid))) / y_sd if np.isfinite(y_sd) and y_sd > 0 else y_valid * 0.0

        mor = global_moran(y_z, w_y, permutations=int(args.permutations))
        global_rows.append(
            {
                "period": spec.period,
                "var": spec.y_col,
                "target": "observed",
                "n": int(len(work)),
                "moran_I": float(mor.I),
                "z_norm": float(mor.z_norm),
                "p_norm": float(mor.p_norm),
                "p_sim": float(mor.p_sim),
            }
        )

        loc = local_moran(y_z, w_y, permutations=int(args.permutations))
        clusters = lisa_cluster_labels(y_z, w_y, loc)
        lisa_df = pd.DataFrame(
            {
                "station_key": work["station_key"].astype(str).tolist(),
                "period": spec.period,
                "var": spec.y_col,
                "z": y_z.tolist(),
                "local_I": loc.Is.tolist(),
                "p_sim": loc.p_sim.tolist(),
                "q": loc.q.tolist(),
                "lisa_cluster": clusters,
            }
        )
        lisa_df.to_csv(out_dir / f"12_spatial__lisa__{spec.period}__{spec.y_col}.csv", index=False)
        lisa_outputs.append(lisa_df)

        # Plot LISA clusters on lon/lat.
        plot_df = work[["station_key", "lng_round", "lat_round"]].copy()
        plot_df["lisa_cluster"] = clusters
        plot_lisa_map(
            plot_df,
            x_col="lng_round",
            y_col="lat_round",
            title=f"LISA clusters ({spec.label_en})",
            out_path=out_dir / f"12_spatial__lisa_map__{spec.period}.png",
        )

        # Baseline residual Moran's I (residuals aligned)
        r = work["baseline_resid"].to_numpy(dtype=float)
        r_sd = float(np.std(r, ddof=0))
        r_z = (r - float(np.mean(r))) / r_sd if np.isfinite(r_sd) and r_sd > 0 else r * 0.0
        mor_r = global_moran(r_z, w_y, permutations=int(args.permutations))
        global_rows.append(
            {
                "period": spec.period,
                "var": spec.y_col,
                "target": "baseline_resid",
                "n": int(len(work)),
                "moran_I": float(mor_r.I),
                "z_norm": float(mor_r.z_norm),
                "p_norm": float(mor_r.p_norm),
                "p_sim": float(mor_r.p_sim),
            }
        )

        # CV comparison for baseline model.
        cv_outputs.append(
            cv_compare(
                df=df0.loc[mask_coords].copy(),
                y_col=spec.y_col,
                x_cols=x_cols,
                cv_folds=int(args.cv_folds),
                blocks=int(args.spatial_blocks),
                seed=int(args.seed),
            )
        )

        # Spatial regression SAR/SEM using same baseline predictors.
        try:
            import statsmodels.api as sm

            base = df0.loc[mask_coords, ["lng_round", "lat_round", spec.y_col] + x_cols].copy()
            base[spec.y_col] = pd.to_numeric(base[spec.y_col], errors="coerce")
            for c in x_cols:
                base[c] = zscore(base[c])
            base = base.dropna(subset=["lng_round", "lat_round", spec.y_col] + x_cols)
            coords2 = base[["lng_round", "lat_round"]].to_numpy(dtype=float)
            w2 = build_knn_weights(coords2, k=int(args.k))
            y2 = base[spec.y_col].to_numpy(dtype=float)
            X2 = sm.add_constant(base[x_cols], has_constant="add").to_numpy(dtype=float)
            sar, sem = fit_sar_sem(y2, X2, w2)
            terms = ["const"] + x_cols
            spreg_outputs.append(extract_spreg_results(sar, "SAR", spec.y_col, terms=terms).assign(period=spec.period))
            spreg_outputs.append(extract_spreg_results(sem, "SEM", spec.y_col, terms=terms).assign(period=spec.period))
        except Exception:
            # If SAR/SEM fails for a period, skip silently but keep other outputs.
            continue

    pd.DataFrame(global_rows).to_csv(out_dir / "12_spatial__moran_global.csv", index=False)
    if coef_outputs:
        pd.concat(coef_outputs, ignore_index=True).to_csv(out_dir / "12_spatial__baseline_coefficients.csv", index=False)
    if cv_outputs:
        pd.concat(cv_outputs, ignore_index=True).to_csv(out_dir / "12_spatial__cv_rmse.csv", index=False)
    if spreg_outputs:
        pd.concat(spreg_outputs, ignore_index=True).to_csv(out_dir / "12_spatial__sar_sem_coefficients.csv", index=False)


if __name__ == "__main__":
    main()
