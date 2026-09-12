from __future__ import annotations

import argparse
import os
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests

# Headless plotting + writable cache
import tempfile
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "mplconfig"))
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

TOPIC_ROOT = Path(__file__).resolve().parents[2]  # repository root
DATA_CLEAN_DIR = TOPIC_ROOT / "intermediate" / "clean"
OUTPUT_DIR = TOPIC_ROOT / "intermediate" / "analysis"


NOISE_OUTCOMES = [
    "laeq_24h_leq_db_a",
    "lday_leq_db_a",
    "lnight_leq_db_a",
    "lpeak_leq_db_a",
    "linter_peak_leq_db_a",
    # Keep legacy/aux metrics for reference.
    "loffpeak_day_leq_db_a",
    "lpeak_db_a",
]


@dataclass(frozen=True)
class KeyVar:
    category: str
    original_pattern: str
    unit: str


# 9 theory-driven categories (proxy NDVI with green land percentage if NDVI not present).
KEY_VARS: list[KeyVar] = [
    KeyVar("road_arterial_density", r"Urban Arterial Road Density", "m/m²"),
    KeyVar("road_expressway_density", r"Elevated/Expressway Density", "m/m²"),
    KeyVar("road_pedestrian_density", r"Pedestrian Road Density", "m/m²"),
    KeyVar("green_land_pct", r"Parks & Green Space Land", "%"),
    KeyVar("commercial_land_pct", r"Commercial Services Land", "%"),
    KeyVar("residential_land_pct", r"Residential Land", "%"),
    KeyVar("water_body_pct", r"Water Body", "%"),
    KeyVar("landuse_mix_shannon", r"Shannon Entropy Index", ""),
    # POI total density will be constructed per scale by summing POI density columns.
    KeyVar("poi_total_density", r"POI Density", "count/km² (?)"),
]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Step 5: multiscale built environment patterns (correlation + partial).")
    parser.add_argument("--station-master", type=str, default=str(DATA_CLEAN_DIR / "station_master.csv"))
    parser.add_argument("--noise", type=str, default=str(DATA_CLEAN_DIR / "noise_metrics_station.csv"))
    parser.add_argument("--colmap", type=str, default=str(DATA_CLEAN_DIR / "landuse_column_map.csv"))
    parser.add_argument("--outdir", type=str, default=str(OUTPUT_DIR))
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--seed", type=int, default=20260205)
    parser.add_argument("--boot", type=int, default=500, help="Bootstrap iterations for Spearman CI.")
    parser.add_argument("--boot-partial", type=int, default=500, help="Bootstrap iterations for partial correlation CI.")
    parser.add_argument("--min-n", type=int, default=30, help="Minimum pairwise N to compute correlation.")
    return parser


def ensure_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing: {path}")


def parse_scale_m(original: str) -> int | None:
    m = re.match(r"^(\d{2,4})m\s+", str(original).strip())
    if m:
        return int(m.group(1))
    # e.g. "... (100m)"
    m2 = re.search(r"\((\d{2,4})m\)", str(original))
    if m2:
        return int(m2.group(1))
    return None


def spearman_r_p(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int]:
    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]
    y = y[ok]
    n = int(len(x))
    if n < 3:
        return np.nan, np.nan, n
    rho, p = stats.spearmanr(x, y)
    return float(rho), float(p), n


def pearson_r_p(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int]:
    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]
    y = y[ok]
    n = int(len(x))
    if n < 3:
        return np.nan, np.nan, n
    if np.nanstd(x) == 0 or np.nanstd(y) == 0:
        return np.nan, np.nan, n
    r, p = stats.pearsonr(x, y)
    return float(r), float(p), n


def bootstrap_ci(
    x: np.ndarray,
    y: np.ndarray,
    fn,
    rng: np.random.Generator,
    b: int,
    alpha: float = 0.05,
) -> tuple[float, float]:
    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]
    y = y[ok]
    n = len(x)
    if n < 5:
        return np.nan, np.nan
    stats_b = []
    for _ in range(b):
        idx = rng.integers(0, n, size=n)
        stats_b.append(fn(x[idx], y[idx]))
    arr = np.array(stats_b, dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) < max(10, b * 0.5):
        return np.nan, np.nan
    lo = float(np.quantile(arr, alpha / 2))
    hi = float(np.quantile(arr, 1 - alpha / 2))
    return lo, hi


def bootstrap_ci_from_residuals(
    rx: np.ndarray,
    ry: np.ndarray,
    rng: np.random.Generator,
    b: int,
    alpha: float = 0.05,
) -> tuple[float, float]:
    ok = np.isfinite(rx) & np.isfinite(ry)
    rx = rx[ok]
    ry = ry[ok]
    n = len(rx)
    if n < 5:
        return np.nan, np.nan
    if np.nanstd(rx) == 0 or np.nanstd(ry) == 0:
        return np.nan, np.nan
    stats_b = []
    for _ in range(b):
        idx = rng.integers(0, n, size=n)
        r, _, _ = pearson_r_p(rx[idx], ry[idx])
        stats_b.append(r)
    arr = np.array(stats_b, dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) < max(10, b * 0.5):
        return np.nan, np.nan
    lo = float(np.quantile(arr, alpha / 2))
    hi = float(np.quantile(arr, 1 - alpha / 2))
    return lo, hi


def design_matrix_for_covariates(df: pd.DataFrame, covariates: list[str]) -> np.ndarray:
    # One-hot encode categorical covariates and add intercept.
    parts = []
    for c in covariates:
        d = pd.get_dummies(df[c].astype("string"), prefix=c, dummy_na=False)
        parts.append(d)
    X = pd.concat(parts, axis=1) if parts else pd.DataFrame(index=df.index)
    X.insert(0, "intercept", 1.0)
    return X.to_numpy(float)


def residualize(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    ok = np.isfinite(y) & np.isfinite(X).all(axis=1)
    yy = y[ok]
    XX = X[ok]
    if len(yy) < 5:
        out = np.full_like(y, np.nan, dtype=float)
        return out
    beta, *_ = np.linalg.lstsq(XX, yy, rcond=None)
    resid = yy - XX @ beta
    out = np.full_like(y, np.nan, dtype=float)
    out[ok] = resid
    return out


def partial_spearman(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    covariates: list[str],
) -> tuple[float, float, int]:
    tmp = df[[x_col, y_col] + covariates].copy()
    tmp[x_col] = pd.to_numeric(tmp[x_col], errors="coerce")
    tmp[y_col] = pd.to_numeric(tmp[y_col], errors="coerce")
    tmp = tmp.dropna(subset=[x_col, y_col] + covariates)
    n = int(len(tmp))
    if n < 10:
        return np.nan, np.nan, n

    x_rank = tmp[x_col].rank(method="average").to_numpy(float)
    y_rank = tmp[y_col].rank(method="average").to_numpy(float)
    X = design_matrix_for_covariates(tmp, covariates=covariates)
    rx = residualize(x_rank, X)
    ry = residualize(y_rank, X)
    r, p, _ = pearson_r_p(rx, ry)
    return r, p, n


def main() -> int:
    args = build_parser().parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    station_master_path = Path(args.station_master)
    noise_path = Path(args.noise)
    colmap_path = Path(args.colmap)
    ensure_exists(station_master_path)
    ensure_exists(noise_path)
    ensure_exists(colmap_path)

    station = pd.read_csv(station_master_path)
    noise = pd.read_csv(noise_path)
    colmap = pd.read_csv(colmap_path)

    # Merge station + noise (prefer noise columns for outcomes).
    overlap = set(station.columns).intersection(set(noise.columns)) - {"station_key"}
    if overlap:
        station = station.drop(columns=sorted(overlap))
    df = station.merge(noise, on="station_key", how="left", validate="one_to_one")

    # Prepare mapping original->clean and scales.
    colmap["scale_m"] = colmap["original"].map(parse_scale_m)
    colmap["original"] = colmap["original"].astype(str)
    colmap["clean_name"] = colmap["clean_name"].astype(str)

    # Identify POI columns per scale from original names.
    poi_map: dict[int, list[str]] = {}
    for _, r in colmap.iterrows():
        orig = r["original"]
        clean = r["clean_name"]
        scale = r["scale_m"]
        if scale is None:
            continue
        if "POI Density" in orig:
            poi_map.setdefault(int(scale), []).append(clean)

    # Build key variable list: for each category and each scale, get clean column.
    key_rows: list[dict] = []
    selected_cols: list[str] = []
    scales = sorted({s for s in colmap["scale_m"].dropna().astype(int).tolist() if s in {100, 300, 500}})

    for kv in KEY_VARS:
        if kv.category == "poi_total_density":
            for s in scales:
                cols = [c for c in poi_map.get(s, []) if c in df.columns]
                if not cols:
                    continue
                new_col = f"poi_total_density_{s}m"
                num = df[cols].apply(pd.to_numeric, errors="coerce")
                df[new_col] = num.sum(axis=1, min_count=1)
                key_rows.append(
                    {
                        "category": kv.category,
                        "scale_m": s,
                        "clean_name": new_col,
                        "original_name": f"SUM(POI Density ({s}m) columns)",
                        "unit": kv.unit,
                        "notes": f"summed {len(cols)} POI density columns",
                    }
                )
                selected_cols.append(new_col)
            continue

        pat = re.compile(kv.original_pattern, flags=re.IGNORECASE)
        for s in scales:
            candidates = colmap.loc[(colmap["scale_m"] == s) & (colmap["original"].map(lambda x: bool(pat.search(x))))].copy()
            if candidates.empty:
                continue
            # Prefer unique match; if multiple, take the first but record.
            cand = candidates.iloc[0]
            clean = str(cand["clean_name"])
            if clean not in df.columns:
                continue
            key_rows.append(
                {
                    "category": kv.category,
                    "scale_m": int(s),
                    "clean_name": clean,
                    "original_name": str(cand["original"]),
                    "unit": kv.unit,
                    "notes": f"matched {len(candidates)} columns" if len(candidates) > 1 else "",
                }
            )
            selected_cols.append(clean)

    # De-fragment after multiple inserts (helps performance).
    df = df.copy()

    keyvars = pd.DataFrame(key_rows).sort_values(["category", "scale_m"], ignore_index=True)
    keyvars.to_csv(outdir / "06_built_env__keyvars.csv", index=False, encoding="utf-8")

    # Correlation analyses
    outcomes = [o for o in NOISE_OUTCOMES if o in df.columns]
    covars = ["line_name", "station_type_code"]
    for c in covars:
        df[c] = df[c].astype("string")

    results: list[dict] = []
    for outcome in outcomes:
        y = pd.to_numeric(df[outcome], errors="coerce").to_numpy(float)
        for _, r in keyvars.iterrows():
            xcol = str(r["clean_name"])
            x = pd.to_numeric(df[xcol], errors="coerce").to_numpy(float)
            rho, p, n = spearman_r_p(x, y)
            if n < args.min_n:
                rho, p = np.nan, np.nan
            # Bootstrap CI for Spearman via rank+Pearson (faster).
            ok = np.isfinite(x) & np.isfinite(y)
            xr = stats.rankdata(x[ok]) if ok.any() else np.array([])
            yr = stats.rankdata(y[ok]) if ok.any() else np.array([])
            ci_lo, ci_hi = bootstrap_ci_from_residuals(xr, yr, rng=rng, b=args.boot)
            results.append(
                {
                    "method": "spearman",
                    "outcome": outcome,
                    "category": r["category"],
                    "scale_m": int(r["scale_m"]),
                    "x": xcol,
                    "rho": rho,
                    "p_value": p,
                    "ci_low": ci_lo,
                    "ci_high": ci_hi,
                    "n": n,
                }
            )

            pr, pp, pn = partial_spearman(df, x_col=xcol, y_col=outcome, covariates=covars)
            if pn < args.min_n:
                pr, pp = np.nan, np.nan
            # Bootstrap CI for partial via fixed residuals (fast).
            tmp = df[[xcol, outcome] + covars].copy()
            tmp[xcol] = pd.to_numeric(tmp[xcol], errors="coerce")
            tmp[outcome] = pd.to_numeric(tmp[outcome], errors="coerce")
            tmp = tmp.dropna(subset=[xcol, outcome] + covars)
            if len(tmp) >= args.min_n:
                x_rank = tmp[xcol].rank(method="average").to_numpy(float)
                y_rank = tmp[outcome].rank(method="average").to_numpy(float)
                X = design_matrix_for_covariates(tmp, covariates=covars)
                rx = residualize(x_rank, X)
                ry = residualize(y_rank, X)
                ci2_lo, ci2_hi = bootstrap_ci_from_residuals(rx, ry, rng=rng, b=args.boot_partial)
            else:
                ci2_lo, ci2_hi = np.nan, np.nan
            results.append(
                {
                    "method": "partial_spearman",
                    "outcome": outcome,
                    "category": r["category"],
                    "scale_m": int(r["scale_m"]),
                    "x": xcol,
                    "rho": pr,
                    "p_value": pp,
                    "ci_low": ci2_lo,
                    "ci_high": ci2_hi,
                    "n": pn,
                }
            )

    res = pd.DataFrame(results)

    # FDR correction within (method, outcome)
    res["q_fdr"] = np.nan
    for (method, outcome), sub_idx in res.groupby(["method", "outcome"]).groups.items():
        pvals = res.loc[sub_idx, "p_value"].to_numpy(float)
        valid = np.isfinite(pvals)
        if valid.any():
            _, qvals, _, _ = multipletests(pvals[valid], method="fdr_bh")
            q_full = np.full_like(pvals, np.nan, dtype=float)
            q_full[valid] = qvals
            res.loc[sub_idx, "q_fdr"] = q_full

    res = res.sort_values(["method", "outcome", "q_fdr", "p_value"], ignore_index=True)
    res.to_csv(outdir / "06_built_env__correlations.csv", index=False, encoding="utf-8")

    # Also export matrices for each (method, outcome): rho and q_fdr.
    for method in ["spearman", "partial_spearman"]:
        for outcome in outcomes:
            sub = res[(res["method"] == method) & (res["outcome"] == outcome)].copy()
            if sub.empty:
                continue
            rho_mat = sub.pivot_table(index="category", columns="scale_m", values="rho", aggfunc="first")
            q_mat = sub.pivot_table(index="category", columns="scale_m", values="q_fdr", aggfunc="first")
            cat_order = [k.category for k in KEY_VARS]
            rho_mat = rho_mat.reindex(cat_order)
            q_mat = q_mat.reindex(cat_order)
            rho_mat.to_csv(outdir / f"06_built_env__matrix_rho__{method}__{outcome}.csv", encoding="utf-8")
            q_mat.to_csv(outdir / f"06_built_env__matrix_q__{method}__{outcome}.csv", encoding="utf-8")

    # Figures: heatmap per outcome, compare scales (rows=category, cols=scale)
    sns.set_theme(style="whitegrid", palette="colorblind")
    for method in ["spearman", "partial_spearman"]:
        for outcome in outcomes:
            sub = res[(res["method"] == method) & (res["outcome"] == outcome)].copy()
            if sub.empty:
                continue
            # Build matrix category x scale with rho.
            mat = sub.pivot_table(index="category", columns="scale_m", values="rho", aggfunc="first")
            # Order categories as in KEY_VARS.
            cat_order = [k.category for k in KEY_VARS]
            mat = mat.reindex(cat_order)
            plt.figure(figsize=(6.5, 4.5))
            ax = sns.heatmap(
                mat,
                cmap="coolwarm",
                center=0.0,
                vmin=-1,
                vmax=1,
                linewidths=0.3,
                linecolor="white",
                cbar_kws={"label": "Spearman rho" if method == "spearman" else "Partial Spearman (residualized)"},
            )
            ax.set_xlabel("Scale (m)")
            ax.set_ylabel("Key built environment category")
            ax.set_title(f"{method}: {outcome} vs built environment (multiscale)")
            plt.tight_layout()
            fig_name = f"06_built_env__{method}__{outcome}__heatmap.png"
            plt.savefig(outdir / fig_name, dpi=args.dpi, bbox_inches="tight")
            plt.close()

    print("Step 5 outputs written to:", outdir)
    print("- 06_built_env__keyvars.csv")
    print("- 06_built_env__correlations.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
