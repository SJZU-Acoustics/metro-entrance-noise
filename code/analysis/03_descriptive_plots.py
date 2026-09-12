from __future__ import annotations

import argparse
import os
from pathlib import Path

import tempfile
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "mplconfig"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

TOPIC_ROOT = Path(__file__).resolve().parents[2]  # repository root
DATA_CLEAN_DIR = TOPIC_ROOT / "intermediate" / "clean"
OUTPUT_DIR = TOPIC_ROOT / "intermediate" / "analysis"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Step 2 descriptive analysis and core figures.")
    parser.add_argument(
        "--station-master",
        type=str,
        default=str(DATA_CLEAN_DIR / "station_master.csv"),
        help="Station-level master table.",
    )
    parser.add_argument(
        "--panel",
        type=str,
        default=str(DATA_CLEAN_DIR / "station_hour_panel.csv"),
        help="Station-hour panel table.",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default=str(OUTPUT_DIR),
        help="Output directory for figures/tables.",
    )
    parser.add_argument("--dpi", type=int, default=300)
    return parser


def ensure_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing: {path}")


def savefig(path: Path, dpi: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close()


def qstats(series: pd.Series, qs: list[float]) -> dict[str, float]:
    s = pd.to_numeric(series, errors="coerce").dropna()
    out: dict[str, float] = {
        "n": float(len(s)),
        "mean": float(s.mean()) if len(s) else float("nan"),
        "median": float(s.median()) if len(s) else float("nan"),
        "min": float(s.min()) if len(s) else float("nan"),
        "max": float(s.max()) if len(s) else float("nan"),
    }
    if len(s):
        qu = s.quantile(qs)
        for q in qs:
            out[f"p{int(q*100):02d}"] = float(qu.loc[q])
    else:
        for q in qs:
            out[f"p{int(q*100):02d}"] = float("nan")
    out["range"] = out["max"] - out["min"] if np.isfinite(out["min"]) and np.isfinite(out["max"]) else float("nan")
    return out


def main() -> int:
    args = build_parser().parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    station_master_path = Path(args.station_master)
    panel_path = Path(args.panel)
    ensure_exists(station_master_path)
    ensure_exists(panel_path)

    sns.set_theme(style="whitegrid", palette="colorblind")

    station = pd.read_csv(station_master_path)
    panel = pd.read_csv(panel_path)

    noise_path = DATA_CLEAN_DIR / "noise_metrics_station.csv"
    noise = pd.read_csv(noise_path) if noise_path.exists() else None
    if noise is not None and "station_key" in noise.columns:
        station = station.merge(noise, on="station_key", how="left", validate="one_to_one", suffixes=("", "_noise"))

    # --- Overall distribution summary (station-level) ---
    key_cols = ["station_key", "station_name_cn", "line_name", "station_type_code", "lng_round", "lat_round"]
    laeq_cols = {
        "laeq_24h_leq_db_a": "LAeq_24h (Leq)",
        "lday_leq_db_a": "Lday (Leq)",
        "lnight_leq_db_a": "Lnight (Leq)",
        "lpeak_leq_db_a": "LPeak (Leq)",
        "linter_peak_leq_db_a": "LInter-peak (Leq)",
        "mp_leq_db_a": "Morning peak (MP) (Leq)",
        "ep_leq_db_a": "Evening peak (EP) (Leq)",
    }
    qs = [0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95]
    rows = []
    for col, label in laeq_cols.items():
        if col not in station.columns:
            continue
        r = qstats(station[col], qs)
        r["metric"] = col
        r["label"] = label
        rows.append(r)
    summary = pd.DataFrame(rows)[
        ["metric", "label", "n", "mean", "median", "min", "max", "range"]
        + [f"p{int(q*100):02d}" for q in qs]
    ]
    summary.to_csv(outdir / "03_descriptive__laeq_station_summary.csv", index=False, encoding="utf-8")

    # Top10 loudest/quietest (using LAeq_24h Leq).
    metric_col = "laeq_24h_leq_db_a"
    top_base = station[key_cols + [metric_col]].copy()
    top_base[metric_col] = pd.to_numeric(top_base[metric_col], errors="coerce")
    top_loud = top_base.sort_values(metric_col, ascending=False, ignore_index=True).head(10)
    top_quiet = top_base.sort_values(metric_col, ascending=True, ignore_index=True).head(10)
    top_loud.to_csv(outdir / "03_descriptive__top10_loudest.csv", index=False, encoding="utf-8")
    top_quiet.to_csv(outdir / "03_descriptive__top10_quietest.csv", index=False, encoding="utf-8")

    # --- Diurnal pattern: hourly mean + quantile bands ---
    panel["hour"] = pd.to_numeric(panel["hour"], errors="coerce").astype("Int64")
    panel["laeq_hour_db_a"] = pd.to_numeric(panel.get("laeq_hour_db_a", panel.get("laeq_db_a")), errors="coerce")
    hour_stats = (
        panel.groupby("hour", as_index=False)["laeq_hour_db_a"]
        .agg(
            mean="mean",
            median="median",
            p10=lambda s: float(pd.to_numeric(s, errors="coerce").dropna().quantile(0.10)) if s.notna().any() else np.nan,
            p25=lambda s: float(pd.to_numeric(s, errors="coerce").dropna().quantile(0.25)) if s.notna().any() else np.nan,
            p75=lambda s: float(pd.to_numeric(s, errors="coerce").dropna().quantile(0.75)) if s.notna().any() else np.nan,
            p90=lambda s: float(pd.to_numeric(s, errors="coerce").dropna().quantile(0.90)) if s.notna().any() else np.nan,
        )
        .sort_values("hour", ignore_index=True)
    )
    hour_stats.to_csv(outdir / "03_descriptive__hourly_curve_stats.csv", index=False, encoding="utf-8")

    plt.figure(figsize=(6.5, 4.0))
    x = hour_stats["hour"].astype(int)
    plt.fill_between(x, hour_stats["p10"], hour_stats["p90"], alpha=0.18, label="P10–P90")
    plt.fill_between(x, hour_stats["p25"], hour_stats["p75"], alpha=0.28, label="P25–P75")
    plt.plot(x, hour_stats["mean"], linewidth=2.0, label="Mean")
    plt.plot(x, hour_stats["median"], linewidth=1.6, linestyle="--", label="Median")
    plt.xticks(range(0, 24, 2))
    plt.xlabel("Hour of day")
    plt.ylabel("LAeq (dB(A))")
    plt.title("Diurnal pattern of hourly LAeq (all stations)")
    plt.legend(frameon=True, fontsize=9)
    savefig(outdir / "03_descriptive__diurnal_curve_laeq.png", dpi=args.dpi)

    # --- Station x hour heatmap ---
    heat = panel[["station_key", "hour", "laeq_hour_db_a"]].copy()
    heat["laeq_hour_db_a"] = pd.to_numeric(heat["laeq_hour_db_a"], errors="coerce")
    order = station[["station_key", metric_col]].copy()
    order[metric_col] = pd.to_numeric(order[metric_col], errors="coerce")
    order = order.sort_values(metric_col, ascending=False, ignore_index=True)["station_key"].astype(str).tolist()
    pivot = heat.pivot_table(index="station_key", columns="hour", values="laeq_hour_db_a", aggfunc="mean")
    pivot = pivot.reindex(order)

    # Shared SPL color scale (used by heatmap + spatial plots) for easier cross-figure reading.
    # Warm colors = louder, cool colors = quieter, consistent with common intuition.
    spl_cmap = "RdYlBu_r"

    def robust_spl_range(values: pd.Series) -> tuple[float, float]:
        s = pd.to_numeric(values, errors="coerce").dropna()
        if len(s) == 0:
            return 0.0, 1.0
        lo = float(s.quantile(0.02))
        hi = float(s.quantile(0.98))
        lo = float(np.floor(lo))
        hi = float(np.ceil(hi))
        if not np.isfinite(lo) or not np.isfinite(hi) or lo >= hi:
            lo, hi = float(s.min()), float(s.max())
        if hi - lo < 5:
            hi = lo + 5.0
        return lo, hi

    # Use station-level + hourly-level SPL values to define the global range.
    spl_sources = [heat["laeq_hour_db_a"]]
    for c in [
        "laeq_24h_leq_db_a",
        "lday_leq_db_a",
        "lnight_leq_db_a",
        "lpeak_leq_db_a",
        "linter_peak_leq_db_a",
        "mp_leq_db_a",
        "ep_leq_db_a",
    ]:
        if c in station.columns:
            spl_sources.append(pd.to_numeric(station[c], errors="coerce"))
    spl_vmin, spl_vmax = robust_spl_range(pd.concat(spl_sources, ignore_index=True))

    plt.figure(figsize=(8.5, 10.0))
    ax = sns.heatmap(
        pivot,
        cmap=spl_cmap,
        vmin=spl_vmin,
        vmax=spl_vmax,
        cbar_kws={"label": "LAeq (dB(A))"},
        linewidths=0.0,
    )
    ax.set_xlabel("Hour of day")
    ax.set_ylabel("Station (sorted by LAeq_24h)")
    ax.set_title("Station × hour heatmap of LAeq")
    savefig(outdir / "03_descriptive__station_hour_heatmap_laeq.png", dpi=args.dpi)

    # --- Spatial scatter plots (rounded coordinates) ---
    geo_cols = [metric_col, "lng_round", "lat_round"]
    for c in [
        "lday_leq_db_a",
        "lnight_leq_db_a",
        "lpeak_leq_db_a",
        "linter_peak_leq_db_a",
        "mp_leq_db_a",
        "ep_leq_db_a",
    ]:
        if c in station.columns:
            geo_cols.append(c)
    geo = station[key_cols + [c for c in geo_cols if c not in key_cols]].copy()
    for c in [metric_col, "lng_round", "lat_round"] + [
        x
        for x in [
            "lday_leq_db_a",
            "lnight_leq_db_a",
            "lpeak_leq_db_a",
            "linter_peak_leq_db_a",
            "mp_leq_db_a",
            "ep_leq_db_a",
        ]
        if x in geo.columns
    ]:
        geo[c] = pd.to_numeric(geo[c], errors="coerce")

    # Use a shared color scale across station-level SPL outcomes for easier comparison.
    shared_value_cols = [metric_col] + [
        c
        for c in [
            "lday_leq_db_a",
            "lnight_leq_db_a",
            "lpeak_leq_db_a",
            "linter_peak_leq_db_a",
            "mp_leq_db_a",
            "ep_leq_db_a",
        ]
        if c in geo.columns
    ]
    stacked = pd.concat([geo[c] for c in shared_value_cols], ignore_index=True)
    stacked = pd.to_numeric(stacked, errors="coerce").dropna()
    # Keep the same global SPL range as the heatmap for cross-figure comparability.
    vmin, vmax = spl_vmin, spl_vmax

    def pca_order_2d(x: np.ndarray, y: np.ndarray) -> np.ndarray:
        coords = np.column_stack([x, y]).astype(float)
        if coords.shape[0] <= 2:
            return np.arange(coords.shape[0])
        centered = coords - np.nanmean(coords, axis=0, keepdims=True)
        centered = np.nan_to_num(centered, nan=0.0, posinf=0.0, neginf=0.0)
        _, _, vh = np.linalg.svd(centered, full_matrices=False)
        pc1 = vh[0]
        proj = centered @ pc1
        return np.argsort(proj, kind="mergesort")

    def mst_edges_2d(x: np.ndarray, y: np.ndarray) -> list[tuple[int, int]]:
        coords = np.column_stack([x, y]).astype(float)
        n = coords.shape[0]
        if n <= 1:
            return []
        start = int(np.lexsort((coords[:, 1], coords[:, 0]))[0])
        visited = np.zeros(n, dtype=bool)
        visited[start] = True

        min_dist = np.full(n, np.inf)
        parent = np.full(n, -1, dtype=int)

        dx = coords[:, 0] - coords[start, 0]
        dy = coords[:, 1] - coords[start, 1]
        min_dist = dx * dx + dy * dy
        min_dist[start] = np.inf
        parent[:] = start
        parent[start] = -1

        edges: list[tuple[int, int]] = []
        for _ in range(n - 1):
            candidates = np.where(~visited)[0]
            if candidates.size == 0:
                break
            idx = candidates[np.argmin(min_dist[candidates])]
            if not np.isfinite(min_dist[idx]):
                break
            visited[idx] = True
            if parent[idx] >= 0:
                edges.append((int(parent[idx]), int(idx)))

            dx = coords[:, 0] - coords[idx, 0]
            dy = coords[:, 1] - coords[idx, 1]
            dist2 = dx * dx + dy * dy
            for j in np.where(~visited)[0]:
                if dist2[j] < min_dist[j]:
                    min_dist[j] = dist2[j]
                    parent[j] = idx
            min_dist[idx] = np.inf
        return edges

    def draw_line_connectors(base: pd.DataFrame, method: str) -> None:
        for _, g in base.groupby("line_name", dropna=True):
            if len(g) < 2:
                continue
            xs = g["lng_round"].to_numpy()
            ys = g["lat_round"].to_numpy()
            if method == "pca":
                order = pca_order_2d(xs, ys)
                xs2 = xs[order]
                ys2 = ys[order]
                plt.plot(xs2, ys2, color="#7a7a7a", linewidth=0.6, alpha=0.35, zorder=1)
            elif method == "mst":
                edges = mst_edges_2d(xs, ys)
                for u, v in edges:
                    plt.plot(
                        [xs[u], xs[v]],
                        [ys[u], ys[v]],
                        color="#7a7a7a",
                        linewidth=0.7,
                        alpha=0.45,
                        zorder=1,
                    )

    def spatial_plot(value_col: str, title: str, outname: str, connect_method: str | None = None) -> None:
        plt.figure(figsize=(6.5, 5.2))
        if connect_method is not None:
            base = geo.dropna(subset=["lng_round", "lat_round", "line_name"]).copy()
            if len(base):
                draw_line_connectors(base, method=connect_method)
        sc = plt.scatter(
            geo["lng_round"],
            geo["lat_round"],
            c=geo[value_col],
            cmap=spl_cmap,
            s=28,
            alpha=0.9,
            edgecolors="none",
            zorder=2,
            vmin=vmin,
            vmax=vmax,
        )
        plt.xlabel("Longitude (rounded)")
        plt.ylabel("Latitude (rounded)")
        plt.title(title)
        cb = plt.colorbar(sc)
        cb.set_label("dB(A)")
        savefig(outdir / outname, dpi=args.dpi)

    spatial_plot(metric_col, "Spatial pattern of LAeq_24h (Leq)", "03_descriptive__spatial_laeq24h.png", connect_method="mst")
    if "lday_leq_db_a" in geo.columns:
        spatial_plot("lday_leq_db_a", "Spatial pattern of Lday (Leq)", "03_descriptive__spatial_lday.png", connect_method="mst")
    if "lnight_leq_db_a" in geo.columns:
        spatial_plot("lnight_leq_db_a", "Spatial pattern of Lnight (Leq)", "03_descriptive__spatial_lnight.png", connect_method="mst")
    if "lpeak_leq_db_a" in geo.columns:
        spatial_plot("lpeak_leq_db_a", "Spatial pattern of LPeak (Leq)", "03_descriptive__spatial_lpeak.png", connect_method="mst")
    if "linter_peak_leq_db_a" in geo.columns:
        spatial_plot(
            "linter_peak_leq_db_a",
            "Spatial pattern of LInter-peak (Leq)",
            "03_descriptive__spatial_linter_peak.png",
            connect_method="mst",
        )
    if "mp_leq_db_a" in geo.columns:
        spatial_plot(
            "mp_leq_db_a",
            "Spatial pattern of morning peak (MP) Leq",
            "03_descriptive__spatial_mp.png",
            connect_method="mst",
        )
    if "ep_leq_db_a" in geo.columns:
        spatial_plot(
            "ep_leq_db_a",
            "Spatial pattern of evening peak (EP) Leq",
            "03_descriptive__spatial_ep.png",
            connect_method="mst",
        )

    print("Step 2 outputs written to:", outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
