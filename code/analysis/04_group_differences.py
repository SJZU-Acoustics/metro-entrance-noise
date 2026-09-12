from __future__ import annotations

import argparse
import itertools
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

TOPIC_ROOT = Path(__file__).resolve().parents[2]  # repository root
DATA_CLEAN_DIR = TOPIC_ROOT / "intermediate" / "clean"
OUTPUT_DIR = TOPIC_ROOT / "intermediate" / "analysis"

STATION_TYPE_CODE_TO_CN = {
    "1": "普通站",
    "2": "换乘站",
    "3": "终点站",
}


DEFAULT_METRICS = [
    "laeq_24h_leq_db_a",
    "lday_leq_db_a",
    "lnight_leq_db_a",
    "lpeak_leq_db_a",
    "linter_peak_leq_db_a",
    "lden_db",
    "laeq_sd_24h_db_a",
    "laeq_peak_trough_diff_db",
]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Step 3: group differences by line and station type (Kruskal–Wallis + Holm posthoc + effect sizes)."
    )
    parser.add_argument(
        "--station-master",
        type=str,
        default=str(DATA_CLEAN_DIR / "station_master.csv"),
        help="Station-level master table.",
    )
    parser.add_argument(
        "--noise",
        type=str,
        default=str(DATA_CLEAN_DIR / "noise_metrics_station.csv"),
        help="Station-level noise metrics table.",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default=str(OUTPUT_DIR),
        help="Output directory.",
    )
    parser.add_argument(
        "--min-group-n",
        type=int,
        default=5,
        help="Minimum N per group to include in tests.",
    )
    parser.add_argument(
        "--metrics",
        type=str,
        default=",".join(DEFAULT_METRICS),
        help="Comma-separated metric columns to test.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="Alpha for reporting significance (used in report; tests always run).",
    )
    return parser


def ensure_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing: {path}")


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    # Efficient O((n+m)log m) using searchsorted on sorted y.
    x = x.astype(float)
    y = y.astype(float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if len(x) == 0 or len(y) == 0:
        return math.nan
    y_sorted = np.sort(y)
    # For each x: count less (y < x) and greater (y > x)
    left = np.searchsorted(y_sorted, x, side="left")
    right = np.searchsorted(y_sorted, x, side="right")
    less = left.sum()
    greater = (len(y_sorted) - right).sum()
    return float((greater - less) / (len(x) * len(y_sorted)))


def cliffs_magnitude(d: float) -> str:
    if not np.isfinite(d):
        return ""
    ad = abs(d)
    if ad < 0.147:
        return "negligible"
    if ad < 0.33:
        return "small"
    if ad < 0.474:
        return "medium"
    return "large"


def epsilon_squared_kw(h: float, k: int, n: int) -> float:
    # epsilon^2 = (H - k + 1) / (n - k)
    if not np.isfinite(h) or n <= k:
        return math.nan
    return float((h - k + 1.0) / (n - k))


def station_type_label_cn(code_like: str) -> str:
    s = str(code_like).strip()
    if s == "" or s.lower() == "nan":
        return s
    try:
        s = str(int(float(s)))
    except Exception:
        pass
    return STATION_TYPE_CODE_TO_CN.get(s, s)


def station_type_codebook_cn_text() -> str:
    return "1=普通站;2=换乘站;3=终点站"


def describe_by_group(df: pd.DataFrame, group_col: str, metric: str) -> pd.DataFrame:
    tmp = df[[group_col, metric]].copy()
    tmp[metric] = pd.to_numeric(tmp[metric], errors="coerce")
    out = (
        tmp.groupby(group_col, as_index=False)[metric]
        .agg(
            n="count",
            mean="mean",
            median="median",
            p25=lambda s: float(s.dropna().quantile(0.25)) if s.notna().any() else math.nan,
            p75=lambda s: float(s.dropna().quantile(0.75)) if s.notna().any() else math.nan,
            min="min",
            max="max",
        )
        .sort_values("n", ascending=False, ignore_index=True)
    )
    out = out.rename(columns={group_col: "group"})
    out["grouping"] = group_col
    out["metric"] = metric
    if group_col == "station_type_code":
        out["group_label_cn"] = out["group"].map(station_type_label_cn)
    else:
        out["group_label_cn"] = out["group"].astype(str)
    out["station_type_codebook_cn"] = station_type_codebook_cn_text()
    return out[
        [
            "grouping",
            "metric",
            "group",
            "group_label_cn",
            "n",
            "mean",
            "median",
            "p25",
            "p75",
            "min",
            "max",
            "station_type_codebook_cn",
        ]
    ]


def run_kw_and_posthoc(df: pd.DataFrame, group_col: str, metric: str, min_group_n: int) -> tuple[dict, pd.DataFrame]:
    tmp = df[[group_col, metric]].copy()
    tmp[metric] = pd.to_numeric(tmp[metric], errors="coerce")
    tmp = tmp.dropna(subset=[group_col, metric])

    group_sizes = tmp.groupby(group_col)[metric].size()
    keep_groups = group_sizes[group_sizes >= min_group_n].index.tolist()
    tmp = tmp[tmp[group_col].isin(keep_groups)].copy()
    k = tmp[group_col].nunique()
    n = len(tmp)

    kw_row = {
        "grouping": group_col,
        "metric": metric,
        "n_total": int(n),
        "k_groups": int(k),
        "h_stat": math.nan,
        "p_value": math.nan,
        "epsilon2": math.nan,
        "min_group_n": int(min_group_n),
        "station_type_codebook_cn": station_type_codebook_cn_text(),
    }

    if k < 2:
        return kw_row, pd.DataFrame(
            columns=[
                "grouping",
                "metric",
                "group_a",
                "group_b",
                "group_a_label_cn",
                "group_b_label_cn",
                "n_a",
                "n_b",
                "u_stat",
                "p_raw",
                "p_holm",
                "cliffs_delta",
                "delta_magnitude",
                "station_type_codebook_cn",
            ]
        )

    arrays = [tmp.loc[tmp[group_col] == g, metric].to_numpy(dtype=float) for g in sorted(keep_groups)]
    h, p = stats.kruskal(*arrays, nan_policy="omit")
    kw_row["h_stat"] = float(h)
    kw_row["p_value"] = float(p)
    kw_row["epsilon2"] = epsilon_squared_kw(float(h), k=int(k), n=int(n))

    # Pairwise posthoc: Mann–Whitney U + Holm correction, with Cliff's delta.
    pairs = list(itertools.combinations(sorted(keep_groups), 2))
    post_rows = []
    pvals = []
    for a, b in pairs:
        xa = tmp.loc[tmp[group_col] == a, metric].to_numpy(dtype=float)
        xb = tmp.loc[tmp[group_col] == b, metric].to_numpy(dtype=float)
        xa = xa[np.isfinite(xa)]
        xb = xb[np.isfinite(xb)]
        if len(xa) == 0 or len(xb) == 0:
            u = math.nan
            pr = math.nan
        else:
            res = stats.mannwhitneyu(xa, xb, alternative="two-sided", method="auto")
            u = float(res.statistic)
            pr = float(res.pvalue)
        pvals.append(pr)
        post_rows.append(
            {
                "grouping": group_col,
                "metric": metric,
                "group_a": str(a),
                "group_b": str(b),
                "group_a_label_cn": station_type_label_cn(a) if group_col == "station_type_code" else str(a),
                "group_b_label_cn": station_type_label_cn(b) if group_col == "station_type_code" else str(b),
                "n_a": int(len(xa)),
                "n_b": int(len(xb)),
                "u_stat": u,
                "p_raw": pr,
                "p_holm": math.nan,
                "cliffs_delta": cliffs_delta(xa, xb),
                "delta_magnitude": "",
                "station_type_codebook_cn": station_type_codebook_cn_text(),
            }
        )

    # Holm correction, ignoring nan p-values.
    pvals_arr = np.array(pvals, dtype=float)
    valid = np.isfinite(pvals_arr)
    if valid.any():
        _, p_adj, _, _ = multipletests(pvals_arr[valid], method="holm")
        j = 0
        for i in range(len(post_rows)):
            if valid[i]:
                post_rows[i]["p_holm"] = float(p_adj[j])
                j += 1

    for r in post_rows:
        r["delta_magnitude"] = cliffs_magnitude(float(r["cliffs_delta"])) if np.isfinite(r["cliffs_delta"]) else ""

    posthoc = pd.DataFrame(post_rows).sort_values(["p_holm", "p_raw"], ascending=[True, True], ignore_index=True)
    return kw_row, posthoc


def main() -> int:
    args = build_parser().parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    station_master_path = Path(args.station_master)
    noise_path = Path(args.noise)
    ensure_exists(station_master_path)
    ensure_exists(noise_path)

    station = pd.read_csv(station_master_path)
    noise = pd.read_csv(noise_path)

    # Prefer noise metrics columns when overlap.
    overlap = set(station.columns).intersection(set(noise.columns)) - {"station_key"}
    if overlap:
        station = station.drop(columns=sorted(overlap))

    df = station.merge(noise, on="station_key", how="left", validate="one_to_one")

    metrics = [m.strip() for m in args.metrics.split(",") if m.strip()]
    metrics = [m for m in metrics if m in df.columns]
    if not metrics:
        raise ValueError("No valid metrics found. Check --metrics and available columns.")

    groupings = [
        ("line_name", "line"),
        ("station_type_code", "station_type"),
    ]

    kw_rows: list[dict] = []
    posthoc_frames: list[pd.DataFrame] = []
    desc_frames: list[pd.DataFrame] = []

    for group_col, slug in groupings:
        if group_col not in df.columns:
            continue
        df[group_col] = df[group_col].astype("string")
        for metric in metrics:
            desc_frames.append(describe_by_group(df, group_col=group_col, metric=metric))
            kw_row, posthoc = run_kw_and_posthoc(
                df=df, group_col=group_col, metric=metric, min_group_n=args.min_group_n
            )
            kw_rows.append(kw_row)
            posthoc_frames.append(posthoc)

    desc_df = pd.concat(desc_frames, ignore_index=True)
    desc_df.to_csv(outdir / "04_groupdiff__descriptives.csv", index=False, encoding="utf-8")
    if "group_label_cn" in desc_df.columns:
        desc_zh = desc_df.copy()
        desc_zh["group_display"] = np.where(
            desc_zh["grouping"] == "station_type_code",
            desc_zh["group"].astype(str) + " (" + desc_zh["group_label_cn"].astype(str) + ")",
            desc_zh["group"].astype(str),
        )
        desc_zh.to_csv(outdir / "04_groupdiff__descriptives_zh.csv", index=False, encoding="utf-8")

    kw_df = pd.DataFrame(kw_rows).sort_values(["grouping", "p_value"], ascending=[True, True], ignore_index=True)
    kw_df.to_csv(outdir / "04_groupdiff__kruskal_wallis.csv", index=False, encoding="utf-8")

    posthoc_df = pd.concat(posthoc_frames, ignore_index=True) if posthoc_frames else pd.DataFrame()
    posthoc_df.to_csv(outdir / "04_groupdiff__posthoc_holm.csv", index=False, encoding="utf-8")
    if not posthoc_df.empty and {"group_a_label_cn", "group_b_label_cn"}.issubset(set(posthoc_df.columns)):
        posthoc_zh = posthoc_df.copy()
        posthoc_zh["group_a_display"] = np.where(
            posthoc_zh["grouping"] == "station_type_code",
            posthoc_zh["group_a"].astype(str) + " (" + posthoc_zh["group_a_label_cn"].astype(str) + ")",
            posthoc_zh["group_a"].astype(str),
        )
        posthoc_zh["group_b_display"] = np.where(
            posthoc_zh["grouping"] == "station_type_code",
            posthoc_zh["group_b"].astype(str) + " (" + posthoc_zh["group_b_label_cn"].astype(str) + ")",
            posthoc_zh["group_b"].astype(str),
        )
        posthoc_zh.to_csv(outdir / "04_groupdiff__posthoc_holm_zh.csv", index=False, encoding="utf-8")

    # Simple console summary.
    sig = kw_df.loc[np.isfinite(kw_df["p_value"]) & (kw_df["p_value"] < args.alpha)]
    print("Group differences summary:")
    print(f"- metrics tested: {len(metrics)}")
    print(f"- significant KW (p<{args.alpha}): {len(sig)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
