from __future__ import annotations

import argparse
import math
import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "mplconfig"))

TOPIC_ROOT = Path(__file__).resolve().parents[2]  # repository root
DATA_CLEAN_DIR = TOPIC_ROOT / "intermediate" / "clean"
OUTPUT_DIR = TOPIC_ROOT / "intermediate" / "analysis"


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Step 7: sound source composition analysis (mechanism / soundscape).")
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
        help="Noise metrics CSV (default: data/clean/noise_metrics_station.csv).",
    )
    p.add_argument(
        "--out-dir",
        type=str,
        default=str(OUTPUT_DIR),
        help="Output directory (default: output/).",
    )
    p.add_argument(
        "--min-dominant-group-size",
        type=int,
        default=5,
        help="Minimum station count for a dominant-event group to be kept (default: 5).",
    )
    p.add_argument(
        "--top-dominant-events",
        type=int,
        default=12,
        help="Top-N dominant events to show in bar plot (default: 12).",
    )
    p.add_argument(
        "--epsilon",
        type=float,
        default=1e-6,
        help="Small epsilon for CLR/ILR-style log transforms when shares include zeros (default: 1e-6).",
    )
    p.add_argument(
        "--include-controls",
        action="store_true",
        help="Include controls (line_name, station_type_code, day_type) in compositional regressions.",
    )
    return p


def _spearmanr(x: pd.Series, y: pd.Series) -> tuple[float, float]:
    xs = pd.to_numeric(x, errors="coerce")
    ys = pd.to_numeric(y, errors="coerce")
    mask = xs.notna() & ys.notna()
    if mask.sum() < 3:
        return math.nan, math.nan
    try:
        from scipy.stats import spearmanr  # type: ignore

        res = spearmanr(xs[mask].to_numpy(), ys[mask].to_numpy())
        return float(res.correlation), float(res.pvalue)
    except Exception:
        xr = xs[mask].rank(method="average")
        yr = ys[mask].rank(method="average")
        rho = float(np.corrcoef(xr, yr)[0, 1])
        return rho, math.nan


def safe_log1p(series: pd.Series) -> pd.Series:
    s = pd.to_numeric(series, errors="coerce")
    if (s < 0).any():
        return s
    return np.log1p(s)


def zscore(series: pd.Series) -> pd.Series:
    s = pd.to_numeric(series, errors="coerce")
    sd = float(s.std(ddof=0))
    if not np.isfinite(sd) or sd == 0:
        return s * 0.0
    return (s - float(s.mean())) / sd


def write_share_summary(df: pd.DataFrame, share_cols: list[str], out_path: Path) -> None:
    rows: list[dict] = []
    for col in share_cols:
        s = pd.to_numeric(df[col], errors="coerce")
        rows.append(
            {
                "share": col,
                "n": int(s.notna().sum()),
                "mean": float(s.mean(skipna=True)),
                "sd": float(s.std(skipna=True, ddof=1)),
                "median": float(s.median(skipna=True)),
                "q25": float(s.quantile(0.25)),
                "q75": float(s.quantile(0.75)),
                "min": float(s.min(skipna=True)),
                "max": float(s.max(skipna=True)),
            }
        )
    pd.DataFrame(rows).to_csv(out_path, index=False)


def make_plots(df: pd.DataFrame, share_cols: list[str], out_dir: Path, top_n: int) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme(style="whitegrid", palette="colorblind")

    # Share distributions.
    long = df[["station_key"] + share_cols].melt(id_vars=["station_key"], var_name="source", value_name="share")
    fig, ax = plt.subplots(figsize=(8.5, 4.2), dpi=300)
    sns.violinplot(data=long, x="source", y="share", inner=None, cut=0, ax=ax)
    sns.boxplot(data=long, x="source", y="share", width=0.25, showfliers=False, ax=ax, color="white")
    ax.set_title("Sound source share distributions (station-level)")
    ax.set_xlabel("Sound source")
    ax.set_ylabel("Share (0–1)")
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(out_dir / "08_sound_source__share_distributions.png", bbox_inches="tight")
    plt.close(fig)

    # Dominant event counts.
    counts_cn = (
        df["dominant_event_cn"]
        .fillna("Unknown")
        .astype(str)
        .value_counts()
        .rename_axis("dominant_event_cn")
        .reset_index(name="station_count")
    )
    counts_cn["dominant_event_en"] = counts_cn["dominant_event_cn"].map(dominant_event_cn_to_en).fillna("Other")
    counts_cn.to_csv(out_dir / "08_sound_source__dominant_event_counts.csv", index=False)

    top = counts_cn.head(top_n).copy()
    fig, ax = plt.subplots(figsize=(9.5, 4.5), dpi=300)
    sns.barplot(data=top, y="dominant_event_en", x="station_count", ax=ax, color="#4C72B0")
    ax.set_title(f"Top dominant sound events (Top {top_n})")
    ax.set_xlabel("Number of stations")
    ax.set_ylabel("Dominant event")
    fig.tight_layout()
    fig.savefig(out_dir / "08_sound_source__dominant_event_counts.png", bbox_inches="tight")
    plt.close(fig)


def correlation_table(df: pd.DataFrame, shares: list[str], outcomes: list[str], out_path: Path, heatmap_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    rows = []
    for y in outcomes:
        for x in shares:
            rho, p = _spearmanr(df[x], df[y])
            rows.append({"outcome": y, "predictor": x, "spearman_rho": rho, "p_value": p})
    out = pd.DataFrame(rows)
    out.to_csv(out_path, index=False)

    # Heatmap (rho only).
    mat = out.pivot(index="outcome", columns="predictor", values="spearman_rho")
    sns.set_theme(style="whitegrid", palette="colorblind")
    fig, ax = plt.subplots(figsize=(7.8, 3.6), dpi=300)
    sns.heatmap(mat, annot=True, fmt=".2f", cmap="vlag", center=0.0, ax=ax)
    ax.set_title("Spearman correlation: noise metrics vs sound source shares")
    ax.set_xlabel("Sound source share")
    ax.set_ylabel("Noise metric")
    fig.tight_layout()
    fig.savefig(heatmap_path, bbox_inches="tight")
    plt.close(fig)


def scatter_noise_vs_shares(df: pd.DataFrame, shares: list[str], y: str, out_dir: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme(style="whitegrid", palette="colorblind")
    for x in shares:
        rho, p = _spearmanr(df[x], df[y])
        fig, ax = plt.subplots(figsize=(5.8, 4.2), dpi=300)
        sns.regplot(
            data=df,
            x=x,
            y=y,
            scatter_kws={"s": 28, "alpha": 0.75},
            line_kws={"color": "#DD8452", "linewidth": 2},
            lowess=True,
            ax=ax,
        )
        ax.set_title(f"{y} vs {x}\nSpearman rho={rho:.2f}, p={p:.3g}" if np.isfinite(p) else f"{y} vs {x}\nSpearman rho={rho:.2f}")
        ax.set_xlabel(x)
        ax.set_ylabel(y)
        fig.tight_layout()
        fig.savefig(out_dir / f"08_sound_source__scatter__{y}__vs__{x}.png", bbox_inches="tight")
        plt.close(fig)


def dominant_event_grouping(series: pd.Series, min_size: int) -> pd.Series:
    s = series.fillna("Unknown").astype(str)
    counts = s.value_counts()
    keep = set(counts.loc[counts >= min_size].index.tolist())
    return s.where(s.isin(keep), other="Other/Rare")


def dominant_event_cn_to_en(label: str) -> str:
    mapping = {
        "汽车行驶声": "Vehicle pass-by",
        "刹车声": "Braking",
        "交谈声": "Conversation",
        "虫鸣声": "Insects",
        "鸟鸣声": "Birds",
        "Unknown": "Unknown",
    }
    return mapping.get(str(label), "Other")



def kw_by_group(df: pd.DataFrame, group_col: str, y_col: str) -> tuple[float, float, int]:
    g = df[[group_col, y_col]].copy()
    g[y_col] = pd.to_numeric(g[y_col], errors="coerce")
    g = g.dropna(subset=[y_col])
    groups = [sub[y_col].to_numpy() for _, sub in g.groupby(group_col)]
    k = len(groups)
    if k < 2:
        return math.nan, math.nan, k
    try:
        from scipy.stats import kruskal  # type: ignore

        res = kruskal(*groups)
        return float(res.statistic), float(res.pvalue), k
    except Exception:
        return math.nan, math.nan, k


def mwu_pairwise_holm(df: pd.DataFrame, group_col: str, y_col: str, out_path: Path) -> None:
    try:
        from scipy.stats import mannwhitneyu  # type: ignore
    except Exception:
        pd.DataFrame([]).to_csv(out_path, index=False)
        return

    work = df[[group_col, y_col]].copy()
    work[y_col] = pd.to_numeric(work[y_col], errors="coerce")
    work = work.dropna(subset=[y_col])
    groups = sorted(work[group_col].astype(str).unique().tolist())
    rows: list[dict] = []
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            a = groups[i]
            b = groups[j]
            xa = work.loc[work[group_col] == a, y_col].to_numpy()
            xb = work.loc[work[group_col] == b, y_col].to_numpy()
            if len(xa) < 3 or len(xb) < 3:
                continue
            res = mannwhitneyu(xa, xb, alternative="two-sided")
            rows.append({"group_a": a, "group_b": b, "u_stat": float(res.statistic), "p_value": float(res.pvalue)})

    if not rows:
        pd.DataFrame([]).to_csv(out_path, index=False)
        return

    out = pd.DataFrame(rows).sort_values("p_value", ignore_index=True)
    # Holm correction.
    m = len(out)
    holm = []
    for rank, p in enumerate(out["p_value"].to_numpy(), start=1):
        holm.append(min(1.0, float(p) * (m - rank + 1)))
    out["p_holm"] = holm
    out.to_csv(out_path, index=False)


def boxplot_by_dominant_event(df: pd.DataFrame, group_col: str, y_col: str, out_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme(style="whitegrid", palette="colorblind")
    work = df[[group_col, y_col]].copy()
    work[y_col] = pd.to_numeric(work[y_col], errors="coerce")
    work = work.dropna(subset=[y_col])
    order = (
        work.groupby(group_col)[y_col]
        .median()
        .sort_values(ascending=False)
        .index.astype(str)
        .tolist()
    )
    fig, ax = plt.subplots(figsize=(10.8, 4.8), dpi=300)
    sns.boxplot(data=work, x=group_col, y=y_col, order=order, showfliers=False, ax=ax)
    sns.stripplot(data=work, x=group_col, y=y_col, order=order, color="black", alpha=0.35, size=3, ax=ax)
    ax.set_title(f"{y_col} by dominant event group")
    ax.set_xlabel("Dominant event group")
    ax.set_ylabel(y_col)
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def compositional_regressions(
    df: pd.DataFrame,
    y_cols: list[str],
    share_cols: list[str],
    controls: list[str],
    eps: float,
    out_dir: Path,
) -> None:
    import statsmodels.api as sm

    # Baseline approach: drop one share as reference to avoid perfect collinearity.
    baseline = "sound_natural_share" if "sound_natural_share" in share_cols else share_cols[-1]
    use_shares = [c for c in share_cols if c != baseline]

    def fit_one(y_col: str, X: pd.DataFrame) -> pd.Series:
        y = pd.to_numeric(df[y_col], errors="coerce")
        Xf = X.copy()
        Xf = sm.add_constant(Xf, has_constant="add").apply(pd.to_numeric, errors="coerce").astype(float)
        mask = y.notna() & Xf.notna().all(axis=1)
        if mask.sum() < max(20, Xf.shape[1] + 5):
            return pd.Series(dtype=float)
        model = sm.OLS(y.loc[mask], Xf.loc[mask]).fit(cov_type="HC3")
        out = model.params.rename("coef").to_frame()
        out["se_hc3"] = model.bse
        out["p_value"] = model.pvalues
        out["outcome"] = y_col
        out["model"] = "baseline_drop_one"
        return out.reset_index(names="term")

    rows: list[pd.DataFrame] = []
    baseX = df[use_shares].apply(pd.to_numeric, errors="coerce")
    if controls:
        cats = df[controls].copy()
        for c in controls:
            cats[c] = cats[c].astype("string")
        dummies = pd.get_dummies(cats, drop_first=True, prefix_sep="=")
        baseX = pd.concat([baseX, dummies], axis=1)
    for y in y_cols:
        rows.append(fit_one(y, baseX))
    out = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame([])
    out.to_csv(out_dir / "08_sound_source__compositional_regression_baseline.csv", index=False)

    # CLR approach: clr(x_i) = log(x_i / gmean(x)) where x are shares.
    share_mat = df[share_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).clip(lower=0.0)
    share_mat = share_mat + float(eps)
    gm = np.exp(np.log(share_mat).mean(axis=1))
    clr = np.log(share_mat.div(gm, axis=0))
    clr = clr.add_prefix("clr__")
    Xclr = clr.copy()
    if controls:
        cats = df[controls].copy()
        for c in controls:
            cats[c] = cats[c].astype("string")
        dummies = pd.get_dummies(cats, drop_first=True, prefix_sep="=")
        Xclr = pd.concat([Xclr, dummies], axis=1)
    rows2: list[pd.DataFrame] = []
    for y in y_cols:
        rows2.append(fit_one(y, Xclr))
    out2 = pd.concat(rows2, ignore_index=True) if rows2 else pd.DataFrame([])
    out2["model"] = "clr"
    out2.to_csv(out_dir / "08_sound_source__compositional_regression_clr.csv", index=False)


def main() -> None:
    args = build_parser().parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    station = pd.read_csv(args.station_master)
    noise = pd.read_csv(args.noise_metrics)
    df = station.merge(noise, on="station_key", how="left", suffixes=("", "__noise"))

    share_cols = [
        "sound_traffic_share",
        "sound_human_share",
        "sound_mechanical_share",
        "sound_natural_share",
    ]
    unknown_col = "sound_unknown_share" if "sound_unknown_share" in df.columns else None
    for c in share_cols + ([unknown_col] if unknown_col else []):
        if c and c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    # Share sum checks.
    df["share_sum_known"] = df[share_cols].sum(axis=1, skipna=False)
    if unknown_col:
        df["share_sum_all"] = df[share_cols + [unknown_col]].sum(axis=1, skipna=False)
    else:
        df["share_sum_all"] = df["share_sum_known"]

    checks = df[["station_key", "share_sum_known", "share_sum_all"]].copy()
    checks["flag_known_not_1"] = (checks["share_sum_known"] - 1.0).abs() > 0.02
    checks["flag_all_not_1"] = (checks["share_sum_all"] - 1.0).abs() > 0.02
    checks.to_csv(out_dir / "08_sound_source__share_sum_checks.csv", index=False)
    df = df.copy()

    # Descriptive summary.
    write_share_summary(df, share_cols + ([unknown_col] if unknown_col else []), out_dir / "08_sound_source__share_summary.csv")

    # Plots: distributions + dominant event frequency.
    if "dominant_event_cn" not in df.columns:
        df["dominant_event_cn"] = pd.NA
    df["dominant_event_en"] = df["dominant_event_cn"].fillna("Unknown").astype(str).map(dominant_event_cn_to_en)
    make_plots(df, share_cols + ([unknown_col] if unknown_col else []), out_dir=out_dir, top_n=args.top_dominant_events)

    # Noise outcomes for correlation/scatter.
    outcomes = [
        c
        for c in [
            "laeq_24h_leq_db_a",
            "lday_leq_db_a",
            "lnight_leq_db_a",
            "lpeak_leq_db_a",
            "linter_peak_leq_db_a",
            "laeq_sd_24h_db_a",
            "laeq_peak_trough_diff_db",
        ]
        if c in df.columns
    ]
    correlation_table(
        df,
        shares=share_cols,
        outcomes=outcomes,
        out_path=out_dir / "08_sound_source__corr_spearman.csv",
        heatmap_path=out_dir / "08_sound_source__corr_heatmap.png",
    )
    if "laeq_24h_leq_db_a" in df.columns:
        scatter_noise_vs_shares(df, shares=share_cols, y="laeq_24h_leq_db_a", out_dir=out_dir)

    # Dominant-event group differences.
    df["dominant_event_group_cn"] = dominant_event_grouping(df["dominant_event_cn"], min_size=args.min_dominant_group_size)
    df["dominant_event_group_en"] = df["dominant_event_group_cn"].astype(str).map(dominant_event_cn_to_en)
    kw_rows: list[dict] = []
    for y in ["laeq_24h_leq_db_a", "lday_leq_db_a", "lnight_leq_db_a", "lpeak_leq_db_a", "linter_peak_leq_db_a"]:
        if y not in df.columns:
            continue
        stat, p, k = kw_by_group(df, "dominant_event_group_en", y)
        kw_rows.append({"outcome": y, "n_groups": k, "kw_stat": stat, "p_value": p})
        boxplot_by_dominant_event(
            df,
            group_col="dominant_event_group_en",
            y_col=y,
            out_path=out_dir / f"08_sound_source__boxplot__{y}__by_dominant_event.png",
        )
        mwu_pairwise_holm(
            df,
            group_col="dominant_event_group_en",
            y_col=y,
            out_path=out_dir / f"08_sound_source__posthoc_holm__{y}.csv",
        )
    pd.DataFrame(kw_rows).to_csv(out_dir / "08_sound_source__kw_dominant_event.csv", index=False)

    # Compositional regression (baseline drop-one + CLR).
    controls = ["line_name", "station_type_code", "day_type"] if args.include_controls else []
    y_cols = [
        c
        for c in ["laeq_24h_leq_db_a", "lday_leq_db_a", "lnight_leq_db_a", "lpeak_leq_db_a", "linter_peak_leq_db_a"]
        if c in df.columns
    ]
    compositional_regressions(
        df=df,
        y_cols=y_cols,
        share_cols=share_cols,
        controls=controls,
        eps=float(args.epsilon),
        out_dir=out_dir,
    )


if __name__ == "__main__":
    main()
